/**
 * This file is part of RRD-SLAM.
 *
 * Copyright (C) 2016 Jae-Hak Kim <jaehak.kim@adelaide.edu.au>
 *
 * RRD-SLAM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * RRD-SLAM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with RRD-SLAM. If not, see <http://www.gnu.org/licenses/>.
 */

//
//  DepthMapForRadialRolling.cpp
//
//  Created by Jae-Hak Kim on 8/06/2016.
//


#include "DepthMapForRadialRolling.h"

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <opencv2/imgproc/imgproc.hpp>

#include "util/settings.h"
#include "DepthEstimation/DepthMapPixelHypothesis.h"
#include "DataStructures/Frame.h"
#include "util/globalFuncs.h"
#include "IOWrapper/ImageDisplay.h"
#include "GlobalMapping/KeyFrameGraph.h"

#include "opencv2/opencv.hpp"

#include "ceres/ceres.h"
#include "Tracking/JetOps.h"

#ifdef __APPLE__
#define isnanf isnan
#endif

namespace lsd_slam
{

//------------------------------------------------------------------------------
std::list<PointEpiCurve>
DepthMapForRadialRolling::
addPtsToEpiCurveForNoRadial(double src_x_undist,
                            double src_y_undist,
                            Sophus::SE3d Xi_ySource_to_rowTarget,
                            Eigen::Matrix<double,3,1> eline,
                            double epx,
                            double epy,
                            double aRow_dist,
                            RunningStats* const stats)
{
    std::list<PointEpiCurve> ptsEpiCurve;
    
    Eigen::Matrix<double,3,1> t = Xi_ySource_to_rowTarget.translation();

    // Upper and lower line
    Eigen::Matrix<double,3,1> upper_eline =
        eline + Eigen::Matrix<double,3,1>(0,0,-2*eline[1]);
    Eigen::Matrix<double,3,1> lower_eline =
        eline + Eigen::Matrix<double,3,1>(0,0,+2*eline[1]);
    

    // Intersection with the current rowTarget
    Eigen::Matrix<double,3,1> pInter
        = eline.cross(Eigen::Matrix<double,3,1>(0, 1, -aRow_dist));
    Eigen::Matrix<double,3,1> upper_pInter
        = upper_eline.cross(Eigen::Matrix<double,3,1>(0, 1, -aRow_dist));
    Eigen::Matrix<double,3,1> lower_pInter
        = lower_eline.cross(Eigen::Matrix<double,3,1>(0, 1, -aRow_dist));
    
    Eigen::Matrix<double,3,1> pInter_h;
    pInter_h << pInter[0]/pInter[2],
    pInter[1]/pInter[2],
    pInter[2]/pInter[2];
    Eigen::Matrix<double,3,1> upper_pInter_h;
    upper_pInter_h << upper_pInter[0]/upper_pInter[2],
    upper_pInter[1]/upper_pInter[2],
    upper_pInter[2]/upper_pInter[2];
    Eigen::Matrix<double,3,1> lower_pInter_h;
    lower_pInter_h << lower_pInter[0]/lower_pInter[2],
    lower_pInter[1]/lower_pInter[2],
    lower_pInter[2]/lower_pInter[2];
    
    // Find a normal to the epipolar line
    double nepx, nepy;
    if (makeAndCheckEPL_forGeneralEpipolarLine((int)round(src_x_undist),
                                               (int)round(src_y_undist),
                                               t,
                                               &nepx, &nepy,
                                               stats) == false)
    {
        
        // No normal of the epipolar line found
        
#if 1 // SKIP IF FAILED IN NORMAL COMPUTATION
        return ptsEpiCurve;
#endif
        
    }
    
    
    // Check points with +/- infinity test
    std::list<PointEpiCurve> partCurve =
    buildCurveForPointAtInfinity(src_x_undist,
                                 src_y_undist,
                                 nepx,
                                 nepy,
                                 Xi_ySource_to_rowTarget,
                                 (int)round(aRow_dist),
                                 upper_pInter,
                                 lower_pInter);
    if (partCurve.size() != 0)
    {
        
        for (std::list<PointEpiCurve>::const_iterator ci
				 = partCurve.begin();
				 ci != partCurve.end();
				 ++ci)
        {
            ptsEpiCurve.push_back(*ci);
        }
        return ptsEpiCurve;
        
    }
    
    // Check boundary
    if (!(pInter_h[0] < 3 || pInter_h[0] > width - 3) &&
       !(pInter_h[1] < 3 || pInter_h[1] > height - 3) &&
       !(upper_pInter_h[0] < 3 || upper_pInter_h[0] > width - 3) &&
       !(upper_pInter_h[1] < 3 || upper_pInter_h[1] > height - 3) &&
       !(lower_pInter_h[0] < 3 || lower_pInter_h[0] > width - 3) &&
       !(lower_pInter_h[1] < 3 || lower_pInter_h[1] > height - 3))
    {
    
        // Add intersections to curve
        pushPointToEpiCurve(src_x_undist, src_y_undist,
                            upper_pInter_h[0], upper_pInter_h[1],
                            nepx, nepy,
                            Xi_ySource_to_rowTarget,
                            &ptsEpiCurve);
        pushPointToEpiCurve(src_x_undist, src_y_undist,
                            lower_pInter_h[0], lower_pInter_h[1],
                            nepx, nepy,
                            Xi_ySource_to_rowTarget,
                            &ptsEpiCurve);
    }
    
    return ptsEpiCurve;
    
}

//------------------------------------------------------------------------------
std::list<PointEpiCurve>
DepthMapForRadialRolling::
addPtsToEpiCurveForRadialRolling(double src_x_undist,
                                 double src_y_undist,
                                 Sophus::SE3d Xi_ySource_to_rowTargetDist,
                                 Eigen::Matrix<double,3,1> eline,
                                 double epx,
                                 double epy,
                                 double tgt_y_dist,
                                 RunningStats* const stats)
{

    std::list<PointEpiCurve> ptsEpiCurve;
    
    Eigen::Matrix<double,3,1> t = Xi_ySource_to_rowTargetDist.translation();
    Eigen::Matrix<double,3,3> R = Xi_ySource_to_rowTargetDist.rotationMatrix();
    
#if 0 // DEBUG
    std::cout << "src_x_undist = " << src_x_undist << std::endl;
    std::cout << "src_y_undist = " << src_y_undist << std::endl;
    std::cout << "eline = " << eline << std::endl;
    std::cout << "epx = " << epx << std::endl;
    std::cout << "epy = " << epy << std::endl;
    std::cout << "tgt_y_dist = " << tgt_y_dist << std::endl;
    std::cout << "t = " << t << std::endl;
    std::cout << "R = " << R << std::endl;
#endif
  
#if 1 // INIT_CHECK (ONE-SIDE)
    
    // Find intersections with the row
    double pInter_x, pInter_y;
    bool found = findInterEpipolar_Optimization(eline,
                                   tgt_y_dist,
                                   src_x_undist,
                                   src_y_undist,
                                   &pInter_x,
                                   &pInter_y);
    if (found == true &&
        checkRescale(src_x_undist,
                     src_y_undist,
                     epx,
                     epy,
                     R,
                     t,
                     0.4,
                     pInter_x,
                     pInter_y) == true)
    {
        
        pushPointToEpiCurve(src_x_undist, src_y_undist,
                            pInter_x, pInter_y,
                            epx, epy,
                            Xi_ySource_to_rowTargetDist,
                            &ptsEpiCurve);
        
    }
    
#else // BOTH-SIDE CHECK
    

    // Find intersections with the row
    // with init_x = source's x-coordinate;
    double pInter_x_left, pInter_y_left;
    bool found_left = findInterEpipolar_Optimization(eline,
                                                tgt_y_dist,
                                                src_x_undist,
                                                src_y_undist,
                                                &pInter_x_left,
                                                &pInter_y_left);
    // Find intersections with the row
    // with init_x = negated source's x-coordinate
    double pInter_x_right, pInter_y_right;
    bool found_right = findInterEpipolar_Optimization(eline,
                                               tgt_y_dist,
                                               width - 1 - src_x_undist,
                                               src_y_undist,
                                               &pInter_x_right,
                                               &pInter_y_right);
    
    if (found_left == true &&
        found_right == true &&
        fabs(pInter_x_left - pInter_x_right) < 1e-2 &&
        fabs(pInter_y_left - pInter_y_right) < 1e-2 &&
        checkRescale(src_x_undist,
                     src_y_undist,
                     epx,
                     epy,
                     R,
                     t,
                     0.4,
                     pInter_x_left,
                     pInter_y_left) == true &&
        checkRescale(src_x_undist,
                     src_y_undist,
                     epx,
                     epy,
                     R,
                     t,
                     0.4,
                     pInter_x_right,
                     pInter_y_right) == true)
    {
        
        pushPointToEpiCurve(src_x_undist, src_y_undist,
                            pInter_x_left, pInter_y_left,
                            epx, epy,
                            Xi_ySource_to_rowTargetDist,
                            &ptsEpiCurve);

    }
    else
    {
    if (found_left == true &&
        checkRescale(src_x_undist,
                     src_y_undist,
                     epx,
                     epy,
                     R,
                     t,
                     0.4,
                     pInter_x_left,
                     pInter_y_left) == true)
    {
        
        pushPointToEpiCurve(src_x_undist, src_y_undist,
                            pInter_x_left, pInter_y_left,
                            epx, epy,
                            Xi_ySource_to_rowTargetDist,
                            &ptsEpiCurve);
        
    }
    
    if (found_right == true &&
        checkRescale(src_x_undist,
                     src_y_undist,
                     epx,
                     epy,
                     R,
                     t,
                     0.4,
                     pInter_x_right,
                     pInter_y_right) == true)
    {
        
        pushPointToEpiCurve(src_x_undist, src_y_undist,
                            pInter_x_right, pInter_y_right,
                            epx, epy,
                            Xi_ySource_to_rowTargetDist,
                            &ptsEpiCurve);
        
    }
    }
    
#endif
    
    return ptsEpiCurve;
    
}

//------------------------------------------------------------------------------
// Extracts image points on the generalised epipolar line band
// for radial-rolling stereo
//------------------------------------------------------------------------------
// It extracts image points on the band of an generalised epipolar line
// from radial-rolling stereo images.
// This version differs from extractPointsOnGeneralEpipolarLine()
// in terms of that input images have both radial and rolling-shutter artifacts.
// It provides a better handling of degenerate cases.
//
//            x : x-coord of point in keyframe (source) contains semi-dense
//                points
//            y : y-coord of point in keyframe (source)
//          ref : reference/currentlyTrackingFrame (target image)
//  ptsEpiCurve : cooridnates of points on epipolar curve
//        stats : Running stats
//        useGT : [bool] Use ground truth motion (Default: false)
//
// returns
//
//      true : success
//     false : fail
//
bool
DepthMapForRadialRolling::
extractPointsOnGeneralEpipolarLineBand(
    const int src_x_undist,
    const int src_y_undist,
    Frame* const ref,
    std::list<PointEpiCurve> *ptsEpiCurve,
    RunningStats* const stats,
    bool useGT)
{
    
    std::list<PointEpiCurve> ptsEpiCurve_front;
    std::list<PointEpiCurve> ptsEpiCurve_back;
    
    Frame *sourceImage = activeKeyFrame;
    Frame *targetImage = ref;
    
    // Distorted pixels for source
    float src_x_dist, src_y_dist;
    sourceImage->undistorter_->distortPoint(src_x_undist, src_y_undist,
                                            &src_x_dist, &src_y_dist);
    
    //--------------------------------------------------------------------------
    // Sweeping rows in the target image - find a plane for the row
    
    // Row y pose in source image
    Eigen::Matrix<double,4,4> world_to_ySource;
    if (useGT == true)
    {
        // Source pose from ground truth
        world_to_ySource =
            sourceImage->getGTPoseAtRow_Distort(src_x_undist, src_y_undist);
        
    }
    else
    {
        
        // Source pose from estimate (radial distortion version)
        world_to_ySource =
            sourceImage->getPoseAtRow_Distort(src_x_undist, src_y_undist);
        
        // Warning if spline is not available
        if (sourceImage->id() != 0 && sourceImage->isSplineValid == false)
        {
            
            printf("Warning: source image %d has no valid spline\n",
                   sourceImage->id());
            
        }
        
    }
    
    // For each row scanline, map the scanline to a curve which is
    // undistortion of the scanline by lens distortion parameter.
    // Then, we sweep each curve in target image to find an epipolar curve
    int searchWindowY =
        (int)round(sourceImage->height_rollingShutter_inputImage*0.1);
    for (float aRow_dist = src_y_dist - searchWindowY;
         aRow_dist <= src_y_dist + searchWindowY;
         aRow_dist = aRow_dist + 1)
//    for (float aRow_dist = 3;
//         aRow_dist <= height - 3;
//         aRow_dist = aRow_dist + 1)
    {
        
        Eigen::Matrix<double,3,Eigen::Dynamic> interPts;
        
        // Check boundary
        if ((aRow_dist < 0) ||
            (aRow_dist > sourceImage->height_rollingShutter_inputImage - 1))
        {
            
            continue;
        
        }
        
        // Get a motion of the row in view 2
        Eigen::Matrix<double,4,4> world_to_rowTargetDist;
        if (useGT == true)
        {
            
            // Target pose from ground truth
            world_to_rowTargetDist =
                targetImage->getGTPoseAtRow((int)round(aRow_dist));
            
        }
        else
        {
            
            // Target pose from estimate as if the targetImage is distorted.
            // We will handle the correct pose later
            world_to_rowTargetDist =
                targetImage->getPoseAtRow(aRow_dist);
            
            // Warning if spline is not available
            if (targetImage->isSplineValid == false)
            {
                
                printf("Warning: target image %d has no valid spline\n",
                       targetImage->id());
                
            }
            
        }
        
        // Pose of the source with respect to the row target in distorted image
        Eigen::Matrix<double,4,4> rowTargetDist_to_ySource
            = world_to_ySource*world_to_rowTargetDist.inverse();
        Sophus::SE3d Xi_rowTargetDist_to_ySource
            = Sophus::SE3d(rowTargetDist_to_ySource);
        Sophus::SE3d Xi_ySource_to_rowTargetDist
            = Xi_rowTargetDist_to_ySource.inverse();
        
        // R and t of the source with respect to the row target in distorted img
        Eigen::Matrix<double,3,1> t = Xi_ySource_to_rowTargetDist.translation();
        Eigen::Matrix<double,3,3> R = Xi_ySource_to_rowTargetDist.rotationMatrix();

        // Fundamental matrix
        Eigen::Matrix<double,3,3> E = skewsym(t)*R;
        Eigen::Matrix<double,3,3> Kd = K.cast<double>();
        Eigen::Matrix<double,3,3> F = Kd.inverse().transpose()*E*Kd.inverse();
        

        // Check normal
        double nepx, nepy;
        if (makeAndCheckEPL_forGeneralEpipolarLine(src_x_undist,
                                                   src_y_undist,
                                                   t,
                                                   &nepx, &nepy,
                                                   stats) == false)
        {
            
            continue;
            
        }
        
        // Global shutter epipolar line
        Eigen::Matrix<double,3,1> eline =
            F*Eigen::Matrix<double,3,1>(src_x_undist,src_y_undist,1);
        
        std::list<PointEpiCurve> partCurve;
        if (targetImage->undistorter_->getInputDistort() == 0.0)
        {
            
            // Find a part of epipolar curve when NO radial distortion exists
            partCurve = addPtsToEpiCurveForNoRadial(src_x_undist,
                                                    src_y_undist,
                                                    Xi_ySource_to_rowTargetDist,
                                                    eline,
                                                    nepx,
                                                    nepy,
                                                    aRow_dist,
                                                    stats);
            
        }
        else
        {
            
            // Find a part of epipolar curve when radial distortion exists
            partCurve = addPtsToEpiCurveForRadialRolling(src_x_undist,
                                                         src_y_undist,
                                                         Xi_ySource_to_rowTargetDist,
                                                         eline,
                                                         nepx,
                                                         nepy,
                                                         aRow_dist,
                                                         stats);
            
        }
        
        // Append the part of curve to form the total epipolar curve
        if (partCurve.size() != 0)
        {
            
            int k = 0;
            for (auto it: partCurve)
            {
                
                if (k % 2 == 0)
                    ptsEpiCurve_front.push_front(it);
                else
                    ptsEpiCurve_back.push_back(it);
                
                k++;
                
            }
            
            continue;
            
        }

        
    } // End for aRow_dist
    
        
    if ((int)ptsEpiCurve_front.size() > 0)
    {

        // Build sample points uniformly distributed on the epi-curve
        build_UniformSamplePoints(&ptsEpiCurve_front, sourceImage, targetImage);
        
    }

    if ((int)ptsEpiCurve_back.size() > 0)
    {
        
        // Build sample points uniformly distributed on the epi-curve
        build_UniformSamplePoints(&ptsEpiCurve_back, sourceImage, targetImage);
        
    }
    
    // Merge
    ptsEpiCurve->insert(ptsEpiCurve->end(),
                        ptsEpiCurve_front.begin(),
                        ptsEpiCurve_front.end());
    ptsEpiCurve->insert(ptsEpiCurve->end(),
                        ptsEpiCurve_back.begin(),
                        ptsEpiCurve_back.end());
    
    if ((int)ptsEpiCurve->size() > 0)
    {
        
        return true;
        
    }
    else
    {
        
        return false;
        
    }
    
}

void
DepthMapForRadialRolling::
build_UniformSamplePoints(std::list<PointEpiCurve> *ptsEpiCurve,
								  Frame *sourceImage,
								  Frame *targetImage)
{

    auto it = ptsEpiCurve->begin(); 
    auto it_prev = it;
    for (;it != ptsEpiCurve->end();it_prev = it, it++)
    {

        // Compare the current and previous target points
        if (it != ptsEpiCurve->begin())
        {

            // Get current target coordinates
            double curr_tgt_x = (*it).xTarget;
            double curr_tgt_y = (*it).yTarget;

            // Get previous target coordinates
            double prev_tgt_x = (*it_prev).xTarget;
            double prev_tgt_y = (*it_prev).yTarget;

            // Distance
            double dist_squared =
                (curr_tgt_x - prev_tgt_x)*(curr_tgt_x - prev_tgt_x) +
                (curr_tgt_y - prev_tgt_y)*(curr_tgt_y - prev_tgt_y);

            // If two points are apart larger than 1.414 pixels
            // insert new samples
            if (dist_squared > 2.0)
            {

                for (double lambda = 0.0;
                     lambda <= 1.0;
                     lambda = lambda + 1.0/(dist_squared/2.0))
                {

                    // Source point
                    double xSource = (*it).xSource;
                    double ySource = (*it).ySource;

                    // Target point and normal by linear combination
                    double xTarget = (1.0 - lambda)*(*it_prev).xTarget +
                                     lambda*(*it).xTarget;
                    double yTarget = (1.0 - lambda)*(*it_prev).yTarget +
                                     lambda*(*it).yTarget;
                    double nepx = (1.0 - lambda)*(*it_prev).nepx +
                                     lambda*(*it).nepx;
                    double nepy = (1.0 - lambda)*(*it_prev).nepy +
                                     lambda*(*it).nepy;
                    
                    // Normalize nepx and nepy
                    double norm_nepx_nepy = sqrt(nepx*nepx + nepy*nepy);
                    float fac = GRADIENT_SAMPLE_DIST / sqrt(norm_nepx_nepy);
                    nepx = nepx*fac;
                    nepy = nepy*fac;
                    
//                    // Normalize nepx and nepy
//                    double len_nep = sqrt(nepx*nepx + nepy+nepy);
//                    if (len_nep > 1e-6)
//                    {
//                        nepx = nepx/len_nep;
//                        nepy = nepy/len_nep;
//                    }
                    
                    // Motion of source wrt to distorted target
                    Sophus::SE3d Xi_world_to_RowTargetDist = 
                        Sophus::SE3d(targetImage->getPoseAtRow_Distort(
                            xTarget, yTarget));
                    Sophus::SE3d Xi_world_to_ySource = 
                        Sophus::SE3d(sourceImage->getPoseAtRow_Distort(
                            xSource, ySource));
                    Sophus::SE3d Xi_ySource_to_rowTargetDist =
                        Xi_world_to_RowTargetDist*Xi_world_to_ySource.inverse();

                    // New point
                    PointEpiCurve pt(xSource,
                                     ySource,
                                     xTarget,
                                     yTarget,
                                     nepx,
                                     nepy,
                                     Xi_ySource_to_rowTargetDist);

                    // Insert the new point before the current point
                    ptsEpiCurve->insert(it, pt);

                } // End for lambda

            } // End If dist

        } // End if it != begin
            
    } // End for auto it

}

bool
DepthMapForRadialRolling::
checkRescale(double u, double v, double nepx, double nepy,
             Eigen::Matrix<double,3,3> R, Eigen::Matrix<double,3,1> t,
             double prior_idepth,
             double u2, double v2)
{

    // Compute point at infinity and point in real
    Eigen::Vector3f KinvP = Eigen::Vector3f(fxi*(float)u + cxi,
                                            fyi*(float)v + cyi,
                                            1.0f);
    Eigen::Vector3f pInf = K*R.cast<float>()*KinvP;
    Eigen::Vector3f pReal = pInf/(float)prior_idepth + K*t.cast<float>();
    
    
    // Check rescale factor
    float rescaleFactor = pReal(2)*(float)prior_idepth;
    if (!(rescaleFactor > 0.7 && rescaleFactor < 1.4))
    {
        
        return false;
        
    }

    // Check normal with rescale factor
    float firstX = (float)(u - 2*nepx*rescaleFactor);
    float firstY = (float)(v - 2*nepy*rescaleFactor);
    float lastX = (float)(u + 2*nepx*rescaleFactor);
    float lastY = (float)(v + 2*nepy*rescaleFactor);
    if (firstX < 0 || firstX > width - 1 ||
        firstY < 0 || firstY > height - 1 ||
        lastX < 0 || lastX > width - 1 ||
        lastY < 0 ||lastY > height - 1)
    {
        // Out of image boundary
        return false;
    
    }

    // Check the boundary in the source
    if (u - 2*nepx*rescaleFactor < 0 ||
        u - 2*nepx*rescaleFactor > width - 1 ||
        u + 2*nepx*rescaleFactor < 0 ||
        u + 2*nepx*rescaleFactor > width - 1 ||
        v - 2*nepy*rescaleFactor < 0 ||
        v - 2*nepy*rescaleFactor > height - 1 ||
        v + 2*nepy*rescaleFactor < 0 ||
        v + 2*nepy*rescaleFactor > height - 1)
    {
        
        // Out of image boundary
        return false;
    
    }

    // Check the boundary in target
    float xTarget = (float)u2;
    float yTarget = (float)v2;
    if (xTarget - 2*nepx*rescaleFactor < 0 ||
        xTarget - 2*nepx*rescaleFactor > width - 1 ||
        xTarget + 2*nepx*rescaleFactor < 0 ||
        xTarget + 2*nepx*rescaleFactor > width - 1 ||
        yTarget - 2*nepy*rescaleFactor < 0 ||
        yTarget - 2*nepy*rescaleFactor > height - 1 ||
        yTarget + 2*nepy*rescaleFactor < 0 ||
        yTarget + 2*nepy*rescaleFactor > height - 1 )
    {
        
        // Out of image boundary
        return false;
        
    }
    
    return true;
    
}

//------------------------------------------------------------------------------
struct Runnested_Y
{
private:
    double w_;
    double yp_;
    double a_;
    double b_;
    double c_;
    double ofx_;
    double ofy_;
    double ocx_;
    double ocy_;
    
public:
    Runnested_Y(double w, double yp, double a, double b, double c,
                double ofx, double ofy, double ocx, double ocy) :
    w_(w), yp_(yp), a_(a), b_(b), c_(c),
    ofx_(ofx), ofy_(ofy),
    ocx_(ocx), ocy_(ocy)
    {
        
    }
    
    template <typename T>
    bool operator() (const T* const param, T* residual) const
    {
        T y = param[0];
        
        // Undistorted x in the image from the epipolar line
        // equation, ax + by + c = 0
        T x = T(-b_/a_)*y - T(c_/a_);
        
        // Normalized undistored x and y
        T x_n = (x - T(ocx_))/T(ofx_);
        T y_n = (y - T(ocy_))/T(ofy_);
        
        // Radius of normalized undistorted
        T r = sqrt(x_n*x_n + y_n*y_n);
        
        // Factor rp/r
        T fac;
        if (fabs(ceres::JetOps<T>::GetScalar(r)) < 1e-8)
        {
            
            fac = T(1.0);
            
        }
        else
        {
            
            fac = (T(w_)*r)/atan(T(2.0)*r*tan(T(w_)/T(2.0)));
            
        }
        
        // Error (to be squared later)
        residual[0] = (fac*T(yp_) - y_n);
                
        return true;
        
    }
    
};

double
runnested_y(double w, double yp, double a, double b, double c, double y0,
            double ofx, double ofy, double ocx, double ocy,
            double *final_cost,
            ceres::TerminationType *termType)
{
    

    
    ceres::Problem problem;

    double y = y0;
    problem.AddResidualBlock(new ceres::AutoDiffCostFunction<Runnested_Y, 1, 1>
                             (new Runnested_Y(w, yp, a, b, c,
                                              ofx, ofy, ocx, ocy)),
                             NULL /* Squared */,
                             &y);
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = false;
    ceres::Solver::Summary summary;
    Solve(options, &problem, &summary);

    *final_cost = summary.final_cost;
    *termType = summary.termination_type;

    return y;
    
}

//------------------------------------------------------------------------------
struct Runnested_X
{
private:
    double w_;
    double yp_;
    double a_;
    double b_;
    double c_;
    double ofx_;
    double ofy_;
    double ocx_;
    double ocy_;
    
public:
    Runnested_X(double w, double yp, double a, double b, double c,
                double ofx, double ofy, double ocx, double ocy) :
    w_(w), yp_(yp), a_(a), b_(b), c_(c),
    ofx_(ofx), ofy_(ofy),
    ocx_(ocx), ocy_(ocy)
    {
        
    }
    
    template <typename T>
    bool operator() (const T* const param, T* residual) const
    {
        T x = param[0];
        
        // Undistorted y in the image from the epipolar line
        // equation, ax + by + c = 0
        T y = T(-a_/b_)*x - T(c_/b_);

        // Normalized undistored x and y
        T x_n = (x - T(ocx_))/T(ofx_);
        T y_n = (y - T(ocy_))/T(ofy_);
        
        // Radius of normalized undistorted
        T r = sqrt(x_n*x_n + y_n*y_n);
        
        // Factor rp/r
        T fac;
        if (fabs(ceres::JetOps<T>::GetScalar(r)) < 1e-8)
        {
            
            fac = T(1.0);
            
        }
        else
        {
            
            fac = (T(w_)*r)/atan(T(2.0)*r*tan(T(w_)/T(2.0)));
                                 
        }
        
        // Error (to be squared later)
        residual[0] = (fac*T(yp_) - y_n);
        
        return true;
        
    }
    
};

double
runnested_x(double w, double yp, double a, double b, double c, double x0,
            double ofx, double ofy, double ocx, double ocy,
            double *final_cost,
            ceres::TerminationType *termType)
{

    ceres::Problem problem;

    double x = x0;
    problem.AddResidualBlock(new ceres::AutoDiffCostFunction<Runnested_X, 1, 1>
                             (new Runnested_X(w, yp, a, b, c,
                                              ofx, ofy, ocx, ocy)),
                             NULL /* Squared */,
                             &x);
    ceres::Solver::Options options;
    options.use_explicit_schur_complement = true;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = false;
    ceres::Solver::Summary summary;
    Solve(options, &problem, &summary);

    *final_cost = summary.final_cost;
    *termType = summary.termination_type;
    
    return x;
    
}

//------------------------------------------------------------------------------
// Find intersections with epipolar line for radial-rolling image
// in the row (y_dist) index in the distorted image
//
//   epline : [double 3x1] Epipolar line induced by a pose at distorted y
//   y_dist : [double] Distorted y coordinate in image
//   x_init : [double] Initial x for optimization
//   y_init : [double] Initial y for optimization
//
// Output
//
//   x_undist : [double] Output of undistorted x coordinate
//   y_undist : [double] Output of undistorted y coordinate
//
bool
DepthMapForRadialRolling::
findInterEpipolar_Optimization(Eigen::Matrix<double,3,1> epline,
                               double y_dist,
                               double x_init,
                               double y_init,
                               double *x_undist,
                               double *y_undist)
{
    
    bool isFound = false;
    Frame *sourceImage = activeKeyFrame;

    cv::Mat K_in = sourceImage->undistorter_->getOriginalK();
    cv::Mat K_out = sourceImage->undistorter_->getK();
    
    // Camera parameters
    // Distorted
    double fy_dist = K_in.at<double>(1, 1);
    double cy_dist = K_in.at<double>(1, 2);
    double fx_undist = K_out.at<double>(0, 0); // Undistorted
    double fy_undist = K_out.at<double>(1, 1);
    double cx_undist = K_out.at<double>(0, 2);
    double cy_undist = K_out.at<double>(1, 2);
    double w = (double)sourceImage->undistorter_->getInputDistort();
    
    // Epipolar line
    double a = epline(0);
    double b = epline(1);
    double c = epline(2);
    
    // Normalized distorted y
    double yp = (y_dist - cy_dist)/fy_dist;
    
    double fval, x_sol, y_sol;
    ceres::TerminationType termType;
    if (fabs(b) < 1e-08 || fabs(a) > fabs(b))
    {
        
        double y0 = y_init; // Initial undistorted y in the image
        y_sol = runnested_y(w, yp, a, b, c, y0,
                                   fx_undist, fy_undist,
                                   cx_undist, cy_undist,
                                   &fval, &termType);
        x_sol = (-b/a)*y_sol - (c/a);
    
    }
    else
    {
        
        double x0 = x_init; // Initial undistorted x in the image
        x_sol = runnested_x(w, yp, a, b, c, x0,
                                   fx_undist, fy_undist,
                                   cx_undist, cy_undist,
                                   &fval, &termType);
        y_sol = (-a/b)*x_sol - (c/b);
    
    }
    
    // Store only good minimum
    double out_height = sourceImage->undistorter_->getOutputHeight();
    double out_width = sourceImage->undistorter_->getOutputWidth();
    if (termType == ceres::CONVERGENCE &&
        0 <= y_sol && y_sol <= out_height - 1 &&
        0 <= x_sol && x_sol <= out_width - 1)
    {

#if 0 // DEBUG
        if (fabs(b) < 2.2e-16 || fabs(a) > fabs(b))
        {
            
        }
        else
        {
            
        
            this->debug_x_sol = x_sol;
            this->debug_y_sol = y_sol;
            this->debug_a = a;
            this->debug_b = b;
            this->debug_c = c;
            this->debug_y_dist = y_dist;
            this->debug_x_init = x_init;
            this->debug_y_init = y_init;
            
            fprintf(stderr, "e = [%f %f %f]'\n",
                    this->debug_a,
                    this->debug_b,
                    this->debug_c);
            fprintf(stderr, "y_dist = %f\n", y_dist);
            fprintf(stderr, "x_init = %f\n", x_init);
            fprintf(stderr, "y_init = %f\n", y_init);
            fprintf(stderr, "x_sol = %f\n", x_sol);
            fprintf(stderr, "y_sol = %f\n", y_sol);
            fprintf(stderr, "----\n");
            
            ;
            
        }
#endif
    
        *x_undist = x_sol;
        *y_undist = y_sol;

        isFound = true;
    
    }
    
    return isFound;
    
}

//------------------------------------------------------------------------------
// Process depth propagation for pixel row in distorted image
// The depth propagated should be rectified later
// This distort version is for speed-up depth propagation
void
DepthMapForRadialRolling::
processPropagationForPixelRow_Distort(
                              int x_undist, int y_undist,
                              const Frame *old_keyframe,
                              Frame *new_keyframe,
                              const DepthMapPixelHypothesis *source,
                              int yMin_dist,
                              int yMax_dist,
                              RunningStats* stats)
{
    
    const bool* trackingWasGood =
    new_keyframe->getTrackingParent() == activeKeyFrame ?
    new_keyframe->refPixelWasGoodNoCreate() : 0;
    const float* activeKFImageData = activeKeyFrame->image(0);
    const float* newKFMaxGrad = new_keyframe->maxGradients(0);
    const float* newKFImageData = new_keyframe->image(0);
    
    
    if (enablePrintDebugInfo)
    {
        
        runningStats.num_prop_attempts++;
        
    }
    
    // A 3D point in the old frame back-projected from 2D image point
    Eigen::Vector3d old2DPoint;
    old2DPoint << x_undist, y_undist, 1.0;
    Eigen::Matrix3d Kinv = old_keyframe->KInv().cast<double>();
    Eigen::Vector4d old3DPoint;
    old3DPoint << Kinv*old2DPoint / source->idepth_smoothed, 1.0;
    
    // Find all 2D point projection in the new frame
    // by checking each row in the new frame
    Sophus::SE3d worldToOld;
    if (useGTmotion_ == true)
    {
        // GT motion
        
        // TODO: THIS NEEDS TO BE IMPLEMENTED FOR RADIAL-ROLLING
        //       CURRENTLY THIS WOULD NOT WORK CORRECTLY
        worldToOld = Sophus::SE3d(old_keyframe->getGTPoseAtRow_Distort(x_undist, y_undist));
        
    }
    else
    {
        // Est motion
        // NOTE: we imagine old_keyframe in distorted image and get a pose
        //       at the row y, but this will be rectified later stage
        worldToOld = Sophus::SE3d(old_keyframe->getPoseAtRow_Distort(x_undist, y_undist));
    }
    
    Eigen::Vector4d new3DPoint;
    for (int rowNew_dist = yMin_dist; rowNew_dist <= yMax_dist; rowNew_dist++)
    {
            
        // Motion of the old to new frame
        Sophus::SE3d worldToNew;
        if (useGTmotion_ == true)
        {
            // GT motion
            
            // TODO: THIS NEEDS TO BE IMPLEMENTED FOR RADIAL-ROLLING
            //       CURRENTLY THIS WOULD NOT WORK CORRECTLY
            
            worldToNew = Sophus::SE3d(new_keyframe->getGTPoseAtRow(rowNew_dist));
            
        }
        else
        {
            worldToNew =
            Sophus::SE3d(new_keyframe->getPoseAtRow(rowNew_dist));
        }
        Sophus::SE3d oldToNew = worldToNew * worldToOld.inverse();
        
        // Transform the 3D point from old to new frame
        new3DPoint = oldToNew.matrix()*old3DPoint;
        
        // Project the 3D point in new frame onto image
        Eigen::Matrix3d K = new_keyframe->K().cast<double>();
        Eigen::Matrix<double,3,1> m = K*new3DPoint.block(0,0,3,1);
        Eigen::Vector3d new2DPoint;
        new2DPoint << m(0)/m(2), m(1)/m(2), 1.0;
        
        // Debug
        this->debugInfo_old2DPoint = old2DPoint;
        this->debugInfo_new2DPoint = new2DPoint;
        this->debugInfo_worldToNew = worldToNew;
        this->debugInfo_worldToOld = worldToOld;
        this->debugInfo_oldToNew = oldToNew;
        this->debugInfo_rowNew = rowNew_dist;
        
        // Set new image point coordinates (undistorted)
        float u_new_undist = (float)new2DPoint[0];
        float v_new_undist = (float)new2DPoint[1];
            
        // Now find its distorted pixel
        float u_new_dist, v_new_dist;
        new_keyframe->undistorter_->distortPoint(u_new_undist,
                                                 v_new_undist,
                                                 &u_new_dist,
                                                 &v_new_dist);
        
        // Check if the row of the 2D point is at the new row
        if (fabs(v_new_dist
                 - rowNew_dist) < 0.5) // A half pixel
        {
            
            // Inverse depth in the new frame
            float new_idepth = 1.0f / (float)new3DPoint[2];
            
            
            // check if still within image, if not: DROP.
            if (!(u_new_undist > 2.1f &&
                  v_new_undist > 2.1f &&
                  u_new_undist < width-3.1f &&
                  v_new_undist < height-3.1f))
            {
                
                if (enablePrintDebugInfo)
                {
                    
                    runningStats.num_prop_removed_out_of_bounds++;
                    
                }
                continue;
                
            }
            
            int newIDX = (int)(u_new_undist+0.5f) +
                            ((int)(v_new_undist+0.5f))*width;
            float destAbsGrad = newKFMaxGrad[newIDX];
            
            if (trackingWasGood != 0)
            {
                
                if (!trackingWasGood[(x_undist >> SE3TRACKING_MIN_LEVEL) +
                                     (width >> SE3TRACKING_MIN_LEVEL) *
                                     (y_undist >> SE3TRACKING_MIN_LEVEL)]
                    || destAbsGrad < MIN_ABS_GRAD_DECREASE)
                {
                    
                    if (enablePrintDebugInfo)
                    {
                        
                        runningStats.num_prop_removed_colorDiff++;
                        
                    }
                    continue;
                    
                }
                
            }
            else
            {
                
                float sourceColor = activeKFImageData[x_undist + y_undist*width];
                float destColor =
                getInterpolatedElement(newKFImageData,
                                       u_new_undist,
                                       v_new_undist,
                                       width);
                
                float residual = destColor - sourceColor;
                
                if (residual*residual /
                    (MAX_DIFF_CONSTANT +
                     MAX_DIFF_GRAD_MULT*destAbsGrad*destAbsGrad) > 1.0f ||
                    destAbsGrad < MIN_ABS_GRAD_DECREASE)
                {
                    
                    if (enablePrintDebugInfo)
                    {
                        
                        runningStats.num_prop_removed_colorDiff++;
                        
                    }
                    continue;
                    
                }
                
            }
            
            DepthMapPixelHypothesis* targetBest = otherDepthMap +  newIDX;
            
            // large idepth = point is near = large increase in variance.
            // small idepth = point is far = small increase in variance.
            float idepth_ratio_4 = new_idepth / source->idepth_smoothed;
            idepth_ratio_4 *= idepth_ratio_4;
            idepth_ratio_4 *= idepth_ratio_4;
            
            float new_var = idepth_ratio_4*source->idepth_var;
            
            // check for occlusion
            if (targetBest->isValid)
            {
                // if they occlude one another, one gets removed.
                float diff = targetBest->idepth - new_idepth;
                if (DIFF_FAC_PROP_MERGE*diff*diff >
                    new_var +
                    targetBest->idepth_var)
                {
                    
                    if (new_idepth < targetBest->idepth)
                    {
                        
                        if (enablePrintDebugInfo)
                        {
                            
                            runningStats.num_prop_occluded++;
                            
                        }
                        continue;
                        
                    }
                    else
                    {
                        
                        if (enablePrintDebugInfo)
                        {
                            
                            runningStats.num_prop_occluded++;
                            
                        }
                        targetBest->isValid = false;
                        
                    }
                    
                }
                
            }
            
            
            if (!targetBest->isValid)
            {
                if (enablePrintDebugInfo)
                {
                    
                    runningStats.num_prop_created++;
                    
                }
                
                *targetBest = DepthMapPixelHypothesis(new_idepth,
                                                      new_var,
                                                      source->validity_counter);
                
            }
            else
            {
                
                if (enablePrintDebugInfo)
                {
                    
                    runningStats.num_prop_merged++;
                    
                }
                
                // merge idepth ekf-style
                float w = new_var / (targetBest->idepth_var + new_var);
                float merged_new_idepth =
                w*targetBest->idepth + (1.0f-w)*new_idepth;
                
                // merge validity
                int merged_validity =
                source->validity_counter + targetBest->validity_counter;
                if(merged_validity >
                   VALIDITY_COUNTER_MAX+(VALIDITY_COUNTER_MAX_VARIABLE))
                {
                    
                    merged_validity = (int)round(VALIDITY_COUNTER_MAX +
                               (VALIDITY_COUNTER_MAX_VARIABLE));
                    
                }
                
                *targetBest =
                DepthMapPixelHypothesis(merged_new_idepth,
                                        1.0f/(1.0f/targetBest->idepth_var
                                              + 1.0f/new_var),
                                        merged_validity);
                
            }
            
        } // end if fabs
        
    } // end for rowNew
    
}

void
DepthMapForRadialRolling::
processPropagationForPixelRow(
                              int x, int y,
                              const Frame *old_keyframe,
                              Frame *new_keyframe,
                              const DepthMapPixelHypothesis *source,
                              int yMin,
                              int yMax,
                              RunningStats* stats)
{
 
    const bool* trackingWasGood =
    new_keyframe->getTrackingParent() == activeKeyFrame ?
    new_keyframe->refPixelWasGoodNoCreate() : 0;
    const float* activeKFImageData = activeKeyFrame->image(0);
    const float* newKFMaxGrad = new_keyframe->maxGradients(0);
    const float* newKFImageData = new_keyframe->image(0);

    
    if (enablePrintDebugInfo)
    {
        
        runningStats.num_prop_attempts++;
        
    }
    
    // A 3D point in the old frame back-projected from 2D image point
    Eigen::Vector3d old2DPoint;
    old2DPoint << x, y, 1.0;
    Eigen::Matrix3d Kinv = old_keyframe->KInv().cast<double>();
    Eigen::Vector4d old3DPoint;
    old3DPoint << Kinv*old2DPoint / source->idepth_smoothed, 1.0;
    
    // Find all 2D point projection in the new frame
    // by checking each row in the new frame
    Sophus::SE3d worldToOld;
    if (useGTmotion_ == true)
    {
        // GT motion
        
        // TODO: THIS NEEDS TO BE IMPLEMENTED FOR RADIAL-ROLLING
        //       CURRENTLY THIS WOULD NOT WORK CORRECTLY
        worldToOld = Sophus::SE3d(old_keyframe->getGTPoseAtRow(y));
        
    }
    else
    {
        // Est motion
        worldToOld = Sophus::SE3d(old_keyframe->getPoseAtRow_Distort(x, y));
    }
    
    bool isFound = false;
    Eigen::Vector4d new3DPoint;
    for (int rowNew = yMin; rowNew <= yMax && isFound == false; rowNew++)
    {
        
        for (int colNew = 0; colNew < width && isFound == false; colNew++)
        {
            
            // Motion of the old to new frame
            Sophus::SE3d worldToNew;
            if (useGTmotion_ == true)
            {
                // GT motion
                
                // TODO: THIS NEEDS TO BE IMPLEMENTED FOR RADIAL-ROLLING
                //       CURRENTLY THIS WOULD NOT WORK CORRECTLY
                
                worldToNew = Sophus::SE3d(new_keyframe->getGTPoseAtRow(rowNew));
                
            }
            else
            {
                worldToNew =
                Sophus::SE3d(new_keyframe->getPoseAtRow_Distort(colNew, rowNew));
            }
            Sophus::SE3d oldToNew = worldToNew * worldToOld.inverse();
            
            // Transform the 3D point from old to new frame
            new3DPoint = oldToNew.matrix()*old3DPoint;
            
            // Project the 3D point in new frame onto image
            Eigen::Matrix3d K = new_keyframe->K().cast<double>();
            Eigen::Matrix<double,3,1> m = K*new3DPoint.block(0,0,3,1);
            Eigen::Vector3d new2DPoint;
            new2DPoint << m(0)/m(2), m(1)/m(2), 1.0;
            
            // Debug
            this->debugInfo_old2DPoint = old2DPoint;
            this->debugInfo_new2DPoint = new2DPoint;
            this->debugInfo_worldToNew = worldToNew;
            this->debugInfo_worldToOld = worldToOld;
            this->debugInfo_oldToNew = oldToNew;
            this->debugInfo_rowNew = rowNew;
            
            // Check if the row of the 2D point is at the row i and
            //       if the column of the 2D point is at the column i
            if (fabs(new2DPoint[1] - rowNew) < 0.5 && // A half pixel
                fabs(new2DPoint[0] - colNew) < 0.5)
            {
                
                // Inverse depth in the new frame
                float new_idepth = 1.0f / (float)new3DPoint[2];
                
                // Set new image point coordinates
                float u_new = (float)new2DPoint[0];
                float v_new = (float)new2DPoint[1];
                
                // check if still within image, if not: DROP.
                if (!(u_new > 2.1f &&
                      v_new > 2.1f &&
                      u_new < width-3.1f &&
                      v_new < height-3.1f))
                {
                    
                    if (enablePrintDebugInfo)
                    {
                        
                        runningStats.num_prop_removed_out_of_bounds++;
                        
                    }
                    continue;
                    
                }
                
                int newIDX = (int)(u_new+0.5f) + ((int)(v_new+0.5f))*width;
                float destAbsGrad = newKFMaxGrad[newIDX];
                
                if (trackingWasGood != 0)
                {
                    
                    if (!trackingWasGood[(x >> SE3TRACKING_MIN_LEVEL) +
                                         (width >> SE3TRACKING_MIN_LEVEL) *
                                         (y >> SE3TRACKING_MIN_LEVEL)]
                        || destAbsGrad < MIN_ABS_GRAD_DECREASE)
                    {
                        
                        if (enablePrintDebugInfo)
                        {
                            
                            runningStats.num_prop_removed_colorDiff++;
                            
                        }
                        continue;
                        
                    }
                    
                }
                else
                {
                    
                    float sourceColor = activeKFImageData[x + y*width];
                    float destColor =
                    getInterpolatedElement(newKFImageData, u_new, v_new, width);
                    
                    float residual = destColor - sourceColor;
                    
                    if (residual*residual /
                        (MAX_DIFF_CONSTANT +
                         MAX_DIFF_GRAD_MULT*destAbsGrad*destAbsGrad) > 1.0f ||
                        destAbsGrad < MIN_ABS_GRAD_DECREASE)
                    {
                        
                        if (enablePrintDebugInfo)
                        {
                            
                            runningStats.num_prop_removed_colorDiff++;
                            
                        }
                        continue;
                        
                    }
                    
                }
                
                DepthMapPixelHypothesis* targetBest = otherDepthMap +  newIDX;
                
                // large idepth = point is near = large increase in variance.
                // small idepth = point is far = small increase in variance.
                float idepth_ratio_4 = new_idepth / source->idepth_smoothed;
                idepth_ratio_4 *= idepth_ratio_4;
                idepth_ratio_4 *= idepth_ratio_4;
                
                float new_var = idepth_ratio_4*source->idepth_var;
                
                // check for occlusion
                if (targetBest->isValid)
                {
                    // if they occlude one another, one gets removed.
                    float diff = targetBest->idepth - new_idepth;
                    if (DIFF_FAC_PROP_MERGE*diff*diff >
                        new_var +
                        targetBest->idepth_var)
                    {
                        
                        if (new_idepth < targetBest->idepth)
                        {
                            
                            if (enablePrintDebugInfo)
                            {
                                
                                runningStats.num_prop_occluded++;
                                
                            }
                            continue;
                            
                        }
                        else
                        {
                            
                            if (enablePrintDebugInfo)
                            {
                                
                                runningStats.num_prop_occluded++;
                                
                            }
                            targetBest->isValid = false;
                            
                        }
                        
                    }
                    
                }
                
                
                if (!targetBest->isValid)
                {
                    if (enablePrintDebugInfo)
                    {
                        
                        runningStats.num_prop_created++;
                        
                    }
                    
                    *targetBest = DepthMapPixelHypothesis(new_idepth,
                                                          new_var,
                                                          source->validity_counter);
                    
                }
                else
                {
                    
                    if (enablePrintDebugInfo)
                    {
                        
                        runningStats.num_prop_merged++;
                        
                    }
                    
                    // merge idepth ekf-style
                    float w = new_var / (targetBest->idepth_var + new_var);
                    float merged_new_idepth =
                    w*targetBest->idepth + (1.0f-w)*new_idepth;
                    
                    // merge validity
                    int merged_validity =
                    source->validity_counter + targetBest->validity_counter;
                    if(merged_validity >
                       VALIDITY_COUNTER_MAX+(VALIDITY_COUNTER_MAX_VARIABLE))
                    {
                        
                        merged_validity =
                        (int)round(VALIDITY_COUNTER_MAX +
                                   (VALIDITY_COUNTER_MAX_VARIABLE));
                        
                    }
                    
                    *targetBest =
                    DepthMapPixelHypothesis(merged_new_idepth,
                                            1.0f/(1.0f/targetBest->idepth_var
                                                  + 1.0f/new_var),
                                            merged_validity);
                    
                    isFound = true;
                    
                }
                
            } // end if fabs
            
        } // end for colNew
        
    } // end for rowNew
    
}

//void
//DepthMapForRadialRolling::
//propagateDepthRow(int yMin, int yMax, RunningStats* stats,
//                    int x, int y,
//                    const Frame *old_keyframe,
//                    const Frame *new_keyframe,
//                    const DepthMapPixelHypothesis *source,
//                    const float* newKFMaxGrad,
//                    const bool *trackingWasGood,
//                    const float* activeKFImageData,
//                    const float* newKFImageData)
//{
//    
//    // Go through all pixels (x, y) of OLD image, propagating forwards.
//    for (int y=yMin;y<yMax;y++)
//    {
//        
//        for (int x=0;x<width;x++)
//        {
//            
//            DepthMapPixelHypothesis* source = currentDepthMap + x + y*width;
//            
//            if (source->isValid == true)
//            {
//                
//                // Process progation for a pixel
//                processPropagationForPixel(x, y,
//                                           old_keyframe,
//                                           new_keyframe,
//                                           source,
//                                           newKFMaxGrad,
//                                           trackingWasGood,
//                                           activeKFImageData,
//                                           newKFImageData);
//                
//            }
//            
//        } // end for x
//        
//    } // end for y
//    
//}

void DepthMapForRadialRolling::propagateDepth(Frame* new_keyframe)
{
    
    runningStats.num_prop_removed_out_of_bounds = 0;
    runningStats.num_prop_removed_colorDiff = 0;
    runningStats.num_prop_removed_validity = 0;
    runningStats.num_prop_grad_decreased = 0;
    runningStats.num_prop_color_decreased = 0;
    runningStats.num_prop_attempts = 0;
    runningStats.num_prop_occluded = 0;
    runningStats.num_prop_created = 0;
    runningStats.num_prop_merged = 0;
    
    // Old keyframe
    Frame *old_keyframe = new_keyframe->getTrackingParent();
    
    // Check if the old keyframe is the same as active keyframe
    if(new_keyframe->getTrackingParent() != activeKeyFrame)
    {
        
        printf("WARNING: propagating depth from frame %d to %d,"
               "which was tracked on a different frame (%d).\n"
               "While this should work, it is not recommended.",
               activeKeyFrame->id(), new_keyframe->id(),
               new_keyframe->getTrackingParent()->id());
        fprintf(stderr, "WARNING: propagating depth from frame %d to %d,"
               "which was tracked on a different frame (%d).\n"
               "While this should work, it is not recommended.",
               activeKeyFrame->id(), new_keyframe->id(),
               new_keyframe->getTrackingParent()->id());
    }
    
    // wipe depthmap
    for(DepthMapPixelHypothesis* pt = otherDepthMap+width*height-1;
        pt >= otherDepthMap; pt--)
    {
        
        pt->isValid = false;
        pt->blacklisted = 0;
        
    }
    
    // re-usable values.
    //    SE3 oldToNew_SE3 =
    //        se3FromSim3(new_keyframe->pose->thisToParent_raw).inverse();
    //    Eigen::Vector3f trafoInv_t = oldToNew_SE3.translation().cast<float>();
    //    Eigen::Matrix3f trafoInv_R =
    //        oldToNew_SE3.rotationMatrix().matrix().cast<float>();
    //  
    
    // Go through all pixels (x, y) of OLD image, propagating forwards.
    for (int y_undist=0;y_undist<height;y_undist++)
    {
        
        for (int x_undist=0;x_undist<width;x_undist++)
        {
            
            int idx = x_undist + y_undist*width;
            DepthMapPixelHypothesis* source = currentDepthMap + idx;
            
            if (source->isValid == true)
            {
                
#if 0 // NO-THREAD

                // Process progation for a pixel
    
#if 0 // FAST-VERSION
                
                RunningStats stats;
                processPropagationForPixelRow_Distort(x, y,
                                                      old_keyframe,
                                                      new_keyframe,
                                                      source,
                                                      3,
                                                      height - 3,
                                                      &stats);
                
#else // SLOW-VERSION
                
                processPropagationForPixel(x, y,
                                          old_keyframe,
                                          new_keyframe,
                                          source,
                                          newKFMaxGrad,
                                          trackingWasGood,
                                          activeKFImageData,
                                          newKFImageData);
#endif
                
#else // THREAD
                
                // Process progation for a pixel
                threadReducer.reduce(boost::bind(
                    &DepthMapForRadialRolling::
                    processPropagationForPixelRow_Distort,
                    this,
                    x_undist,
                    y_undist,
                    old_keyframe,
                    new_keyframe,
                    source,
                    _1, _2, _3),
                3, // start row
                new_keyframe->height_rollingShutter_inputImage - 3, // end row
                10); // step

#endif
            }
            
        } // end for x
        
    } // end for y
    
    

    
    // swap!
    std::swap(currentDepthMap, otherDepthMap);
    
    if (enablePrintDebugInfo && printPropagationStatistics)
    {
        
        printf("PROPAGATE: %d: %d drop (%d oob, %d color); %d created;"
               " %d merged; %d occluded. %d col-dec, %d grad-dec.\n",
               runningStats.num_prop_attempts,
               runningStats.num_prop_removed_validity +
               runningStats.num_prop_removed_out_of_bounds +
               runningStats.num_prop_removed_colorDiff,
               runningStats.num_prop_removed_out_of_bounds,
               runningStats.num_prop_removed_colorDiff,
               runningStats.num_prop_created,
               runningStats.num_prop_merged,
               runningStats.num_prop_occluded,
               runningStats.num_prop_color_decreased,
               runningStats.num_prop_grad_decreased);
        
    }
    
}


    
//------------------------------------------------------------------------------
//      x : [int] Pixel x-coordinate of an image point in the old frame
//      y : [int] Pixel y-coordinate of an image point in the old frame
// source : [DepthMapPixelHypothesis] Depth map of the old frame
void DepthMapForRadialRolling::
processPropagationForPixel(int x, int y,
                           const Frame *old_keyframe,
                           const Frame *new_keyframe,
                           const DepthMapPixelHypothesis *source,
                           const float* newKFMaxGrad,
                           const bool *trackingWasGood,
                           const float* activeKFImageData,
                           const float* newKFImageData)
{
    
    if (enablePrintDebugInfo)
    {
        
        runningStats.num_prop_attempts++;
        
    }
    
    // A 3D point in the old frame back-projected from 2D image point
    Eigen::Vector3d old2DPoint;
    old2DPoint << x, y, 1.0;
    Eigen::Matrix3d Kinv = old_keyframe->KInv().cast<double>();
    Eigen::Vector4d old3DPoint;
    old3DPoint << Kinv*old2DPoint / source->idepth_smoothed, 1.0;
    
    // Find all 2D point projection in the new frame
    // by checking each row in the new frame
    Sophus::SE3d worldToOld;
    if (useGTmotion_ == true)
    {
        // GT motion
        
        // TODO: THIS NEEDS TO BE IMPLEMENTED FOR RADIAL-ROLLING
        //       CURRENTLY THIS WOULD NOT WORK CORRECTLY
        worldToOld = Sophus::SE3d(old_keyframe->getGTPoseAtRow(y));
        
    }
    else
    {
        // Est motion
        worldToOld = Sophus::SE3d(old_keyframe->getPoseAtRow_Distort(x, y));
    }
    
    bool isFound = false;
    Eigen::Vector4d new3DPoint;
    for (int rowNew = 0; rowNew < height && isFound == false; rowNew++)
    {
        
        for (int colNew = 0; colNew < width && isFound == false; colNew++)
        {

            // Motion of the old to new frame
            Sophus::SE3d worldToNew;
            if (useGTmotion_ == true)
            {
                // GT motion

                // TODO: THIS NEEDS TO BE IMPLEMENTED FOR RADIAL-ROLLING
                //       CURRENTLY THIS WOULD NOT WORK CORRECTLY

                worldToNew = Sophus::SE3d(new_keyframe->getGTPoseAtRow(rowNew));
                
            }
            else
            {
                worldToNew =
                Sophus::SE3d(new_keyframe->getPoseAtRow_Distort(colNew, rowNew));
            }
            Sophus::SE3d oldToNew = worldToNew * worldToOld.inverse();
            
            // Transform the 3D point from old to new frame
            new3DPoint = oldToNew.matrix()*old3DPoint;
            
            // Project the 3D point in new frame onto image
            Eigen::Matrix3d K = new_keyframe->K().cast<double>();
            Eigen::Matrix<double,3,1> m = K*new3DPoint.block(0,0,3,1);
            Eigen::Vector3d new2DPoint;
            new2DPoint << m(0)/m(2), m(1)/m(2), 1.0;
            
            // Debug
            this->debugInfo_old2DPoint = old2DPoint;
            this->debugInfo_new2DPoint = new2DPoint;
            this->debugInfo_worldToNew = worldToNew;
            this->debugInfo_worldToOld = worldToOld;
            this->debugInfo_oldToNew = oldToNew;
            this->debugInfo_rowNew = rowNew;
            
            // Check if the row of the 2D point is at the row i and
            //       if the column of the 2D point is at the column i
            if (fabs(new2DPoint[1] - rowNew) < 0.5 && // A half pixel
                fabs(new2DPoint[0] - colNew) < 0.5)
            {
                
                // Inverse depth in the new frame
                float new_idepth = 1.0f / (float)new3DPoint[2];
                
                // Set new image point coordinates
                float u_new = (float)new2DPoint[0];
                float v_new = (float)new2DPoint[1];
                
                // check if still within image, if not: DROP.
                if (!(u_new > 2.1f &&
                      v_new > 2.1f &&
                      u_new < width-3.1f &&
                      v_new < height-3.1f))
                {
                    
                    if (enablePrintDebugInfo)
                    {
                        
                        runningStats.num_prop_removed_out_of_bounds++;
                        
                    }
                    continue;
                    
                }
                
                int newIDX = (int)(u_new+0.5f) + ((int)(v_new+0.5f))*width;
                float destAbsGrad = newKFMaxGrad[newIDX];
                
                if (trackingWasGood != 0)
                {
                    
                    if (!trackingWasGood[(x >> SE3TRACKING_MIN_LEVEL) +
                                         (width >> SE3TRACKING_MIN_LEVEL) *
                                         (y >> SE3TRACKING_MIN_LEVEL)]
                        || destAbsGrad < MIN_ABS_GRAD_DECREASE)
                    {
                        
                        if (enablePrintDebugInfo)
                        {
                            
                            runningStats.num_prop_removed_colorDiff++;
                            
                        }
                        continue;
                        
                    }
                    
                }
                else
                {
                    
                    float sourceColor = activeKFImageData[x + y*width];
                    float destColor =
                    getInterpolatedElement(newKFImageData, u_new, v_new, width);
                    
                    float residual = destColor - sourceColor;
                    
                    if (residual*residual /
                        (MAX_DIFF_CONSTANT +
                         MAX_DIFF_GRAD_MULT*destAbsGrad*destAbsGrad) > 1.0f ||
                        destAbsGrad < MIN_ABS_GRAD_DECREASE)
                    {
                        
                        if (enablePrintDebugInfo)
                        {
                            
                            runningStats.num_prop_removed_colorDiff++;
                            
                        }
                        continue;
                        
                    }
                    
                }
                
                DepthMapPixelHypothesis* targetBest = otherDepthMap +  newIDX;
                
                // large idepth = point is near = large increase in variance.
                // small idepth = point is far = small increase in variance.
                float idepth_ratio_4 = new_idepth / source->idepth_smoothed;
                idepth_ratio_4 *= idepth_ratio_4;
                idepth_ratio_4 *= idepth_ratio_4;
                
                float new_var = idepth_ratio_4*source->idepth_var;
                
                // check for occlusion
                if (targetBest->isValid)
                {
                    // if they occlude one another, one gets removed.
                    float diff = targetBest->idepth - new_idepth;
                    if (DIFF_FAC_PROP_MERGE*diff*diff >
                        new_var +
                        targetBest->idepth_var)
                    {
                        
                        if (new_idepth < targetBest->idepth)
                        {
                            
                            if (enablePrintDebugInfo)
                            {
                                
                                runningStats.num_prop_occluded++;
                                
                            }
                            continue;
                            
                        }
                        else
                        {
                            
                            if (enablePrintDebugInfo)
                            {
                                
                                runningStats.num_prop_occluded++;
                                
                            }
                            targetBest->isValid = false;
                            
                        }
                        
                    }
                    
                }
                
                
                if (!targetBest->isValid)
                {
                    if (enablePrintDebugInfo)
                    {
                        
                        runningStats.num_prop_created++;
                        
                    }
                    
                    *targetBest = DepthMapPixelHypothesis(new_idepth,
                                                          new_var,
                                                          source->validity_counter);
                    
                }
                else
                {
                    
                    if (enablePrintDebugInfo)
                    {
                        
                        runningStats.num_prop_merged++;
                        
                    }
                    
                    // merge idepth ekf-style
                    float w = new_var / (targetBest->idepth_var + new_var);
                    float merged_new_idepth =
                    w*targetBest->idepth + (1.0f-w)*new_idepth;
                    
                    // merge validity
                    int merged_validity =
                    source->validity_counter + targetBest->validity_counter;
                    if(merged_validity >
                       VALIDITY_COUNTER_MAX+(VALIDITY_COUNTER_MAX_VARIABLE))
                    {
                        
                        merged_validity =
                        VALIDITY_COUNTER_MAX+(VALIDITY_COUNTER_MAX_VARIABLE);
                        
                    }
                    
                    *targetBest =
                    DepthMapPixelHypothesis(merged_new_idepth,
                                            1.0f/(1.0f/targetBest->idepth_var
                                                  + 1.0f/new_var),
                                            merged_validity);
                    
                    isFound = true;
                    
                }
                
            } // end if fabs
            
        } // end for colNew
        
    } // end for rowNew
    
}
    
void DepthMapForRadialRolling::initializeFromGTDepth(Frame* new_frame)
{
    assert(new_frame->hasIDepthBeenSet());
    
    activeKeyFramelock = new_frame->getActiveLock();
    activeKeyFrame = new_frame;
    activeKeyFrameImageData = activeKeyFrame->image(0);
    activeKeyFrameIsReactivated = false;
    
    const float* idepth = new_frame->idepth();
    
    
    float averageGTIDepthSum = 0;
    int averageGTIDepthNum = 0;
    for(int y=0;y<height;y++)
    {
        for(int x=0;x<width;x++)
        {
            float idepthValue = idepth[x+y*width];
            if(!isnanf(idepthValue) && idepthValue > 0)
            {
                averageGTIDepthSum += idepthValue;
                averageGTIDepthNum ++;
            }
        }
    }
    
    
    for(int y=0;y<height;y++)
    {
        for(int x=0;x<width;x++)
        {
            float idepthValue = idepth[x+y*width];
            
            if(!isnan(idepthValue) && idepthValue > 0)
            {
                currentDepthMap[x+y*width] =
                    DepthMapPixelHypothesis(
                    idepthValue,
                    idepthValue,
                    VAR_GT_INIT_INITIAL,
                    VAR_GT_INIT_INITIAL,
                    VALIDITY_COUNTER_INITIAL_OBSERVE);
            }
            else
            {
                currentDepthMap[x+y*width].isValid = false;
                currentDepthMap[x+y*width].blacklisted = 0;
            }
        }
    }
    
    
    activeKeyFrame->setDepth(currentDepthMap);
}

} // end of namespace lsd_slam
