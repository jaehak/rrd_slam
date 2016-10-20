//
//  DepthMapForRollingShutter.cpp
//
//  Created by Jae-Hak Kim on 14/05/2015.
//
//
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

#include "DepthMapForRollingShutter.h"

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

#ifdef __APPLE__
#define isnanf isnan
#endif

namespace lsd_slam
{
 
cv::Mat
DepthMapForRollingShutter::
plotGeneralStereoImages(Frame* const sourceFrame,
                        Frame* const targetFrame,
                        const std::list<PointEpiCurve> *ptsEpiCurve,
                        PointEpiCurve *best_ptEC)
{

    if (activeKeyFrame == 0)
    {
        
        return cv::Mat(height, width*2, CV_8UC3);
        
    }
    
    //--------------------------------------------------------------------------
    // Prepare for image structures
    cv::Mat srcImgContainer(sourceFrame->height(),
                   sourceFrame->width(),
                   CV_32F,
                   const_cast<float*>(sourceFrame->image(0)));
    cv::Mat srcImg;
    srcImgContainer.convertTo(srcImg, CV_8UC1);
    cv::cvtColor(srcImg, srcImg, CV_GRAY2RGB);
    
    cv::Mat trgImgContainer(targetFrame->height(),
                   targetFrame->width(),
                   CV_32F,
                   const_cast<float*>(targetFrame->image(0)));
    cv::Mat trgImg;
    trgImgContainer.convertTo(trgImg, CV_8UC1);
    cv::cvtColor(trgImg, trgImg, CV_GRAY2RGB);
    
    //--------------------------------------------------------------------------
    // Draw epipoar lines and normals in the target
    // for (int i=0; i<(int)ptsEpiCurve->size() - 1; i++)
    std::list<PointEpiCurve>::const_iterator ci = ptsEpiCurve->begin();
    std::list<PointEpiCurve>::const_iterator ci2 = ci;
    std::advance(ci2, 1);
    for (;
			ci2 != ptsEpiCurve->end();
			++ci, ++ci2)
    {
        
        float x = (*ci).xTarget;
        float y = (*ci).yTarget;
        trgImg.at<cv::Vec3b>(y, x) = cv::Vec3b(0,0,255);
        
        // Draw epipolar points in the target image
        cv::circle(trgImg, cv::Point(x, y), 1, cv::Scalar(0,255,0));
        float x2 = (*ci2).xTarget;
        float y2 = (*ci2).yTarget;
        cv::line(trgImg, cv::Point(x, y),
                 cv::Point(x2, y2),
                 cv::Scalar(0,255,0), 1); // Green for epipolar curve
        
    }
   
#if 0
    // Draw normals
     for (int i=0; i<(int)ptsEpiCurve->size() - 1 - 10; i = i + 10)
    {
        
        float x = (*ci).xTarget;
        float y = (*ci).yTarget;
        trgImg.at<cv::Vec3b>(y, x) = cv::Vec3b(0,0,255);
        float nepx = (*ci).nepx;
        float nepy = (*ci).nepy;

        // Draw normal on epipolar curve
        cv::line(trgImg, cv::Point(x, y),
                 cv::Point(x + 5*nepx, y + 5*nepy),
                 cv::Scalar(255,255,255), 1); // White for normal
        
    }
#endif
    
    // Draw the best point in the target
    if ((ptsEpiCurve->size() > 0) && (best_ptEC != NULL))
    {
        
        float x = (*best_ptEC).xTarget;
        float y = (*best_ptEC).yTarget;
        float nepx = (*best_ptEC).nepx;
        float nepy = (*best_ptEC).nepy;
        float rescaleFactor = (*best_ptEC).rescaleFactor;
        
        cv::circle(trgImg, cv::Point(x, y), 7, cv::Scalar(0,50,255), 3);
        cv::line(trgImg, cv::Point(x-5,y-5), cv::Point(x+5,y+5),
                 cv::Scalar(0,50,255), 1);
        cv::line(trgImg, cv::Point(x+5,y-5), cv::Point(x-5,y+5),
                 cv::Scalar(0,50,255), 1);
        
#if 0 // DRAW 3x3 window

        cv::line(trgImg, cv::Point(x, y),
                 cv::Point(x - 0.0*len_incxincy,
                           y - 3.0*len_incxincy),
                 cv::Scalar(0,80,127), 2);
        cv::line(trgImg, cv::Point(x, y),
                 cv::Point(x - 3.0*len_incxincy,
                           y - 0.0*len_incxincy),
                 cv::Scalar(0,80,127), 2);
        cv::line(trgImg, cv::Point(x, y),
                 cv::Point(x + 3.0*len_incxincy,
                           y + 0.0*len_incxincy),
                 cv::Scalar(0,80,127), 2);
        cv::line(trgImg, cv::Point(x, y),
                 cv::Point(x + 0.0*len_incxincy,
                           y + 3.0*len_incxincy),
                 cv::Scalar(0,80,127), 2);
        
#else // DRAW 5x1 window

        
        cv::line(trgImg, cv::Point(x, y),
                 cv::Point(x - 4.0*nepx*rescaleFactor,
                           y - 4.0*nepy*rescaleFactor),
                 cv::Scalar(0,80,127), 2);
        cv::line(trgImg, cv::Point(x, y),
                 cv::Point(x + 4.0*nepx*rescaleFactor,
                           y + 4.0*nepy*rescaleFactor),
                 cv::Scalar(0,80,127), 2);
        
#endif
    }
    
    // Draw a point in the source
    if (ptsEpiCurve->size() > 0)
    {
       
		  ci = ptsEpiCurve->begin(); 
        float x = (*ci).xSource;
        float y = (*ci).ySource;
        cv::circle(srcImg, cv::Point(x, y), 7, cv::Scalar(0,255,0), 2);
        cv::line(srcImg, cv::Point(x-5,y-5), cv::Point(x+5,y+5),
                 cv::Scalar(0,255,0), 1);
        cv::line(srcImg, cv::Point(x+5,y-5), cv::Point(x-5,y+5),
                 cv::Scalar(0,255,0), 1);
        
    }
    
    // Attach two images in horizontal
    cv::Mat outImg(sourceFrame->height(),
                   sourceFrame->width()*2,
                   CV_8UC3);
    for (int i=0; i<sourceFrame->height(); i++)
    {
        for (int j=0; j<sourceFrame->width(); j++)
        {
            
            outImg.at<cv::Vec3b>(i, j) = srcImg.at<cv::Vec3b>(i, j);
            outImg.at<cv::Vec3b>(i, j + sourceFrame->width())
                = trgImg.at<cv::Vec3b>(i, j);
            
        }
    }
    
    return outImg;
    
}
    
cv::Mat
DepthMapForRollingShutter::
plotGeneralStereoImagesWithGT(Frame* const sourceFrame,
                        Frame* const targetFrame,
                        const std::list<PointEpiCurve> *ptsEpiCurve,
                        const std::list<PointEpiCurve> *GT_ptsEpiCurve,
                        PointEpiCurve *best_ptEC)
{
    
    if (activeKeyFrame == 0)
    {
        
        return cv::Mat(height, width*2, CV_8UC3);
        
    }
    
    //--------------------------------------------------------------------------
    // Prepare for image structures
    cv::Mat srcImgContainer(sourceFrame->height(),
                            sourceFrame->width(),
                            CV_32F,
                            const_cast<float*>(sourceFrame->image(0)));
    cv::Mat srcImg;
    srcImgContainer.convertTo(srcImg, CV_8UC1);
    cv::cvtColor(srcImg, srcImg, CV_GRAY2RGB);
    
    cv::Mat trgImgContainer(targetFrame->height(),
                            targetFrame->width(),
                            CV_32F,
                            const_cast<float*>(targetFrame->image(0)));
    cv::Mat trgImg;
    trgImgContainer.convertTo(trgImg, CV_8UC1);
    cv::cvtColor(trgImg, trgImg, CV_GRAY2RGB);
    
    //--------------------------------------------------------------------------
    // Draw epipoar lines and normals in the target
    // for (int i=0; i<(int)ptsEpiCurve->size() - 1; i++)
    std::list<PointEpiCurve>::const_iterator ci = ptsEpiCurve->begin();
    std::list<PointEpiCurve>::const_iterator ci2 = ci;
	 std::advance(ci2, 1);
    for (;
			ci2 != ptsEpiCurve->end();
			++ci, ++ci2)
    {
        
        float x = (*ci).xTarget;
        float y = (*ci).yTarget;
        trgImg.at<cv::Vec3b>(y, x) = cv::Vec3b(0,0,255);
        float nepx = (*ci).nepx;
        float nepy = (*ci).nepy;
        float rescaleFactor = (*ci).rescaleFactor;
        float incx = nepx*rescaleFactor;
        float incy = nepy*rescaleFactor;
        
        // Draw epipolar points in the target image in green (Estimate)
        cv::circle(trgImg, cv::Point(x, y), 1, cv::Scalar(0,255,0));
        float x2 = (*ci2).xTarget;
        float y2 = (*ci2).yTarget;
        cv::line(trgImg, cv::Point(x, y),
                 cv::Point(x2, y2),
                 cv::Scalar(0,255,0), 1); // Estimate green line
        cv::line(trgImg, cv::Point(x, y),
                 cv::Point(x + 2*incx, y),
                 cv::Scalar(255,0,0), 1); // Blue
        cv::line(trgImg,
                 cv::Point(x, y),
                 cv::Point(x, y + 2*incy),
                 cv::Scalar(0,0,255), 1); // Red
        
    }
    
    //--------------------------------------------------------------------------
    // Draw GT epipoar lines and normals in the target
    // for (int i=0; i<(int)GT_ptsEpiCurve->size() - 1; i++)
    ci = GT_ptsEpiCurve->begin();
	 ci2 = ci;
	 std::advance(ci2, 1);
    for (;
			ci2 != GT_ptsEpiCurve->end();
			++ci, ++ci2)
    {
        
        float x = (*ci).xTarget;
        float y = (*ci).yTarget;
        trgImg.at<cv::Vec3b>(y, x) = cv::Vec3b(0,0,255);
        float nepx = (*ci).nepx;
        float nepy = (*ci).nepy;
        float rescaleFactor = (*ci).rescaleFactor;
        float incx = nepx*rescaleFactor;
        float incy = nepy*rescaleFactor;
        
        // Draw epipolar points in the target image in red (Ground truth)
        cv::circle(trgImg, cv::Point(x, y), 1, cv::Scalar(0,0,255));
        float x2 = (*ci2).xTarget;
        float y2 = (*ci2).yTarget;
        cv::line(trgImg, cv::Point(x, y),
                 cv::Point(x2, y2),
                 cv::Scalar(0,0,255), 1); // GT Red line
        cv::line(trgImg, cv::Point(x, y),
                 cv::Point(x + 2*incx, y),
                 cv::Scalar(255,0,0), 1); // Blue
        cv::line(trgImg,
                 cv::Point(x, y),
                 cv::Point(x, y + 2*incy),
                 cv::Scalar(0,0,255), 1); // Red
        
    }
    
    // Draw the best point in the target
    if ((ptsEpiCurve->size() > 0) && (best_ptEC != NULL))
    {
        
        float x = (*best_ptEC).xTarget;
        float y = (*best_ptEC).yTarget;
        cv::circle(trgImg, cv::Point(x, y), 50, cv::Scalar(0,255,255));
        
    }
    
    // Draw a point in the source
    if (ptsEpiCurve->size() > 0)
    {
        
		  ci = ptsEpiCurve->begin(); 
        float x = (*ci).xSource;
        float y = (*ci).ySource;
        cv::circle(srcImg, cv::Point(x, y), 3, cv::Scalar(0,255,0));
        
    }
    
    // Attach two images in horizontal
    cv::Mat outImg(sourceFrame->height(),
                   sourceFrame->width()*2,
                   CV_8UC3);
    for (int i=0; i<sourceFrame->height(); i++)
    {
        for (int j=0; j<sourceFrame->width(); j++)
        {
            
            outImg.at<cv::Vec3b>(i, j) = srcImg.at<cv::Vec3b>(i, j);
            outImg.at<cv::Vec3b>(i, j + sourceFrame->width())
            = trgImg.at<cv::Vec3b>(i, j);
            
        }
    }
    
    return outImg;
    
}
    
Eigen::Matrix<double,3,3> skewsym(Eigen::Matrix<double,3,1> v)
{
    
    Eigen::Matrix<double,3,3> M;
    
    M <<     0,  -v(2),   v(1),
          v(2),      0,  -v(0),
         -v(1),   v(0),      0;
    
    return M;
    
}
    
//------------------------------------------------------------------------------
// Extracts image points on the generalised epipolar line band
//------------------------------------------------------------------------------
// It extracts image points on the band of an generalised epipolar line.
// This version differs from extractPointsOnGeneralEpipolarLine()
// and it provides a better handling of degenerate cases.
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
DepthMapForRollingShutter::
extractPointsOnGeneralEpipolarLineBand(const int x,
                                   const int y,
                                   Frame* const ref,
                                   std::list<PointEpiCurve> *ptsEpiCurve,
                                   RunningStats* const stats,
                                   bool useGT)
{
    
    Frame *sourceImage = activeKeyFrame;
    Frame *targetImage = ref;
    
    //--------------------------------------------------------------------------
    // Sweeping rows in the target image - find a plane for the row
    
    // Row y pose in source image
    Eigen::Matrix<double,4,4> world_to_ySource;
    if (useGT == true)
    {
        // Source pose from ground truth
        world_to_ySource = sourceImage->getGTPoseAtRow(y);
        
    }
    else
    {
        
        // Source pose from estimate
        world_to_ySource = sourceImage->getPoseAtRow(y);
        
        // Warning if spline is not available
        if (sourceImage->id() != 0 && sourceImage->isSplineValid == false)
        {
            
            printf("Warning: source image %d has no valid spline\n",
                   sourceImage->id());
            
        }
        
    }
    
    // Sweeping each row in target image
    for (int rowTarget = 3; rowTarget <= height - 3; rowTarget++)
    {
        
        //----------------------------------------------------------------------
        // Find an intersection with the ray of the given point
        // in the source image
        // Project the intersection onto the target image row
        //
        // Alternative: find an intersection of the row and the normal
        //              epipolar line
        
        // Row pose in target image with respect to row y in source image
        Eigen::Matrix<double,4,4> world_to_rowTarget;
        if (useGT == true)
        {
            
            // Target pose from ground truth
            world_to_rowTarget = targetImage->getGTPoseAtRow(rowTarget);
            
        }
        else
        {
            
            // Target pose from estimate
            world_to_rowTarget = targetImage->getPoseAtRow(rowTarget);
            
            // Warning if spline is not available
            if (targetImage->isSplineValid == false)
            {
                
                printf("Warning: target image %d has no valid spline\n",
                       targetImage->id());
                
            }
            
        }
        
        Eigen::Matrix<double,4,4> rowTarget_to_ySource
            = world_to_ySource*world_to_rowTarget.inverse();

        Sophus::SE3d Xi_rowTarget_to_ySource
            = Sophus::SE3d(rowTarget_to_ySource);
        
        // Motion from source to target
        Sophus::SE3d Xi_ySource_to_rowTarget
            = Xi_rowTarget_to_ySource.inverse();
        
#if 0 //TEST_THIS_MOTION TARGET with respect to SOURCE
        Eigen::Matrix<double,3,1> t = Xi_rowTarget_to_ySource.translation();
        Eigen::Matrix<double,3,3> R = Xi_rowTarget_to_ySource.rotationMatrix();
#else
        Eigen::Matrix<double,3,1> t = Xi_ySource_to_rowTarget.translation();
        Eigen::Matrix<double,3,3> R = Xi_ySource_to_rowTarget.rotationMatrix();
#endif
        // Fundamental matrix
        Eigen::Matrix<double,3,3> E = skewsym(t)*R;
        Eigen::Matrix<double,3,3> Kd = K.cast<double>();
        Eigen::Matrix<double,3,3> F = Kd.inverse().transpose()*E*Kd.inverse();
        
        // Normal epipolar line in target
        Eigen::Matrix<double,3,1> eline = F*Eigen::Matrix<double,3,1>(x,y,1);
        
        // Band (top/bottom) of epipolar line
        // Given ax + by + c = 0, the addition of delta in y-direction
        // gives ax + by + c - b*delta = 0.
        Eigen::Matrix<double,3,1> upper_eline
            = eline + Eigen::Matrix<double,3,1>(0,0,-1*eline[1]);
        Eigen::Matrix<double,3,1> lower_eline
            = eline + Eigen::Matrix<double,3,1>(0,0,+1*eline[1]);
        
        // Intersection with the current rowTarget
        Eigen::Matrix<double,3,1> pInter
        = eline.cross(Eigen::Matrix<double,3,1>(0, 1, -rowTarget));
        Eigen::Matrix<double,3,1> upper_pInter
        = upper_eline.cross(Eigen::Matrix<double,3,1>(0, 1, -rowTarget));
        Eigen::Matrix<double,3,1> lower_pInter
        = lower_eline.cross(Eigen::Matrix<double,3,1>(0, 1, -rowTarget));

#if 0 // USE HCOORD
        
        // Homogeneous coords
        Eigen::Matrix<double,3,1> pInter_h = hcoord(pInter);
        Eigen::Matrix<double,3,1> upper_pInter_h = hcoord(upper_pInter);
        Eigen::Matrix<double,3,1> lower_pInter_h = hcoord(lower_pInter);
#else // USE EXPLICIT DIVISION
        
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
#endif
        
        // Find a normal to the epipolar line
        double nepx, nepy;
#if 0 // CONSTANT nepx and nepy
        nepx = 0.7071;
        nepy = 0.7071;
#else
        if (makeAndCheckEPL_forGeneralEpipolarLine(x, y, t,
                                                   &nepx, &nepy,
                                                   stats) == false)
        {
            
            // No normal of the epipolar line found
            
#if 1 // SKIP IF FAILED IN NORMAL COMPUTATION
            continue;
#endif
            
        }
#endif


        // Check points with +/- infinity test
        std::list<PointEpiCurve> partCurve =
            buildCurveForPointAtInfinity(x, y, nepx, nepy,
                                         Xi_ySource_to_rowTarget,
                                         rowTarget,
                                         upper_pInter,
                                         lower_pInter);
        if (partCurve.size() != 0)
        {
        
            for (std::list<PointEpiCurve>
                     ::const_iterator ci = partCurve.begin();
                     ci != partCurve.end(); ++ci)
            {
                ptsEpiCurve->push_back(*ci);
            }
            continue;
        
        }
        
        // Check boundary
        if (pInter_h[0] < 3 || pInter_h[0] > width - 3)
            continue;
        if (pInter_h[1] < 3 || pInter_h[1] > height - 3)
            continue;
        if (upper_pInter_h[0] < 3 || upper_pInter_h[0] > width - 3)
            continue;
        if (upper_pInter_h[1] < 3 || upper_pInter_h[1] > height - 3)
            continue;
        if (lower_pInter_h[0] < 3 || lower_pInter_h[0] > width - 3)
            continue;
        if (lower_pInter_h[1] < 3 || lower_pInter_h[1] > height - 3)
            continue;

        // Add intersections to curve
        pushPointToEpiCurve(x, y,
                            upper_pInter_h[0], upper_pInter_h[1],
                            nepx, nepy,
                            Xi_ySource_to_rowTarget,
                            ptsEpiCurve);
        pushPointToEpiCurve(x, y,
                            lower_pInter_h[0], lower_pInter_h[1],
                            nepx, nepy,
                            Xi_ySource_to_rowTarget,
                            ptsEpiCurve);

        
        
        
    } // end of for
    
    if (ptsEpiCurve->size() > 0)
    {
        
        // Sort by the inverse depth
        //sortPointsOnGeneralEpipolarLineByIDepth(ptsEpiCurve);
        
        return true;
        
    }
    else
    {
        
        return false;
        
    }
    
}
    
// xs : x source
// ys : y source
// xt : x target
// yt : y target
// nepx : normal of x in the epipolar line
// nepy : normal of y in the epipolar line
// motion : motion from source to target
void DepthMapForRollingShutter::
    pushPointToEpiCurve(double xs, double ys,
                        double xt, double yt,
                        double nepx, double nepy,
                        Sophus::SE3d motion,
                        std::list<PointEpiCurve> *ptsEpiCurve)
{
    
    PointEpiCurve pt(xs, ys, xt, yt, nepx, nepy, motion);

#if 0 // COMPUTE INV DEPTH
    // Set an inverse depth of the point in the source
    pt.computeInvDepth(fxi, fyi, cxi, cyi);
    
#if 0 // POSITIVE and minimum depth ONLY
    
    // Save to the vector only if the inverse depth is positive
    if ((pt.invDepth > 1.0e-6) && (pt.invDepth < 1.0e+6))
    {
        
        ptsEpiCurve->push_back(pt);
        
    }
#else // ALL (POSITIVE + NEGATIVE)
    
    ptsEpiCurve->push_back(pt);
    
#endif
    
#else // DO NOT COMPUTE INV DEPTH HERE but PUSH BACK
    
    ptsEpiCurve->push_back(pt);
    
#endif
    
}

//------------------------------------------------------------------------------
// Extracts image points on the generalised epipolar line
//------------------------------------------------------------------------------
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
DepthMapForRollingShutter::
extractPointsOnGeneralEpipolarLine(const int x,
                    const int y,
                    Frame* const ref,
                    std::vector<PointEpiCurve> *ptsEpiCurve,
                    RunningStats* const stats,
                    bool useGT)
{
    
    Frame *sourceImage = activeKeyFrame;
    Frame *targetImage = ref;
    
    //--------------------------------------------------------------------------
    // Sweeping rows in the target image - find a plane for the row
    
    // Row y pose in source image
    Eigen::Matrix<double,4,4> world_to_ySource, GT_world_to_ySource;
    if (useGT == true)
    {
        // Source pose from ground truth
        world_to_ySource = sourceImage->getGTPoseAtRow(y);
        
    }
    else
    {
        
        // Source pose from estimate
        world_to_ySource = sourceImage->getPoseAtRow(y);
        
        // GT source pose
        GT_world_to_ySource = sourceImage->getGTPoseAtRow(y);
        
        // Compare with GT
        Sophus::SE3d diffPose =
            Sophus::SE3d(world_to_ySource.inverse()*GT_world_to_ySource);
        if (diffPose.log().sum() > 5.0)
        {
            std::cout << "  pose:"
                      << Sophus::SE3d(world_to_ySource).log()
                      << std::endl;
            std::cout << "GTpose:"
                      << Sophus::SE3d(GT_world_to_ySource).log()
                      << std::endl;
        }
        
    }
    
    // Sweeping each row in target image
    for (int rowTarget = 3; rowTarget <= height - 3; rowTarget++)
    {
        
        //----------------------------------------------------------------------
        // Find an intersection with the ray of the given point
        // in the source image
        // Project the intersection onto the target image row
        //
        // Alternative: find an intersection of the row and the normal
        //              epipolar line
        
        // Row pose in target image with respect to row y in source image
        Eigen::Matrix<double,4,4> world_to_rowTarget, GT_world_to_rowTarget;
        if (useGT == true)
        {
        
            // Target pose from ground truth
            world_to_rowTarget = targetImage->getGTPoseAtRow(rowTarget);
        
        }
        else
        {
        
            // Target pose from estimate
            world_to_rowTarget = targetImage->getPoseAtRow(rowTarget);
            
            // GT source pose
            GT_world_to_rowTarget = targetImage->getGTPoseAtRow(rowTarget);
            
            // Compare with GT
            Sophus::SE3d diffPose =
            Sophus::SE3d(world_to_rowTarget.inverse()*GT_world_to_rowTarget);
            if (diffPose.log().sum() > 5.0)
            {
                std::cout << "  pose:"
                << Sophus::SE3d(world_to_rowTarget).log()
                << std::endl;
                std::cout << "GTpose:"
                << Sophus::SE3d(GT_world_to_rowTarget).log()
                << std::endl;
            }
            
        }
        
        Sophus::SE3d Xi_world_to_rowTarget = Sophus::SE3d(world_to_rowTarget);
        if (fabs(Xi_world_to_rowTarget.unit_quaternion().w()) < 1.e-6)
        {
            std::cout << world_to_rowTarget << std::endl;
        }
        
        Eigen::Matrix<double,4,4> rowTarget_to_ySource
            = world_to_ySource*world_to_rowTarget.inverse();
        
        Sophus::SE3d Xi2_rowTarget_to_ySource =
            Sophus::SE3d(rowTarget_to_ySource);

        // Motion from source to target
        Sophus::SE3d Xi_ySource_to_rowTarget
            = Xi2_rowTarget_to_ySource.inverse();
        
#if 0 //TEST_THIS_MOTION TARGET with respect to SOURCE
        Eigen::Matrix<double,3,1> t = Xi_rowTarget_to_ySource.translation();
        Eigen::Matrix<double,3,3> R = Xi_rowTarget_to_ySource.rotationMatrix();
#else
        Eigen::Matrix<double,3,1> t = Xi_ySource_to_rowTarget.translation();
        Eigen::Matrix<double,3,3> R = Xi_ySource_to_rowTarget.rotationMatrix();
#endif
        // Fundamental matrix
        Eigen::Matrix<double,3,3> E = skewsym(t)*R;
        Eigen::Matrix<double,3,3> Kd = K.cast<double>();
        Eigen::Matrix<double,3,3> F = Kd.inverse().transpose()*E*Kd.inverse();

        // Normal epipolar line in target
        Eigen::Matrix<double,3,1> epLineTarget =
            F*Eigen::Matrix<double,3,1>(x,y,1);

        // Intersection with the current rowTarget
        Eigen::Matrix<double,3,1> rowPointOnEpCurve
            = epLineTarget.cross(Eigen::Matrix<double,3,1>(0, 1, -rowTarget));
        
        // Check the intersection exist
        double xTarget_double = rowPointOnEpCurve[0]/rowPointOnEpCurve[2];
        if (isnan(xTarget_double) == true)
        {
            
            // The intersection point does not exist
            // since the ray is parallel to the plane.
            continue;
            
        }
        
        //----------------------------------------------------------------------
        // Store the projection into buffer
        int xTarget = (int)floor(xTarget_double);
        
        // Check xTarget value
        if ((xTarget < 3) || (xTarget > width - 3))
        {
            
            // Out of the width
            continue;
        
        }
        
        //----------------------------------------------------------------------
        // Find a normal to the epipolar line
        double nepx, nepy;
        if (makeAndCheckEPL_forGeneralEpipolarLine(x, y, t,
                                                   &nepx, &nepy,
                                                   stats) == false)
        {
            
            // No normal of the epipolar line found
            continue;
            
        }
        
        //----------------------------------------------------------------------
        // Store points into a vector of PointEpiCurve
        PointEpiCurve pt(x, y,
                         xTarget_double, (double)rowTarget,
                         nepx, nepy,
                         Xi_ySource_to_rowTarget);
        
        // Set an inverse depth of the point in the source
        // THIS IS INVALID SINCE NO INCX and INCY available
        // pt.computeInvDepth(fxi, fyi, cxi, cyi);
        
#if 0 // POSITIVE and minimum depth ONLY
        
        // Save to the vector only if the inverse depth is positive
        if ((pt.invDepth > 1.0e-6) && (pt.invDepth < 1.0e+6))
        {
            
            
            ptsEpiCurve->push_back(pt);
            
        }
#else // ALL (POSITIVE + NEGATIVE)
            ptsEpiCurve->push_back(pt);
#endif
        
    }
    
    if (ptsEpiCurve->size() > 0)
    {
        
        // Sort by the inverse depth
        sortPointsOnGeneralEpipolarLineByIDepth(ptsEpiCurve);
        
        return true;
        
    }
    else
    {
        
        return false;
        
    }
    
}

//------------------------------------------------------------------------------
// Sort points on epipolar curve by depth
// Since a vector of ptsEpiCurve structures extracted from the method
// extractPointsOnGeneralEpipolarLine() is in the order of rows in the target
// image but not in the inverse depth of a ray in the source camera,
// here we sort elements in the vector in order of the inverse depth
//------------------------------------------------------------------------------
//      ptsEpiCurve : coordinates of points on generalised epipolar line
//                    in the target image.
//------------------------------------------------------------------------------
void DepthMapForRollingShutter::
sortPointsOnGeneralEpipolarLineByIDepth(std::vector<PointEpiCurve> *ptsEpiCurve)
{

    std::sort(ptsEpiCurve->begin(), ptsEpiCurve->end());
    
}
    
 
//------------------------------------------------------------------------------
// Copy image bands to buffer
//------------------------------------------------------------------------------
//     prior_idepth : Prior of inverse depth
//      targetImage : Target image
//      ptsEpiCurve : coordinates of points on generalised epipolar line
//                    in the target image.
//                    Intensity values will be stored into a member variable
//                    of this structure for image bands along
//                    the generalised epipolar line.
//                    The image band has in total 5 pixels which are
//                    from 2 equi-distance point on the generalised epipolar
//                    line.
//            stats : status return
//
// returns
//
//            [int] : number of stored values in the image band (> 0)
//
//------------------------------------------------------------------------------
int
DepthMapForRollingShutter::
copyImageBandsToBuffer(const float min_idepth,
                       const float prior_idepth,
                       float max_idepth,
                       const float *targetImage,
                       std::list<PointEpiCurve> *ptsEpiCurve,
                       RunningStats* const stats)
{
    
    int numStoredImageBand = 0;
    
    // For each point on the epipolar curve or generalised epipolar line
    for (std::list<PointEpiCurve>::iterator it = ptsEpiCurve->begin();
         it != ptsEpiCurve->end(); )
    {

        // Point on epipolar curve
        //PointEpiCurve *ptEC = &((*ptsEpiCurve)[i]);
        PointEpiCurve *ptEC = &(*it);

        // Motion of target to source for the point on epipolar curve
        const Sophus::SE3d Xi = ptEC->Xi;
        
        // Source point
        float u = (float)ptEC->xSource;
        float v = (float)ptEC->ySource;

        // Calculate the point at infinity on the normal epipolar line
        Eigen::Vector3f KinvP = Eigen::Vector3f(fxi*u + cxi,
                                                fyi*v + cyi,
                                                1.0f);
        Eigen::Vector3f pInf = K*Xi.rotationMatrix().cast<float>()*KinvP;
        Eigen::Vector3f Kt = K*Xi.translation().cast<float>();
        Eigen::Vector3f pReal = pInf/prior_idepth + Kt;

        // Rescale factor
        float rescaleFactor = pReal[2]*prior_idepth;
#if 0 // CONSTANT_RESCALE_FACTOR
        rescaleFactor = 1.0;
#endif
        
        if(!(rescaleFactor > 0.7f && rescaleFactor < 1.4f))
        {
        
            if (enablePrintDebugInfo)
            {
                
                stats->num_stereo_rescale_oob++;
                
            }
            
            // Out of rescale factor bounds
            it = ptsEpiCurve->erase(it);
            continue;
        
        }

        // Normal of the epipolar curve
        float normalEpx = (float)ptEC->nepx;
        float normalEpy = (float)ptEC->nepy;
        
        // Prepare firstX, firstY, lastX, lastY
        float firstX = u - 2*normalEpx*rescaleFactor;
        float firstY = v - 2*normalEpy*rescaleFactor;
        float lastX = u + 2*normalEpx*rescaleFactor;
        float lastY = v + 2*normalEpy*rescaleFactor;
        if (firstX <= 0 ||
            firstX >= width - 2 ||
            firstY <= 0 ||
            firstY >= height - 2 ||
            lastX <= 0 ||
            lastX >= width - 2 ||
            lastY <= 0 ||
            lastY >= height - 2)
        {
            
            // Out of image boundary
            it = ptsEpiCurve->erase(it);
            continue;
        
        }
        
        //----------------------------------------------------------------------
        // Check other conditions
        Eigen::Vector3f pClose = pInf + Kt*max_idepth;
        // if the assumed close-point lies behind the
        // image, have to change that.
        if(pClose[2] < 0.001f)
        {
            max_idepth = (0.001f-pInf[2]) / Kt[2];
            pClose = pInf + Kt*max_idepth;
        }
        pClose = pClose / pClose[2]; // pos in new image of point (xy), assuming max_idepth
        
        Eigen::Vector3f pFar = pInf + Kt*min_idepth;
        // if the assumed far-point lies behind the image or closter than the near-point,
        // we moved past the Point it and should stop.
        if(pFar[2] < 0.001f || max_idepth < min_idepth)
        {
            if(enablePrintDebugInfo) stats->num_stereo_inf_oob++;
            it = ptsEpiCurve->erase(it);
            continue;
        }
        pFar = pFar / pFar[2]; // pos in new image of point (xy), assuming min_idepth
        
        
        // check for nan due to eg division by zero.
        if(isnanf((float)(pFar[0]+pClose[0])))
        {
            // Out of image boundary
            it = ptsEpiCurve->erase(it);
            continue;
        }
        
        // calculate increments in which we will step through the epipolar line.
        // they are sampleDist (or half sample dist) long
        float incx = pClose[0] - pFar[0];
        float incy = pClose[1] - pFar[1];
        float eplLength = sqrtf(incx*incx+incy*incy);
        float incx_epl = incx;
        float incy_epl = incy;
        
        if(!(eplLength > 0) || std::isinf(eplLength))
        {
        
            it = ptsEpiCurve->erase(it);
            continue;
        
        }
    
        if(eplLength > MAX_EPL_LENGTH_CROP)
        {
            pClose[0] = pFar[0] + incx*MAX_EPL_LENGTH_CROP/eplLength;
            pClose[1] = pFar[1] + incy*MAX_EPL_LENGTH_CROP/eplLength;
        }
        
        incx *= GRADIENT_SAMPLE_DIST/eplLength;
        incy *= GRADIENT_SAMPLE_DIST/eplLength;
        
        // extend one sample_dist to left & right.
        pFar[0] -= incx;
        pFar[1] -= incy;
        pClose[0] += incx;
        pClose[1] += incy;
        
        
        // make epl long enough (pad a little bit).
        if(eplLength < MIN_EPL_LENGTH_CROP)
        {
            float pad = (MIN_EPL_LENGTH_CROP - (eplLength)) / 2.0f;
            pFar[0] -= incx*pad;
            pFar[1] -= incy*pad;
            
            pClose[0] += incx*pad;
            pClose[1] += incy*pad;
        }
        
        // if inf point is outside of image: skip pixel.
        if(
           pFar[0] <= SAMPLE_POINT_TO_BORDER ||
           pFar[0] >= width-SAMPLE_POINT_TO_BORDER ||
           pFar[1] <= SAMPLE_POINT_TO_BORDER ||
           pFar[1] >= height-SAMPLE_POINT_TO_BORDER)
        {
            if(enablePrintDebugInfo) stats->num_stereo_inf_oob++;
            it = ptsEpiCurve->erase(it);
            continue;
        }
        
        
        
        // if near point is outside: move inside, and test length again.
        if(
           pClose[0] <= SAMPLE_POINT_TO_BORDER ||
           pClose[0] >= width-SAMPLE_POINT_TO_BORDER ||
           pClose[1] <= SAMPLE_POINT_TO_BORDER ||
           pClose[1] >= height-SAMPLE_POINT_TO_BORDER)
        {
            if(pClose[0] <= SAMPLE_POINT_TO_BORDER)
            {
                float toAdd = (SAMPLE_POINT_TO_BORDER - pClose[0]) / incx;
                pClose[0] += toAdd * incx;
                pClose[1] += toAdd * incy;
            }
            else if(pClose[0] >= width-SAMPLE_POINT_TO_BORDER)
            {
                float toAdd = (width-SAMPLE_POINT_TO_BORDER - pClose[0]) / incx;
                pClose[0] += toAdd * incx;
                pClose[1] += toAdd * incy;
            }
            
            if(pClose[1] <= SAMPLE_POINT_TO_BORDER)
            {
                float toAdd = (SAMPLE_POINT_TO_BORDER - pClose[1]) / incy;
                pClose[0] += toAdd * incx;
                pClose[1] += toAdd * incy;
            }
            else if(pClose[1] >= height-SAMPLE_POINT_TO_BORDER)
            {
                float toAdd = (height-SAMPLE_POINT_TO_BORDER - pClose[1]) / incy;
                pClose[0] += toAdd * incx;
                pClose[1] += toAdd * incy;
            }
            
            // get new epl length
            float fincx = pClose[0] - pFar[0];
            float fincy = pClose[1] - pFar[1];
            float newEplLength = sqrtf(fincx*fincx+fincy*fincy);
            
            // test again
            if(
               pClose[0] <= SAMPLE_POINT_TO_BORDER ||
               pClose[0] >= width-SAMPLE_POINT_TO_BORDER ||
               pClose[1] <= SAMPLE_POINT_TO_BORDER ||
               pClose[1] >= height-SAMPLE_POINT_TO_BORDER ||
               newEplLength < 8.0f
               )
            {
                if(enablePrintDebugInfo) stats->num_stereo_near_oob++;
                it = ptsEpiCurve->erase(it);
                continue;
            }
            
            
        } // end if pClose[0] ...


        //----------------------------------------------------------------------
        // Calculate values in source to search for
        
        // Check the boundary
        if (u - 2*normalEpx*rescaleFactor <= 2 ||
            u - 2*normalEpx*rescaleFactor >= width - 3||
            u + 2*normalEpx*rescaleFactor <= 2 ||
            u + 2*normalEpx*rescaleFactor >= width - 3||
            v - 2*normalEpy*rescaleFactor <= 2 ||
            v - 2*normalEpy*rescaleFactor >= height - 3||
            v + 2*normalEpy*rescaleFactor <= 2 ||
            v + 2*normalEpy*rescaleFactor >= height - 3)
        {
            
            // Out of image boundary
            it = ptsEpiCurve->erase(it);
            continue;
            
        }
        
        //----------------------------------------------------------------------
        // Get an image band in target
        float xTarget = (float)ptEC->xTarget;
        float yTarget = (float)ptEC->yTarget;
        float abs_incx = (float)fabs(incx);
        float abs_incy = (float)fabs(incy);
        float abs_nepx = (float)fabs(normalEpx);
        float abs_nepy = (float)fabs(normalEpy);
        
        // Check the boundary
        if (xTarget - 4*abs_nepx*rescaleFactor <= 3 ||
            xTarget - 4*abs_nepx*rescaleFactor >= width - 4 ||
            xTarget + 4*abs_nepx*rescaleFactor <= 3 ||
            xTarget + 4*abs_nepx*rescaleFactor >= width - 4 ||
            yTarget - 4*abs_nepy*rescaleFactor <= 3 ||
            yTarget - 4*abs_nepy*rescaleFactor >= height - 4 ||
            yTarget + 4*abs_nepy*rescaleFactor <= 3 ||
            yTarget + 4*abs_nepy*rescaleFactor >= height - 4 )
        {
            
            // Out of image boundary
            it = ptsEpiCurve->erase(it);
            continue;
            
        }

        // Check the boundary
        if (xTarget - 4*abs_incx <= 3 ||
            xTarget - 4*abs_incx >= width - 4 ||
            xTarget + 4*abs_incx <= 3 ||
            xTarget + 4*abs_incx >= width - 4 ||
            yTarget - 4*abs_incy <= 3 ||
            yTarget - 4*abs_incy >= height - 4 ||
            yTarget + 4*abs_incy <= 3 ||
            yTarget + 4*abs_incy >= height - 4 )
        {

            // Out of image boundary
            it = ptsEpiCurve->erase(it);
            continue;

        }

#if 1 // PIXELS ALONG THE EPIPOLAR LINE
        
        // Get intensity
//
//        float trgVal_m2 = getInterpolatedElement(targetImage,
//                                         xTarget - 2.0f*abs_incx,
//                                         yTarget - 2.0f*abs_incy,
//                                         width);
//        float trgVal_m1 = getInterpolatedElement(targetImage,
//                                         xTarget - 1.0f*abs_incx,
//                                         yTarget - 1.0f*abs_incy,
//                                         width);
//        float trgVal = getInterpolatedElement(targetImage,
//                                              xTarget,
//                                              yTarget,
//                                              width);
//        float trgVal_p1 = getInterpolatedElement(targetImage,
//                                         xTarget + 1.0f*abs_incx,
//                                         yTarget + 1.0f*abs_incy,
//                                         width);
//        float trgVal_p2 = getInterpolatedElement(targetImage,
//                                         xTarget + 2.0f*abs_incx,
//                                         yTarget + 2.0f*abs_incy,
//                                         width);
        
        float nepx = normalEpx;
        float nepy = normalEpy;
        
        if (xTarget <= 2 || xTarget >= width - 2 ||
            yTarget <= 2 || yTarget >= height - 2 ||
            xTarget - 1.0f*nepx*rescaleFactor <= 2 ||
            xTarget - 1.0f*nepx*rescaleFactor >= width - 2 ||
            yTarget - 1.0f*nepy*rescaleFactor <= 2 ||
            yTarget - 1.0f*nepy*rescaleFactor >= height - 2 ||
            xTarget + 1.0f*nepx*rescaleFactor <= 2 ||
            xTarget + 1.0f*nepx*rescaleFactor >= width - 2 ||
            yTarget + 1.0f*nepy*rescaleFactor <= 2 ||
            yTarget + 1.0f*nepy*rescaleFactor >= height - 2 ||
            xTarget - 2.0f*nepx*rescaleFactor <= 2 ||
            xTarget - 2.0f*nepx*rescaleFactor >= width - 2 ||
            yTarget - 2.0f*nepy*rescaleFactor <= 2 ||
            yTarget - 2.0f*nepy*rescaleFactor >= height - 2 ||
            xTarget + 2.0f*nepx*rescaleFactor <= 2 ||
            xTarget + 2.0f*nepx*rescaleFactor >= width - 2 ||
            yTarget + 2.0f*nepy*rescaleFactor <= 2 ||
            yTarget + 2.0f*nepy*rescaleFactor >= height - 2)
        {
            
            // Out of image boundary
            it = ptsEpiCurve->erase(it);
            continue;
            
        }

        float trgVal_m2 = getInterpolatedElement(targetImage,
                                                 xTarget - 2.0f*nepx*rescaleFactor,
                                                 yTarget - 2.0f*nepy*rescaleFactor,
                                                 width);
        float trgVal_m1 = getInterpolatedElement(targetImage,
                                                 xTarget - 1.0f*nepx*rescaleFactor,
                                                 yTarget - 1.0f*nepy*rescaleFactor,
                                                 width);
        float trgVal = getInterpolatedElement(targetImage,
                                              xTarget,
                                              yTarget,
                                              width);
        float trgVal_p1 = getInterpolatedElement(targetImage,
                                                 xTarget + 1.0f*nepx*rescaleFactor,
                                                 yTarget + 1.0f*nepy*rescaleFactor,
                                                 width);
        float trgVal_p2 = getInterpolatedElement(targetImage,
                                                 xTarget + 2.0f*nepx*rescaleFactor,
                                                 yTarget + 2.0f*nepy*rescaleFactor,
                                                 width);
        
#else
#if 0 // CROSS by normal and rescale
        // PIXELS ON A 3x3 WINDOW
        //         (0, -1)
        // (-1,0)  (0,  0)   (1,0)
        //         (0,  1)
        // Get intensity
        float trgVal_m2 = getInterpolatedElement(targetImage,
                                                 xTarget - 0*normalEpx*rescaleFactor,
                                                 yTarget - 1*normalEpy*rescaleFactor,
                                                 width);
        float trgVal_m1 = getInterpolatedElement(targetImage,
                                                 xTarget - 1*normalEpx*rescaleFactor,
                                                 yTarget - 0*normalEpy*rescaleFactor,
                                                 width);
        float trgVal = getInterpolatedElement(targetImage,
                                              xTarget,
                                              yTarget,
                                              width);
        float trgVal_p1 = getInterpolatedElement(targetImage,
                                                 xTarget + 1*normalEpx*rescaleFactor,
                                                 yTarget + 0*normalEpy*rescaleFactor,
                                                 width);
        float trgVal_p2 = getInterpolatedElement(targetImage,
                                                 xTarget + 0*normalEpx*rescaleFactor,
                                                 yTarget + 1*normalEpy*rescaleFactor,
                                                 width);
#else   // CROSS without normal and rescale
        // PIXELS ON A 3x3 WINDOW
        //         (0, -1)
        // (-1,0)  (0,  0)   (1,0)
        //         (0,  1)
        // Get intensity
        float abs_incx = (float)fabs(incx);
        float abs_incy = (float)fabs(incy);
        float len_incxincy = sqrt(incx*incx + incy*incy);
        float trgVal_m2 = getInterpolatedElement(targetImage,
                                                 xTarget - 0.0*len_incxincy,
                                                 yTarget - 3.0*len_incxincy,
                                                 width);
        float trgVal_m1 = getInterpolatedElement(targetImage,
                                                 xTarget - 3.0*len_incxincy,
                                                 yTarget - 0.0*len_incxincy,
                                                 width);
        float trgVal = getInterpolatedElement(targetImage,
                                              xTarget,
                                              yTarget,
                                              width);
        float trgVal_p1 = getInterpolatedElement(targetImage,
                                                 xTarget + 3.0*len_incxincy,
                                                 yTarget + 0.0*len_incxincy,
                                                 width);
        float trgVal_p2 = getInterpolatedElement(targetImage,
                                                 xTarget + 0.0*len_incxincy,
                                                 yTarget + 3.0*len_incxincy,
                                                 width);

#endif
#endif
        
        //----------------------------------------------------------------------
        // Set incx and incy parameter and compute inverse depth for checking
        ptEC->setIncXandY(incx, incy);
#if 1 // CHECK MIN_MAX_IDEPTH
        // Check min_idepth and max_idepth (Dependency incx and incy)
        ptEC->computeInvDepth(fxi, fyi, cxi, cyi);
        if (ptEC->invDepth < min_idepth ||
            ptEC->invDepth > max_idepth)
        {
            
            it = ptsEpiCurve->erase(it);
            continue;
            
        }
#endif
        
        //----------------------------------------------------------------------
        // Copy the image band and set intensity values
        Eigen::Matrix<float,5,1> vals;
        vals << trgVal_m2, trgVal_m1, trgVal, trgVal_p1, trgVal_p2;
        ptEC->setVals(vals);

        // Set extra
        ptEC->setRescaleFactor(rescaleFactor);
        ptEC->setIncXandYForEpl(incx_epl, incy_epl);
        ptEC->setEplLength(eplLength);

        numStoredImageBand++;
        
        // Increase iterator
        it++;

    }
    
//    //--------------------------------------------------------------------------
//    // Remove invalid data
//    for (std::vector<PointEpiCurve>::iterator it = ptsEpiCurve->begin();
//         it != ptsEpiCurve->end();
//         ++it)
//    {
//
//        if (it->isValsAvailable == false)
//        {
//            ptsEpiCurve->erase(it);
//        }
//        
//    }
    
#ifdef DEBUG
    
    if (numStoredImageBand > 0)
    {
        
        printf("numStoredImageBand = %d\n", numStoredImageBand);
        
    }
    
#endif
    
    return numStoredImageBand;
    
}
    
inline Eigen::Matrix<float,5,1>
DepthMapForRollingShutter::
getSourceValsFromEpiPtCurve(const PointEpiCurve* ptEC)
{
    
    // Rescale factor
    float rescaleFactor = ptEC->rescaleFactor;
    
    // A set of image intensities in the source image to be searched for
    float u = (float)ptEC->xSource;
    float v = (float)ptEC->ySource;
    float nepx = (float)ptEC->nepx;
    float nepy = (float)ptEC->nepy;
    
#if 1 // PIXELS ALONG THE EPIPOLAR LINE
    if (u <= 2 || u >= width - 2 ||
        v <= 2 || v >= height - 2 ||
        u - 1.0f*nepx*rescaleFactor <= 2 ||
        u - 1.0f*nepx*rescaleFactor >= width - 2 ||
        v - 1.0f*nepy*rescaleFactor <= 2 ||
        v - 1.0f*nepy*rescaleFactor >= height - 2 ||
        u + 1.0f*nepx*rescaleFactor <= 2 ||
        u + 1.0f*nepx*rescaleFactor >= width - 2 ||
        v + 1.0f*nepy*rescaleFactor <= 2 ||
        v + 1.0f*nepy*rescaleFactor >= height - 2 ||
        u - 2.0f*nepx*rescaleFactor <= 2 ||
        u - 2.0f*nepx*rescaleFactor >= width - 2 ||
        v - 2.0f*nepy*rescaleFactor <= 2 ||
        v - 2.0f*nepy*rescaleFactor >= height - 2 ||
        u + 2.0f*nepx*rescaleFactor <= 2 ||
        u + 2.0f*nepx*rescaleFactor >= width - 2 ||
        v + 2.0f*nepy*rescaleFactor <= 2 ||
        v + 2.0f*nepy*rescaleFactor >= height - 2)
    {
        
        Eigen::Matrix<float,5,1> sourceVals;
        sourceVals << -1, -1, -1, -1, -1;
        return sourceVals;
        
    }
        
    float realVal_mmxy2 = getInterpolatedElement(activeKeyFrameImageData,
                                              u - 2.0f*nepx*rescaleFactor,
                                              v - 2.0f*nepy*rescaleFactor,
                                              width);
    float realVal_mmxy1 = getInterpolatedElement(activeKeyFrameImageData,
                                              u - 1.0f*nepx*rescaleFactor,
                                              v - 1.0f*nepy*rescaleFactor,
                                              width);
    float realVal = getInterpolatedElement(activeKeyFrameImageData,
                                           u,
                                           v,
                                           width);
    float realVal_ppxy1 = getInterpolatedElement(activeKeyFrameImageData,
                                              u + 1.0f*nepx*rescaleFactor,
                                              v + 1.0f*nepy*rescaleFactor,
                                              width);
    float realVal_ppxy2 = getInterpolatedElement(activeKeyFrameImageData,
                                              u + 2.0f*nepx*rescaleFactor,
                                              v + 2.0f*nepy*rescaleFactor,
                                              width);
//    float realVal_mpyx1 = getInterpolatedElement(activeKeyFrameImageData,
//                                              u - 1.0*abs_nepy*rescaleFactor,
//                                              v + 1.0*abs_nepx*rescaleFactor,
//                                              width);
//    float realVal_mpyx2 = getInterpolatedElement(activeKeyFrameImageData,
//                                              u - 2.0*abs_nepy*rescaleFactor,
//                                              v + 2.0*abs_nepx*rescaleFactor,
//                                              width);
//    float realVal_ppyx1 = getInterpolatedElement(activeKeyFrameImageData,
//                                              u + 1.0*abs_nepy*rescaleFactor,
//                                              v + 1.0*abs_nepx*rescaleFactor,
//                                              width);
//    float realVal_ppyx2 = getInterpolatedElement(activeKeyFrameImageData,
//                                              u + 2.0*abs_nepy*rescaleFactor,
//                                              v + 2.0*abs_nepx*rescaleFactor,
//                                              width);
#else
#if 0 // CROSS with normal and rescale
    
    // PIXELS ON A 3x3 WINDOW
    //         (0, -1)
    // (-1,0)  (0,  0)   (1,0)
    //         (0,  1)
    float realVal_m2 = getInterpolatedElement(activeKeyFrameImageData,
                                              u - 0*normalEpx*rescaleFactor,
                                              v - 1*normalEpy*rescaleFactor,
                                              width);
    float realVal_m1 = getInterpolatedElement(activeKeyFrameImageData,
                                              u - 1*normalEpx*rescaleFactor,
                                              v - 0*normalEpy*rescaleFactor,
                                              width);
    float realVal = getInterpolatedElement(activeKeyFrameImageData,
                                           u,
                                           v,
                                           width);
    float realVal_p1 = getInterpolatedElement(activeKeyFrameImageData,
                                              u + 1*normalEpx*rescaleFactor,
                                              v + 0*normalEpy*rescaleFactor,
                                              width);
    float realVal_p2 = getInterpolatedElement(activeKeyFrameImageData,
                                              u + 0*normalEpx*rescaleFactor,
                                              v + 1*normalEpy*rescaleFactor,
                                              width);
#else // CROSS without normal and rescale
    
    // PIXELS ON A 3x3 WINDOW
    //         (0, -1)
    // (-1,0)  (0,  0)   (1,0)
    //         (0,  1)
    float incx = ptEC->incx;
    float incy = ptEC->incy;
    float abs_incx = (float)fabs(incx);
    float abs_incy = (float)fabs(incy);
    float len_incxincy = sqrt(incx*incx + incy*incy);
    if (u - 3.0f*len_incxincy <= 0 ||
        u - 3.0f*len_incxincy >= width - 2 ||
        v - 3.0f*len_incxincy <= 0 ||
        v - 3.0f*len_incxincy >= height - 2 ||
        u + 3.0f*len_incxincy <= 0 ||
        u + 3.0f*len_incxincy >= width - 2 ||
        v + 3.0f*len_incxincy <= 0 ||
        v + 3.0f*len_incxincy >= height - 2)
    {
        
        Eigen::Matrix<float,5,1> sourceVals;
        sourceVals << -1, -1, -1, -1, -1;
        return sourceVals;
        
    }
    float realVal_m2 = getInterpolatedElement(activeKeyFrameImageData,
                                              u - 0.0*len_incxincy,
                                              v - 3.0*len_incxincy,
                                              width);
    float realVal_m1 = getInterpolatedElement(activeKeyFrameImageData,
                                              u - 3.0*len_incxincy,
                                              v - 0.0*len_incxincy,
                                              width);
    float realVal = getInterpolatedElement(activeKeyFrameImageData,
                                           u,
                                           v,
                                           width);
    float realVal_p1 = getInterpolatedElement(activeKeyFrameImageData,
                                              u + 3.0*len_incxincy,
                                              v + 0.0*len_incxincy,
                                              width);
    float realVal_p2 = getInterpolatedElement(activeKeyFrameImageData,
                                              u + 0.0*len_incxincy,
                                              v + 3.0*len_incxincy,
                                              width);
#endif
#endif
    
#if 0
    
    // Nine values from the source
    Eigen::Matrix<float,9,1> sourceVals;
    sourceVals << realVal_mmxy2, realVal_mmxy1, // Minus Minus X Y
    realVal,                      // Origin
    realVal_ppxy1, realVal_ppxy2, // Plus Plus X Y
    realVal_mpyx1, realVal_mpyx2, // Minus Plus Y X
    realVal_ppyx1, realVal_ppyx2; // Plus Plus Y X
    
    
#else // (ORIGINAL) FIVE VALUES
    
    // Five values from the source
    Eigen::Matrix<float,5,1> sourceVals;
    sourceVals << realVal_mmxy2, realVal_mmxy1, // Minus Minus X Y
                  realVal,                      // Origin
                  realVal_ppxy1, realVal_ppxy2;
    
#endif
    
    return sourceVals;
    
}

//------------------------------------------------------------------------------
// Generalised epipolar line stereo, returns a best matching error score.
//------------------------------------------------------------------------------
//
//             ptsEpiCurve : A vector of PointEpiCurve structures
//                           that is sorted by the inverse depth
//                           in the ascending order
//              min_idepth : Minimum inverse depth
//              max_idepth : Maximum inverse depth
//          referenceFrame : refFrame (target)
//     referenceFrameImage : tracking image
//           result_idepth : Result inverse depth
//              result_var : Result variance
//        result_eplLength : Result epipolar length
//                   stats : Return status
//
// If any error happens, following error code returns
//   -1 : out of bounds
//   -2 : not found
//
inline float
DepthMapForRollingShutter::
doGeneralEpipolarLineStereo(
    const std::list<PointEpiCurve> *ptsEpiCurve,
    const float min_idepth,
    const float prior_idepth,
    float max_idepth,
    lsd_slam::Frame *const referenceFrame,
    const float *referenceFrameImage,
    float &result_idepth,
    float &result_var,
    float &result_eplLength,
    lsd_slam::StereoDebugInfo &stereoDebugInfo,
    PointEpiCurve &best_ptEC,
    lsd_slam::RunningStats *const stats)
{
    
    //--------------------------------------------------------------------------
    // Variables
    float best_match_x = -1;
    float best_match_y = -1;
    float best_match_err = FLT_MAX;
    float second_best_match_err = FLT_MAX;
    int loopCBest = -1, loopCSecond = -1;
    float e1A = NAN, e1B = NAN, e2A = NAN, e2B = NAN, e3A = NAN;
    float e3B = NAN, e4A = NAN, e4B = NAN, e5A = NAN, e5B = NAN;
//    float e6A = NAN, e6B = NAN, e7A = NAN, e7B = NAN, e8A = NAN;
//    float e8B = NAN, e9A = NAN, e9B = NAN;
    float eeLast = -1;
    float best_match_errPre = NAN, best_match_errPost = NAN;
    float best_match_DiffErrPre = NAN, best_match_DiffErrPost = NAN;
    bool bestWasLastLoop = false;
    Sophus::SE3d best_motion; // Motion found by the stereo match
    //const PointEpiCurve *best_ptEC; // Best point on EpiCurve found in matching
    
    // For each element in the image band from the buffer
    std::list<PointEpiCurve>::const_iterator ci = ptsEpiCurve->begin();
    Eigen::Matrix<float,5,1> sourceVals =
        getSourceValsFromEpiPtCurve(&(*ci));
    if (sourceVals[0] == -1) // No source vals available
    {
        return -1;
    }

    for (int loopCounter=0;
        	ci != ptsEpiCurve->end(); ++ci, loopCounter++)
    {
        
        // Point structure on the epioplar curve
        const PointEpiCurve *ptEC = &(*ci);
        
        // Five values from the buffer (target)
        Eigen::Matrix<float,5,1> targetVals = ptEC->vals;
        
        if (ptEC->isValsAvailable == false)
        {
            
            continue;
            
        }

        //----------------------------------------------------------------------
        // Search a match for disparity
        
        // Compute error
        // Calcuate errors and accumulate sums
        float ee = 0;
#if 1 // USE ALL FIVE EQUIDISTANT PIXELS
        if (loopCounter%2 == 0)
        {
            // calc error and accumulate sums.
            e1A = targetVals[0] - sourceVals[0]; ee += e1A*e1A;
            e2A = targetVals[1] - sourceVals[1]; ee += e2A*e2A;
            e3A = targetVals[2] - sourceVals[2]; ee += e3A*e3A;
            e4A = targetVals[3] - sourceVals[3]; ee += e4A*e4A;
            e5A = targetVals[4] - sourceVals[4]; ee += e5A*e5A;
        
//            e6A = targetVals[5] - sourceVals[5]; ee += e6A*e6A;
//            e7A = targetVals[6] - sourceVals[6]; ee += e7A*e7A;
//            e8A = targetVals[7] - sourceVals[7]; ee += e8A*e8A;
//            e9A = targetVals[8] - sourceVals[8]; ee += e9A*e9A;
        }
        else
        {
            // calc error and accumulate sums.
            e1B = targetVals[0] - sourceVals[0]; ee += e1B*e1B;
            e2B = targetVals[1] - sourceVals[1]; ee += e2B*e2B;
            e3B = targetVals[2] - sourceVals[2]; ee += e3B*e3B;
            e4B = targetVals[3] - sourceVals[3]; ee += e4B*e4B;
            e5B = targetVals[4] - sourceVals[4]; ee += e5B*e5B;

//            e6B = targetVals[5] - sourceVals[5]; ee += e6B*e6B;
//            e7B = targetVals[6] - sourceVals[6]; ee += e7B*e7B;
//            e8B = targetVals[7] - sourceVals[7]; ee += e8B*e8B;
//            e9B = targetVals[8] - sourceVals[8]; ee += e9B*e9B;
        }
#else

        e3A = targetVals[2] - sourceVals[2]; ee += e3A*e3A;
        
#endif
        
        // Do I have a new winner?
        // If so, set
        if (ee < best_match_err)
        {
            
            // Put to second-best
            second_best_match_err = best_match_err;
            loopCSecond = loopCBest;
            
            // Set best
            best_match_err = ee;
            loopCBest = loopCounter;
            
            best_match_errPre = eeLast;
            best_match_DiffErrPre = e1A*e1B + e2A*e2B + e3A*e3B +
                                    e4A*e4B + e5A*e5B;
//            best_match_DiffErrPre += e6A*e6B + e7A*e7B +
//                                    e8A*e8B + e9A*e9B;
            best_match_errPost = -1;
            best_match_DiffErrPost = -1;
            
            best_match_x = ptEC->xTarget;
            best_match_y = ptEC->yTarget;
            bestWasLastLoop = true;
            
            best_motion = ptEC->Xi;
            best_ptEC = *ptEC;
            
        }
        // Otherwise: the last might be the current winner,
        // in which case i have to save these values.
        else
        {
            
            if (bestWasLastLoop)
            {
                
                best_match_errPost = ee;
                best_match_DiffErrPost = e1A*e1B + e2A*e2B + e3A*e3B +
                                            e4A*e4B + e5A*e5B;
//                best_match_DiffErrPost += e6A*e6B + e7A*e7B +
//                                            e8A*e8B + e9A*e9B;
                bestWasLastLoop = false;
                
            }
            
            // Collect second-best:
            // just take the best of all that are NOT equal to current best.
            if (ee < second_best_match_err)
            {
                
                second_best_match_err = ee;
                loopCSecond = loopCounter;
                
            }
            
        }
        
        // shift everything one further.
        eeLast = ee;
        
        if (enablePrintDebugInfo) stats->num_stereo_comparisons++;
        
    }
    
    // if error too big, will return -3, otherwise -2.
    if(best_match_err > 4.0f*(float)MAX_ERROR_STEREO)
    {
        if (enablePrintDebugInfo) stats->num_stereo_invalid_bigErr++;
        return -3;
    }
    
    // check if clear enough winner
    if (abs(loopCBest - loopCSecond) > 1.0f &&
        MIN_DISTANCE_ERROR_STEREO * best_match_err > second_best_match_err)
    {
        if(enablePrintDebugInfo) stats->num_stereo_invalid_unclear_winner++;
        return -2;
    }
    
    bool didSubpixel = false;
#if 0// DEACTIVATE SUBPIXEL
    useSubpixelStereo = false;
#endif
    if(useSubpixelStereo)
    {
        // ================== compute exact match =========================
        // compute gradients (they are actually only half the real gradient)
        float gradPre_pre = -(best_match_errPre - best_match_DiffErrPre);
        float gradPre_this = +(best_match_err - best_match_DiffErrPre);
        float gradPost_this = -(best_match_err - best_match_DiffErrPost);
        float gradPost_post = +(best_match_errPost - best_match_DiffErrPost);
        
        // final decisions here.
        bool interpPost = false;
        bool interpPre = false;
        
        // if one is oob: return false.
        if(enablePrintDebugInfo && (best_match_errPre < 0 || best_match_errPost < 0))
        {
            stats->num_stereo_invalid_atEnd++;
        }
        
        
        // - if zero-crossing occurs exactly in between (gradient Inconsistent),
        else if((gradPost_this < 0) ^ (gradPre_this < 0))
        {
            // return exact pos, if both central gradients are small compared to their counterpart.
            if(enablePrintDebugInfo && (gradPost_this*gradPost_this > 0.1f*0.1f*gradPost_post*gradPost_post ||
                                        gradPre_this*gradPre_this > 0.1f*0.1f*gradPre_pre*gradPre_pre))
                stats->num_stereo_invalid_inexistantCrossing++;
        }
        
        // if pre has zero-crossing
        else if((gradPre_pre < 0) ^ (gradPre_this < 0))
        {
            // if post has zero-crossing
            if((gradPost_post < 0) ^ (gradPost_this < 0))
            {
                if(enablePrintDebugInfo) stats->num_stereo_invalid_twoCrossing++;
            }
            else
                interpPre = true;
        }
        
        // if post has zero-crossing
        else if((gradPost_post < 0) ^ (gradPost_this < 0))
        {
            interpPost = true;
        }
        
        // if none has zero-crossing
        else
        {
            if(enablePrintDebugInfo) stats->num_stereo_invalid_noCrossing++;
        }
        
        
        // DO interpolation!
        // minimum occurs at zero-crossing of gradient, which is a straight line => easy to compute.
        // the error at that point is also computed by just integrating.
        if(interpPre)
        {
            float d = gradPre_this / (gradPre_this - gradPre_pre);
            float incx = best_ptEC.incx;
            float incy = best_ptEC.incy;
            best_match_x -= d*incx;
            best_match_y -= d*incy;
            best_match_err = best_match_err - 2*d*gradPre_this - (gradPre_pre - gradPre_this)*d*d;
            if(enablePrintDebugInfo) stats->num_stereo_interpPre++;
            didSubpixel = true;
            
        }
        else if(interpPost)
        {
            float d = gradPost_this / (gradPost_this - gradPost_post);
            float incx = best_ptEC.incx;
            float incy = best_ptEC.incy;
            best_match_x += d*incx;
            best_match_y += d*incy;
            best_match_err = best_match_err + 2*d*gradPost_this + (gradPost_post - gradPost_this)*d*d;
            if(enablePrintDebugInfo) stats->num_stereo_interpPost++;
            didSubpixel = true;
        }
        else
        {
            if(enablePrintDebugInfo) stats->num_stereo_interpNone++;
        }
    }
    
    
    //----------------------------------------------------------------------
    // Convert disparity to inverse depth
    // Propagate variance of the inverse depth
    
    // sampleDist is the distance in pixel at which the realVal's were sampled
    float rescaleFactor = best_ptEC.rescaleFactor;
    float sampleDist = GRADIENT_SAMPLE_DIST*rescaleFactor;
    float realVal_m2 = sourceVals[0];
    float realVal_m1 = sourceVals[1];
    float realVal = sourceVals[2];
    float realVal_p1 = sourceVals[3];
    float realVal_p2 = sourceVals[4];
    
    float gradAlongLine = 0;
    float tmp = realVal_p2 - realVal_p1;  gradAlongLine+=tmp*tmp;
    tmp = realVal_p1 - realVal;  gradAlongLine+=tmp*tmp;
    tmp = realVal - realVal_m1;  gradAlongLine+=tmp*tmp;
    tmp = realVal_m1 - realVal_m2;  gradAlongLine+=tmp*tmp;
    
    gradAlongLine /= sampleDist*sampleDist;
    
    // check if interpolated error is OK. use evil hack to allow more error if there is a lot of gradient.
    if(best_match_err > (float)MAX_ERROR_STEREO + sqrtf( gradAlongLine)*20)
    {
        if(enablePrintDebugInfo) stats->num_stereo_invalid_bigErr++;
        return -3;
    }
    
    // ================= calc depth (in KF) ====================
    // * KinvP = Kinv * (x,y,1); where x,y are pixel coordinates of point we search for, in the KF.
    // * best_match_x = x-coordinate of found correspondence in the reference frame.
    float u = best_ptEC.xSource;
    float v = best_ptEC.ySource;
    
    float idnew_best_match;	// depth in the new image
    float alpha; // d(idnew_best_match) / d(disparity in pixel) == conputed inverse depth derived by the pixel-disparity.

#if 0 // COMPUTE DEPTH
    // Inverse depth compuatation in x-direction
    //if(incx*incx>incy*incy)
    {
        Sophus::SE3d otherToThis = best_motion;
        Eigen::Matrix<float,3,1> t = otherToThis.translation().cast<float>();
        Eigen::Matrix<float,3,3> R = otherToThis.rotationMatrix().cast<float>();
        float oldX = fxi*best_match_x+cxi;
        float nominator = (oldX*t[2] - t[0]);
        float dot0 = KinvP.dot(R.row(0));
        float dot2 = KinvP.dot(R.row(2));
        
        idnew_best_match = (dot0 - oldX*dot2) / nominator;
        alpha = 1.0*fxi*(dot0*t[2] - dot2*t[0]) / (nominator*nominator);
        
    }
#else // GET DEPTH FROM ptEC structure
    
    // Compute inverse depth and alpha
    if (didSubpixel == true)
    {
        // Update from subpixel match
        best_ptEC.xTarget = best_match_x;
        best_ptEC.yTarget = best_match_y;
        
    }
    best_ptEC.computeInvDepth(fxi, fyi, cxi, cyi);
    idnew_best_match = best_ptEC.invDepth;
    alpha = best_ptEC.alpha;
    
#endif
    
//    else
//    {
//        float oldY = fyi*best_match_y+cyi;
//        
//        float nominator = (oldY*referenceFrame->otherToThis_t[2] - referenceFrame->otherToThis_t[1]);
//        float dot1 = KinvP.dot(referenceFrame->otherToThis_R_row1);
//        float dot2 = KinvP.dot(referenceFrame->otherToThis_R_row2);
//        
//        idnew_best_match = (dot1 - oldY*dot2) / nominator;
//        alpha = incy*fyi*(dot1*referenceFrame->otherToThis_t[2] - dot2*referenceFrame->otherToThis_t[1]) / (nominator*nominator);
//        
//    }

    if (idnew_best_match < 0)
    {
        if(enablePrintDebugInfo) stats->num_stereo_negative++;
        if(!allowNegativeIdepths)
            return -2;
    }
    
    if (enablePrintDebugInfo) stats->num_stereo_successfull++;
    
    // ================= calc var (in NEW image) ====================
    
    // calculate error from photometric noise
    float photoDispError = 4.0f * cameraPixelNoise2 /
                            (gradAlongLine + DIVISION_EPS);
    
    float trackingErrorFac =
            0.25f*(1.0f+referenceFrame->initialTrackedResidual);
    
    // calculate error from geometric noise (wrong camera pose / calibration)
    Eigen::Vector2f gradsInterp =
        getInterpolatedElement42(activeKeyFrame->gradients(0), u, v, width);
    float epxn = best_ptEC.nepx;
    float epyn = best_ptEC.nepy;
    float geoDispError_denom = (gradsInterp[0]*epxn + gradsInterp[1]*epyn) +
                            DIVISION_EPS;
    this->debugInfo_geoDispError_first = geoDispError_denom;

     
    float geoDispError =
        trackingErrorFac*trackingErrorFac*(gradsInterp[0]*gradsInterp[0] +
                                           gradsInterp[1]*gradsInterp[1]) /
        (geoDispError_denom*geoDispError_denom);
    
#if 1 // ORIGINAL result_var
    
    result_var = alpha*alpha*((didSubpixel ? 0.05f : 0.5f)*
                              sampleDist*sampleDist +  geoDispError +
                              photoDispError);	// square to make variance
    
#else // ENFORCE constant variance for testing
    
    result_var = 0.001;
    
#endif

    stereoDebugInfo.alpha = alpha;
    stereoDebugInfo.sampleDist = sampleDist;
    stereoDebugInfo.geoDispError = geoDispError;
    stereoDebugInfo.photoDispError = photoDispError;
    stereoDebugInfo.rescaleFactor = rescaleFactor;
    stereoDebugInfo.gradAlongLine = gradAlongLine;
    stereoDebugInfo.trackingErrorFac = trackingErrorFac;
    stereoDebugInfo.geoDispError_denom = geoDispError_denom;
    stereoDebugInfo.gradsInterpX = gradsInterp[0];
    stereoDebugInfo.gradsInterpY = gradsInterp[1];
    stereoDebugInfo.initialTrackedResidual =
        referenceFrame->initialTrackedResidual;
    
    // Debug info
    this->debugInfo_rescaleFactor = rescaleFactor;
    this->debugInfo_epxn = epxn;
    this->debugInfo_epyn = epyn;
    this->debugInfo_alpha = alpha;
    this->debugInfo_sampleDist = sampleDist;
    this->debugInfo_geoDispError = geoDispError;
    this->debugInfo_photoDispError = photoDispError;
    this->debugInfo_result_var = result_var;
    this->debugInfo_trackingErrorFac = trackingErrorFac;
    this->debugInfo_gradAlongLine = gradAlongLine;
    this->debugInfo_gradsInterp_0 = gradsInterp[0];
    this->debugInfo_gradsInterp_1 = gradsInterp[1];
    this->debugInfo_referenceFrame_initialTrackedResidual =
        referenceFrame->initialTrackedResidual;
    this->debugInfo_xSource = u;
    this->debugInfo_ySource = v;
    this->debugInfo_xTarget = best_ptEC.xTarget;
    this->debugInfo_yTarget = best_ptEC.yTarget;
    this->debugInfo_best_match_err = best_match_err;
    
    // Obtain pFar and pClose assumming that the vector of ptsEpiCurve
    // structures is sorted
    Eigen::Vector3f pFar((*ptsEpiCurve).begin()->xTarget,
                           (*ptsEpiCurve).begin()->yTarget,
                           1.0f);
    Eigen::Vector3f pClose((*ptsEpiCurve).end()->xTarget,
                           (*ptsEpiCurve).end()->yTarget,
                           1.0f);
    
    if (plotStereoImages)
    {
        
        if (rand()%5==0)
        {

            float fac = best_match_err /
                ((float)MAX_ERROR_STEREO + sqrtf( gradAlongLine)*20);
            
            cv::Scalar color = cv::Scalar(255*fac, 255-255*fac, 0);// bw
            
            cv::line(debugImageStereoLines,
                     cv::Point2f(pClose[0], pClose[1]),
                     cv::Point2f(pFar[0], pFar[1]),
                     color, 1, 8, 0);
            
        }
        
    }

    // Output the inverse depth
    result_idepth = idnew_best_match;
    
#if 1
    // MORE CARE NEEDED for eplLength and this is approxmiate.
    // ...
//    float incx = pClose[0] - pFar[0];
//    float incy = pClose[1] - pFar[1];
//    float eplLength = sqrt(incx*incx+incy*incy);
    result_eplLength = best_ptEC.eplLength;
    
#endif
    
    return best_match_err;
    
}
    
// find pixel in image (do stereo along epipolar line).
// mat: NEW image
// KinvP: point in OLD image (Kinv * (u_old, v_old, 1)), projected
// trafo: x_old = trafo * x_new; (from new to old image)
// realVal: descriptor in OLD image.
// returns: result_idepth : point depth in new camera's coordinate system
// returns: result_u/v : point's coordinates in new camera's coordinate system
// returns: idepth_var: (approximated) measurement variance of inverse depth of result_point_NEW
// returns error if sucessful; -1 if out of bounds, -2 if not found.
inline float DepthMapForRollingShutter::doLineStereo(
    const float u, const float v, const float epxn, const float epyn,
    const float min_idepth, const float prior_idepth, float max_idepth,
    Frame* const referenceFrame, const float* referenceFrameImage,
    float &result_idepth, float &result_var, float &result_eplLength,
    RunningStats* stats)
// u : x-coord of the semi-dense point in the source/key image
// v : y-coord of the semi-dense point in the source/key image
// referenceFrame : target/ref/currentTrackingFrame image
{

    if(enablePrintDebugInfo) stats->num_stereo_calls++;
    
    // calculate epipolar line start and end point in old image
    Eigen::Vector3f KinvP = Eigen::Vector3f(fxi*u+cxi,fyi*v+cyi,1.0f);
    Eigen::Vector3f pInf = referenceFrame->K_otherToThis_R * KinvP;
    Eigen::Vector3f pReal = pInf / prior_idepth + referenceFrame->K_otherToThis_t;
    
    float rescaleFactor = pReal[2] * prior_idepth;
    
    float firstX = u - 2*epxn*rescaleFactor;
    float firstY = v - 2*epyn*rescaleFactor;
    float lastX = u + 2*epxn*rescaleFactor;
    float lastY = v + 2*epyn*rescaleFactor;
    // width - 2 and height - 2 comes from the one-sided gradient calculation at the bottom
    if (firstX <= 0 || firstX >= width - 2
        || firstY <= 0 || firstY >= height - 2
        || lastX <= 0 || lastX >= width - 2
        || lastY <= 0 || lastY >= height - 2) {
        return -1;
    }
    
    if(!(rescaleFactor > 0.7f && rescaleFactor < 1.4f))
    {
        if(enablePrintDebugInfo) stats->num_stereo_rescale_oob++;
        return -1;
    }
    
    // calculate values to search for
    float realVal_p1 = getInterpolatedElement(activeKeyFrameImageData,u + epxn*rescaleFactor, v + epyn*rescaleFactor, width);
    float realVal_m1 = getInterpolatedElement(activeKeyFrameImageData,u - epxn*rescaleFactor, v - epyn*rescaleFactor, width);
    float realVal = getInterpolatedElement(activeKeyFrameImageData,u, v, width);
    float realVal_m2 = getInterpolatedElement(activeKeyFrameImageData,u - 2*epxn*rescaleFactor, v - 2*epyn*rescaleFactor, width);
    float realVal_p2 = getInterpolatedElement(activeKeyFrameImageData,u + 2*epxn*rescaleFactor, v + 2*epyn*rescaleFactor, width);
    
    
    
    //	if(referenceFrame->K_otherToThis_t[2] * max_idepth + pInf[2] < 0.01)
    
    
    Eigen::Vector3f pClose = pInf + referenceFrame->K_otherToThis_t*max_idepth;
    // if the assumed close-point lies behind the
    // image, have to change that.
    if(pClose[2] < 0.001f)
    {
        max_idepth = (0.001f-pInf[2]) / referenceFrame->K_otherToThis_t[2];
        pClose = pInf + referenceFrame->K_otherToThis_t*max_idepth;
    }
    pClose = pClose / pClose[2]; // pos in new image of point (xy), assuming max_idepth
    
    Eigen::Vector3f pFar = pInf + referenceFrame->K_otherToThis_t*min_idepth;
    // if the assumed far-point lies behind the image or closter than the near-point,
    // we moved past the Point it and should stop.
    if(pFar[2] < 0.001f || max_idepth < min_idepth)
    {
        if(enablePrintDebugInfo) stats->num_stereo_inf_oob++;
        return -1;
    }
    pFar = pFar / pFar[2]; // pos in new image of point (xy), assuming min_idepth
    
    
    // check for nan due to eg division by zero.
    if(isnanf((float)(pFar[0]+pClose[0])))
        return -4;
    
    // calculate increments in which we will step through the epipolar line.
    // they are sampleDist (or half sample dist) long
    float incx = pClose[0] - pFar[0];
    float incy = pClose[1] - pFar[1];
    float eplLength = sqrt(incx*incx+incy*incy);
    if(!(eplLength > 0) || std::isinf(eplLength)) return -4;
    
    if(eplLength > MAX_EPL_LENGTH_CROP)
    {
        pClose[0] = pFar[0] + incx*MAX_EPL_LENGTH_CROP/eplLength;
        pClose[1] = pFar[1] + incy*MAX_EPL_LENGTH_CROP/eplLength;
    }
    
    incx *= GRADIENT_SAMPLE_DIST/eplLength;
    incy *= GRADIENT_SAMPLE_DIST/eplLength;
    
    
    // extend one sample_dist to left & right.
    pFar[0] -= incx;
    pFar[1] -= incy;
    pClose[0] += incx;
    pClose[1] += incy;
    
    
    // make epl long enough (pad a little bit).
    if(eplLength < MIN_EPL_LENGTH_CROP)
    {
        float pad = (MIN_EPL_LENGTH_CROP - (eplLength)) / 2.0f;
        pFar[0] -= incx*pad;
        pFar[1] -= incy*pad;
        
        pClose[0] += incx*pad;
        pClose[1] += incy*pad;
    }
    
    // if inf point is outside of image: skip pixel.
    if(
       pFar[0] <= SAMPLE_POINT_TO_BORDER ||
       pFar[0] >= width-SAMPLE_POINT_TO_BORDER ||
       pFar[1] <= SAMPLE_POINT_TO_BORDER ||
       pFar[1] >= height-SAMPLE_POINT_TO_BORDER)
    {
        if(enablePrintDebugInfo) stats->num_stereo_inf_oob++;
        return -1;
    }
    
    
    
    // if near point is outside: move inside, and test length again.
    if(
       pClose[0] <= SAMPLE_POINT_TO_BORDER ||
       pClose[0] >= width-SAMPLE_POINT_TO_BORDER ||
       pClose[1] <= SAMPLE_POINT_TO_BORDER ||
       pClose[1] >= height-SAMPLE_POINT_TO_BORDER)
    {
        if(pClose[0] <= SAMPLE_POINT_TO_BORDER)
        {
            float toAdd = (SAMPLE_POINT_TO_BORDER - pClose[0]) / incx;
            pClose[0] += toAdd * incx;
            pClose[1] += toAdd * incy;
        }
        else if(pClose[0] >= width-SAMPLE_POINT_TO_BORDER)
        {
            float toAdd = (width-SAMPLE_POINT_TO_BORDER - pClose[0]) / incx;
            pClose[0] += toAdd * incx;
            pClose[1] += toAdd * incy;
        }
        
        if(pClose[1] <= SAMPLE_POINT_TO_BORDER)
        {
            float toAdd = (SAMPLE_POINT_TO_BORDER - pClose[1]) / incy;
            pClose[0] += toAdd * incx;
            pClose[1] += toAdd * incy;
        }
        else if(pClose[1] >= height-SAMPLE_POINT_TO_BORDER)
        {
            float toAdd = (height-SAMPLE_POINT_TO_BORDER - pClose[1]) / incy;
            pClose[0] += toAdd * incx;
            pClose[1] += toAdd * incy;
        }
        
        // get new epl length
        float fincx = pClose[0] - pFar[0];
        float fincy = pClose[1] - pFar[1];
        float newEplLength = sqrt(fincx*fincx+fincy*fincy);
        
        // test again
        if(
           pClose[0] <= SAMPLE_POINT_TO_BORDER ||
           pClose[0] >= width-SAMPLE_POINT_TO_BORDER ||
           pClose[1] <= SAMPLE_POINT_TO_BORDER ||
           pClose[1] >= height-SAMPLE_POINT_TO_BORDER ||
           newEplLength < 8.0f
           )
        {
            if(enablePrintDebugInfo) stats->num_stereo_near_oob++;
            return -1;
        }
        
        
    }
    
    
    // from here on:
    // - pInf: search start-point
    // - p0: search end-point
    // - incx, incy: search steps in pixel
    // - eplLength, min_idepth, max_idepth: determines search-resolution, i.e. the result's variance.
    
    
    float cpx = pFar[0];
    float cpy =  pFar[1];
    
    float val_cp_m2 = getInterpolatedElement(referenceFrameImage,cpx-2.0f*incx, cpy-2.0f*incy, width);
    float val_cp_m1 = getInterpolatedElement(referenceFrameImage,cpx-incx, cpy-incy, width);
    float val_cp = getInterpolatedElement(referenceFrameImage,cpx, cpy, width);
    float val_cp_p1 = getInterpolatedElement(referenceFrameImage,cpx+incx, cpy+incy, width);
    float val_cp_p2;
    
    
    
    /*
     * Subsequent exact minimum is found the following way:
     * - assuming lin. interpolation, the gradient of Error at p1 (towards p2) is given by
     *   dE1 = -2sum(e1*e1 - e1*e2)
     *   where e1 and e2 are summed over, and are the residuals (not squared).
     *
     * - the gradient at p2 (coming from p1) is given by
     * 	 dE2 = +2sum(e2*e2 - e1*e2)
     *
     * - linear interpolation => gradient changes linearely; zero-crossing is hence given by
     *   p1 + d*(p2-p1) with d = -dE1 / (-dE1 + dE2).
     *
     *
     *
     * => I for later exact min calculation, I need sum(e_i*e_i),sum(e_{i-1}*e_{i-1}),sum(e_{i+1}*e_{i+1})
     *    and sum(e_i * e_{i-1}) and sum(e_i * e_{i+1}),
     *    where i is the respective winning index.
     */
    
    
    // walk in equally sized steps, starting at depth=infinity.
    int loopCounter = 0;
    float best_match_x = -1;
    float best_match_y = -1;
    float best_match_err = 1e50;
    float second_best_match_err = 1e50;
    
    // best pre and post errors.
    float best_match_errPre=NAN, best_match_errPost=NAN, best_match_DiffErrPre=NAN, best_match_DiffErrPost=NAN;
    bool bestWasLastLoop = false;
    
    float eeLast = -1; // final error of last comp.
    
    // alternating intermediate vars
    float e1A=NAN, e1B=NAN, e2A=NAN, e2B=NAN, e3A=NAN, e3B=NAN, e4A=NAN, e4B=NAN, e5A=NAN, e5B=NAN;
    
    int loopCBest=-1, loopCSecond =-1;
    while(((incx < 0) == (cpx > pClose[0]) && (incy < 0) == (cpy > pClose[1])) || loopCounter == 0)
    {
        // interpolate one new point
        val_cp_p2 = getInterpolatedElement(referenceFrameImage,cpx+2*incx, cpy+2*incy, width);
        
        
        // hacky but fast way to get error and differential error: switch buffer variables for last loop.
        float ee = 0;
        if(loopCounter%2==0)
        {
            // calc error and accumulate sums.
            e1A = val_cp_p2 - realVal_p2;ee += e1A*e1A;
            e2A = val_cp_p1 - realVal_p1;ee += e2A*e2A;
            e3A = val_cp - realVal;      ee += e3A*e3A;
            e4A = val_cp_m1 - realVal_m1;ee += e4A*e4A;
            e5A = val_cp_m2 - realVal_m2;ee += e5A*e5A;
        }
        else
        {
            // calc error and accumulate sums.
            e1B = val_cp_p2 - realVal_p2;ee += e1B*e1B;
            e2B = val_cp_p1 - realVal_p1;ee += e2B*e2B;
            e3B = val_cp - realVal;      ee += e3B*e3B;
            e4B = val_cp_m1 - realVal_m1;ee += e4B*e4B;
            e5B = val_cp_m2 - realVal_m2;ee += e5B*e5B;
        }
        
        
        // do I have a new winner??
        // if so: set.
        if(ee < best_match_err)
        {
            // put to second-best
            second_best_match_err=best_match_err;
            loopCSecond = loopCBest;
            
            // set best.
            best_match_err = ee;
            loopCBest = loopCounter;
            
            best_match_errPre = eeLast;
            best_match_DiffErrPre = e1A*e1B + e2A*e2B + e3A*e3B + e4A*e4B + e5A*e5B;
            best_match_errPost = -1;
            best_match_DiffErrPost = -1;
            
            best_match_x = cpx;
            best_match_y = cpy;
            bestWasLastLoop = true;
        }
        // otherwise: the last might be the current winner, in which case i have to save these values.
        else
        {
            if(bestWasLastLoop)
            {
                best_match_errPost = ee;
                best_match_DiffErrPost = e1A*e1B + e2A*e2B + e3A*e3B + e4A*e4B + e5A*e5B;
                bestWasLastLoop = false;
            }
            
            // collect second-best:
            // just take the best of all that are NOT equal to current best.
            if(ee < second_best_match_err)
            {
                second_best_match_err=ee;
                loopCSecond = loopCounter;
            }
        }
        
        
        // shift everything one further.
        eeLast = ee;
        val_cp_m2 = val_cp_m1; val_cp_m1 = val_cp; val_cp = val_cp_p1; val_cp_p1 = val_cp_p2;
        
        if(enablePrintDebugInfo) stats->num_stereo_comparisons++;
        
        cpx += incx;
        cpy += incy;
        
        loopCounter++;
    }
    
    // if error too big, will return -3, otherwise -2.
    if(best_match_err > 4.0f*(float)MAX_ERROR_STEREO)
    {
        if(enablePrintDebugInfo) stats->num_stereo_invalid_bigErr++;
        return -3;
    }
    
    
    // check if clear enough winner
    if(abs(loopCBest - loopCSecond) > 1.0f && MIN_DISTANCE_ERROR_STEREO * best_match_err > second_best_match_err)
    {
        if(enablePrintDebugInfo) stats->num_stereo_invalid_unclear_winner++;
        return -2;
    }
    
    bool didSubpixel = false;
    if(useSubpixelStereo)
    {
        // ================== compute exact match =========================
        // compute gradients (they are actually only half the real gradient)
        float gradPre_pre = -(best_match_errPre - best_match_DiffErrPre);
        float gradPre_this = +(best_match_err - best_match_DiffErrPre);
        float gradPost_this = -(best_match_err - best_match_DiffErrPost);
        float gradPost_post = +(best_match_errPost - best_match_DiffErrPost);
        
        // final decisions here.
        bool interpPost = false;
        bool interpPre = false;
        
        // if one is oob: return false.
        if(enablePrintDebugInfo && (best_match_errPre < 0 || best_match_errPost < 0))
        {
            stats->num_stereo_invalid_atEnd++;
        }
        
        
        // - if zero-crossing occurs exactly in between (gradient Inconsistent),
        else if((gradPost_this < 0) ^ (gradPre_this < 0))
        {
            // return exact pos, if both central gradients are small compared to their counterpart.
            if(enablePrintDebugInfo && (gradPost_this*gradPost_this > 0.1f*0.1f*gradPost_post*gradPost_post ||
                                        gradPre_this*gradPre_this > 0.1f*0.1f*gradPre_pre*gradPre_pre))
                stats->num_stereo_invalid_inexistantCrossing++;
        }
        
        // if pre has zero-crossing
        else if((gradPre_pre < 0) ^ (gradPre_this < 0))
        {
            // if post has zero-crossing
            if((gradPost_post < 0) ^ (gradPost_this < 0))
            {
                if(enablePrintDebugInfo) stats->num_stereo_invalid_twoCrossing++;
            }
            else
                interpPre = true;
        }
        
        // if post has zero-crossing
        else if((gradPost_post < 0) ^ (gradPost_this < 0))
        {
            interpPost = true;
        }
        
        // if none has zero-crossing
        else
        {
            if(enablePrintDebugInfo) stats->num_stereo_invalid_noCrossing++;
        }
        
        
        // DO interpolation!
        // minimum occurs at zero-crossing of gradient, which is a straight line => easy to compute.
        // the error at that point is also computed by just integrating.
        if(interpPre)
        {
            float d = gradPre_this / (gradPre_this - gradPre_pre);
            best_match_x -= d*incx;
            best_match_y -= d*incy;
            best_match_err = best_match_err - 2*d*gradPre_this - (gradPre_pre - gradPre_this)*d*d;
            if(enablePrintDebugInfo) stats->num_stereo_interpPre++;
            didSubpixel = true;
            
        }
        else if(interpPost)
        {
            float d = gradPost_this / (gradPost_this - gradPost_post);
            best_match_x += d*incx;
            best_match_y += d*incy;
            best_match_err = best_match_err + 2*d*gradPost_this + (gradPost_post - gradPost_this)*d*d;
            if(enablePrintDebugInfo) stats->num_stereo_interpPost++;
            didSubpixel = true;
        }
        else
        {
            if(enablePrintDebugInfo) stats->num_stereo_interpNone++;
        }
    }
    
    
    // sampleDist is the distance in pixel at which the realVal's were sampled
    float sampleDist = GRADIENT_SAMPLE_DIST*rescaleFactor;
    
    float gradAlongLine = 0;
    float tmp = realVal_p2 - realVal_p1;  gradAlongLine+=tmp*tmp;
    tmp = realVal_p1 - realVal;  gradAlongLine+=tmp*tmp;
    tmp = realVal - realVal_m1;  gradAlongLine+=tmp*tmp;
    tmp = realVal_m1 - realVal_m2;  gradAlongLine+=tmp*tmp;
    
    gradAlongLine /= sampleDist*sampleDist;
    
    // check if interpolated error is OK. use evil hack to allow more error if there is a lot of gradient.
    if(best_match_err > (float)MAX_ERROR_STEREO + sqrtf( gradAlongLine)*20)
    {
        if(enablePrintDebugInfo) stats->num_stereo_invalid_bigErr++;
        return -3;
    }
    
    
    // ================= calc depth (in KF) ====================
    // * KinvP = Kinv * (x,y,1); where x,y are pixel coordinates of point we search for, in the KF.
    // * best_match_x = x-coordinate of found correspondence in the reference frame.
    
    float idnew_best_match;	// depth in the new image
    float alpha; // d(idnew_best_match) / d(disparity in pixel) == conputed inverse depth derived by the pixel-disparity.
    if(incx*incx>incy*incy)
    {
        float oldX = fxi*best_match_x+cxi;
        float nominator = (oldX*referenceFrame->otherToThis_t[2] - referenceFrame->otherToThis_t[0]);
        float dot0 = KinvP.dot(referenceFrame->otherToThis_R_row0);
        float dot2 = KinvP.dot(referenceFrame->otherToThis_R_row2);
        
        idnew_best_match = (dot0 - oldX*dot2) / nominator;
        alpha = incx*fxi*(dot0*referenceFrame->otherToThis_t[2] - dot2*referenceFrame->otherToThis_t[0]) / (nominator*nominator);
        
    }
    else
    {
        float oldY = fyi*best_match_y+cyi;
        
        float nominator = (oldY*referenceFrame->otherToThis_t[2] - referenceFrame->otherToThis_t[1]);
        float dot1 = KinvP.dot(referenceFrame->otherToThis_R_row1);
        float dot2 = KinvP.dot(referenceFrame->otherToThis_R_row2);
        
        idnew_best_match = (dot1 - oldY*dot2) / nominator;
        alpha = incy*fyi*(dot1*referenceFrame->otherToThis_t[2] - dot2*referenceFrame->otherToThis_t[1]) / (nominator*nominator);
        
    }
    
    if(idnew_best_match < 0)
    {
        if(enablePrintDebugInfo) stats->num_stereo_negative++;
        if(!allowNegativeIdepths)
            return -2;
    }
    
    if(enablePrintDebugInfo) stats->num_stereo_successfull++;
    
    // ================= calc var (in NEW image) ====================
    
    // calculate error from photometric noise
    float photoDispError = 4.0f * cameraPixelNoise2 / (gradAlongLine + DIVISION_EPS);
    
    float trackingErrorFac = 0.25f*(1.0f+referenceFrame->initialTrackedResidual);
    
    // calculate error from geometric noise (wrong camera pose / calibration)
    Eigen::Vector2f gradsInterp = getInterpolatedElement42(activeKeyFrame->gradients(0), u, v, width);
    float geoDispError = (gradsInterp[0]*epxn + gradsInterp[1]*epyn) + DIVISION_EPS;
    geoDispError = trackingErrorFac*trackingErrorFac*(gradsInterp[0]*gradsInterp[0] + gradsInterp[1]*gradsInterp[1]) / (geoDispError*geoDispError);
    
    
    //geoDispError *= (0.5 + 0.5 *result_idepth) * (0.5 + 0.5 *result_idepth);
    
    // final error consists of a small constant part (discretization error),
    // geometric and photometric error.
    result_var = alpha*alpha*((didSubpixel ? 0.05f : 0.5f)*sampleDist*sampleDist +  geoDispError + photoDispError);	// square to make variance
    
    if(plotStereoImages)
    {
        if(rand()%5==0)
        {
            //if(rand()%500 == 0)
            //	printf("geo: %f, photo: %f, alpha: %f\n", sqrt(geoDispError), sqrt(photoDispError), alpha, sqrt(result_var));
            
            
            //int idDiff = (keyFrame->pyramidID - referenceFrame->id);
            //cv::Scalar color = cv::Scalar(0,0, 2*idDiff);// bw
            
            //cv::Scalar color = cv::Scalar(sqrt(result_var)*2000, 255-sqrt(result_var)*2000, 0);// bw
            
            //			float eplLengthF = std::min((float)MIN_EPL_LENGTH_CROP,(float)eplLength);
            //			eplLengthF = std::max((float)MAX_EPL_LENGTH_CROP,(float)eplLengthF);
            //
            //			float pixelDistFound = sqrtf((float)((pReal[0]/pReal[2] - best_match_x)*(pReal[0]/pReal[2] - best_match_x)
            //					+ (pReal[1]/pReal[2] - best_match_y)*(pReal[1]/pReal[2] - best_match_y)));
            //
            float fac = best_match_err / ((float)MAX_ERROR_STEREO + sqrtf( gradAlongLine)*20);
            
            cv::Scalar color = cv::Scalar(255*fac, 255-255*fac, 0);// bw
            
            
            /*
             if(rescaleFactor > 1)
             color = cv::Scalar(500*(rescaleFactor-1),0,0);
             else
             color = cv::Scalar(0,500*(1-rescaleFactor),500*(1-rescaleFactor));
             */
            
            cv::line(debugImageStereoLines,cv::Point2f(pClose[0], pClose[1]),cv::Point2f(pFar[0], pFar[1]),color,1,8,0);
        }
    }
    
    result_idepth = idnew_best_match;
    
    result_eplLength = eplLength;
    
    return best_match_err;
    
}
    
// find pixel in image (do stereo along epipolar line).
// mat: NEW image
// KinvP: point in OLD image (Kinv * (u_old, v_old, 1)), projected
// trafo: x_old = trafo * x_new; (from new to old image)
// realVal: descriptor in OLD image.
// returns: result_idepth : point depth in new camera's coordinate system
// returns: result_u/v : point's coordinates in new camera's coordinate system
// returns: idepth_var: (approximated) measurement variance of inverse depth of result_point_NEW
// returns error if sucessful; -1 if out of bounds, -2 if not found.
inline float DepthMapForRollingShutter::doLineStereo_ver2(
     const float u, const float v, const float epxn, const float epyn,
     const float min_idepth, const float prior_idepth, float max_idepth,
     Frame* const referenceFrame, const float* referenceFrameImage,
     Sophus::Sim3d Xi_worldToKF,
     float &result_idepth, float &result_var, float &result_eplLength,
     RunningStats* stats)
{
    
    //--------------------------------------------------------------------------
    // Prepare stereo for a pixel with the row pose
    //--------------------------------------------------------------------------
    // Note: here, the "referenceFrame" is an unmapped tracking frame.
    //       Therefore, "refToKf" is same as "FrameToRef" in trackFrame().
    Sim3 refToKf = Xi_worldToKF*referenceFrame->getScaledCamToWorld();
    
    referenceFrame->prepareForRollingShutterStereo(activeKeyFrame,
                                                   refToKf, v, K, 0);
    //--------------------------------------------------------------------------
    
    if(enablePrintDebugInfo) stats->num_stereo_calls++;
    
    // calculate epipolar line start and end point in old image
    Eigen::Vector3f KinvP = Eigen::Vector3f(fxi*u+cxi,fyi*v+cyi,1.0f);
    Eigen::Vector3f pInf = referenceFrame->K_otherToThis_R * KinvP;
    Eigen::Vector3f pReal = pInf / prior_idepth + referenceFrame->K_otherToThis_t;
    
    float rescaleFactor = pReal[2] * prior_idepth;
    
    float firstX = u - 2*epxn*rescaleFactor;
    float firstY = v - 2*epyn*rescaleFactor;
    float lastX = u + 2*epxn*rescaleFactor;
    float lastY = v + 2*epyn*rescaleFactor;
    // width - 2 and height - 2 comes from the one-sided gradient calculation at the bottom
    if (firstX <= 0 || firstX >= width - 2
        || firstY <= 0 || firstY >= height - 2
        || lastX <= 0 || lastX >= width - 2
        || lastY <= 0 || lastY >= height - 2) {
        return -1;
    }
    
    if(!(rescaleFactor > 0.7f && rescaleFactor < 1.4f))
    {
        if(enablePrintDebugInfo) stats->num_stereo_rescale_oob++;
        return -1;
    }
    
    // calculate values to search for
    float realVal_p1 = getInterpolatedElement(activeKeyFrameImageData,u + epxn*rescaleFactor, v + epyn*rescaleFactor, width);
    float realVal_m1 = getInterpolatedElement(activeKeyFrameImageData,u - epxn*rescaleFactor, v - epyn*rescaleFactor, width);
    float realVal = getInterpolatedElement(activeKeyFrameImageData,u, v, width);
    float realVal_m2 = getInterpolatedElement(activeKeyFrameImageData,u - 2*epxn*rescaleFactor, v - 2*epyn*rescaleFactor, width);
    float realVal_p2 = getInterpolatedElement(activeKeyFrameImageData,u + 2*epxn*rescaleFactor, v + 2*epyn*rescaleFactor, width);
    
    
    
    //	if(referenceFrame->K_otherToThis_t[2] * max_idepth + pInf[2] < 0.01)
    
    
    Eigen::Vector3f pClose = pInf + referenceFrame->K_otherToThis_t*max_idepth;
    // if the assumed close-point lies behind the
    // image, have to change that.
    if(pClose[2] < 0.001f)
    {
        max_idepth = (0.001f-pInf[2]) / referenceFrame->K_otherToThis_t[2];
        pClose = pInf + referenceFrame->K_otherToThis_t*max_idepth;
    }
    pClose = pClose / pClose[2]; // pos in new image of point (xy), assuming max_idepth
    
    Eigen::Vector3f pFar = pInf + referenceFrame->K_otherToThis_t*min_idepth;
    // if the assumed far-point lies behind the image or closter than the near-point,
    // we moved past the Point it and should stop.
    if(pFar[2] < 0.001f || max_idepth < min_idepth)
    {
        if(enablePrintDebugInfo) stats->num_stereo_inf_oob++;
        return -1;
    }
    pFar = pFar / pFar[2]; // pos in new image of point (xy), assuming min_idepth
    
    
    // check for nan due to eg division by zero.
    if(isnanf((float)(pFar[0]+pClose[0])))
        return -4;
    
    // calculate increments in which we will step through the epipolar line.
    // they are sampleDist (or half sample dist) long
    float incx = pClose[0] - pFar[0];
    float incy = pClose[1] - pFar[1];
    float eplLength = sqrt(incx*incx+incy*incy);
    if(!(eplLength > 0) || std::isinf(eplLength)) return -4;
    
    if(eplLength > MAX_EPL_LENGTH_CROP)
    {
        pClose[0] = pFar[0] + incx*MAX_EPL_LENGTH_CROP/eplLength;
        pClose[1] = pFar[1] + incy*MAX_EPL_LENGTH_CROP/eplLength;
    }
    
    incx *= GRADIENT_SAMPLE_DIST/eplLength;
    incy *= GRADIENT_SAMPLE_DIST/eplLength;
    
    
    // extend one sample_dist to left & right.
    pFar[0] -= incx;
    pFar[1] -= incy;
    pClose[0] += incx;
    pClose[1] += incy;
    
    
    // make epl long enough (pad a little bit).
    if(eplLength < MIN_EPL_LENGTH_CROP)
    {
        float pad = (MIN_EPL_LENGTH_CROP - (eplLength)) / 2.0f;
        pFar[0] -= incx*pad;
        pFar[1] -= incy*pad;
        
        pClose[0] += incx*pad;
        pClose[1] += incy*pad;
    }
    
    // if inf point is outside of image: skip pixel.
    if(
       pFar[0] <= SAMPLE_POINT_TO_BORDER ||
       pFar[0] >= width-SAMPLE_POINT_TO_BORDER ||
       pFar[1] <= SAMPLE_POINT_TO_BORDER ||
       pFar[1] >= height-SAMPLE_POINT_TO_BORDER)
    {
        if(enablePrintDebugInfo) stats->num_stereo_inf_oob++;
        return -1;
    }
    
    
    
    // if near point is outside: move inside, and test length again.
    if(
       pClose[0] <= SAMPLE_POINT_TO_BORDER ||
       pClose[0] >= width-SAMPLE_POINT_TO_BORDER ||
       pClose[1] <= SAMPLE_POINT_TO_BORDER ||
       pClose[1] >= height-SAMPLE_POINT_TO_BORDER)
    {
        if(pClose[0] <= SAMPLE_POINT_TO_BORDER)
        {
            float toAdd = (SAMPLE_POINT_TO_BORDER - pClose[0]) / incx;
            pClose[0] += toAdd * incx;
            pClose[1] += toAdd * incy;
        }
        else if(pClose[0] >= width-SAMPLE_POINT_TO_BORDER)
        {
            float toAdd = (width-SAMPLE_POINT_TO_BORDER - pClose[0]) / incx;
            pClose[0] += toAdd * incx;
            pClose[1] += toAdd * incy;
        }
        
        if(pClose[1] <= SAMPLE_POINT_TO_BORDER)
        {
            float toAdd = (SAMPLE_POINT_TO_BORDER - pClose[1]) / incy;
            pClose[0] += toAdd * incx;
            pClose[1] += toAdd * incy;
        }
        else if(pClose[1] >= height-SAMPLE_POINT_TO_BORDER)
        {
            float toAdd = (height-SAMPLE_POINT_TO_BORDER - pClose[1]) / incy;
            pClose[0] += toAdd * incx;
            pClose[1] += toAdd * incy;
        }
        
        // get new epl length
        float fincx = pClose[0] - pFar[0];
        float fincy = pClose[1] - pFar[1];
        float newEplLength = sqrt(fincx*fincx+fincy*fincy);
        
        // test again
        if(
           pClose[0] <= SAMPLE_POINT_TO_BORDER ||
           pClose[0] >= width-SAMPLE_POINT_TO_BORDER ||
           pClose[1] <= SAMPLE_POINT_TO_BORDER ||
           pClose[1] >= height-SAMPLE_POINT_TO_BORDER ||
           newEplLength < 8.0f
           )
        {
            if(enablePrintDebugInfo) stats->num_stereo_near_oob++;
            return -1;
        }
        
        
    }
    
    
    // from here on:
    // - pInf: search start-point
    // - p0: search end-point
    // - incx, incy: search steps in pixel
    // - eplLength, min_idepth, max_idepth: determines search-resolution, i.e. the result's variance.
    
    
    float cpx = pFar[0];
    float cpy =  pFar[1];
    
    float val_cp_m2 = getInterpolatedElement(referenceFrameImage,cpx-2.0f*incx, cpy-2.0f*incy, width);
    float val_cp_m1 = getInterpolatedElement(referenceFrameImage,cpx-incx, cpy-incy, width);
    float val_cp = getInterpolatedElement(referenceFrameImage,cpx, cpy, width);
    float val_cp_p1 = getInterpolatedElement(referenceFrameImage,cpx+incx, cpy+incy, width);
    float val_cp_p2;
    
    
    
    /*
     * Subsequent exact minimum is found the following way:
     * - assuming lin. interpolation, the gradient of Error at p1 (towards p2) is given by
     *   dE1 = -2sum(e1*e1 - e1*e2)
     *   where e1 and e2 are summed over, and are the residuals (not squared).
     *
     * - the gradient at p2 (coming from p1) is given by
     * 	 dE2 = +2sum(e2*e2 - e1*e2)
     *
     * - linear interpolation => gradient changes linearely; zero-crossing is hence given by
     *   p1 + d*(p2-p1) with d = -dE1 / (-dE1 + dE2).
     *
     *
     *
     * => I for later exact min calculation, I need sum(e_i*e_i),sum(e_{i-1}*e_{i-1}),sum(e_{i+1}*e_{i+1})
     *    and sum(e_i * e_{i-1}) and sum(e_i * e_{i+1}),
     *    where i is the respective winning index.
     */
    
    
    // walk in equally sized steps, starting at depth=infinity.
    int loopCounter = 0;
    float best_match_x = -1;
    float best_match_y = -1;
    float best_match_err = 1e50;
    float second_best_match_err = 1e50;
    
    // best pre and post errors.
    float best_match_errPre=NAN, best_match_errPost=NAN, best_match_DiffErrPre=NAN, best_match_DiffErrPost=NAN;
    bool bestWasLastLoop = false;
    
    float eeLast = -1; // final error of last comp.
    
    // alternating intermediate vars
    float e1A=NAN, e1B=NAN, e2A=NAN, e2B=NAN, e3A=NAN, e3B=NAN, e4A=NAN, e4B=NAN, e5A=NAN, e5B=NAN;
    
    int loopCBest=-1, loopCSecond =-1;
    while(((incx < 0) == (cpx > pClose[0]) && (incy < 0) == (cpy > pClose[1])) || loopCounter == 0)
    {
        // interpolate one new point
        val_cp_p2 = getInterpolatedElement(referenceFrameImage,cpx+2*incx, cpy+2*incy, width);
        
        
        // hacky but fast way to get error and differential error: switch buffer variables for last loop.
        float ee = 0;
        if(loopCounter%2==0)
        {
            // calc error and accumulate sums.
            e1A = val_cp_p2 - realVal_p2;ee += e1A*e1A;
            e2A = val_cp_p1 - realVal_p1;ee += e2A*e2A;
            e3A = val_cp - realVal;      ee += e3A*e3A;
            e4A = val_cp_m1 - realVal_m1;ee += e4A*e4A;
            e5A = val_cp_m2 - realVal_m2;ee += e5A*e5A;
        }
        else
        {
            // calc error and accumulate sums.
            e1B = val_cp_p2 - realVal_p2;ee += e1B*e1B;
            e2B = val_cp_p1 - realVal_p1;ee += e2B*e2B;
            e3B = val_cp - realVal;      ee += e3B*e3B;
            e4B = val_cp_m1 - realVal_m1;ee += e4B*e4B;
            e5B = val_cp_m2 - realVal_m2;ee += e5B*e5B;
        }
        
        
        // do I have a new winner??
        // if so: set.
        if(ee < best_match_err)
        {
            // put to second-best
            second_best_match_err=best_match_err;
            loopCSecond = loopCBest;
            
            // set best.
            best_match_err = ee;
            loopCBest = loopCounter;
            
            best_match_errPre = eeLast;
            best_match_DiffErrPre = e1A*e1B + e2A*e2B + e3A*e3B + e4A*e4B + e5A*e5B;
            best_match_errPost = -1;
            best_match_DiffErrPost = -1;
            
            best_match_x = cpx;
            best_match_y = cpy;
            bestWasLastLoop = true;
        }
        // otherwise: the last might be the current winner, in which case i have to save these values.
        else
        {
            if(bestWasLastLoop)
            {
                best_match_errPost = ee;
                best_match_DiffErrPost = e1A*e1B + e2A*e2B + e3A*e3B + e4A*e4B + e5A*e5B;
                bestWasLastLoop = false;
            }
            
            // collect second-best:
            // just take the best of all that are NOT equal to current best.
            if(ee < second_best_match_err)
            {
                second_best_match_err=ee;
                loopCSecond = loopCounter;
            }
        }
        
        
        // shift everything one further.
        eeLast = ee;
        val_cp_m2 = val_cp_m1; val_cp_m1 = val_cp; val_cp = val_cp_p1; val_cp_p1 = val_cp_p2;
        
        if(enablePrintDebugInfo) stats->num_stereo_comparisons++;
        
        cpx += incx;
        cpy += incy;
        
        loopCounter++;
    }
    
    // if error too big, will return -3, otherwise -2.
    if(best_match_err > 4.0f*(float)MAX_ERROR_STEREO)
    {
        if(enablePrintDebugInfo) stats->num_stereo_invalid_bigErr++;
        return -3;
    }
    
    
    // check if clear enough winner
    if(abs(loopCBest - loopCSecond) > 1.0f && MIN_DISTANCE_ERROR_STEREO * best_match_err > second_best_match_err)
    {
        if(enablePrintDebugInfo) stats->num_stereo_invalid_unclear_winner++;
        return -2;
    }
    
    bool didSubpixel = false;
    if(useSubpixelStereo)
    {
        // ================== compute exact match =========================
        // compute gradients (they are actually only half the real gradient)
        float gradPre_pre = -(best_match_errPre - best_match_DiffErrPre);
        float gradPre_this = +(best_match_err - best_match_DiffErrPre);
        float gradPost_this = -(best_match_err - best_match_DiffErrPost);
        float gradPost_post = +(best_match_errPost - best_match_DiffErrPost);
        
        // final decisions here.
        bool interpPost = false;
        bool interpPre = false;
        
        // if one is oob: return false.
        if(enablePrintDebugInfo && (best_match_errPre < 0 || best_match_errPost < 0))
        {
            stats->num_stereo_invalid_atEnd++;
        }
        
        
        // - if zero-crossing occurs exactly in between (gradient Inconsistent),
        else if((gradPost_this < 0) ^ (gradPre_this < 0))
        {
            // return exact pos, if both central gradients are small compared to their counterpart.
            if(enablePrintDebugInfo && (gradPost_this*gradPost_this > 0.1f*0.1f*gradPost_post*gradPost_post ||
                                        gradPre_this*gradPre_this > 0.1f*0.1f*gradPre_pre*gradPre_pre))
                stats->num_stereo_invalid_inexistantCrossing++;
        }
        
        // if pre has zero-crossing
        else if((gradPre_pre < 0) ^ (gradPre_this < 0))
        {
            // if post has zero-crossing
            if((gradPost_post < 0) ^ (gradPost_this < 0))
            {
                if(enablePrintDebugInfo) stats->num_stereo_invalid_twoCrossing++;
            }
            else
                interpPre = true;
        }
        
        // if post has zero-crossing
        else if((gradPost_post < 0) ^ (gradPost_this < 0))
        {
            interpPost = true;
        }
        
        // if none has zero-crossing
        else
        {
            if(enablePrintDebugInfo) stats->num_stereo_invalid_noCrossing++;
        }
        
        
        // DO interpolation!
        // minimum occurs at zero-crossing of gradient, which is a straight line => easy to compute.
        // the error at that point is also computed by just integrating.
        if(interpPre)
        {
            float d = gradPre_this / (gradPre_this - gradPre_pre);
            best_match_x -= d*incx;
            best_match_y -= d*incy;
            best_match_err = best_match_err - 2*d*gradPre_this - (gradPre_pre - gradPre_this)*d*d;
            if(enablePrintDebugInfo) stats->num_stereo_interpPre++;
            didSubpixel = true;
            
        }
        else if(interpPost)
        {
            float d = gradPost_this / (gradPost_this - gradPost_post);
            best_match_x += d*incx;
            best_match_y += d*incy;
            best_match_err = best_match_err + 2*d*gradPost_this + (gradPost_post - gradPost_this)*d*d;
            if(enablePrintDebugInfo) stats->num_stereo_interpPost++;
            didSubpixel = true;
        }
        else
        {
            if(enablePrintDebugInfo) stats->num_stereo_interpNone++;
        }
    }
    
    
    // sampleDist is the distance in pixel at which the realVal's were sampled
    float sampleDist = GRADIENT_SAMPLE_DIST*rescaleFactor;
    
    float gradAlongLine = 0;
    float tmp = realVal_p2 - realVal_p1;  gradAlongLine+=tmp*tmp;
    tmp = realVal_p1 - realVal;  gradAlongLine+=tmp*tmp;
    tmp = realVal - realVal_m1;  gradAlongLine+=tmp*tmp;
    tmp = realVal_m1 - realVal_m2;  gradAlongLine+=tmp*tmp;
    
    gradAlongLine /= sampleDist*sampleDist;
    
    // check if interpolated error is OK. use evil hack to allow more error if there is a lot of gradient.
    if(best_match_err > (float)MAX_ERROR_STEREO + sqrtf( gradAlongLine)*20)
    {
        if(enablePrintDebugInfo) stats->num_stereo_invalid_bigErr++;
        return -3;
    }
    
    
    // ================= calc depth (in KF) ====================
    // * KinvP = Kinv * (x,y,1); where x,y are pixel coordinates of point we search for, in the KF.
    // * best_match_x = x-coordinate of found correspondence in the reference frame.
    
    float idnew_best_match;	// depth in the new image
    float alpha; // d(idnew_best_match) / d(disparity in pixel) == conputed inverse depth derived by the pixel-disparity.
    if(incx*incx>incy*incy)
    {
        float oldX = fxi*best_match_x+cxi;
        float nominator = (oldX*referenceFrame->otherToThis_t[2] - referenceFrame->otherToThis_t[0]);
        float dot0 = KinvP.dot(referenceFrame->otherToThis_R_row0);
        float dot2 = KinvP.dot(referenceFrame->otherToThis_R_row2);
        
        idnew_best_match = (dot0 - oldX*dot2) / nominator;
        alpha = incx*fxi*(dot0*referenceFrame->otherToThis_t[2] - dot2*referenceFrame->otherToThis_t[0]) / (nominator*nominator);
        
    }
    else
    {
        float oldY = fyi*best_match_y+cyi;
        
        float nominator = (oldY*referenceFrame->otherToThis_t[2] - referenceFrame->otherToThis_t[1]);
        float dot1 = KinvP.dot(referenceFrame->otherToThis_R_row1);
        float dot2 = KinvP.dot(referenceFrame->otherToThis_R_row2);
        
        idnew_best_match = (dot1 - oldY*dot2) / nominator;
        alpha = incy*fyi*(dot1*referenceFrame->otherToThis_t[2] - dot2*referenceFrame->otherToThis_t[1]) / (nominator*nominator);
        
    }
    
    
    
    
    
    if(idnew_best_match < 0)
    {
        if(enablePrintDebugInfo) stats->num_stereo_negative++;
        if(!allowNegativeIdepths)
            return -2;
    }
    
    if(enablePrintDebugInfo) stats->num_stereo_successfull++;
    
    // ================= calc var (in NEW image) ====================
    
    // calculate error from photometric noise
    float photoDispError = 4.0f * cameraPixelNoise2 / (gradAlongLine + DIVISION_EPS);
    
    float trackingErrorFac = 0.25f*(1.0f+referenceFrame->initialTrackedResidual);
    
    // calculate error from geometric noise (wrong camera pose / calibration)
    Eigen::Vector2f gradsInterp = getInterpolatedElement42(activeKeyFrame->gradients(0), u, v, width);
    float geoDispError = (gradsInterp[0]*epxn + gradsInterp[1]*epyn) + DIVISION_EPS;
    geoDispError = trackingErrorFac*trackingErrorFac*(gradsInterp[0]*gradsInterp[0] + gradsInterp[1]*gradsInterp[1]) / (geoDispError*geoDispError);
    
    
    //geoDispError *= (0.5 + 0.5 *result_idepth) * (0.5 + 0.5 *result_idepth);
    
    // final error consists of a small constant part (discretization error),
    // geometric and photometric error.
    result_var = alpha*alpha*((didSubpixel ? 0.05f : 0.5f)*sampleDist*sampleDist +  geoDispError + photoDispError);	// square to make variance
    
    if(plotStereoImages)
    {
        if(rand()%5==0)
        {
            //if(rand()%500 == 0)
            //	printf("geo: %f, photo: %f, alpha: %f\n", sqrt(geoDispError), sqrt(photoDispError), alpha, sqrt(result_var));
            
            
            //int idDiff = (keyFrame->pyramidID - referenceFrame->id);
            //cv::Scalar color = cv::Scalar(0,0, 2*idDiff);// bw
            
            //cv::Scalar color = cv::Scalar(sqrt(result_var)*2000, 255-sqrt(result_var)*2000, 0);// bw
            
            //			float eplLengthF = std::min((float)MIN_EPL_LENGTH_CROP,(float)eplLength);
            //			eplLengthF = std::max((float)MAX_EPL_LENGTH_CROP,(float)eplLengthF);
            //
            //			float pixelDistFound = sqrtf((float)((pReal[0]/pReal[2] - best_match_x)*(pReal[0]/pReal[2] - best_match_x)
            //					+ (pReal[1]/pReal[2] - best_match_y)*(pReal[1]/pReal[2] - best_match_y)));
            //
            float fac = best_match_err / ((float)MAX_ERROR_STEREO + sqrtf( gradAlongLine)*20);
            
            cv::Scalar color = cv::Scalar(255*fac, 255-255*fac, 0);// bw
            
            
            /*
             if(rescaleFactor > 1)
             color = cv::Scalar(500*(rescaleFactor-1),0,0);
             else
             color = cv::Scalar(0,500*(1-rescaleFactor),500*(1-rescaleFactor));
             */
            
            cv::line(debugImageStereoLines,cv::Point2f(pClose[0], pClose[1]),cv::Point2f(pFar[0], pFar[1]),color,1,8,0);
        }
    }
    
    result_idepth = idnew_best_match;
    
    result_eplLength = eplLength;
    
    return best_match_err;
    
}
    
bool DepthMapForRollingShutter::observeDepthCreate(const int &x, const int &y,
                                                   const int &idx,
                                                   RunningStats* const &stats,
                                                   bool debugPlot)
//  x : x-coord of the semi-dense point in the source/key image
//  y : y-coord of the semi-dense point in the source/key image
//  refFrame : target/ref/currentTrackingFrame image
{
    DepthMapPixelHypothesis* target = currentDepthMap+idx;
    
    Frame* refFrame = activeKeyFrameIsReactivated ? newest_referenceFrame : oldest_referenceFrame;
    
    if(refFrame->getTrackingParent() == activeKeyFrame)
    {
        bool* wasGoodDuringTracking = refFrame->refPixelWasGoodNoCreate();
        if(wasGoodDuringTracking != 0 && !wasGoodDuringTracking[(x >> SE3TRACKING_MIN_LEVEL) + (width >> SE3TRACKING_MIN_LEVEL)*(y >> SE3TRACKING_MIN_LEVEL)])
        {
            if(plotStereoImages)
                debugImageHypothesisHandling.at<cv::Vec3b>(y, x) = cv::Vec3b(255,0,0); // BLUE for SKIPPED NOT GOOD TRACKED
            return false;
        }
    }
    
#if 0// ORIGINAL (GLOBAL SHUTTER)
    
    float epx, epy; // epipolar line
    bool isGood = makeAndCheckEPL(x, y, refFrame, &epx, &epy, stats);
    if(!isGood) return false;
    
    if(enablePrintDebugInfo) stats->num_observe_create_attempted++;
    
    float new_u = x;
    float new_v = y;
    float result_idepth, result_var, result_eplLength;
    float error = doLineStereo(
                               new_u,new_v,epx,epy,
                               0.0f, 1.0f, 1.0f/MIN_DEPTH,
                               refFrame, refFrame->image(0),
                               result_idepth, result_var, result_eplLength, stats);
    
    if(error == -3 || error == -2)
    {
        target->blacklisted--;
        if(enablePrintDebugInfo) stats->num_observe_blacklisted++;
    }
    
    if(error < 0 || result_var > MAX_VAR)
        return false;

#else // ROLLING SHUTTER
    
    std::list<PointEpiCurve> buffer_ptsEpiCurve;
    if (useGTmotion_ == true)
    {
        // GT motion

        // Extract points on generalised epipolar line
        extractPointsOnGeneralEpipolarLineBand(x, y, refFrame,
                                               &buffer_ptsEpiCurve,
                                               stats, true);

    }
    else
    {

        // Original
        // Extract points on generalised epipolar line
        extractPointsOnGeneralEpipolarLineBand(x, y, refFrame,
                                               &buffer_ptsEpiCurve,
                                               stats);
    }
    
    // Copy image bands along the extracted points into a buffer of
    // PtsEpiCurve
    // which exact point to track, and where from.
    float min_idepth = 0.0f;
    float prior_idepth = 1.0f;
    float max_idepth = 1.0f/MIN_DEPTH;
    int numBuffer = copyImageBandsToBuffer(min_idepth, prior_idepth, max_idepth,
                           refFrame->image(0), &buffer_ptsEpiCurve, stats);
    
    // Search for a match in the buffer
//    const float min_idepth = 0.0f;
//    const float prior_idepth = 1.0f;
//    const float max_idepth = 1.0f/MIN_DEPTH;
    float result_idepth = 0, result_var = 0, result_eplLength = 0;
    PointEpiCurve result_ptEC;
    float error = -1;
    lsd_slam::StereoDebugInfo stereoDebugInfo;
    
    if (buffer_ptsEpiCurve.size() > 0 && numBuffer > 0)
    {
        error = doGeneralEpipolarLineStereo(&buffer_ptsEpiCurve,
                                            min_idepth,
                                            prior_idepth,
                                            max_idepth,
                                            refFrame,
                                            refFrame->image(0),
                                            result_idepth,
                                            result_var,
                                            result_eplLength = 0,
                                            stereoDebugInfo,
                                            result_ptEC,
                                            stats);
       
    }
    else
    {
        // No epipolar curve available
//        target->blacklisted--;
        return false;
        
    }
    
#if 0
    // Too close depth (Closer than 0.5m)
    if (result_idepth > 1.0/0.5)
    {
//        target->blacklisted--;
        return false;
    }
    
//    // Epl length to short
//    if (result_eplLength < 100)
//    {
//        target->blacklisted--;
//        return false;
//    }
#endif
    
    if ((result_idepth < 0.0) &&
        (lsd_slam::allowNegativeIdepths == false))
    {
        
//        target->blacklisted--;
        return false;
    
    }
    
#if 1 // AVOID NAN INVERSE DEPTH
    // Any nan value inverse depth should be avoided
    if (isnan(result_idepth) == true)
    {
        
//        target->blacklisted--;
        return false;
        
    }
#endif
    
#if 0 // AVOID NEAR ZERO INVERSE DEPTH (Farther than 20m)
    // near zero inverse depth should be avoided
    if ((result_idepth < 1.0/20.0) &&
        (result_idepth >= 0.0))
    {
        
//        target->blacklisted--;
        return false;
        
    }
#endif
    
    if (error == -3 || error == -2)
    {
        
        target->blacklisted--;
        if (enablePrintDebugInfo) stats->num_observe_blacklisted++;
        
    }
    
    if (error < 0 || result_var > MAX_VAR)
        return false;
    
    //--------------------------------------------------------------------------
#endif
    
    result_idepth = UNZERO(result_idepth);
    
    // add hypothesis
    *target = DepthMapPixelHypothesis(
                                      result_idepth,
                                      result_var,
                                      VALIDITY_COUNTER_INITIAL_OBSERVE);
    
    if(plotStereoImages)
        debugImageHypothesisHandling.at<cv::Vec3b>(y, x) = cv::Vec3b(255,255,255); // white for GOT CREATED
    
    if(enablePrintDebugInfo) stats->num_observe_created++;
    
    // Plot stereo with general epipolar curve
    if ((buffer_ptsEpiCurve.size() > 0) && (error > 0) &&
        (debugPlot == true))
    {
        
        cv::Mat tmpImage(height, width*2, CV_8UC3);
        tmpImage = plotGeneralStereoImages(activeKeyFrame,
                                           refFrame,
                                           &buffer_ptsEpiCurve,
                                           &result_ptEC);
        //Util::displayImage("plot general epipolar curve", tmpImage, false);
        
        for(int x=0;x<tmpImage.cols;x++)
            for(int y=0; y<80;y++)
                tmpImage.at<cv::Vec3b>(y,x) *= 0.5;

        char text[255];
        sprintf(text, "crt: idepth: %f (depth = %f), var: %f, eplL: %f, err: %f",
                result_idepth, 1.0/result_idepth,
                result_var, result_eplLength, error);
        cv::putText(tmpImage, text, cv::Point(10,10),
                    cv::FONT_HERSHEY_SIMPLEX, 0.4, cv::Scalar(0,255,0), 1, 8);
        sprintf(text, "(%f, %f) => (%f, %f)\n",
                result_ptEC.xSource, result_ptEC.ySource,
                result_ptEC.xTarget, result_ptEC.yTarget);
        cv::putText(tmpImage, text, cv::Point(10,23),
                    cv::FONT_HERSHEY_SIMPLEX, 0.4, cv::Scalar(0,255,0), 1, 8);
        sprintf(text, "incx, incy = %f %f;  rescaleFactor = %f",
                result_ptEC.incx,
                result_ptEC.incy,
                stereoDebugInfo.rescaleFactor);
        cv::putText(tmpImage, text, cv::Point(10,36),
                    cv::FONT_HERSHEY_SIMPLEX, 0.4, cv::Scalar(0,255,0), 1, 8);
        
        sprintf(text,
                "     geoDispErr: %f, photoDispErr: %f",
                stereoDebugInfo.geoDispError,
                stereoDebugInfo.photoDispError
                );
        cv::putText(tmpImage, text, cv::Point(10,49),
                    cv::FONT_HERSHEY_SIMPLEX, 0.4,
                    cv::Scalar(127,255,127), 1, 8);
        sprintf(text,
                "     alpha: %f, sampleDist: %f, gradAlongLine: %f",
                stereoDebugInfo.alpha,
                stereoDebugInfo.sampleDist,
                stereoDebugInfo.gradAlongLine
                );
        cv::putText(tmpImage, text, cv::Point(10,62),
                    cv::FONT_HERSHEY_SIMPLEX, 0.4,
                    cv::Scalar(127,255,127), 1, 8);
        
        sprintf(text,
                "[result_ptEC] incx: %f"
                ", incy: %f"
                ", incx_epl: %f"
                ", incy_epl: %f"
                ", nepx: %f"
                ", nepy: %f"
                ", rescaleFactor: %f",
                result_ptEC.incx,
                result_ptEC.incy,
                result_ptEC.incx_epl,
                result_ptEC.incy_epl,
                result_ptEC.nepx,
                result_ptEC.nepy,
                result_ptEC.rescaleFactor);
        cv::putText(tmpImage, text, cv::Point(10,75),
                    cv::FONT_HERSHEY_SIMPLEX, 0.4,
                    cv::Scalar(127,255,127), 1, 8);
        
        char buf[500];
        sprintf(buf, "create_key%04d_cur%04d_%04d.jpg",
                activeKeyFrame->id(), refFrame->id(), idx);
        cv::imwrite(buf, tmpImage);
        printf("Save: create_key%04d_cur%04d_%04d.jpg",
               activeKeyFrame->id(), refFrame->id(), idx);
        
    }
    
    return true;
}
    
bool DepthMapForRollingShutter::observeDepthUpdate(const int &x, const int &y, const int &idx, const float* keyFrameMaxGradBuf, RunningStats* const &stats, bool debugPlot)
{
    DepthMapPixelHypothesis* target = currentDepthMap+idx;
    Frame* refFrame;
    
    
    if(!activeKeyFrameIsReactivated)
    {
        if((int)target->nextStereoFrameMinID - referenceFrameByID_offset >= (int)referenceFrameByID.size())
        {
            if(plotStereoImages)
                debugImageHypothesisHandling.at<cv::Vec3b>(y, x) = cv::Vec3b(0,255,0);	// GREEN FOR skip
            
            if(enablePrintDebugInfo) stats->num_observe_skip_alreadyGood++;
            return false;
        }
        
        if((int)target->nextStereoFrameMinID - referenceFrameByID_offset < 0)
            refFrame = oldest_referenceFrame;
        else
            refFrame = referenceFrameByID[(int)target->nextStereoFrameMinID - referenceFrameByID_offset];
    }
    else
        refFrame = newest_referenceFrame;
    
    
    if(refFrame->getTrackingParent() == activeKeyFrame)
    {
        bool* wasGoodDuringTracking = refFrame->refPixelWasGoodNoCreate();
        if(wasGoodDuringTracking != 0 && !wasGoodDuringTracking[(x >> SE3TRACKING_MIN_LEVEL) + (width >> SE3TRACKING_MIN_LEVEL)*(y >> SE3TRACKING_MIN_LEVEL)])
        {
            if(plotStereoImages)
                debugImageHypothesisHandling.at<cv::Vec3b>(y, x) = cv::Vec3b(255,0,0); // BLUE for SKIPPED NOT GOOD TRACKED
            return false;
        }
    }

#if 0 // ORIGINAL (GLOBAL SHUTTER)
    
    float epx, epy;
    bool isGood = makeAndCheckEPL(x, y, refFrame, &epx, &epy, stats);
    if(!isGood) return false;
    
    // which exact point to track, and where from.
    float sv = sqrt(target->idepth_var_smoothed);
    float min_idepth = target->idepth_smoothed - sv*STEREO_EPL_VAR_FAC;
    float max_idepth = target->idepth_smoothed + sv*STEREO_EPL_VAR_FAC;

    if(min_idepth < 0) min_idepth = 0;
    if(max_idepth > 1/MIN_DEPTH) max_idepth = 1/MIN_DEPTH;
    
    stats->num_observe_update_attempted++;
    
    float result_idepth, result_var, result_eplLength;
    
    float error = doLineStereo(
                               x,y,epx,epy,
                               min_idepth, target->idepth_smoothed ,max_idepth,
                               refFrame, refFrame->image(0),
                               result_idepth, result_var, result_eplLength, stats);
    
#endif
    // ROLLING SHUTTER
    
    
    // which exact point to track, and where from.
    float sv = sqrt(target->idepth_var_smoothed);
#if 1 // ORIGINAL MIN and MAX DEPTH

    float min_idepth = target->idepth_smoothed - sv*STEREO_EPL_VAR_FAC;
    float max_idepth = target->idepth_smoothed + sv*STEREO_EPL_VAR_FAC;

#else // MIN and MAX DEPTH by non-inverse (ordinary) depth variance
    
    float min_idepth = target->idepth_smoothed - 12.0f*sv*STEREO_EPL_VAR_FAC;
    float max_idepth = target->idepth_smoothed + 12.0f*sv*STEREO_EPL_VAR_FAC;
    
#endif
    
    if(min_idepth < 0) min_idepth = 0;
    if(max_idepth > 1/MIN_DEPTH) max_idepth = 1/MIN_DEPTH;

    
    std::list<PointEpiCurve> GT_buffer_ptsEpiCurve;
    std::list<PointEpiCurve> buffer_ptsEpiCurve;
    if (useGTmotion_ == true)
    {
    
        // GT motion
        // Extract points on generalised epipolar line
        extractPointsOnGeneralEpipolarLineBand(x, y, refFrame,
                                               &GT_buffer_ptsEpiCurve,
                                               stats, true);

        // Copy image bands along the extracted points into a buffer of
        // PtsEpiCurve
        copyImageBandsToBuffer(min_idepth, target->idepth_smoothed, max_idepth,
                               refFrame->image(0),
                               &GT_buffer_ptsEpiCurve, stats);

    }
    else
    {

        // Estimation
        // Extract points on generalised epipolar line
        extractPointsOnGeneralEpipolarLineBand(x, y, refFrame,
                                               &buffer_ptsEpiCurve,
                                               stats, false);

        // Copy image bands along the extracted points into a buffer of
        // PtsEpiCurve
        copyImageBandsToBuffer(min_idepth, target->idepth_smoothed, max_idepth,
                               refFrame->image(0), &buffer_ptsEpiCurve, stats);

    }
    
    // Search for a match in the buffer
    const float prior_idepth = target->idepth_smoothed;
    float result_idepth = 0, result_var = 0, result_eplLength = 0;
    PointEpiCurve result_ptEC;
    float error = -1;
    StereoDebugInfo stereoDebugInfo;
    
    if (useGTmotion_ == true)
    {

        if (GT_buffer_ptsEpiCurve.size() > 0)
        {
            // GT motion
            error = doGeneralEpipolarLineStereo(&GT_buffer_ptsEpiCurve,
                                                min_idepth,
                                                prior_idepth,
                                                max_idepth,
                                                refFrame,
                                                refFrame->image(0),
                                                result_idepth,
                                                result_var,
                                                result_eplLength,
                                                stereoDebugInfo,
                                                result_ptEC,
                                                stats);
        }
        else
        {

            // No epipolar curve available
            target->blacklisted--;
            return false;

        }

    }
    else
    {
        if (buffer_ptsEpiCurve.size() > 0)
        {

            // Estimate
            error = doGeneralEpipolarLineStereo(&buffer_ptsEpiCurve,
                                                min_idepth,
                                                prior_idepth,
                                                max_idepth,
                                                refFrame,
                                                refFrame->image(0),
                                                result_idepth,
                                                result_var,
                                                result_eplLength,
                                                stereoDebugInfo,
                                                result_ptEC,
                                                stats);

        }

        else
        {

            // No epipolar curve available
//            target->blacklisted--;
            return false;

        }

    }

#if 0
    // Too close depth (Closer than 0.5m)
    if (result_idepth > 1.0/0.5)
    {
//        target->blacklisted--;
        return false;
    }
    
//    // Epl length to short
//    if (result_eplLength < 100)
//    {
//        target->blacklisted--;
//        return false;
//    }

#endif

    if ((result_idepth < 0.0) &&
        (lsd_slam::allowNegativeIdepths == false))
    {

//        target->blacklisted--;
        return false;

    }
    
#if 1 // AVOID NAN INVERSE DEPTH
    // Any nan value inverse depth should be avoided
    if (isnan(result_idepth) == true)
    {
        
//        target->blacklisted--;
        return false;
        
    }
#endif
    
#if 0 // AVOID NEAR ZERO INVERSE DEPTH (Farther than 10m)
    // near zero inverse depth should be avoided
    if ((result_idepth < 1.0/20.0) &&
        (result_idepth >= 0.0))
    {
        
//        target->blacklisted--;
        return false;
        
    }
#endif

    float diff = result_idepth - target->idepth_smoothed;
    
    // if oob: (really out of bounds)
    if(error == -1)
    {
        // do nothing, pixel got oob, but is still in bounds in original. I will want to try again.
        if(enablePrintDebugInfo) stats->num_observe_skip_oob++;
        
        if(plotStereoImages)
            debugImageHypothesisHandling.at<cv::Vec3b>(y, x) = cv::Vec3b(0,0,255);	// RED FOR OOB
        return false;
    }
    
    // if just not good for stereo (e.g. some inf / nan occured; has inconsistent minimum; ..)
    else if(error == -2)
    {
        if(enablePrintDebugInfo) stats->num_observe_skip_fail++;
        
        if(plotStereoImages)
            debugImageHypothesisHandling.at<cv::Vec3b>(y, x) = cv::Vec3b(255,0,255);	// PURPLE FOR NON-GOOD
        
        
        target->validity_counter -= VALIDITY_COUNTER_DEC;
        if(target->validity_counter < 0) target->validity_counter = 0;
        
        
        target->nextStereoFrameMinID = 0;
        
        target->idepth_var *= FAIL_VAR_INC_FAC;
        if(target->idepth_var > MAX_VAR)
        {
            target->isValid = false;
            target->blacklisted--;
        }
        return false;
    }
    
    // if not found (error too high)
    else if(error == -3)
    {
        if(enablePrintDebugInfo) stats->num_observe_notfound++;
        if(plotStereoImages)
            debugImageHypothesisHandling.at<cv::Vec3b>(y, x) = cv::Vec3b(0,0,0);	// BLACK FOR big not-found
        
        
        return false;
    }
    
    else if(error == -4)
    {
        if(plotStereoImages)
            debugImageHypothesisHandling.at<cv::Vec3b>(y, x) = cv::Vec3b(0,0,0);	// BLACK FOR big arithmetic error
        
        return false;
    }
    
    // if inconsistent
    else if(DIFF_FAC_OBSERVE*diff*diff > result_var + target->idepth_var_smoothed)
    {
        if(enablePrintDebugInfo) stats->num_observe_inconsistent++;
        if(plotStereoImages)
            debugImageHypothesisHandling.at<cv::Vec3b>(y, x) = cv::Vec3b(255,255,0);	// Turkoise FOR big inconsistent
        
        target->idepth_var *= FAIL_VAR_INC_FAC;
        if(target->idepth_var > MAX_VAR) target->isValid = false;
        
        return false;
    }
    else
    {
        // one more successful observation!
        if(enablePrintDebugInfo) stats->num_observe_good++;
        
        if(enablePrintDebugInfo) stats->num_observe_updated++;
        
        // Store old and update Info
        float old_idepth = target->idepth;
        float old_var = target->idepth_var;
        int old_validity_counter = target->validity_counter;
        float upd_idepth = result_idepth;
        float upd_var = result_var;
        
        // do textbook ekf update:
        // increase var by a little (prediction-uncertainty)
        float id_var = target->idepth_var*SUCC_VAR_INC_FAC;
        
        // update var with observation
        float w = result_var / (result_var + id_var);
        float new_idepth = (1-w)*result_idepth + w*target->idepth;
        target->idepth = UNZERO(new_idepth);
        
        if (isnan(target->idepth) == true)
        {
            
            printf("TARGET idepth = %f \n", target->idepth);
            printf("       new_idepth = %f\n", new_idepth);
            printf("       result_idepth = %f\n", result_idepth);
            printf("       result_var = %f\n", result_var);
            printf("       id_var = %f\n", id_var);
            printf("       validty_counter = %d\n", target->validity_counter);
            
        }
        
        // variance can only decrease from observation; never increase.
        id_var = id_var * w;
        if(id_var < target->idepth_var)
            target->idepth_var = id_var;
        
        // increase validity!
        target->validity_counter += VALIDITY_COUNTER_INC;
        float absGrad = keyFrameMaxGradBuf[idx];
        if(target->validity_counter > VALIDITY_COUNTER_MAX+absGrad*(VALIDITY_COUNTER_MAX_VARIABLE)/255.0f)
            target->validity_counter = VALIDITY_COUNTER_MAX+absGrad*(VALIDITY_COUNTER_MAX_VARIABLE)/255.0f;
        
        // increase Skip!
        if(result_eplLength < MIN_EPL_LENGTH_CROP)
        {
            float inc = activeKeyFrame->numFramesTrackedOnThis / (float)(activeKeyFrame->numMappedOnThis+5);
            if(inc < 3) inc = 3;
            
            inc +=  ((int)(result_eplLength*10000)%2);
            
            if(enablePrintDebugInfo) stats->num_observe_addSkip++;
            
            if(result_eplLength < 0.5*MIN_EPL_LENGTH_CROP)
                inc *= 3;
            
            
            target->nextStereoFrameMinID = refFrame->id() + inc;
        }
        
        // Store new Info
        new_idepth = target->idepth;
        float new_var = target->idepth_var;
        int new_validity_counter = target->validity_counter;
        
        if(plotStereoImages)
            debugImageHypothesisHandling.at<cv::Vec3b>(y, x) = cv::Vec3b(0,255,255); // yellow for GOT UPDATED
        
        
        // GT motion and save debug image
//        if (useGTmotion_ == true)
//        {

            // Plot stereo with general epipolar curve
            if ((debugPlot == true) && (error > 0))
            {

                cv::Mat tmpImage(height, width*2, CV_8UC3);
                if (useGTmotion_ == true)
                {
                    // Plot for GT motion
                    tmpImage = plotGeneralStereoImages(activeKeyFrame,
                                                       refFrame,
                                                       &GT_buffer_ptsEpiCurve,
                                                       &result_ptEC);
                
                }
                else
                {
                    // Plot for w/o GT motion
                    tmpImage = plotGeneralStereoImages(activeKeyFrame,
                                                       refFrame,
                                                       &buffer_ptsEpiCurve,
                                                       &result_ptEC);
                
                }
                //Util::displayImage("plot general epipolar curve", tmpImage, false);

                for(int x=0;x<tmpImage.cols;x++)
                    for(int y=0; y<90;y++)
                        tmpImage.at<cv::Vec3b>(y,x) *= 0.5;
                
                // Put info text on image
                char text[255];
                sprintf(text, "Left px (%f, %f) --> Right px (%f, %f)",
                        result_ptEC.xSource, result_ptEC.ySource,
                        result_ptEC.xTarget, result_ptEC.yTarget
                        );
                cv::putText(tmpImage, text, cv::Point(10,10),
                            cv::FONT_HERSHEY_SIMPLEX, 0.4,
                            cv::Scalar(165,255,255), 1, 8);

                sprintf(text, "old; idepth: %f (depth = %f), "
                        "var: %f, vldcnt: %d",
                        old_idepth, 1.0/old_idepth,
                        old_var, old_validity_counter
                        );
                cv::putText(tmpImage, text, cv::Point(10,23),
                            cv::FONT_HERSHEY_SIMPLEX, 0.4,
                            cv::Scalar(255,255,255), 1, 8);

                sprintf(text,
                        "upd; idepth: %f (depth = %f), "
                        "var: %f, (eplL: %f, err: %f)",
                        upd_idepth, 1.0/upd_idepth,
                        upd_var, result_eplLength, error
                        );
                cv::putText(tmpImage, text, cv::Point(10,36),
                            cv::FONT_HERSHEY_SIMPLEX, 0.4,
                            cv::Scalar(127,255,127), 1, 8);

                sprintf(text,
                        "     geoDispErr: %f, photoDispErr: %f"
                        ", target.initialTrackedResidual: %f",
                        stereoDebugInfo.geoDispError,
                        stereoDebugInfo.photoDispError,
                        stereoDebugInfo.initialTrackedResidual
                        );
                cv::putText(tmpImage, text, cv::Point(10,49),
                            cv::FONT_HERSHEY_SIMPLEX, 0.4,
                            cv::Scalar(127,255,127), 1, 8);
                sprintf(text,
                        "     alpha: %f, sampleDist: %f, gradAlongLine: %f"
                        ", trackingErrorFac: %f"
                        ", geoDispError_denom: %f"
                        ", gradsInterpX: %f"
                        ", gradsInterpY: %f",
                        stereoDebugInfo.alpha,
                        stereoDebugInfo.sampleDist,
                        stereoDebugInfo.gradAlongLine,
                        stereoDebugInfo.trackingErrorFac,
                        stereoDebugInfo.geoDispError_denom,
                        stereoDebugInfo.gradsInterpX,
                        stereoDebugInfo.gradsInterpY
                        );
                cv::putText(tmpImage, text, cv::Point(10,62),
                            cv::FONT_HERSHEY_SIMPLEX, 0.4,
                            cv::Scalar(127,255,127), 1, 8);
                
                sprintf(text,
                        "[result_ptEC] incx: %f"
                        ", incy: %f"
                        ", incx_epl: %f"
                        ", incy_epl: %f"
                        ", nepx: %f"
                        ", nepy: %f"
                        ", rescaleFactor: %f",
                        result_ptEC.incx,
                        result_ptEC.incy,
                        result_ptEC.incx_epl,
                        result_ptEC.incy_epl,
                        result_ptEC.nepx,
                        result_ptEC.nepy,
                        result_ptEC.rescaleFactor);
                cv::putText(tmpImage, text, cv::Point(10,75),
                            cv::FONT_HERSHEY_SIMPLEX, 0.4,
                            cv::Scalar(127,255,127), 1, 8);

                sprintf(text,
                        "new; idepth: %f (depth = %f), "
                        "var: %f, vldcnt: %d",
                        new_idepth, 1.0/new_idepth,
                        new_var, new_validity_counter
                        );
                cv::putText(tmpImage, text, cv::Point(10,88),
                            cv::FONT_HERSHEY_SIMPLEX, 0.4,
                            cv::Scalar(255,255,165), 1, 8);

                // Save the image
                char buf[500];
                sprintf(buf, "update_key%04d_cur%04d_%04d.jpg",
                        activeKeyFrame->id(), refFrame->id(), idx);
                cv::imwrite(buf, tmpImage);

            }

//        } // end if useGTmotion
        
        return true;
    }
}
    
bool DepthMapForRollingShutter::makeAndCheckEPL(
        const int x, // x-coord of semi-dense points in the source/key image
        const int y, // y-coord of semi-dense points in the source/key image
        const Frame* const ref, // the target/ref/currentTrackingFrame image
        float* pepx,
        float* pepy,
        RunningStats* const stats)
{
    
    int idx = x+y*width;
    
    // ======= make epl ========
    // calculate the plane spanned by the two camera centers and the point (x,y,1)
    // intersect it with the keyframe's image plane (at depth=1)
    //--------------------------------------------------------------------------
    // Note: Generally, e' = K*t when P = K[I|0] and P' = K[R|t]
    //       because of P'*[0 0 0 1]'.
    //
    //       Wee define P and P' are the target and source, respectively.
    //       Now, we assume P is a target image (ref) and P' is a source
    //       image (keyframe). The point x' = (x, y, 1) is defined in the source
    //       image. Also, note that P = K[I|0] and P' = K[R|t].
    //       Therefore, the epipole e' = K*t is an epipole in the source image.
    //       Then, the epipolar line in the source image found by the cross
    //       product of x' and e' as follows:
    //
    //          cross(e', x') = [cy*t2 + fy*t1 - t2*y
    //                           t2*x - fx*t0 - cx*t2
    //                           y*(cx*t2 + fx*t0) - x*(cy*t2 + fy*t1)];
    //
    //          cross(x', e') = [t2*y - fy*t1 - cy*t2
    //                           cx*t2 + fx*t0 - t2*x
    //                           x*(cy*t2 + fy*t1) - y*(cx*t2 + fx*t0)];
    //
    //
    float epx = - fx * ref->thisToOther_t[0] + ref->thisToOther_t[2]*(x - cx);
    float epy = - fy * ref->thisToOther_t[1] + ref->thisToOther_t[2]*(y - cy);
    //--------------------------------------------------------------------------
    
    if(isnanf(epx+epy))
        return false;
    
    
    // ======== check epl length =========
    float eplLengthSquared = epx*epx+epy*epy;
    if(eplLengthSquared < MIN_EPL_LENGTH_SQUARED)
    {
        if(enablePrintDebugInfo) stats->num_observe_skipped_small_epl++;
        return false;
    }
    
    
    // ===== check epl-grad magnitude ======
    float gx = activeKeyFrameImageData[idx+1] - activeKeyFrameImageData[idx-1];
    float gy = activeKeyFrameImageData[idx+width] - activeKeyFrameImageData[idx-width];
    float eplGradSquared = gx * epx + gy * epy;
    eplGradSquared = eplGradSquared*eplGradSquared / eplLengthSquared;	// square and norm with epl-length
    
    if(eplGradSquared < MIN_EPL_GRAD_SQUARED)
    {
        if(enablePrintDebugInfo) stats->num_observe_skipped_small_epl_grad++;
        return false;
    }
    
    
    // ===== check epl-grad angle ======
    if(eplGradSquared / (gx*gx+gy*gy) < MIN_EPL_ANGLE_SQUARED)
    {
        if(enablePrintDebugInfo) stats->num_observe_skipped_small_epl_angle++;
        return false;
    }
    
    
    // ===== DONE - return "normalized" epl =====
    float fac = GRADIENT_SAMPLE_DIST / sqrt(eplLengthSquared);
    *pepx = epx * fac;
    *pepy = epy * fac;
    
    return true;
}

bool DepthMapForRollingShutter::makeAndCheckEPL_ver2(
    const int x, // x-coord of semi-dense points in the source/key image
    const int y, // y-coord of semi-dense points in the source/key image
    const Frame* const ref, // the target/ref/currentTrackingFrame image
    Sophus::Sim3d Xi_worldToKF,
    float* pepx,
    float* pepy,
    RunningStats* const stats)
{
    
    //--------------------------------------------------------------------------
    // Prepare stereo for a pixel with the row pose
    //--------------------------------------------------------------------------
    // Note: here, the "referenceFrame" is an unmapped tracking frame.
    //       Therefore, "refToKf" is same as "FrameToRef" in trackFrame().
    Frame ref_copy = *ref;
    Sim3 refToKf = Xi_worldToKF*ref_copy.pose->getCamToWorld();
    
    ref_copy.prepareForRollingShutterStereo(activeKeyFrame,
                                            refToKf, y, K, 0);
    //--------------------------------------------------------------------------
    
    
    int idx = x+y*width;
    
    // ======= make epl ========
    // calculate the plane spanned by the two camera centers and the point (x,y,1)
    // intersect it with the keyframe's image plane (at depth=1)
    float epx = - fx * ref_copy.thisToOther_t[0] + ref_copy.thisToOther_t[2]*(x - cx);
    float epy = - fy * ref_copy.thisToOther_t[1] + ref_copy.thisToOther_t[2]*(y - cy);
    
    if(isnanf(epx+epy))
        return false;
    
    
    // ======== check epl length =========
    float eplLengthSquared = epx*epx+epy*epy;
    if(eplLengthSquared < MIN_EPL_LENGTH_SQUARED)
    {
        if(enablePrintDebugInfo) stats->num_observe_skipped_small_epl++;
        return false;
    }
    
    
    // ===== check epl-grad magnitude ======
    float gx = activeKeyFrameImageData[idx+1] - activeKeyFrameImageData[idx-1];
    float gy = activeKeyFrameImageData[idx+width] - activeKeyFrameImageData[idx-width];
    float eplGradSquared = gx * epx + gy * epy;
    eplGradSquared = eplGradSquared*eplGradSquared / eplLengthSquared;	// square and norm with epl-length
    
    if(eplGradSquared < MIN_EPL_GRAD_SQUARED)
    {
        if(enablePrintDebugInfo) stats->num_observe_skipped_small_epl_grad++;
        return false;
    }
    
    
    // ===== check epl-grad angle ======
    if(eplGradSquared / (gx*gx+gy*gy) < MIN_EPL_ANGLE_SQUARED)
    {
        if(enablePrintDebugInfo) stats->num_observe_skipped_small_epl_angle++;
        return false;
    }
    
    
    // ===== DONE - return "normalized" epl =====
    float fac = GRADIENT_SAMPLE_DIST / sqrtf(eplLengthSquared);
    *pepx = epx * fac;
    *pepy = epy * fac;
    
    return true;
}
   
bool DepthMapForRollingShutter::makeAndCheckEPL_forGeneralEpipolarLine(
    const int x, // x-coord of semi-dense points in the source/key image
    const int y, // y-coord of semi-dense points in the source/key image
    Eigen::Matrix<double,3,1> thisToOther_t, // translation this to other
    double* pepx, // x-coord of the normal of epipolar line
    double* pepy, // y-coord of the normal of epipolar line
    RunningStats* const stats)
{
    
    int idx = x+y*width;
    
    // ======= make epl ========
    // calculate the plane spanned by the two camera centers and the point (x,y,1)
    // intersect it with the keyframe's image plane (at depth=1)
    //--------------------------------------------------------------------------
    // Note: Generally, e' = K*t when P = K[I|0] and P' = K[R|t]
    //       because of P'*[0 0 0 1]'.
    //
    //       We define P and P' are the target and source, respectively.
    //       Now, we assume P is a target image (ref) and P' is a source
    //       image (keyframe). The point x' = (x, y, 1) is defined in the source
    //       image. Also, note that P = K[I|0] and P' = K[R|t].
    //       Therefore, the epipole e' = K*t is an epipole in the source image.
    //       Then, the epipolar line in the source image found by the cross
    //       product of x' and e' as follows:
    //
    //          cross(e', x') = [cy*t2 + fy*t1 - t2*y
    //                           t2*x - fx*t0 - cx*t2
    //                           y*(cx*t2 + fx*t0) - x*(cy*t2 + fy*t1)];
    //
    //          cross(x', e') = [t2*y - fy*t1 - cy*t2
    //                           cx*t2 + fx*t0 - t2*x
    //                           x*(cy*t2 + fy*t1) - y*(cx*t2 + fx*t0)];
    //
    //
    float epx = -fx*(float)thisToOther_t[0] + (float)thisToOther_t[2]*(x - cx);
    float epy = -fy*(float)thisToOther_t[1] + (float)thisToOther_t[2]*(y - cy);
    //--------------------------------------------------------------------------
    
    if(isnanf(epx+epy))
        return false;
    
    
    // ======== check epl length =========
    float eplLengthSquared = epx*epx+epy*epy;
    if(eplLengthSquared < MIN_EPL_LENGTH_SQUARED)
    {
        if(enablePrintDebugInfo) stats->num_observe_skipped_small_epl++;
        return false;
    }
    
    
    // ===== check epl-grad magnitude ======
    float gx = activeKeyFrameImageData[idx+1] - activeKeyFrameImageData[idx-1];
    float gy = activeKeyFrameImageData[idx+width] - activeKeyFrameImageData[idx-width];
    float eplGradSquared = gx * epx + gy * epy;
    eplGradSquared = eplGradSquared*eplGradSquared / eplLengthSquared;	// square and norm with epl-length
    
#if 1 // DO NOT CHECK FOR PARALLEL EPIPOLAR LINE
#else // ORIGINAL
    if(eplGradSquared < MIN_EPL_GRAD_SQUARED)
    {
        if(enablePrintDebugInfo) stats->num_observe_skipped_small_epl_grad++;
        return false;
    }
    
    
    // ===== check epl-grad angle ======
    if(eplGradSquared / (gx*gx+gy*gy) < MIN_EPL_ANGLE_SQUARED)
    {
        if(enablePrintDebugInfo) stats->num_observe_skipped_small_epl_angle++;
        return false;
    }
#endif
    
    // ===== DONE - return "normalized" epl =====
    float fac = GRADIENT_SAMPLE_DIST / sqrtf(eplLengthSquared);
    *pepx = epx * fac;
    *pepy = epy * fac;
    
    return true;
}

void DepthMapForRollingShutter::observeDepthRow(int yMin, int yMax, RunningStats* stats)
{
    
    const float* keyFrameMaxGradBuf = activeKeyFrame->maxGradients(0);
    
    int successes = 0;
    
    int yRand = rand()%(yMax - yMin) + yMin;
    
    for(int y=yMin;y<yMax; y++)
        for(int x=3;x<width-3;x++)
        {
            int idx = x+y*width;
            DepthMapPixelHypothesis* target = currentDepthMap+idx;
            bool hasHypothesis = target->isValid;
            
            // ======== 1. check absolute grad =========
            if(hasHypothesis && keyFrameMaxGradBuf[idx] < MIN_ABS_GRAD_DECREASE)
            {
                target->isValid = false;
                continue;
            }
            
            if(keyFrameMaxGradBuf[idx] < MIN_ABS_GRAD_CREATE || target->blacklisted < MIN_BLACKLIST)
                continue;
            
            
            bool success;
            if(!hasHypothesis)
            {
#if 1 // ENABLE CREATE DEPTH (Origin)
                if (y == yRand && fmod(x,150) == 0)
                {
                    success = observeDepthCreate(x, y, idx, stats, true);
                }
                else
                {
                    success = observeDepthCreate(x, y, idx, stats, false);
                }
#else // DISABLE CREATE DEPTH
                success = false;
#endif
            }
            else
            {
#if 1 // ENABLE UPDATE DEPTH (Origin)
                if (y == yRand && fmod(x,150) == 0)
                {
                    success = observeDepthUpdate(x, y, idx, keyFrameMaxGradBuf, stats, true);
                }
                else
                {
                    success = observeDepthUpdate(x, y, idx, keyFrameMaxGradBuf, stats, false);
                }
#else // DISABLE UPDATE DEPTH
                success = false;
#endif
            }
            
            if(success)
                successes++;
        }
    
}
    
void DepthMapForRollingShutter::updateKeyframe(std::deque< std::shared_ptr<Frame> > referenceFrames)
{
    assert(isValid());
    
    struct timeval tv_start_all, tv_end_all;
    gettimeofday(&tv_start_all, NULL);
    
    oldest_referenceFrame = referenceFrames.front().get();
    newest_referenceFrame = referenceFrames.back().get();
    referenceFrameByID.clear();
    referenceFrameByID_offset = oldest_referenceFrame->id();
    
    for(std::shared_ptr<Frame> frame : referenceFrames)
    {
        assert(frame->hasTrackingParent());
        
        if(frame->getTrackingParent()->id() != activeKeyFrame->id())
        {
            printf("WARNING: updating frame %d (with %d, "
                   "which was tracked on a different frame (%d).\n"
                   "While this should work, it is not recommended.\n",
                   activeKeyFrame->id(), frame->id(),
                   frame->getTrackingParent()->id());
            fprintf(stderr, "WARNING: updating frame %d with %d, "
                   "which was tracked on a different frame (%d).\n"
                   "While this should work, it is not recommended.\n",
                   activeKeyFrame->id(), frame->id(),
                   frame->getTrackingParent()->id());
            printf("  activeKetframe %d Parent %d\n",
                   activeKeyFrame->id(),
                   activeKeyFrame->getTrackingParent()->id());
            printf("  frame %d Parent %d\n",
                   frame->id(),
                   frame->getTrackingParent()->id());
            fprintf(stderr, "  activeKetframe %d Parent %d\n",
                   activeKeyFrame->id(),
                   activeKeyFrame->getTrackingParent()->id());
            fprintf(stderr, "  frame %d Parent %d\n",
                   frame->id(),
                   frame->getTrackingParent()->id());

#if 0 // DROP
            fprintf(stderr, "  WARNING: No depth updating frame %d\n",
                    frame->id());
            return; // DO NOTHING
#endif
            
        }
        
        Sim3 refToKf;
        if(frame->pose->trackingParent->frameID == activeKeyFrame->id())
        {
            
            refToKf = frame->pose->thisToParent_raw;
#if 1//DEBUG
#else
            
            printf("refToKf = \n");
            std::cout << refToKf.log().transpose() << std::endl;
            
#endif
            
        }
        else
        {
            refToKf = activeKeyFrame->getScaledCamToWorld().inverse() *  frame->getScaledCamToWorld();
#if 1//DEBUG
            
            printf("activeKeyFrame.pose = \n");
            std::cout << activeKeyFrame->getScaledCamToWorld().log().transpose()
            << std::endl;
            
            printf("frame.pose = \n");
            std::cout << frame->getScaledCamToWorld().log().transpose()
            << std::endl;
            
            printf("refToKf = \n");
            std::cout << refToKf.log().transpose() << std::endl;
            
#endif
        }
        
        
        frame->prepareForStereoWith(activeKeyFrame, refToKf, K, 0);
        
        while((int)referenceFrameByID.size() + referenceFrameByID_offset <= frame->id())
            referenceFrameByID.push_back(frame.get());
    }
    
    resetCounters();
    
    
    if(plotStereoImages)
    {
        cv::Mat keyFrameImage(activeKeyFrame->height(), activeKeyFrame->width(), CV_32F, const_cast<float*>(activeKeyFrameImageData));
        keyFrameImage.convertTo(debugImageHypothesisHandling, CV_8UC1);
        cv::cvtColor(debugImageHypothesisHandling, debugImageHypothesisHandling, CV_GRAY2RGB);
        
        cv::Mat oldest_refImage(oldest_referenceFrame->height(), oldest_referenceFrame->width(), CV_32F, const_cast<float*>(oldest_referenceFrame->image(0)));
        cv::Mat newest_refImage(newest_referenceFrame->height(), newest_referenceFrame->width(), CV_32F, const_cast<float*>(newest_referenceFrame->image(0)));
        cv::Mat rfimg = 0.5f*oldest_refImage + 0.5f*newest_refImage;
        rfimg.convertTo(debugImageStereoLines, CV_8UC1);
        cv::cvtColor(debugImageStereoLines, debugImageStereoLines, CV_GRAY2RGB);
    }
    
    struct timeval tv_start, tv_end;
    

    fprintf(stderr, "  Depth mapping started frame %d\n",
            referenceFrames.front().get()->id());
    gettimeofday(&tv_start, NULL);

#if 0 //DISABLE depth update or create
#else //ORIGINAL

    observeDepth();
#endif
    
    gettimeofday(&tv_end, NULL);
    
    fprintf(stderr, "  Depth mapping completed frame %d (%5.4f sec)\n",
            referenceFrames.front().get()->id(),
            (double)(tv_end.tv_sec - tv_start.tv_sec) +
            (double)(tv_end.tv_usec - tv_start.tv_usec)*1e-06);
    msObserve = 0.9*msObserve + 0.1*((tv_end.tv_sec-tv_start.tv_sec)*1000.0f + (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
    nObserve++;
    

    
    //if(rand()%10==0)
    {
        gettimeofday(&tv_start, NULL);
        regularizeDepthMapFillHoles();
        gettimeofday(&tv_end, NULL);
        msFillHoles = 0.9*msFillHoles + 0.1*((tv_end.tv_sec-tv_start.tv_sec)*1000.0f + (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
        nFillHoles++;
    }
    
    
    gettimeofday(&tv_start, NULL);
    regularizeDepthMap(false, VAL_SUM_MIN_FOR_KEEP);
    gettimeofday(&tv_end, NULL);
    msRegularize = 0.9*msRegularize + 0.1*((tv_end.tv_sec-tv_start.tv_sec)*1000.0f + (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
    nRegularize++;
    
    
    // Update depth in keyframe
    if(!activeKeyFrame->depthHasBeenUpdatedFlag)
    {
        
#if 0// Update the current depth from GT
        
        
        setDepthFromGT();
        
#else // ORIGINAL - do not update every keyframe depth from GT (ground truth)
      // The depth will be updated or created by stereo matching
        
#endif
        
        gettimeofday(&tv_start, NULL);
        activeKeyFrame->setDepth(currentDepthMap);
        gettimeofday(&tv_end, NULL);
        msSetDepth = 0.9*msSetDepth + 0.1*((tv_end.tv_sec-tv_start.tv_sec)*1000.0f + (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
        nSetDepth++;

    }
    
    
    gettimeofday(&tv_end_all, NULL);
    msUpdate = 0.9*msUpdate + 0.1*((tv_end_all.tv_sec-tv_start_all.tv_sec)*1000.0f + (tv_end_all.tv_usec-tv_start_all.tv_usec)/1000.0f);
    nUpdate++;
    
    
    activeKeyFrame->numMappedOnThis++;
    activeKeyFrame->numMappedOnThisTotal++;
    
    
    if(plotStereoImages)
    {
        Util::displayImage( "Stereo Key Frame", debugImageHypothesisHandling, false );
        Util::displayImage( "Stereo Reference Frame", debugImageStereoLines, false );
        
        // Save stereo plot images as a file
        char buf[255];
        sprintf(buf, "stereo_key_frame_%06d.png", activeKeyFrame->id());
        cv::imwrite(buf, debugImageHypothesisHandling);
        sprintf(buf, "stereo_ref_frame_%06d.png", activeKeyFrame->id());
        cv::imwrite(buf, debugImageStereoLines);
        
    }
    
    
    
    if(enablePrintDebugInfo && printLineStereoStatistics)
    {
        printf("ST: calls %6d, comp %6d, int %7d; good %6d (%.0f%%), neg %6d (%.0f%%); interp %6d / %6d / %6d\n",
               runningStats.num_stereo_calls,
               runningStats.num_stereo_comparisons,
               runningStats.num_pixelInterpolations,
               runningStats.num_stereo_successfull,
               100*runningStats.num_stereo_successfull / (float) runningStats.num_stereo_calls,
               runningStats.num_stereo_negative,
               100*runningStats.num_stereo_negative / (float) runningStats.num_stereo_successfull,
               runningStats.num_stereo_interpPre,
               runningStats.num_stereo_interpNone,
               runningStats.num_stereo_interpPost);
    }
    if(enablePrintDebugInfo && printLineStereoFails)
    {
        printf("ST-ERR: oob %d (scale %d, inf %d, near %d); err %d (%d uncl; %d end; zro: %d btw, %d no, %d two; %d big)\n",
               runningStats.num_stereo_rescale_oob+
               runningStats.num_stereo_inf_oob+
               runningStats.num_stereo_near_oob,
               runningStats.num_stereo_rescale_oob,
               runningStats.num_stereo_inf_oob,
               runningStats.num_stereo_near_oob,
               runningStats.num_stereo_invalid_unclear_winner+
               runningStats.num_stereo_invalid_atEnd+
               runningStats.num_stereo_invalid_inexistantCrossing+
               runningStats.num_stereo_invalid_noCrossing+
               runningStats.num_stereo_invalid_twoCrossing+
               runningStats.num_stereo_invalid_bigErr,
               runningStats.num_stereo_invalid_unclear_winner,
               runningStats.num_stereo_invalid_atEnd,
               runningStats.num_stereo_invalid_inexistantCrossing,
               runningStats.num_stereo_invalid_noCrossing,
               runningStats.num_stereo_invalid_twoCrossing,
               runningStats.num_stereo_invalid_bigErr);
    }
}
    
void DepthMapForRollingShutter::observeDepth()
{
    
#if 1 // ORIGINAL
    
    threadReducer.reduce(boost::bind(&DepthMapForRollingShutter::
                                     observeDepthRow, this, _1, _2, _3),
                         3,
                         height-3,
                         10 // Number of threads (ORIGINAL)
                         );

#else // NO THREAD

    RunningStats stats;
    observeDepthRow(3, height - 3, &stats);

#endif
    
    if(enablePrintDebugInfo && printObserveStatistics)
    {
        printf("OBSERVE (%d): %d / %d created; %d / %d updated; %d skipped; %d init-blacklisted\n",
               activeKeyFrame->id(),
               runningStats.num_observe_created,
               runningStats.num_observe_create_attempted,
               runningStats.num_observe_updated,
               runningStats.num_observe_update_attempted,
               runningStats.num_observe_skip_alreadyGood,
               runningStats.num_observe_blacklisted
               );
    }
    
    
    if(enablePrintDebugInfo && printObservePurgeStatistics)
    {
        printf("OBS-PRG (%d): Good: %d; inconsistent: %d; notfound: %d; oob: %d; failed: %d; addSkip: %d;\n",
               activeKeyFrame->id(),
               runningStats.num_observe_good,
               runningStats.num_observe_inconsistent,
               runningStats.num_observe_notfound,
               runningStats.num_observe_skip_oob,
               runningStats.num_observe_skip_fail,
               runningStats.num_observe_addSkip
               );
    }
}

void DepthMapForRollingShutter::regularizeDepthMapFillHoles()
{
    
    buildRegIntegralBuffer();
    
    runningStats.num_reg_created=0;
    
    memcpy(otherDepthMap,currentDepthMap,width*height*sizeof(DepthMapPixelHypothesis));
    threadReducer.reduce(boost::bind(&DepthMapForRollingShutter::regularizeDepthMapFillHolesRow, this, _1, _2, _3), 3, height-2, 10);
    if(enablePrintDebugInfo && printFillHolesStatistics)
        printf("FillHoles (discreteDepth): %d created\n",
               runningStats.num_reg_created);
    
}
    
void DepthMapForRollingShutter::buildRegIntegralBuffer()
{
    threadReducer.reduce(boost::bind(&DepthMapForRollingShutter::buildRegIntegralBufferRow1, this, _1, _2,_3), 0, height);
    
    int* validityIntegralBufferPT = validityIntegralBuffer;
    int* validityIntegralBufferPT_T = validityIntegralBuffer+width;
    
    int wh = height*width;
    for(int idx=width;idx<wh;idx++)
        *(validityIntegralBufferPT_T++) += *(validityIntegralBufferPT++);
    
}

void DepthMapForRollingShutter::regularizeDepthMapFillHolesRow(int yMin, int yMax, RunningStats* stats)
{
    // =========== regularize fill holes
    const float* keyFrameMaxGradBuf = activeKeyFrame->maxGradients(0);
    
    for(int y=yMin; y<yMax; y++)
    {
        for(int x=3;x<width-2;x++)
        {
            int idx = x+y*width;
            DepthMapPixelHypothesis* dest = otherDepthMap + idx;
            if(dest->isValid) continue;
            if(keyFrameMaxGradBuf[idx]<MIN_ABS_GRAD_DECREASE) continue;
            
            int* io = validityIntegralBuffer + idx;
            int val = io[2+2*width] - io[2-3*width] - io[-3+2*width] + io[-3-3*width];
            
            
            if((dest->blacklisted >= MIN_BLACKLIST && val > VAL_SUM_MIN_FOR_CREATE) || val > VAL_SUM_MIN_FOR_UNBLACKLIST)
            {
                float sumIdepthObs = 0, sumIVarObs = 0;
                int num = 0;
                
                DepthMapPixelHypothesis* s1max = otherDepthMap + (x-2) + (y+3)*width;
                for (DepthMapPixelHypothesis* s1 = otherDepthMap + (x-2) + (y-2)*width; s1 < s1max; s1+=width)
                    for(DepthMapPixelHypothesis* source = s1; source < s1+5; source++)
                    {
                        if(!source->isValid) continue;
                        
                        sumIdepthObs += source->idepth /source->idepth_var;
                        sumIVarObs += 1.0f/source->idepth_var;
                        num++;
                    }
                
                float idepthObs = sumIdepthObs / sumIVarObs;
                idepthObs = UNZERO(idepthObs);
                
                currentDepthMap[idx] =
                DepthMapPixelHypothesis(
                                        idepthObs,
                                        VAR_RANDOM_INIT_INITIAL,
                                        0);
                
                if(enablePrintDebugInfo) stats->num_reg_created++;
            }
        }
    }
}

void DepthMapForRollingShutter::buildRegIntegralBufferRow1(int yMin, int yMax, RunningStats* stats)
{
    // ============ build inegral buffers
    int* validityIntegralBufferPT = validityIntegralBuffer+yMin*width;
    DepthMapPixelHypothesis* ptSrc = currentDepthMap+yMin*width;
    for(int y=yMin;y<yMax;y++)
    {
        int validityIntegralBufferSUM = 0;
        
        for(int x=0;x<width;x++)
        {
            if(ptSrc->isValid)
                validityIntegralBufferSUM += ptSrc->validity_counter;
            
            *(validityIntegralBufferPT++) = validityIntegralBufferSUM;
            ptSrc++;
        }
    }
}

void DepthMapForRollingShutter::setDepthFromGT()
{
//    threadReducer.reduce(
//        boost::bind(&DepthMapForRollingShutter::setDepthFromGTRow,
//                         this, _1, _2, _3), 3, height - 3, 10);
    RunningStats stats;
    setDepthFromGTRow(3, height - 3, &stats);

}

void DepthMapForRollingShutter::setDepthFromGTRow(
        int yMin, int yMax, RunningStats* stats)
{

    if (activeKeyFrame != 0)
    {
        
        const float* keyFrameMaxGradBuf = activeKeyFrame->maxGradients(0);
        
        // Read from GT file
        std::string filename = activeKeyFrame->depthFilename;
        std::shared_ptr<float> depthGT = readDepthFile(filename);
        
        for (int y = yMin; y<yMax; y++)
        {
            
            for (int x=3; x<width - 3; x++)
            {
                
                int idx = x + y*width;
                DepthMapPixelHypothesis* target = currentDepthMap + idx;
                bool hasHypothesis = target->isValid;
                
                if (hasHypothesis &&
                    keyFrameMaxGradBuf[idx] < MIN_ABS_GRAD_DECREASE)
                {
                    
                    target->isValid = false;
                    continue;
                    
                }
                
                if (keyFrameMaxGradBuf[idx] < MIN_ABS_GRAD_CREATE
                    || target->blacklisted < MIN_BLACKLIST)
                    continue;
#if 1 // GT DEPTH for ALL
#else // GT DEPTH ONLY for existing
                if (hasHypothesis)
#endif
                {
                    
                    //----------------------------------------------------------
                    // Set from the depth GT loaded from a file
                    float depthValue = depthGT.get()[idx];
                    
                    // Set from GT
                    target->isValid = true;
                    target->blacklisted = 0;
                    target->nextStereoFrameMinID = 0;
                    target->validity_counter = 20;
                    target->idepth = 1.0f / depthValue; // Set from the GT file
                    target->idepth_smoothed = 1.0f / depthValue; // Set from the GT
                    target->idepth_var = VAR_GT_INIT_INITIAL;
                    target->idepth_var_smoothed = VAR_GT_INIT_INITIAL;
                    //----------------------------------------------------------
                    
                }
                
            } // end for (int x ...
            
        } // end for (int y ...
        
    } // end if (activeKeyFrame ...
    
}
                
std::shared_ptr<float>
DepthMapForRollingShutter::readDepthFile(std::string filename)
{
    
    //float *depthInput = (float*)malloc(sizeof(float)*width*height);
    std::shared_ptr<float> depthInput(new float[width*height]);
    printf("---------------------------------------------\n");
    printf("Reading %s\n", filename.c_str());
    
    FILE* fid = fopen(filename.c_str(),"r");
    float *ptrDepth = depthInput.get();
    float depthVal;
    int u = 0;
    int v = 0;
    while (fscanf(fid, "%f", &depthVal)!=EOF)
    {
        
        if (u == width) // u for column index
        {
            v++; // v for row index
            u = 0;
        }
        
        // Set depth value from file via calibration information
        float u_u0_by_fx = ((float)u - cx)/fx;
        float v_v0_by_fy = ((float)v - cy)/fy;
        float u_u0_by_fx_sq = u_u0_by_fx*u_u0_by_fx;
        float v_v0_by_fy_sq = v_v0_by_fy*v_v0_by_fy;
        // Refer to compute3Dpositions.m from ICL-NUIM dataset
        float zVal = depthVal /
        sqrt(u_u0_by_fx_sq + v_v0_by_fy_sq + 1.0);
        
        //printf("%d %d %f\n", u, v, zVal);
        // Inverse depth
        //float invDepthVal = 1.0/depthVal;
        
        *ptrDepth = zVal;
        ptrDepth++;
        
        u++;
        
    }
    
    fclose(fid);
    
    return depthInput;
    
}
    
void DepthMapForRollingShutter::propagateDepth(Frame* new_keyframe)
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
    const bool* trackingWasGood =
        new_keyframe->getTrackingParent() == activeKeyFrame ?
        new_keyframe->refPixelWasGoodNoCreate() : 0;
    
    const float* activeKFImageData = activeKeyFrame->image(0);
    const float* newKFMaxGrad = new_keyframe->maxGradients(0);
    const float* newKFImageData = new_keyframe->image(0);
    
    // Go through all pixels (x, y) of OLD image, propagating forwards.
    for (int y=0;y<height;y++)
    {
        
        for (int x=0;x<width;x++)
        {
            
            DepthMapPixelHypothesis* source = currentDepthMap + x + y*width;
            
            if (source->isValid == true)
            {
                
                // Process progation for a pixel
                processPropagationForPixel(x, y,
                                           old_keyframe,
                                           new_keyframe,
                                           source,
                                           newKFMaxGrad,
                                           trackingWasGood,
                                           activeKFImageData,
                                           newKFImageData);
                
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
void DepthMapForRollingShutter::
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
        worldToOld = Sophus::SE3d(old_keyframe->getGTPoseAtRow(y));

    }
    else
    {
        // Est motion
        worldToOld = Sophus::SE3d(old_keyframe->getPoseAtRow(y));
    }

    Eigen::Vector4d new3DPoint;
    for (int rowNew = 0; rowNew < height; rowNew++)
    {
        
        // Motion of the old to new frame
        Sophus::SE3d worldToNew;
        if (useGTmotion_ == true)
        {
            // GT motion
            worldToNew = Sophus::SE3d(new_keyframe->getGTPoseAtRow(rowNew));

        }
        else
        {
            worldToNew = Sophus::SE3d(new_keyframe->getPoseAtRow(rowNew));
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

        // Check if the row of the 2D point is at the row i
        if (fabs(new2DPoint[1] - rowNew) < 0.5) // A half pixel
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
                
            }
        
        } // end if fabs
        
    } // end for rowNew
    
}
    
void DepthMapForRollingShutter::createKeyFrame(Frame* new_keyframe)
{
    assert(isValid());
    assert(new_keyframe != nullptr);
    assert(new_keyframe->hasTrackingParent());
    
    //boost::shared_lock<boost::shared_mutex> lock = activeKeyFrame->getActiveLock();
    boost::shared_lock<boost::shared_mutex> lock2 = new_keyframe->getActiveLock();
    
    struct timeval tv_start_all, tv_end_all;
    gettimeofday(&tv_start_all, NULL);
    
    
    resetCounters();
    
    if(plotStereoImages)
    {
        cv::Mat keyFrameImage(new_keyframe->height(), new_keyframe->width(), CV_32F, const_cast<float*>(new_keyframe->image(0)));
        keyFrameImage.convertTo(debugImageHypothesisPropagation, CV_8UC1);
        cv::cvtColor(debugImageHypothesisPropagation, debugImageHypothesisPropagation, CV_GRAY2RGB);
    }
    
    
    
    SE3 oldToNew_SE3 = se3FromSim3(new_keyframe->pose->thisToParent_raw).inverse();
    
    struct timeval tv_start, tv_end;
    gettimeofday(&tv_start, NULL);
    propagateDepth(new_keyframe);
    gettimeofday(&tv_end, NULL);
    msPropagate = 0.9*msPropagate + 0.1*((tv_end.tv_sec-tv_start.tv_sec)*1000.0f + (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
    nPropagate++;
    
    activeKeyFrame = new_keyframe;
    activeKeyFramelock = activeKeyFrame->getActiveLock();
    activeKeyFrame->setUndistorter(new_keyframe->undistorter_); // Set undistorter
    activeKeyFrame->setSpline(new_keyframe->spline_);
    activeKeyFrameImageData = new_keyframe->image(0);
    activeKeyFrameIsReactivated = false;
    
    
    
    gettimeofday(&tv_start, NULL);
    regularizeDepthMap(true, VAL_SUM_MIN_FOR_KEEP);
    gettimeofday(&tv_end, NULL);
    msRegularize = 0.9*msRegularize + 0.1*((tv_end.tv_sec-tv_start.tv_sec)*1000.0f + (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
    nRegularize++;
    
    
    gettimeofday(&tv_start, NULL);
    regularizeDepthMapFillHoles();
    gettimeofday(&tv_end, NULL);
    msFillHoles = 0.9*msFillHoles + 0.1*((tv_end.tv_sec-tv_start.tv_sec)*1000.0f + (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
    nFillHoles++;
    
    
    gettimeofday(&tv_start, NULL);
    regularizeDepthMap(false, VAL_SUM_MIN_FOR_KEEP);
    gettimeofday(&tv_end, NULL);
    msRegularize = 0.9*msRegularize + 0.1*((tv_end.tv_sec-tv_start.tv_sec)*1000.0f + (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
    nRegularize++;
    
    
    
#if 1// TURN OFF MAKING MEAN INVERSE DEPTH TO BE ONE (Jae-Hak Kim)
    activeKeyFrame->pose->thisToParent_raw = sim3FromSE3(oldToNew_SE3.inverse(), 1.0);
    activeKeyFrame->pose->invalidateCache();
#else // TURN ON MAKING MEAN INVERSE DEPTH TO BE ONE (Original)
    // make mean inverse depth be one.
    float sumIdepth=0, numIdepth=0;
    for(DepthMapPixelHypothesis* source = currentDepthMap; source < currentDepthMap+width*height; source++)
    {
        if(!source->isValid)
            continue;
        sumIdepth += source->idepth_smoothed;
        numIdepth++;
    }
    float rescaleFactor = numIdepth / sumIdepth;
    float rescaleFactor2 = rescaleFactor*rescaleFactor;
    for(DepthMapPixelHypothesis* source = currentDepthMap; source < currentDepthMap+width*height; source++)
    {
        if(!source->isValid)
            continue;
        source->idepth *= rescaleFactor;
        source->idepth_smoothed *= rescaleFactor;
        source->idepth_var *= rescaleFactor2;
        source->idepth_var_smoothed *= rescaleFactor2;
    }
    activeKeyFrame->pose->thisToParent_raw = sim3FromSE3(oldToNew_SE3.inverse(), rescaleFactor);
    activeKeyFrame->pose->invalidateCache();
#endif
    
    // Update depth in keyframe
    
    gettimeofday(&tv_start, NULL);
#if 0//TURN OFF DEPTH MAP UPDATE (Jae-Hak Kim)
#else // TURN ON DEPTH MAP UPDATE (Original)
    activeKeyFrame->setDepth(currentDepthMap);
#endif
    gettimeofday(&tv_end, NULL);
    msSetDepth = 0.9*msSetDepth + 0.1*((tv_end.tv_sec-tv_start.tv_sec)*1000.0f + (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
    nSetDepth++;
    
    gettimeofday(&tv_end_all, NULL);
    msCreate = 0.9*msCreate + 0.1*((tv_end_all.tv_sec-tv_start_all.tv_sec)*1000.0f + (tv_end_all.tv_usec-tv_start_all.tv_usec)/1000.0f);
    nCreate++;
    
    
    
    if(plotStereoImages)
    {
        //Util::displayImage( "KeyFramePropagation", debugImageHypothesisPropagation );
    }
    
}
    
//------------------------------------------------------------------------------
// Save the depth map pixel hypothesis as a text file with the following
// format of each row consists of x, y, z coordinates and the variance
// on the inverse depth (z).
//
//      filename: [string] output text filename
//
void DepthMapForRollingShutter::saveDepthAsFile(char *filename)
{

    FILE *fid = fopen(filename, "w");
    
    for (int y=3; y<height - 3; y++)
    {
        
        for (int x=3; x<width - 3; x++)
        {
            
            // Depth for the pixel x and y
            int idx = x + y*width;
            DepthMapPixelHypothesis* target = currentDepthMap + idx;
            
            // Save a valid depth
            if (target->isValid)
            {

//                // Set depth value from file via calibration information
//                float u_u0_by_fx = ((float)x - cx)/fx;
//                float v_v0_by_fy = ((float)y - cy)/fy;
//                float u_u0_by_fx_sq = u_u0_by_fx*u_u0_by_fx;
//                float v_v0_by_fy_sq = v_v0_by_fy*v_v0_by_fy;
//
//                float zVal = target->idepth /
//                    sqrt(u_u0_by_fx_sq + v_v0_by_fy_sq + 1.0);

                float zVal = 1.0/target->idepth;
                // Get x, y and z coordinates in the frame
                Eigen::Vector3f KinvP =
                    Eigen::Vector3f(fxi*x + cxi, fyi*y + cyi, 1.0f);
                float zcoord = zVal;
                float xcoord = KinvP[0]*zVal;
                float ycoord = KinvP[1]*zVal;
#if 0 // DEBUG
                // Save pixel coordinates with the inverse depth
                fprintf(fid, "xcoord=%f ycoord=%f zcoord=%f idepth_var=%f "
                             "idepth=%f fxi=%f fyi=%f cxi=%f cyi=%f"
                             "KinvP0=%f KinvP1=%f\n",
                        xcoord, ycoord, zcoord, target->idepth_var,
                        target->idepth, fxi, fyi, cxi, cyi,
                        KinvP[0], KinvP[1]);
#else
                // Save pixel coordinates with the inverse depth
                fprintf(fid, "%f %f %f %f\n",
                        xcoord, ycoord, zcoord, target->idepth_var);
#endif
            }
            
        }
        
    }
                
    fclose(fid);
    
}
    
void DepthMapForRollingShutter::regularizeDepthMap(bool removeOcclusions,
                                                   int validityTH)
{
    runningStats.num_reg_smeared=0;
    runningStats.num_reg_total=0;
    runningStats.num_reg_deleted_secondary=0;
    runningStats.num_reg_deleted_occluded=0;
    runningStats.num_reg_blacklisted=0;
    runningStats.num_reg_setBlacklisted=0;
    
    memcpy(otherDepthMap,
           currentDepthMap,
           width*height*sizeof(DepthMapPixelHypothesis));
    
    
    if(removeOcclusions)
        threadReducer.reduce(boost::bind(&DepthMapForRollingShutter::
                                         regularizeDepthMapRow<true>,
                                         this,
                                         validityTH,
                                         _1, _2, _3),
                             2, height-2, 10);
    else
        threadReducer.reduce(boost::bind(&DepthMapForRollingShutter::
                                         regularizeDepthMapRow<false>,
                                         this,
                                         validityTH,
                                         _1, _2, _3),
                             2, height-2, 10);
    
    
    if(enablePrintDebugInfo && printRegularizeStatistics)
        printf("REGULARIZE (%d): %d smeared; %d blacklisted /%d new); "
               "%d deleted; %d occluded; %d filled\n",
               activeKeyFrame->id(),
               runningStats.num_reg_smeared,
               runningStats.num_reg_blacklisted,
               runningStats.num_reg_setBlacklisted,
               runningStats.num_reg_deleted_secondary,
               runningStats.num_reg_deleted_occluded,
               runningStats.num_reg_created);
}

template<bool removeOcclusions>
void DepthMapForRollingShutter::regularizeDepthMapRow(int validityTH,
                                                      int yMin,
                                                      int yMax,
                                                      RunningStats* stats)
{
    const int regularize_radius = 2;
    
    const float regDistVar = REG_DIST_VAR;
    
    for(int y=yMin;y<yMax;y++)
    {
        for(int x=regularize_radius;x<width-regularize_radius;x++)
        {
            DepthMapPixelHypothesis* dest = currentDepthMap + x + y*width;
            DepthMapPixelHypothesis* destRead = otherDepthMap + x + y*width;
            
            // if isValid need to do better examination and then update.
            
            if(enablePrintDebugInfo && destRead->blacklisted < MIN_BLACKLIST)
                stats->num_reg_blacklisted++;
            
            if(!destRead->isValid)
                continue;
            
            float sum=0, val_sum=0, sumIvar=0;//, min_varObs = 1e20;
            int numOccluding = 0, numNotOccluding = 0;
            
            for(int dx=-regularize_radius; dx<=regularize_radius;dx++)
                for(int dy=-regularize_radius; dy<=regularize_radius;dy++)
                {
                    DepthMapPixelHypothesis* source = destRead + dx + dy*width;
                    
                    if(!source->isValid) continue;
                    //					stats->num_reg_total++;
                    
                    float diff =source->idepth - destRead->idepth;
                    if(DIFF_FAC_SMOOTHING*diff*diff >
                       source->idepth_var + destRead->idepth_var)
                    {
                        if(removeOcclusions)
                        {
                            if(source->idepth > destRead->idepth)
                                numOccluding++;
                        }
                        continue;
                    }
                    
                    val_sum += source->validity_counter;
                    
                    if(removeOcclusions)
                        numNotOccluding++;
                    
                    float distFac = (float)(dx*dx+dy*dy)*regDistVar;
                    float ivar = 1.0f/(source->idepth_var + distFac);
                    
                    sum += source->idepth * ivar;
                    sumIvar += ivar;
                    
                    
                }
            
            if(val_sum < validityTH)
            {
                dest->isValid = false;
                if(enablePrintDebugInfo) stats->num_reg_deleted_secondary++;
                dest->blacklisted--;
                
                if(enablePrintDebugInfo) stats->num_reg_setBlacklisted++;
                continue;
            }
            
            
            if(removeOcclusions)
            {
                if(numOccluding > numNotOccluding)
                {
                    dest->isValid = false;
                    if(enablePrintDebugInfo) stats->num_reg_deleted_occluded++;
                    
                    continue;
                }
            }
            
            sum = sum / sumIvar;
            sum = UNZERO(sum);
            
            
            // update!
            dest->idepth_smoothed = sum;
            dest->idepth_var_smoothed = 1.0f/sumIvar;
            
            if(enablePrintDebugInfo) stats->num_reg_smeared++;
        }
    }
}
template void DepthMapForRollingShutter::
    regularizeDepthMapRow<true>(int validityTH, int yMin,
                                int yMax, RunningStats* stats);
template void DepthMapForRollingShutter::
    regularizeDepthMapRow<false>(int validityTH, int yMin,
                                 int yMax, RunningStats* stats);


void DepthMapForRollingShutter::initializeFromGTDepthWithNoise(Frame* new_frame)
{
    assert(new_frame->hasIDepthBeenSet());
    
    activeKeyFramelock = new_frame->getActiveLock();
    activeKeyFrame = new_frame;
    activeKeyFrame->setUndistorter(new_frame->undistorter_);
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
    
    float averageGTIDepth = averageGTIDepthSum / averageGTIDepthNum;
    float varGTIDepth_sum = 0;
    for(int y=0;y<height;y++)
    {
        for(int x=0;x<width;x++)
        {
            float idepthValue = idepth[x+y*width];
            if(!isnanf(idepthValue) && idepthValue > 0)
            {
                varGTIDepth_sum = varGTIDepth_sum +
                    (idepthValue - averageGTIDepth)*
                    (idepthValue - averageGTIDepth);
            }
        }
    }
    float varGTIDepth = varGTIDepth_sum / averageGTIDepthNum;
    float stdDevGTIDepth = sqrt(varGTIDepth);
    
    
    for(int y=0;y<height;y++)
    {
        for(int x=0;x<width;x++)
        {
            float idepthValue = idepth[x+y*width];
            
            if(!isnan(idepthValue) && idepthValue > 0)
            {
                currentDepthMap[x+y*width] =
                    DepthMapPixelHypothesis(idepthValue,
                                            idepthValue,
                                            2.0*stdDevGTIDepth,
                                            2.0*stdDevGTIDepth,
                                            20);
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
