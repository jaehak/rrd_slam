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
//==============================================================================
// BN3
//
// B-spline with support of
// Neighbour 7 (Seven) segments by 4 Control Points
// Multiple keyframes
// Rolling shutter keyframe
// Version 3
//
// Author: Jae-Hak Kim <jaehak.kim@adelaide.edu.au>
// Copyright@2015, University of Adelaide
//==============================================================================
#include "Tracking/SE3TrackerRollingShutter.h"
#include <opencv2/highgui/highgui.hpp>
#include "DataStructures/Frame.h"
#include "Tracking/TrackingReference.h"
#include "util/globalFuncs.h"
#include "IOWrapper/ImageDisplay.h"
#include "Tracking/least_squares.h"
#include <math.h>

#include "ceres/ceres.h"
#include "Tracking/Spline.h"

#include "Tracking/WarpResidual_BN3.h"

#include <unsupported/Eigen/MatrixFunctions>

namespace lsd_slam
{
    
#define CERES_THREADS 8

SE3 SE3TrackerRollingShutter::trackFrame_Ceres_BN3(
        TrackingReference* reference,
        std::deque<std::shared_ptr<TrackingReference> > *referenceList,
        std::deque<std::shared_ptr<Frame> > *frameList, // for tracking a frame
        const SE3& frameToReference_initialEstimate,
        std::deque<int> startRowImageList,
        std::deque<int> startRowSegList,
        double rowPerSegment,
        double framePerSegment,
        double scale,
        std::vector<std::string> files)
{

    // Copy file list only once
    if (files_.size() == 0)
    {
        files_ = files;
    }

    printf("SPLINE_K = %d\n", SPLINE_K);

    // Set frame number
    frameNumberToDisplay = frameList->front()->id();

    for (int i=0; i<(int)frameList->size(); i++)
    {

        boost::shared_lock<boost::shared_mutex> lock =
                (*frameList)[i]->getActiveLock();

    }
    diverged = false;
    trackingWasGood = true;
    affineEstimation_a = 1; affineEstimation_b = 0;

    if (saveAllTrackingStages)
    {

        saveAllTrackingStages = false;
        saveAllTrackingStagesInternal = true;

    }

    if (plotTrackingIterationInfo)
    {

        const float* frameImage = frameList->back()->image();
        for (int row = 0; row < height; ++ row)
        {

            for (int col = 0; col < width; ++ col)
            {

                setPixelInCvMat(&debugImageSecondFrameRollingShutter,
                                getGrayCvPixel(frameImage[col+row*width]),
                                col, row, 1);
            }

        }

    }

    // ============ track frame ============
#if 1// USE inverse of SE3
    Sophus::SE3d referenceToFrame =
            frameToReference_initialEstimate.inverse();
#else // USE negative of Lie algebra elements for the inverse
      // This should be more stable in numerically
    Sophus::SE3d referenceToFrame =
            Sophus::SE3d::exp(
                -1.0*frameToReference_initialEstimate.log());
#endif


#if 1//PRINT REF
    Eigen::Matrix<float,6,1> tmp2
            = frameToReference_initialEstimate.log().cast<float>();
    std::cout << "[REF] frameToReference_initEst = "
              << tmp2.transpose() << std::endl;
    std::cout << "[REF] refererenceToFrame = ";
    Eigen::Matrix<float,6,1> tmp = referenceToFrame.log().cast<float>();
    std::cout << tmp.transpose() << std::endl;
#endif

    int numCalcResidualCalls[PYRAMID_LEVELS];
    int numCalcWarpUpdateCalls[PYRAMID_LEVELS];

    //--------------------------------------------------------------------------
    // Set all initial 5 control poses (WorldToControl)
    if (this->controlPoses.size() == 0)
    {

        for (int i=0; i<INIT_NUM_CTRLPTS; i++)
        {
#if 1 // SET ALL ZERO INITIALS
            Sophus::SE3d::Tangent zeroVec6;
            zeroVec6 << 0, 0, 0, 0, 0, 0;
            Sophus::SE3d origin = Sophus::SE3d::exp(zeroVec6);
            controlPoses.push_back(origin);
#else

#if 1 // SET GT CTRL PTS for starting frame 38
//            0                                   0                   0
//            0.181850000000000                   0                   0
//            0.363580000000000                   0                   0
//            0.545480000000000                   0                   0
//            0.727280000000000                   0                   0
//            0.907380000000000                   0  -0.001690000000000
//            1.069480000000000                   0  -0.021490000000000
//            1.189380000000000                   0  -0.083330000000000
//            1.251180000000000                   0  -0.203320000000000
//            1.270980000000000                   0  -0.365390000000000
//            1.272680000000000                   0  -0.545450000000000
            
            if (i==0)
            {
                Sophus::SE3d::Tangent zeroVec6;
                zeroVec6 << 0, 0, 0, 0, 0, 0; // t and omega
                Sophus::SE3d origin = Sophus::SE3d::exp(zeroVec6);
                controlPoses.push_back(origin);
            }
            if (i==1)
            {
                Sophus::SE3d::Tangent zeroVec6;
                zeroVec6 << 0, 0, 0, 0, 0, 0; // t and omega
                Sophus::SE3d origin = Sophus::SE3d::exp(zeroVec6);
                controlPoses.push_back(origin);
            }
            if (i==2)
            {
                float noise = 0.0f*(rand() % 100001)/100000.0f;
                Sophus::SE3d::Tangent zeroVec6;
                zeroVec6 << -0.18185 + noise, 0, 0, 0, 0, 0; // t and omega
                //zeroVec6 << -0.18185, 0, 0, 0, 0, 0; // t and omega
                Sophus::SE3d origin = Sophus::SE3d::exp(zeroVec6);
                controlPoses.push_back(origin);
            }
            if (i==3)
            {
                Sophus::SE3d::Tangent zeroVec6;
                //zeroVec6 << -0.36358, 0, 0, 0, 0, 0; // t and omega
                //float noise = 1.0f*(rand() % 100001)/100000.0f;
                float noise = 0.0f;
                zeroVec6 << -0.36358 + noise, 0, 0, 0, 0, 0; // t and omega
                Sophus::SE3d origin = Sophus::SE3d::exp(zeroVec6);
                controlPoses.push_back(origin);
            }
#endif

#if 0 // SET GT CTRL PTS for starting frame 45
//            %    0.0618         0    -0.1200
//            %    0.0816         0    -0.2821
//            %    0.0833         0    -0.4621
//            %    0.0833         0    -0.6427
//            %    0.0833         0    -0.8185
//            %    0.0833         0    -0.9857
//            %    0.0833         0    -1.1443
//            %    0.0833         0    -1.2973
//            %    0.0833         0    -1.4479
//            %    0.0833         0    -1.5980
//            %    0.0833         0    -1.7483
            if (i==0)
            {
                Sophus::SE3d::Tangent zeroVec6;
                zeroVec6 << 0, 0, 0, 0, 0, 0; // t and omega
                Sophus::SE3d origin = Sophus::SE3d::exp(zeroVec6);
                controlPoses.push_back(origin);
            }
            if (i==1)
            {
                Sophus::SE3d::Tangent zeroVec6;
                zeroVec6 << 0, 0, 0, 0, 0, 0; // t and omega
                Sophus::SE3d origin = Sophus::SE3d::exp(zeroVec6);
                controlPoses.push_back(origin);
            }
            if (i==2)
            {
                //float noise = 0.0f*(rand() % 100001)/100000.0f;
                Sophus::SE3d::Tangent zeroVec6;
                zeroVec6 << -0.0618, 0, 0.1200, 0, 0, 0; // t and omega
                Sophus::SE3d origin = Sophus::SE3d::exp(zeroVec6);
                controlPoses.push_back(origin);
            }
            if (i==3)
            {
                Sophus::SE3d::Tangent zeroVec6;
                //float noise = 1.0f*(rand() % 100001)/100000.0f;
                //float noise = 0.1f;
                zeroVec6 << -0.0816, 0, 0.2821, 0, 0, 0; // t and omega
                Sophus::SE3d origin = Sophus::SE3d::exp(zeroVec6);
                controlPoses.push_back(origin);
            }

#endif

#if 0 // SET GT CTRL PTS for starting frame 52
//            %          0         0         0
//            %          0         0   -0.1530
//            %          0         0   -0.3036
//            %          0         0   -0.4537
//            %          0         0   -0.6041
            if (i==0)
            {
                Sophus::SE3d::Tangent zeroVec6;
                zeroVec6 << 0, 0, 0, 0, 0, 0; // t and omega
                Sophus::SE3d origin = Sophus::SE3d::exp(zeroVec6);
                controlPoses.push_back(origin);
            }
            if (i==1)
            {
                Sophus::SE3d::Tangent zeroVec6;
                zeroVec6 << 0, 0, 0, 0, 0, 0; // t and omega
                Sophus::SE3d origin = Sophus::SE3d::exp(zeroVec6);
                controlPoses.push_back(origin);
            }
            if (i==2)
            {
                //float noise = 0.0f*(rand() % 100001)/100000.0f;
                Sophus::SE3d::Tangent zeroVec6;
                zeroVec6 << 0, 0, 0.1530, 0, 0, 0; // t and omega
                Sophus::SE3d origin = Sophus::SE3d::exp(zeroVec6);
                controlPoses.push_back(origin);
            }
            if (i==3)
            {
                Sophus::SE3d::Tangent zeroVec6;
                //float noise = 1.0f*(rand() % 100001)/100000.0f;
                //float noise = 0.1f;
                zeroVec6 << 0, 0, 0.3036, 0, 0, 0; // t and omega
                Sophus::SE3d origin = Sophus::SE3d::exp(zeroVec6);
                controlPoses.push_back(origin);
            }

#endif
#endif
        }

    }
    else
    {

#if 0// SEGHZ_EQ_CAMHZ
        // Add a new control pose
        printf("Adding a new control point\n");
        Sophus::SE3d copyLast = controlPoses.back();
        controlPoses.push_back(copyLast);
#endif
        
#if 0 // (ORIGINAL CONDITION ADDING A NEW CONTROL POINT)
        // Add a new control pose if startRow crosses rowPerSegment
        printf("fmod(%d, %f)=%f\n",
               startRowImageList.back(), rowPerSegment,
               fmod((double)startRowImageList.back(), rowPerSegment));
        if (fmod((double)startRowImageList.back(), rowPerSegment) == 0)
#else
        // NEW CONDITION ADDING A NEW CONTROL POINT
        // Add a new control point if any start row crosses rowPerSegment
        bool needNewControlPoint = false;
        for (int i=0;i<(int)frameList->size();i++)
        {
            
            // If any frame has an index larger than the given size
            // of the control points, we tick that new control point needed
            int idxSCP = (*frameList)[i]->indexSplineCtrlPoint();
            if (idxSCP + (SPLINE_K - 1) >= (int)controlPoses.size())
            {
                
                needNewControlPoint = true;
                printf("New control point needed at startRowImg %d;"
                       "idxSCP %d; frameList %d\n",
                       startRowImageList.back(),
                       idxSCP, i);
                
            }
            
        }
        if (needNewControlPoint == true)
            
#endif
        {

#if 1//COPY LAST
            
            Sophus::SE3d copyLast = controlPoses.back();
            controlPoses.push_back(copyLast);
            
#else // EXTRAPOLATE
            
            // Instead of copying the last, use an extrapolation to infer
            // the location of this new control point
            if (controlPoses.size() >= 3)
            {

                // Note: they are "World to Control"
                Sophus::SE3d cp2 = controlPoses.back();
                Sophus::SE3d cp1 = controlPoses[controlPoses.size() - 2];
                Sophus::SE3d cp0 = controlPoses[controlPoses.size() - 3];
                Sophus::SE3d cp1ToCp2 = cp2*cp1.inverse();
                Sophus::SE3d cp0ToCp1 = cp1*cp0.inverse();
                //Sophus::SE3d newCP = cp2*cp1ToCp2;
                //Sophus::SE3d cp1toCp2 = cp2.inverse()*cp1;
                //Sophus::SE3d newCP = cp1ToCp2*cp2; // Correct one
                Sophus::SE3d newCP = cp0ToCp1*cp2; // More stable one
                controlPoses.push_back(newCP);

            }
            else
            {

                Sophus::SE3d copyLast = controlPoses.back();
                controlPoses.push_back(copyLast);

            }
            
#endif

        }

    }

    // Run optimization to refine the 4 control poses
    bool isOptimized = optimizeUsingCeres_BN3(reference,
                                             referenceList,
                                             frameList,
                                             &referenceToFrame,
                                             &controlPoses, // WorldToControl
                                             startRowImageList,
                                             startRowSegList,
                                             rowPerSegment,
                                             scale);

    if (isOptimized == false)
    {

        diverged = true;
        trackingWasGood = false;

        printf("Optimization failed.. DIVERGED.\n");

        return SE3();

    }

    //--------------------------------------------------------------------------
    // Print control poses
    for (unsigned int i=0; i<controlPoses.size(); i++)
    {

        printf("CTRLPOSE_Frame%08d_Ctrl%04d ", frameList->front()->id(), i);

        Eigen::Matrix<double,4,4> M = controlPoses[i].matrix();
        for (int j=0;j<4;j++)
        {
            for (int k=0;k<4;k++)
            {

                printf("%f ", M(j,k));

            }
        }
        printf("\n");

    }

#ifdef DEBUG
    //--------------------------------------------------------------------------
    // Print segments by the control poses

    Spline<double> spline(controlPoses, SPLINE_K);
    for (unsigned int i=0; i<controlPoses.size() - 3; i++)
    {

        for (float u=0.0; u<1.0; u=u+0.01)
        {

            printf("POSE_Frame%08d_Ctrl%04d ", frameList->front()->id(), i);

            Eigen::Matrix<double,4,4> M = spline.getPoseOnSpline(u, i);
            for (int j=0;j<4;j++)
            {
                for (int k=0;k<4;k++)
                {

                    printf("%f ", M(j,k));

                }
            }
            printf("\n");


        }

    }
    
#endif

    //--------------------------------------------------------------------------

    if (plotTracking)
    {

        Util::displayImage("TrackingResidual",
                           debugImageResidualsRollingShutter, false);

        if (isWindowMoved2 == false)
        {

            isWindowMoved2 = true;
            cv::moveWindow("TrackingResidual",  650*2 + 70, 520*2 + 30);

        }

    }

    if (enablePrintDebugInfo && printTrackingIterationInfo)
    {

        printf("Tracking: ");
        for(int lvl=PYRAMID_LEVELS-1;lvl >= 0;lvl--)
        {
            printf("lvl %d: %d (%d); ",
                   lvl,
                   numCalcResidualCalls[lvl],
                   numCalcWarpUpdateCalls[lvl]);
        }

        printf("\n");

    }

    saveAllTrackingStagesInternal = true;

    //lastResidual = last_residual;
    float lastResidual = (float)last_residual_;

    trackingWasGood = !diverged &&
            lastGoodCount / (frameList->front()->width(SE3TRACKING_MIN_LEVEL) *
                             frameList->front()->height(SE3TRACKING_MIN_LEVEL)
                             *lsd_slam::pixelDensityForTracking)
            > MIN_GOODPERALL_PIXEL &&
            lastGoodCount / (lastGoodCount + lastBadCount) >
            MIN_GOODPERGOODBAD_PIXEL;

    if (!trackingWasGood)
    {

        printf("Tracking was NOT good\n");
        printf("    lastGoodCount = %f\n", lastGoodCount);
        printf("    (frameList.front()->width(SE3TRACKING_MIN_LEVEL) * "
               "frameList.front()->height(SE3TRACKING_MIN_LEVEL)) = %d\n",
                frameList->front()->width(SE3TRACKING_MIN_LEVEL) *
                frameList->front()->height(SE3TRACKING_MIN_LEVEL));
        printf("    MIN_GOODPERALL_PIXEL = %f\n", MIN_GOODPERALL_PIXEL);
        printf("    (lastGoodCount + lastBadCount) = %f\n",
               (lastGoodCount + lastBadCount));
        printf("    MIN_GOODPERGOODBAD_PIXEL = %f\n",
               MIN_GOODPERGOODBAD_PIXEL);

    }

#if 0 // USE NUM FRAMES TRACKED on BACK
    if (trackingWasGood)
    {

        referenceList->back()->keyframe->numFramesTrackedOnThis++;
        reference->keyframe->numFramesTrackedOnThis++;

    }
#else
    // USE NUM FRAMES TRACKED on FRONT
    if (trackingWasGood)
    {
        
        //double framePerSegment = rowPerSegment / (height + delayRows);
        int numToPop = (int)framePerSegment; // Number of pop from the list
        
        for (int i=0;i<numToPop;i++)
        {
        
            (*referenceList)[i].get()->keyframe->numFramesTrackedOnThis++;
        
        }
        reference->keyframe->numFramesTrackedOnThis++;
        
    }
#endif
    
#if 0 // SET INFO and SPLINE for ALL FRAME LIST
    
    for (int i=0;i<(int)frameList->size();i++)
    {
        
        // Set newTrackingFrame information
        (*frameList)[i].get()->initialTrackedResidual = lastResidual / pointUsage;

        // Debug info
        this->debugInfo_lastResidual = lastResidual;
        this->debugInfo_pointUsage = pointUsage;
        this->debugInfo_initialTrackedResidual = (*frameList)[i].get()->initialTrackedResidual;

        (*frameList)[i].get()->pose->thisToParent_raw =
            sim3FromSE3(toSophus(referenceToFrame.inverse()),1);
        (*frameList)[i].get()->pose->trackingParent =
            (*referenceList)[i].get()->keyframe->pose;

        std::cout << "referenceToFrame = " << std::endl;
        std::cout << referenceToFrame.matrix() << std::endl;

        // Set spline pointer to the newTrackingFrame for rolling shutter image
        std::shared_ptr<Spline<double> > newSpline
        (new Spline<double>(controlPoses, SPLINE_K));
        
        // Set new spline
        (*frameList)[i].get()->setSpline(newSpline);
        (*frameList)[i].get()->setUndistorter(this->undistorter);
        
    }
#else
    
    #if 1// SET INFO and SPLINE for BACK of FRAME LIST (FIRST IMPLEMENTATION)
    
        #if 0 // NEW SPLINE IN EVERY UPDATE
    
    // Set spline pointer to the newTrackingFrame for rolling shutter image
    std::shared_ptr<Spline<double> > newSpline
    (new Spline<double>(controlPoses, SPLINE_K));
    
        #else // UPDATE INTERNAL SPLINE

    internalSpline->updateControlPoses(controlPoses);
    
        #endif
    
    
    //double framePerSegment = rowPerSegment / height;
    int numToPop = (int)framePerSegment; // Number of pop from the list
    
    // Set newTrackingFrame information to the last back blocks
    for (int i=numToPop; i<(int)frameList->size(); i++)
    {
        
        (*frameList)[i].get()->initialTrackedResidual
            = lastResidual / pointUsage;
        (*frameList)[i].get()->pose->thisToParent_raw
            = sim3FromSE3(toSophus(referenceToFrame.inverse()),1);
        (*frameList)[i].get()->pose->trackingParent
            = (*referenceList)[i].get()->keyframe->pose;

        std::cout << "referenceToFrame = " << std::endl;
        std::cout << referenceToFrame.matrix() << std::endl;

        #if 0 // NEW SPLINE AT EVERY UPDATE
        
        (*frameList)[i].get()->setSpline(newSpline);
        (*frameList)[i].get()->setUndistorter(this->undistorter);
        
        #else // UPDATE INTERNAL SPLINE
        
        (*frameList)[i].get()->setSpline(internalSpline);
        (*frameList)[i].get()->setUndistorter(this->undistorter);
        
        #endif
        
    }

    #else //// SET INFO and SPLINE for FRONT of FRAME LIST

    // Set newTrackingFrame information
    frameList->front()->initialTrackedResidual = lastResidual / pointUsage;
    
    // Debug info
    this->debugInfo_lastResidual = lastResidual;
    this->debugInfo_pointUsage = pointUsage;
    this->debugInfo_initialTrackedResidual = frameList->front()->initialTrackedResidual;
    
    frameList->front()->pose->thisToParent_raw =
        sim3FromSE3(toSophus(referenceToFrame.inverse()),1);
    frameList->front()->pose->trackingParent =
        referenceList->front()->keyframe->pose;
    
    std::cout << "referenceToFrame = " << std::endl;
    std::cout << referenceToFrame.matrix() << std::endl;
    
    // Set spline pointer to the newTrackingFrame for rolling shutter image
    std::shared_ptr<Spline<double> > newSpline
    (new Spline<double>(controlPoses, SPLINE_K));
    
    // Set only the front
    frameList->front()->setSpline(newSpline);
    frameList->front()->setUndistorter(this->undistorter);
    
        #if 0
    // Set all with new spline estimate
    for (std::deque<std::shared_ptr<Frame> >::iterator
         iter = frameList->begin();
         iter != frameList->end();
         iter++)

             
    {
        
        Frame *frame = iter->get();
        frame->setSpline(newSpline);
        frame->setUndistorter(this->undistorter);
        
    }
        #endif
    
    #endif
    
#endif

    return toSophus(referenceToFrame.inverse());

}

//------------------------------------------------------------------------------
// Given reference/key (source) image and new (target) segment image,
// it finds control poses describing the trajectory as B-Spline
//------------------------------------------------------------------------------
// Input:
//           reference : [TrackingReference] reference key as tracking source
//               frame : [Frame] input image frame
//    referenceToFrame : [nx1 SE3d] a vector of control poses SE3d
//
// Output:
//
//    referenceToFrame : [nx1 SE3d] a vector of control poses SE3d
//
//------------------------------------------------------------------------------
bool
SE3TrackerRollingShutter::optimizeUsingCeres_BN3(
        TrackingReference *reference,
        std::deque<std::shared_ptr<TrackingReference> > *referenceList,
        std::deque<std::shared_ptr<Frame> > *frameList,
        Sophus::SE3d *referenceToFrame,
        std::vector<Sophus::SE3d> *controlPoses, // WorldToControl
        std::deque<int> startRowImageList,
        std::deque<int> startRowSegList,
        double rowPerSegment,
        double scale)
{

    // Check the number of control points to be estimated
    const int numCtrlPtsToBeEst = SPLINE_K;
    if ((int)controlPoses->size() < numCtrlPtsToBeEst)
    {

        printf("Error: Number of control points to be estimated must be"
               " larger than and equal to %d.\n", numCtrlPtsToBeEst);
        exit(1);

    }

    //--------------------------------------------------------------------------
    // Set the initial control poses into an array
    double *pose = (double *) malloc(sizeof(double)*6*numCtrlPtsToBeEst); // 24

    //--------------------------------------------------------------------------
    // Index for spline control points
    double uaccRow = startRowImageList.back() / rowPerSegment;
#if 1//ONE-START
    int idxSplineCtrlPts = (int)floor(uaccRow) + NUM_EXTRA_CTRLPTS; // One-start
#else // ZERO-START

    int idxSplineCtrlPts = (int)floor(uaccRow); // Zero-start

#endif



    printf("idxSplineCtrlPts = %d\n", idxSplineCtrlPts);

    //--------------------------------------------------------------------------
    // Set uppper and lower bound of u time values which depend on start row
    double lowerBound_uValue = (double)
                               (0 // first row
                                + startRowImageList.back()
                                - startRowSegList.back())
                               / rowPerSegment;
    double upperBound_uValue =
    (double) (frameList->back()->height_rollingShutter_inputImage  // last row
                                + startRowImageList.back()
                                - startRowSegList.back())
                               / rowPerSegment;
    printf("Lower uValue %f, Upper uValue %f\n",
           lowerBound_uValue, upperBound_uValue);

#if 0 // USE WORLD TO CONTROL
    
    //--------------------------------------------------------------------------
    // Copy k control points to pose (World To Control)
    for (int i=idxSplineCtrlPts, k = 0;
         i<(int)controlPoses->size(); i++, k++)
    {

        for (int j=0; j<6; j++)
        {
            // World To Control
            pose[6*k + j] = (*controlPoses)[i].log()[j];

        }

    }

#else // USE CONTROL TO FRONT KEYFRAME
    
    //--------------------------------------------------------------------------
    // Copy k control points to pose (Control to frontkeyframe)

    // Index of control point for front keyframe
    int idxSC_frontKF
            = (int)floor(referenceList->front().get()->keyframe->startRowImage
            / referenceList->front().get()->keyframe->rowPerSegment)
            + NUM_EXTRA_CTRLPTS;
    
    // World to Frontkeyframe
    Sophus::SE3d worldToFrontkeyframe = (*controlPoses)[idxSC_frontKF];
    
#ifdef DEBUG
    std::cout << "worldToFrontkeyframe " << worldToFrontkeyframe.matrix() << std::endl;
    std::cout << worldToFrontkeyframe.log() << std::endl;
#endif
    
    // Assingn Control to Frontkeyframe to pose
    for (int i=idxSplineCtrlPts, k = 0;
         i<(int)controlPoses->size(); i++, k++)
    {
        
        // i-th Control to Frontkeyframe
        Sophus::SE3d worldToControl = (*controlPoses)[i];
        Sophus::SE3d controlToWorld = worldToControl.inverse();
        Sophus::SE3d controlToFrontkeyframe
            = worldToFrontkeyframe * controlToWorld;
        
#ifdef DEBUG
        std::cout << "worldToControl:" << worldToControl.log() << std::endl;
        std::cout << "controlToWorld:" << controlToWorld.log() << std::endl;
#endif
        
        for (int j=0; j<6; j++)
        {
            // Copy the control to frontkeyframe to pose
            
#ifdef DEBUG
            printf("pose[%d](%f) <= %f\n",
                   6*k + j,
                   pose[6*k + j],
                   controlToFrontkeyframe.log()[j]);
#endif
            pose[6*k + j] = controlToFrontkeyframe.log()[j];
            
        }
        
    }
    
#endif
    
    
    
   // WORK BACK (CODE Is broken)
//    //--------------------------------------------------------------------------
//    // Copy k control points to pose (First row keyframe to Control)
//    for (int i=idxSplineCtrlPts, k = 0;
//         i<(int)controlPoses->size(); i++, k++)
//    {
//        
//        // Transform as FirstRowKeyFrame to control
//        // RefToControl = WorldToControl*FirstRowKeyFrameToWorld
//        (*controlPoses)[i].log() * firstRowKeyFrameToWorld;
//        for (int j=0; j<6; j++)
//        {
//            // World To Control
//            pose[6*k + j] = (*controlPoses)[i].log()[j];
//            
//        }
//        
//    }

    

    // Print initial pose
    for (int i=0; i<(int)controlPoses->size(); i++)
    {

        printf("CTRL Initial control pose to frame %d:  ", i);
        //std::cout << (*controlPoses)[i].log().transpose() << std::endl;
        for (int j=0; j<6; j++)
            printf("%.6f ", (*controlPoses)[i].log()[j]);
        printf("\n");

    }

    //--------------------------------------------------------------------------
    // For each pyramid level, image warping is computed. The the result pose
    // is estimated and reused for the next pyramid level
    bool isSolved = false;
    for (int pyramidLevel = SE3TRACKING_MAX_LEVEL - 1;
         //int pyramidLevel = 4;
         //pyramidLevel >= 4;
         //pyramidLevel >= 3;
         //pyramidLevel >= 2;
         pyramidLevel >= SE3TRACKING_MIN_LEVEL;
         //pyramidLevel >= 0;
         pyramidLevel--)
    {

        // Compute pose for image warping for the given pyramid level
        isSolved = optimizeUsingCeres_BN3_PyramidLevelAt(
                    reference,
                    pyramidLevel,
#if 0 // USE WORLD TO CONTROL (Original)
                    pose, // WorldToControl (OUTPUT)
#else // USE CONTROL TO FRONT KEYFRAME
                    pose, // Control to front keyframe (OUTPUT)
#endif
                    referenceList,
                    frameList, // A queue of Frame objects
                    numCtrlPtsToBeEst,
                    startRowImageList, // A queue of starting rows
                    startRowSegList,
                    rowPerSegment,
                    lowerBound_uValue,
                    upperBound_uValue,
                    *controlPoses, // Pass-by-value (no change in this function)
                    idxSplineCtrlPts,
                    scale);

        //----------------------------------------------------------------------
        // Store statistics
        // Total number of warped pixel points from source to target image
        int totalWarpedPixels =
                WarpResidual_BN3::lastGoodCount_ +
                WarpResidual_BN3::lastBadCount_;

        // Store point usage
#if 1 // ORIGINAL
        pointUsage = WarpResidual_BN3::usageCount_ /
                (float) referenceList->front()->numData[pyramidLevel];

        this->debugInfo_pointUsage = pointUsage;
        this->debugInfo_usageCount = WarpResidual_BN3::usageCount_;
        this->debugInfo_numData =  (float) referenceList->front()->numData[pyramidLevel];

#else // NEW usage Count
        pointUsage = WarpResidual_BN3::usageCount_ /
        (float) referenceList->front()->keyframe->permaRefNumPts;
        
        this->debugInfo_pointUsage = pointUsage;
        this->debugInfo_usageCount = WarpResidual_BN3::usageCount_;
        this->debugInfo_permaRefNumPts = (float)referenceList->front()->keyframe->permaRefNumPts;
#endif
        // Check tracking valid
        if (totalWarpedPixels <
                MIN_GOODPERALL_PIXEL_ABSMIN *
                (frameList->front()->width(pyramidLevel)*
                 (frameList->front()->height(pyramidLevel)*
                  lsd_slam::pixelDensityForTracking)))
        {

            // Too less number of pixels is warped
            diverged = true;
            trackingWasGood = false;

            printf("\nDIVERGED -- Set divered = TRUE\n\n");

            std::cerr << "DIVERGED -- Set diverged = TRUE" << std::endl;
            exit(1);

            isSolved = false;
            return isSolved;

        }

        //----------------------------------------------------------------------
        // Output pose is being store into controlPoses
        //----------------------------------------------------------------------
        //----------------------------------------------------------------------
        // Copy back the result
#if 0// COPY ONLY ONE
        for (int i=0; i<1; i++)
#else // COPY ALL
        
        for (int i=0; i<SPLINE_K; i++)
#endif
        {

#if 0 // USE WORLD TO CONTROL (Original)

            Eigen::Matrix<double, 6, 1> poseOut;
            for (int j=0; j<6; j++)
            {

                poseOut[j] = pose[6*i + j]; // World to Control
                
            }
            Sophus::SE3d worldToControl = Sophus::SE3d::exp(poseOut);
            (*controlPoses)[idxSplineCtrlPts + i] = worldToControl;
            
#else // USE CONTROL TO FRONT KEYFRAME
           
            Eigen::Matrix<double, 6, 1> poseOut;
            for (int j=0; j<6; j++)
            {
                
                poseOut[j] = pose[6*i + j]; // Control to Front Keyframe
                
            }
            Sophus::SE3d controlToFrontKeyframe = Sophus::SE3d::exp(poseOut);
            
            // World to front keyframe
            Eigen::Matrix<double,6,1> worldToFrontKeyframe_log;
            for (int j=0; j<6; j++)
            {
                
                // ControlPoses_[idxSC_frontKF] is a control point related to
                // the front keyframe
                worldToFrontKeyframe_log[j]
                = (*controlPoses)[idxSC_frontKF].log()[j];
                
            }
            Sophus::SE3d worldToFrontKeyframe
                = Sophus::SE3d::exp(worldToFrontKeyframe_log);
            
            // World to Control
            Sophus::SE3d worldToControl
                = controlToFrontKeyframe.inverse() * worldToFrontKeyframe;
            (*controlPoses)[idxSplineCtrlPts + i] = worldToControl;
            
#endif
            
#if 1//PRINT OUT

            printf("\n\nCTRL_LEVEL%d [Segment %d; Control %d] "
                   "Final control to frame = "
                   "%f %f %f %f %f %f ===> ",
                   pyramidLevel,
                   frameList->back()->id(), i,
                   poseOut[0], poseOut[1], poseOut[2],
                   poseOut[3], poseOut[4], poseOut[5]);

#if 1// WORLD_TO_CONTROL
            std::cout << "CTRL_LEVEL" << pyramidLevel
                      << " World to Control " << i << " " <<
                worldToControl.log().transpose() << std::endl;
#endif

#endif

        }
        printf("\n");

    } // end of for - pyramid

/*
 *
 *
 * SET refrenceToFrame from controlPoses (worldToControl)
 *
 *
 *
 */

#if  0// ASSIGN MID ROW
    Spline<double> spline(*controlPoses, SPLINE_K);
    Eigen::Matrix4d poseMidRow =
            spline.getPoseOnSpline((upperBound_uValue + lowerBound_uValue)/2.0,
            (int)idxSplineCtrlPts);

    Sophus::SE3d worldToPoseMidRow = Sophus::SE3d(poseMidRow);
    Sophus::SE3d refToWorld =
        se3FromSim3(reference->keyframe->pose->getCamToWorld());
    // WorldToPosLastRow @ ReferenceToWorld = ReferenceToPoseLastRow
    *referenceToFrame = worldToPoseMidRow * refToWorld;
#endif
#if 0//ASSIGN_LAST_ROW_FOR_REFERENCE_OUTPUT

    Spline<double> spline(*controlPoses, SPLINE_K);
    Eigen::Matrix4d poseLastRow =
            spline.getPoseOnSpline(upperBound_uValue, (int)idxSplineCtrlPts);

    Sophus::SE3d worldToPoseLastRow = Sophus::SE3d(poseLastRow);
    Sophus::SE3d refToWorld =
        se3FromSim3(reference->keyframe->pose->getCamToWorld());
    // WorldToPosLastRow @ ReferenceToWorld = ReferenceToPoseLastRow
    *referenceToFrame = worldToPoseLastRow * refToWorld;

#endif
#if 1// ASSIGN_FIRST_ROW_FOR_REFERENCE_OUTPUT

    Spline<double> spline(*controlPoses, SPLINE_K);
    Eigen::Matrix4d poseFirstRow =
            spline.getPoseOnSpline(lowerBound_uValue, (int)idxSplineCtrlPts);

    // Pose for the first row in the image will be assigned for the reference
    // to frame
    Sophus::SE3d worldToPoseFirstRow = Sophus::SE3d(poseFirstRow);
    Sophus::SE3d refToWorld =
        se3FromSim3(referenceList->front()->keyframe->pose->getCamToWorld());
    // WorldToPosLastRow @ ReferenceToWorld = ReferenceToPoseFirstRow
#if 0 // TEST new ref to frame
    *referenceToFrame = worldToPoseFirstRow * refToWorld.inverse();
#else
    *referenceToFrame = worldToPoseFirstRow * refToWorld;
#endif

#endif

    lastGoodCount = WarpResidual_BN3::lastGoodCount_;
    lastBadCount = WarpResidual_BN3::lastBadCount_;

    delete pose;

    return isSolved;

}

//------------------------------------------------------------------------------
// optimizeUsingCeres_BN3_PyramidLevelAt()
//------------------------------------------------------------------------------
//          pose : An array of 4 control points in SE(3) (24 elements = 4*6)
//     frameList : A list of frames. Its index is shared with the startRowList.
//                 Each element is a Frame object which is an image as a part
//                 of poses on a segment. (Index shared with other List)
//  startRowList : A list of starting row. Its index is shared with
//                 the frameList, i.e. the same number of elements as the number
//                 of frames. Each element is a row index accumulated from the
//                 starting time stamp. (Index shared with other List)
// referenceList : A list of tracking references. Its index is shared with
//                 frameList, startRowImageList and startRowSegList.
//
//------------------------------------------------------------------------------
bool
SE3TrackerRollingShutter::optimizeUsingCeres_BN3_PyramidLevelAt(
        TrackingReference *reference,
        int level,
#if 0 // USE WORLD TO CONTROL pose (ORIGINAL)
        double *pose, // World To Control
#else // USE REF TO CONTROL pose
                                                                
        // An array of k control points as 6-vector in se(3)
        // representing transformations
        // from a control point for the current tracking frame
        // to control point for the front keyframe
        double *pose_of_kCtrlPts_to_frontKeyframe, // Control to FrontKeyframe
                                                                
                                                                
#endif
        std::deque<std::shared_ptr<TrackingReference> > *referenceList,
        std::deque<std::shared_ptr<Frame> > *frameList, // A list of frames
        const int numCtrlPtsToBeEst,
        std::deque<int> startRowImageList, // A list of starting rows
        std::deque<int> startRowSegList,
        double rowPerSegment,
        double lowerBound_uValue,
        double upperBound_uValue,
        const std::vector<Sophus::SE3d> controlPoses, // WorldToControl
        int idxSplineCtrlPts,
        double scale)
{

    // Check number of control poses
    if (numCtrlPtsToBeEst != SPLINE_K)
    {

        printf("Error: Number of control points to be estimated "
               "should be %d.\n", SPLINE_K);
        exit(1);

    }

    // Debug start
    calcResidualAndBuffers_debugStart();

    if (lsd_slam::plotTracking)
    {

        debugImageResidualsRollingShutter.setTo(0);

    }

    for (int i=0;i<(int)referenceList->size();i++)
    {

        (*referenceList)[i]->makePointCloud(level);

    }

    // Check pyramid gradients
    const Eigen::Vector4f* frame_gradients =
            frameList->back()->gradients(level);

    // Gradient X image
    if (lsd_slam::plotTracking)
    {

        cv::Mat *gradImageX = &debugGradientXImage0_RollingShutter;
        cv::Mat *gradImageY = &debugGradientXImage1_RollingShutter;

        gradImageX->setTo(0);
        gradImageY->setTo(0);

        for (int x = 0; x<frameList->back()->width(level); x++)
        {

            for (int y = 0; y<frameList->back()->height(level); y++)
            {

                // gradVal[0]: gradX
                // gradVal[1]: gradY
                // gradVal[2]: intesity
                const Eigen::Vector4f* gradVal =
                        frame_gradients + x +
                        y*frameList->back()->width(level);

                setPixelInCvMat(gradImageX,
                        getGrayCvPixel((float)(*gradVal)[0] + 100),
                        x, y,
                        frameList->back()->width(0)/
                        frameList->back()->width(level));
                setPixelInCvMat(gradImageY,
                        getGrayCvPixel((float)(*gradVal)[1] + 100),
                        x, y,
                        frameList->back()->width(0)/
                        frameList->back()->width(level));

            }

        }

    }

    // Display tracking referene
    if (lsd_slam::plotTracking)
    {

        if (referenceList->size() >= 1)
        {
        cv::Mat *refImage0 = &debugRefImage0;
        refImage0->setTo(0);
        int width = (*referenceList)[0]->keyframe->width(level);
        int height = (*referenceList)[0]->keyframe->height(level);
        for (int x=0;x<width;x++)
        {
            for (int y=0;y<height;y++)
            {

                // gradVal[0]: gradX
                // gradVal[1]: gradY
                // gradVal[2]: intesity
                const Eigen::Vector4f* gradVal =
                        (*referenceList)[0]->keyframe->gradients(level) +
                        x + y*width;

                setPixelInCvMat(refImage0,
                        getGrayCvPixel((float)(*gradVal)[1] + 100),
                        x, y,
                        (*referenceList)[0]->keyframe->width(0)/
                        width);

            }
        }
        }

        if (referenceList->size() >= 2)
        {
        cv::Mat *refImage1 = &debugRefImage1;
        refImage1->setTo(0);
        int width = (*referenceList)[1]->keyframe->width(level);
        int height = (*referenceList)[1]->keyframe->height(level);
        for (int x=0;x<width;x++)
        {
            for (int y=0;y<height;y++)
            {

                // gradVal[0]: gradX
                // gradVal[1]: gradY
                // gradVal[2]: intesity
                const Eigen::Vector4f* gradVal =
                        (*referenceList)[1]->keyframe->gradients(level) +
                        x + y*width;

                setPixelInCvMat(refImage1,
                        getGrayCvPixel((float)(*gradVal)[1] + 100),
                        x, y,
                        (*referenceList)[1]->keyframe->width(0)/
                        width);

            }
        }
        }

        if (referenceList->size() >= 3)
        {
        cv::Mat *refImage2 = &debugRefImage2;
        refImage2->setTo(0);
        int width = (*referenceList)[2]->keyframe->width(level);
        int height = (*referenceList)[2]->keyframe->height(level);
        for (int x=0;x<width;x++)
        {
            for (int y=0;y<height;y++)
            {

                // gradVal[0]: gradX
                // gradVal[1]: gradY
                // gradVal[2]: intesity
                const Eigen::Vector4f* gradVal =
                        (*referenceList)[2]->keyframe->gradients(level) +
                        x + y*width;

                setPixelInCvMat(refImage2,
                        getGrayCvPixel((float)(*gradVal)[1] + 100),
                        x, y,
                        (*referenceList)[2]->keyframe->width(0)/
                        width);

            }
        }
        }
    }

#if 0// Sophus TEST
    printf("Sophus test\n");
    Sophus::SE3Group<double>::Tangent s;
    s << -1, -2, -3, 1, 1, 0;
    Eigen::Matrix<double,4,4> M = Sophus::SE3Group<double>::exp(s).matrix();
    std::cout << "M=" << std::endl;
    std::cout << M << std::endl;
    Eigen::Matrix<double,4,4> M2 = Sophus::SE3Group<double>::hat(s).exp();
    std::cout << "M2=" << std::endl;
    std::cout << M2 << std::endl;
    Sophus::SE3Group<double>::Tangent s2 =
            Sophus::SE3Group<double>::exp(s).log();
    std::cout << "s2=" << std::endl;
    std::cout << s2 << std::endl;
#endif

    //--------------------------------------------------------------------------
    // Set the optimization using Ceres-solver
    //--------------------------------------------------------------------------
    Problem problem;

    // For each pixel in the reference frame (keyframe), a residual block
    // computing the error of image warping is computed
    // NOTE: the image is from pyramid, so their size varies on each level

    // Pointer to first semi-dense point for the given pyramid level
    // It consists of (x, y, z) for x-coordinate, y-coordinate and
    // z inverse depth. So, the refPoint is a 3D point in space with respect
    // to the reference camera frame.
    Eigen::Vector3f* refPoint;

    // Pointer to last semi-dense point in the reference image
    Eigen::Vector3f* refPoint_max;

    // Pointer to first semi-dense point containing colour and variance
    // information in the given pyramid level
    Eigen::Vector2f* refColourVariance;

    // Allocate memory for ycoordinate for warped points
    //      The size depends on the level and the number of source images
    //const int numSemiPtsInLevel = referenceList->back()->numData[level];
    //double *qycoord  = new double[numSemiPtsInLevel];

    TrackingReference *referenceAt;
    
    // Index of control point for front keyframe
    int idxSC_frontKF
        = (int)floor(referenceList->front().get()->keyframe->startRowImage
                 / referenceList->front().get()->keyframe->rowPerSegment)
                 + NUM_EXTRA_CTRLPTS;
    
    //--------------------------------------------------------------------------
    // Prepare qycoordList
    std::vector<std::shared_ptr<double> >qycoordList;
    for (int i=0;i<(int)referenceList->size();i++)
    {

        int numSemiPtsInLevel = (*referenceList)[i]->numData[level];
        std::shared_ptr<double> qycoordPtr(new double[numSemiPtsInLevel]);
        qycoordList.push_back(std::move(qycoordPtr));
        printf("Dim of qycoordList[%d] = %d;", i, numSemiPtsInLevel);
        printf("Ref.frameID = %d for Curr.frameID %d;",
               (*referenceList)[i]->frameID,
               (*frameList)[i].get()->id());
        printf("Ref.pose: [");
        std::cout <<
        se3FromSim3((*referenceList)[i]->keyframe
                    ->getScaledCamToWorld()).log().transpose();
        printf("] Ref.parent: [");
        if ((*referenceList)[i]->keyframe->hasTrackingParent())
        {

            std::cout <<
            se3FromSim3((*referenceList)[i]->keyframe
                                     ->getTrackingParent()
                                     ->getScaledCamToWorld()).log().transpose();

        }
        printf("] ");
        printf("Curr.pose: [");
        std::cout <<
        se3FromSim3((*frameList)[i]
                    ->getScaledCamToWorld()).log().transpose();
        printf("] Curr.parent: [");
        if ((*frameList)[i]->hasTrackingParent())
        {

            std::cout <<
            se3FromSim3((*frameList)[i]
                        ->getTrackingParent()
                        ->getScaledCamToWorld()).log().transpose();

        }
        printf("]\n");

    }

    printf("K[level %d] = \n", level);
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {

            printf("%f ", frameList->back()->K(level)(i,j));

        }
        printf("\n");
    }

    printf("referenceList->back()->numData[level] = %d\n",
           referenceList->back()->numData[level]);

    //--------------------------------------------------------------------------
    // Update rules
    ceres::LocalParameterization* local_parameterization =
            new ceres::AutoDiffLocalParameterization<se3PlusBN3,
                        6*SPLINE_K, 6*SPLINE_K>;

    // Mininum level bin
    std::vector<int> bin_frontLevel;
    
    //--------------------------------------------------------------------------
    // For each avaialble frame
    int numPixUsed = 0;
    for (unsigned int idxFrameList = 0;
         idxFrameList < frameList->size();
         idxFrameList++)
    {

        // Get a Frame and starting row
        Frame *frame = (*frameList)[idxFrameList].get();
        int startRowImage = startRowImageList[idxFrameList];
        int startRowSeg = startRowSegList[idxFrameList];

#if 1//TEMP - check it is working
        // Get a tracking reference from the reference list
        referenceAt = (*referenceList)[idxFrameList].get();
#else
        // Put the last reference always (so it makes all the same reference)
        TrackingReference *reference = referenceList->back();

#endif
        // Reset refPoint and refColourVariance
        refPoint = referenceAt->posData[level];
        refColourVariance = referenceAt->colorAndVarData[level];
        refPoint_max = refPoint + referenceAt->numData[level];
        
#if 0 // FORCE only ONE POINT WARPING FOR DEBUGGING
        
        refPoint_max = refPoint + 2;
        
#endif

        // An array for ycoordinate of the warped point q
        double *qycoord = qycoordList[idxFrameList].get();

        std::vector<int> bin;
        if (lsd_slam::sparseTracking == false)
        {
            
            // Entire semi-dense points
            for(int i=0;i<(int)referenceAt->numData[level];i++)
            {
            
                bin.push_back(i);
            
            }
            
        }
        else
        {
                
            //----------------------------------------------------------------------
            // Select sparse and confident points from all semi-dense points
            
            // Collect all variance
            std::vector<float> varianceList;
            for (;refPoint < refPoint_max; refPoint++, refColourVariance++)
            {
                
                varianceList.push_back((*refColourVariance)[1]*rand()*0.001f);
                
            }

            // Reset refPoint and refColourVariance
            refPoint = referenceAt->posData[level];
            refColourVariance = referenceAt->colorAndVarData[level];
            refPoint_max = refPoint + referenceAt->numData[level];
            int numSemiPtsInLevel = referenceAt->numData[level];

            // Sort index by the inverse depth variance
            std::vector<int> sort_idx(varianceList.size());
            iota(sort_idx.begin(), sort_idx.end(), 0);
            sort(sort_idx.begin(),
                 sort_idx.end(),
                 [&varianceList](const int i1,
                    const int i2) -> bool
                 {
                     return varianceList[i1] < varianceList[i2];
                 });
        
            // Set 9 bins
            std::vector<int> bin11;
            std::vector<int> bin12;
            std::vector<int> bin13;
            std::vector<int> bin21;
            std::vector<int> bin22;
            std::vector<int> bin23;
            std::vector<int> bin31;
            std::vector<int> bin32;
            std::vector<int> bin33;
            int MAX_BIN_SIZE = 50*(int)pow(2,level);
            
            // Nine bins for 2D image coordinates
            for (auto it: sort_idx)
            {
                // Point at the index i
                Eigen::Vector3f* point = refPoint + it;
                
                // X and Y coordiantes
                float x = ((*point)[0]/(*point)[2])*frame->K(level)(0,0) +
                                frame->K(level)(0,2);
                float y = ((*point)[1]/(*point)[2])*frame->K(level)(1,1) +
                                frame->K(level)(1,2);
                
                // Drop samples to the bins
                if ((x >= 0) &&
                    (x < frame->width(level)/3.0) &&
                    (y >= 0) &&
                    (y < frame->height(level)/3.0))
                {
                    // Bin (1, 1)
                    if ((int)bin11.size() < MAX_BIN_SIZE)
                    {
                        
                        bin11.push_back(it);
                    
                    }
                    
                }
                else if ((x >= 0) &&
                    (x < frame->width(level)/3.0) &&
                    (y >= frame->height(level)/3.0) &&
                    (y < frame->height(level)*2.0/3.0))
                {
                    // Bin (1, 2)
                    if ((int)bin12.size() < MAX_BIN_SIZE)
                    {
                        
                        bin12.push_back(it);
                        
                    }
                    
                }
                else if ((x >= 0) &&
                    (x < frame->width(level)/3.0) &&
                    (y >= frame->height(level)*2.0/3.0) &&
                    (y < frame->height(level)*3.0/3.0))
                {
                    // Bin (1, 3)
                    if ((int)bin13.size() < MAX_BIN_SIZE)
                    {
                        
                        bin13.push_back(it);
                        
                    }
                    
                }
                else if ((x >= frame->width(level)/3.0) &&
                    (x < frame->width(level)*2.0/3.0) &&
                    (y >= 0) &&
                    (y < frame->height(level)/3.0))
                {
                    // Bin (2,1)
                    if ((int)bin21.size() < MAX_BIN_SIZE)
                    {
                        
                        bin21.push_back(it);
                        
                    }
                    
                }
                else if ((x >= frame->width(level)/3.0) &&
                    (x < frame->width(level)*2.0/3.0) &&
                    (y >= frame->height(level)/3.0) &&
                    (y < frame->height(level)*2.0/3.0))
                {
                    // Bin (2,2)
                    if ((int)bin22.size() < MAX_BIN_SIZE)
                    {
                        
                        bin22.push_back(it);
                        
                    }
                    
                }
                else if ((x >= frame->width(level)/3.0) &&
                    (x < frame->width(level)*2.0/3.0) &&
                    (y >= frame->height(level)*2.0/3.0) &&
                    (y < frame->height(level)*3.0/3.0))
                {
                    // Bin (2,3)
                    if ((int)bin23.size() < MAX_BIN_SIZE)
                    {
                        
                        bin23.push_back(it);
                        
                    }
                    
                }
                else if ((x >= frame->width(level)*2.0/3.0) &&
                    (x < frame->width(level)*3.0/3.0) &&
                    (y >= 0) &&
                    (y < frame->height(level)/3.0))
                {
                    // Bin (3,1)
                    if ((int)bin31.size() < MAX_BIN_SIZE)
                    {
                        
                        bin31.push_back(it);
                        
                    }
                    
                }
                else if ((x >= frame->width(level)*2.0/3.0) &&
                    (x < frame->width(level)*3.0/3.0) &&
                    (y >= frame->height(level)/3.0) &&
                    (y < frame->height(level)*2.0/3.0))
                {
                    // Bin (3,2)
                    if ((int)bin32.size() < MAX_BIN_SIZE)
                    {
                        
                        bin32.push_back(it);
                        
                    }
                    
                }
                else if ((x >= frame->width(level)*2.0/3.0) &&
                    (x < frame->width(level)*3.0/3.0) &&
                    (y >= frame->height(level)*2.0/3.0) &&
                    (y < frame->height(level)*3.0/3.0))
                {
                    // Bin (3,3)
                    if ((int)bin33.size() < MAX_BIN_SIZE)
                    {
                        
                        bin33.push_back(it);
                        
                    }
                    
                }
                
            } // End for sort_idx
            
            //----------------------------------------------------------------------
            // Merge all bins into a single bin
            
            int numSparsePts = bin11.size() + bin12.size() + bin13.size() +
                    bin21.size() + bin22.size() + bin23.size() +
                    bin31.size() + bin32.size() + bin33.size();
            if (numSparsePts < MIN_GOODPERALL_PIXEL_ABSMIN*
                frame->width(level)*
                frame->height(level)*
                lsd_slam::pixelDensityForTracking)
            {
                
                // use entire
                fprintf(stderr, "numSparsePts reset %d with %d < %f\n",
                        numSparsePts, numSemiPtsInLevel,
                        frame->width(level)*
                        frame->height(level)*
                        lsd_slam::pixelDensityForTracking);
                bin = sort_idx;
                
                
            }
            else
            {
            
                // preallocate memory
                bin.reserve(bin11.size() + bin12.size() + bin13.size() +
                            bin21.size() + bin22.size() + bin23.size() +
                            bin31.size() + bin32.size() + bin33.size());
                bin.insert( bin.end(), bin11.begin(), bin11.end() );
                bin.insert( bin.end(), bin12.begin(), bin12.end() );
                bin.insert( bin.end(), bin13.begin(), bin13.end() );
                bin.insert( bin.end(), bin21.begin(), bin21.end() );
                bin.insert( bin.end(), bin22.begin(), bin22.end() );
                bin.insert( bin.end(), bin23.begin(), bin23.end() );
                bin.insert( bin.end(), bin31.begin(), bin31.end() );
                bin.insert( bin.end(), bin32.begin(), bin32.end() );
                bin.insert( bin.end(), bin33.begin(), bin33.end() );
            
            }
            
        } // End if sparse_tracking
        
        // Store the minimum level bin if it is the front frame
        if (idxFrameList == 0)
        {
            
            bin_frontLevel = bin; // This will be used in BN3_callback
            
        }

        //----------------------------------------------------------------------
        // For each available semi-dense point, a residual block is built

#if 0 // USE_SPARSE_STEP

        // Determine sparse skip step
        int SPARSE_STEP = 1;
        if (lsd_slam::sparseTracking == false)
        {
            SPARSE_STEP = 1; // No-sparse tracking
        }
        else
        {
            SPARSE_STEP = 5; // Sparse skip pixels
        }
        
        // Use all semi dense pts if too sparse
        int numSparsePts = numSemiPtsInLevel/SPARSE_STEP;
        if (numSparsePts <
            MIN_GOODPERALL_PIXEL_ABSMIN*
            frame->width(level)*
            frame->height(level)*
            lsd_slam::pixelDensityForTracking)
        {
            
            SPARSE_STEP = 1;
            fprintf(stderr, "numSparsePts reset %d with %d < %f\n",
                    numSparsePts, numSemiPtsInLevel,
                    frame->width(level)*
                    frame->height(level)*
                    lsd_slam::pixelDensityForTracking);
            numSparsePts = numSemiPtsInLevel;
            
        }

        // For avaiable pixels,
        for (int idx = 0;
             refPoint < refPoint_max - SPARSE_STEP;
             idx+=SPARSE_STEP,
             refPoint+=SPARSE_STEP,
             refColourVariance+=SPARSE_STEP)
        {
            
#else // USE_BIN_SORTED_MERGED
 
        int idx = -1;
        for (auto it: bin)
        {
            
            // Init pointer to point and color with variance
            refPoint = referenceAt->posData[level];
            refColourVariance = referenceAt->colorAndVarData[level];

            // Pointing the sample out of the bin
            idx++;
            refPoint = refPoint + it;
            refColourVariance = refColourVariance + it;
            
#endif
            
            // Initialize from the point in the reference frame
            const int indexQyCoord = idx;

            // Scale the reference point by a scalar (e.g. 1.0)
            Eigen::Vector3f scaled_refPoint = (*refPoint)*(float)scale;
            
            // Check if the refPoint contains nan
            if (isnan(scaled_refPoint[0]) == true ||
                isnan(scaled_refPoint[1]) == true ||
                isnan(scaled_refPoint[2]) == true)
            {
                
                // 3D point contains a NaN value therefore do not consider
                // in optimization for image warping
                continue;
                
            }

#if 0//// new_refPoint is a reference point wrt world

            //------------------------------------------------------------------
            // RefPoint transformed by its reference pose
            // new_refPoint (3D point is with respect to the world???)
            // Then, how it projects to the reference camera correctly??
            // I guess getScaledCamToWorld actullay gives us
            // a pose of the world with respect to the reference camera.
            // then, this new_refPoint is actually a 3D point
            // with respect to the camera reference. then
            // the original scaled_refPoint was actually
            // a 3D point with respect to the world system
            // Let me test..
            Eigen::Matrix<double,3,3> rotRef =
                    se3FromSim3(referenceAt->keyframe
                               ->getScaledCamToWorld()).rotationMatrix();
            Eigen::Matrix<double,3,1> transRef =
                    se3FromSim3(referenceAt->keyframe
                                ->getScaledCamToWorld()).translation();

            Eigen::Vector3f new_refPoint = // This should be 3D point wrt World
                    rotRef.cast<float>()*scaled_refPoint
                    + transRef.cast<float>();

#else // new_refPoint is a reference point wrt reference

            Eigen::Vector3f new_refPoint = scaled_refPoint;

#endif
            Eigen::Vector2f new_refColourVariance = *refColourVariance;

#if 0
            Sophus::SE3d rowPose_RefPoint; // TEMP TO BE DONE
#else
            //------------------------------------------------------------------
            // A row pose of an image point in the reference image
            Spline<double> spline(controlPoses, SPLINE_K); // World to control
            double xRef = ((new_refPoint)[0]/(new_refPoint)[2]) *
                            frame->K(level)(0,0) +
                            frame->K(level)(0,2);
            double yRef = ((new_refPoint)[1]/(new_refPoint)[2]) *
                            frame->K(level)(1,1) +
                            frame->K(level)(1,2);
            double yRef_update = yRef;
            if (undistorter->isNoRectification() == false &&
                methodRollingShutter == 1)
            {
                
                // Update yRef via undistortion
                float x_dist, y_dist;
                undistorter->distortPointLevel((float)xRef, (float)yRef,
                                               &x_dist, &y_dist,
                                               level);
                yRef_update = (double)y_dist;
                
            }
            double uTime_RefPoint =
                    (yRef_update + 0
                     + ((int) referenceAt->keyframe->startRowImage >> level)
                     - ((int)referenceAt->keyframe->startRowSeg >> level))
                    / ((int)referenceAt->keyframe->rowPerSegment >> level);
            int idxSC_RefPoint =
                    (int)floor(referenceAt->keyframe->startRowImage
                               / referenceAt->keyframe->rowPerSegment)
                               + NUM_EXTRA_CTRLPTS;
            Sophus::SE3d rowPose_RefPoint = // World to RowRef
                    Sophus::SE3d(spline.getPoseOnSpline(uTime_RefPoint,
                                                        idxSC_RefPoint));
#endif
            //------------------------------------------------------------------
            // Set a reasonable initial estimate rathern than
            // a constant value (e.g. 0.5) between 0 and 1
            //
            // qycoord[indexQyCoord] = 0.5; // Any value between 0 and 1
            //
            // Use the refPoint ycoordinate as initial
            double x = ((new_refPoint)[0]/(new_refPoint)[2]) *
                    frame->K(level)(0,0) +
                    frame->K(level)(0,2);
            double y = ((new_refPoint)[1]/(new_refPoint)[2]) *
                    frame->K(level)(1,1) +
                    frame->K(level)(1,2);
            double y_new = y;
            if (undistorter->isNoRectification() == false &&
                methodRollingShutter == 1)
            {
                
                // Update yRef via undistortion
                float x_dist, y_dist;
                undistorter->distortPointLevel((float)x, (float)y,
                                               &x_dist, &y_dist, level);
                y_new = (double)y_dist;
                
            }
            double uacc = (y_new + 0 + ((int)startRowImage >> level)
                             - ((int)startRowSeg >> level))
                          / ((int)rowPerSegment >> level);
            qycoord[indexQyCoord] = uacc;
            
            
            // Set uppper and lower bound of u time values which
            // depend on start row
            double lowerBound_uValue = (double)
                (0 + 0 // first row
                 + ((int)startRowImage >> level)
                 - ((int)startRowSeg >> level))
                / ((int)rowPerSegment >> level);
            double upperBound_uValue = (double)
                ((frame->height_rollingShutter_inputImage >> level) + 0 // last row
                 + ((int)startRowImage >> level)
                 - ((int)startRowSeg >> level))
                / ((int)rowPerSegment >> level);

            //------------------------------------------------------------------
            // Add a residual block for each semi-dense pixel
            // computing warping error and returning the estimate of pose
            // Parameter dimensions:
            //    2 for residual
            //    24 for 4 control poses
            //    1 for ycoordinate of the warped pixel q
            //      but be careful this ycoordinate should be different
            //      depending on the source image (idxSrcSegment)
            //
            //------------------------------------------------------------------
            // NOTE: Any values passed by reference (i.e. pointers) need
            //       to be alive until the time of Solver being called.
            //       E.g. referenceAt, frame and undistorter
            //            should not be a dangling reference when
            //            the solver is called.
            //------------------------------------------------------------------
            // Print-out
            if (idx == 0)
            {
                printf("AddResidualBlock@idx=%d,idxFrameList=%d:\n"
                       "    level = %d\n"
                       "    new_refPoint = %f %f %f\n"
                       "    new_refColourVariance = %f %f\n"
                       "    referenceAt frameID = %d\n"
                       "    frame ID = %d\n"
                       "    startRowImage = %d\n"
                       "    startRowSeg = %d\n"
                       "    rowPerSegment = %f\n"
                       "    controlPoses->size() = %lu\n"
                       "    idxSC_frontKF = %d\n"
                       "    idxSC_RefPoint = %d\n"
                       "    idxSplineCtrlPts = %d\n"
                       "    scale = %f\n"
                       "    indexQyCoord = %d\n"
                       "    uTime_RefPoint = %f\n"
                       "    yRef_update = %f\n"
                       "    lowerBound_uValue = %f\n"
                       "    upperBound_uValue = %f\n",
                       idx, idxFrameList,
                       level,
                       new_refPoint[0], new_refPoint[1], new_refPoint[2],
                       new_refColourVariance[0], new_refColourVariance[1],
                       referenceAt->frameID,
                       frame->id(),
                       startRowImage,
                       startRowSeg,
                       rowPerSegment,
                       controlPoses.size(),
                       idxSC_frontKF,
                       idxSC_RefPoint,
                       idxSplineCtrlPts,
                       scale,
                       indexQyCoord,
                       uTime_RefPoint,
                       yRef_update,
                       lowerBound_uValue,
                       upperBound_uValue);
            }
            
            numPixUsed++;
            ceres::ResidualBlockId rbid = problem.AddResidualBlock(
                    new AutoDiffCostFunction<WarpResidual_BN3,
                    2, 6*SPLINE_K, 1>
                    // 1 residual, 24 for k control poses, 1 ycoord
                    (new WarpResidual_BN3(
                         level,
                         new_refPoint, // 3D point wrt reference
                         new_refColourVariance,
                         referenceAt,
                         rowPose_RefPoint, // World to rowRef
                         frame,
                         undistorter,
                         startRowImage,
                         startRowSeg,
                         rowPerSegment,
                         controlPoses, // world to control
                         idxSC_frontKF, // idx of ctrl pt for front keyframe
                         idxSC_RefPoint, // idx of ctrl pt for each keyframe
                         idxSplineCtrlPts, // idx of ctrl pt to be estimated
                         scale,
                         idx,
                         settings,
                         files_,
                         methodRollingShutter)),
                    //new ceres::HuberLoss(
                    //        settings.huber_d*settings.huber_d/4.0),
                    new ceres::HuberLoss(1e-9),
                    //new ceres::TolerantLoss(255, 0.1),
                    //NULL,
                    //new lsd_slam::JHK_PseudoHuberLoss(1e-4),

#if 0 // USE WORLD_TO_CONTROL pose (ORIGINAL IMPLEMENTATION)

                    //pose, // WorldToControl

#else // USE REF_TO_CONTROL pose

                    // An array of k control points as 6-vector in se(3)
                    // representing transformations
                    // from a control point for the current tracking frame
                    // to control point for the front keyframe
                    pose_of_kCtrlPts_to_frontKeyframe, // Ctrl to FrontKF

#endif
                    qycoord + indexQyCoord);


#if 1 // BOUND
            // Set qycoord lower and upper bound in estimation
            problem.SetParameterLowerBound(qycoord + indexQyCoord,
                                            0, lowerBound_uValue); // POSSIBLE IMPROVE?
            problem.SetParameterUpperBound(qycoord + indexQyCoord,
                                            0, upperBound_uValue); // USE uValue upper/lower?
#endif

#if 1 // UPDATING RULE
            // Set local parametrization for rule updates in Lie algebra
#if 0 // USE WORLD_TO_CONTROL pose (ORIGINAL IMPLEMENTATION)
            problem.SetParameterization(pose, local_parameterization);
#else // USE REF_TO_CONTROL pose
            problem.SetParameterization(pose_of_kCtrlPts_to_frontKeyframe,
                                        local_parameterization);
#endif
#endif


        } // end for refPoint
        
    } // end for frameList
    
    // Print-out
    fprintf(stderr, "Number of pixels used"
            " in tracking [frame %d, level %d] = %d\n",
            (*frameList)[0].get()->id(), level , numPixUsed);




    //--------------------------------------------------------------------------
    // Make rotation parameter constant
    for (int i=0;i<SPLINE_K;i++)
    {
        // The first 3 out of 6 parameters is for translation
        // The last 3 out of 6 is for rotation
        {

#if 0 //SMALL ROTATION
            float delta = M_1_PI/6.0f; // 30 degree
            problem.SetParameterLowerBound(pose_of_kCtrlPts_to_frontKeyframe,
                                           6*i + 3, 0.0 - delta);
            problem.SetParameterUpperBound(pose_of_kCtrlPts_to_frontKeyframe,
                                           6*i + 3, 0.0 + delta);
            problem.SetParameterLowerBound(pose_of_kCtrlPts_to_frontKeyframe,
                                           6*i + 4, 0.0 - delta);
            problem.SetParameterUpperBound(pose_of_kCtrlPts_to_frontKeyframe,
                                           6*i + 4, 0.0 + delta);
            problem.SetParameterLowerBound(pose_of_kCtrlPts_to_frontKeyframe,
                                           6*i + 5, 0.0 - delta);
            problem.SetParameterUpperBound(pose_of_kCtrlPts_to_frontKeyframe,
                                           6*i + 5, 0.0 + delta);
#endif
            // DEBUG Purpose
            // Fix z-axis of translation as zeros
            //problem.SetParameterLowerBound(pose, 6*i + 2, 0.0 - delta);
            //problem.SetParameterUpperBound(pose, 6*i + 2, 0.0 + delta);

#if 0 // SMALL_TRANSLATION +/-100
            // Small motion to translation x, y and z directions
            float smallMotion = 0.5; // 0.5 metre
            problem.SetParameterLowerBound(pose_of_kCtrlPts_to_frontKeyframe,
                                           6*i + 0, 0.0 - smallMotion);
            problem.SetParameterUpperBound(pose_of_kCtrlPts_to_frontKeyframe,
                                           6*i + 0, 0.0 + smallMotion);
            problem.SetParameterLowerBound(pose_of_kCtrlPts_to_frontKeyframe,
                                           6*i + 1, 0.0 - smallMotion);
            problem.SetParameterUpperBound(pose_of_kCtrlPts_to_frontKeyframe,
                                           6*i + 1, 0.0 + smallMotion);
            problem.SetParameterLowerBound(pose_of_kCtrlPts_to_frontKeyframe,
                                           6*i + 2, 0.0 - smallMotion);
            problem.SetParameterUpperBound(pose_of_kCtrlPts_to_frontKeyframe,
                                           6*i + 2, 0.0 + smallMotion);
#endif
        }

    }

#if 0 // FIX CONTROL POINTS
    //--------------------------------------------------------------------------
    // Constraint of fixing the control points
    if ((fmod((double)startRowList.back(), rowPerSegment) == 0)
            && (startRowList.back() > 0))
    {
        for (int i=6*0; i<6*(SPLINE_K - 1); i++)
        {

            float delta = 1.0e-12;
            problem.SetParameterLowerBound(pose, i, pose[i] - delta);
            problem.SetParameterUpperBound(pose, i, pose[i] + delta);

        }
    }
#endif

    //--------------------------------------------------------------------------
    // Store problem into problems
    //problemList.push_back(&problem);

    //--------------------------------------------------------------------------
    // Solve using Ceres-solver

    // Set options for solver
    ceres::Solver::Options options;

    // Set a maximum number of iterations for each pyramid level
    //const int maxIterations[6] = {40, 40, 80, 100, 200, 400};
    //const int maxIterations[6] = {80, 80, 160, 200, 400, 800};
    //const int maxIterations[6] = {10, 10, 20, 40, 80, 160};
    //const int maxIterations[6] = {5, 5, 10, 20, 40, 80};
    const int m = 10; // 8 // Default: 8
    //const int maxIterations[6] = {m*5, m*5, m*10, m*20, m*40, m*80};
    const int maxIterations[6] = {m*5, m*20, m*50, m*100, m*100, m*100};
    //const int maxIterations[6] = {1, 1, 1, 1, 1, 1};
    options.max_num_iterations = maxIterations[level];

    if (CERES_THREADS > 1)
    {
        options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY; // If SuiteSparse installed
        //options.linear_solver_type = ceres::SPARSE_SCHUR;
        //options.linear_solver_type = ceres::ITERATIVE_SCHUR;
        //options.preconditioner_type = ceres::SCHUR_JACOBI;
        options.num_threads = CERES_THREADS; // Default: 2; // Maximum number of threads on your system available
        options.num_linear_solver_threads = CERES_THREADS; // Default: 2; // Maximum number of threads
    }
    else
    {
        //options.linear_solver_type = ceres::DENSE_QR;
        options.linear_solver_type = ceres::DENSE_SCHUR;
    }

    options.minimizer_progress_to_stdout = true;
    //options.use_explicit_schur_complement = true; // Boost performance for medium size problem

    //options.use_approximate_eigenvalue_bfgs_scaling = true;

    // Set the parameter tolerance less than the default (1e-8)
    // for more iterations of optimization
    //options.parameter_tolerance = 1e-12;

    // Set the function tolerance less than the default (1e-6)
    //options.function_tolerance = 1e-10;

    // Summary
    ceres::Solver::Summary summary;
    
//    // Evaluate
//    ceres::CRSMatrix jacobian;
//    double cost[3];
//    std::vector<double> initial_residuals;
//    problem.Evaluate(Problem::EvaluateOptions(),
//                     cost,
//                     &initial_residuals,
//                     NULL, // No gradients
//                     &jacobian);
    
//    // Remove outlier residual blocks
//    for (unsigned int i=0;i<initial_residuals.size();i++)
//    {
//        ceres::ResidualBlockID rbid = i;
//        problem.RemoveResidualBlock(i);
//        
//    }
    
#if 0
    // Initial evaluation and print out
    std::cout << std::endl << "Run init evaluation... " << std::endl;
    ceres::CRSMatrix jacobian;
    std::vector<double> residual;
    std::vector<double> gradient;
    double cost[2];
    problem.Evaluate(Problem::EvaluateOptions(),
                     cost, &residual, &gradient, &jacobian);
    std::cout << "Done" << std::endl << std::endl;
    
    // Print init cost
    for (int i=0; i<2; i++)
        printf("cost[%d] = %24.18f\n", i, cost[i]);
    
    printf("Ceres Jacobian = \n\n");
    Eigen::Matrix<double,2,6*SPLINE_K + 1> J;
    for (int i=0; i<(int)jacobian.values.size();)
    {
        for (int p=0;p<J.rows();p++)
            for (int q=0;q<J.cols();q++)
            {
                J(p,q) = jacobian.values[i];
                i++;
            }
    }
    
    for (int j=0;j<J.cols(); j++)
    {
        for (int i=0;i<J.rows(); i++)
        {
            
            printf("%24.18f  ", J(i,j));
            
        }
        printf("\n");
    }
//    printf("Ceres gradients = \n\n");
//    for (int i=0;i<(int)gradient.size();i++)
//    {
//        printf("%03d: %24.18f\n", i, gradient[i]);
//    }
#endif


    // Solve
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.FullReport() << "\n";
    std::cout << summary.message << std::endl;
    final_cost_ = summary.final_cost;
    num_residual_blocks_ = summary.num_residual_blocks;
    
#if 0
    // Final evaluation and print out
    std::cout << std::endl << "Run final evaluation... " << std::endl;
    ceres::CRSMatrix jacobian;
    std::vector<double> gradient;
    std::vector<double> residual;
    gradient.clear();
    double cost[2];
    problem.Evaluate(Problem::EvaluateOptions(),
                     cost, &residual, &gradient, &jacobian);
    std::cout << "Done" << std::endl << std::endl;
    
    // Print final cost
    for (int i=0; i<2; i++)
        printf("cost[%d] = %24.18f\n", i, cost[i]);
    
    printf("Ceres Jacobian = \n\n");
        Eigen::Matrix<double,2,6*SPLINE_K + 1> J;
    for (int i=0; i<(int)jacobian.values.size();)
    {
        for (int p=0;p<J.rows();p++)
            for (int q=0;q<J.cols();q++)
            {
                J(p,q) = jacobian.values[i];
                i++;
            }
    }
    
    for (int j=0;j<J.cols(); j++)
    {
        for (int i=0;i<J.rows(); i++)
        {
            
            printf("%24.18f  ", J(i,j));
            
        }
        printf("\n");
    }
#endif
    
//    printf("Ceres gradients = \n\n");
//    for (int i=0;i<(int)gradient.size();i++)
//    {
//        printf("%03d: %24.18f\n", i, gradient[i]);
//    }

#if 1 // Use the first element for callback2

    // Update counters using callback
    BN3_UserCallback callback2(level,
                               referenceList->front().get(),
                               frameList->front().get(),
                               controlPoses, // World to Control
#if 0 // USE WORLD_TO_CONTROL pose (ORIGINAL IMPLEMENTATION)
                               pose, // WorldToControl
#else // USE REF_TO_CONTROL pose
                               pose_of_kCtrlPts_to_frontKeyframe, // Reference to Control
#endif
                               idxSC_frontKF,
                               idxSplineCtrlPts,
                               qycoordList.front().get(),
                               rowPerSegment,
                               startRowImageList.front(),
                               startRowSegList.front(),
                               scale,
                               settings,
                               &debugImageWeightsRollingShutter,
                               &debugImageResidualsRollingShutter,
                               &debugImageSecondFrameRollingShutter,
                               &debugImageOldImageWarpedRollingShutter,
                               &debugImageOldImageSourceRollingShutter,
                               undistorter,
                               methodRollingShutter,
                               bin_frontLevel,
                               true);

#else // Use the last element for callback2

    // Update counters using callback
    BN3_UserCallback callback2(level,
                              referenceList->back().get(),
                              frameList->back().get(),
                              controlPoses, // World to Control
                              pose, // WorldToControl
                              idxSplineCtrlPts,
                              qycoordList.back().get(),
                              rowPerSegment,
                              startRowImageList.back(),
                              startRowSegList.back(),
                              scale,
                              settings,
                              &debugImageWeightsRollingShutter,
                              &debugImageResidualsRollingShutter,
                              &debugImageSecondFrameRollingShutter,
                              &debugImageOldImageWarpedRollingShutter,
                              &debugImageOldImageSourceRollingShutter,
                              undistorter,
                              methodRollingShutter,
                              true);
#endif

    // Accumulate counters for each source segment
    int sum_lastGoodCount = 0;
    int sum_lastBadCount = 0;
    float sum_lastUsageCount = 0.0;

    int tmpGoodCount = 0;
    int tmpBadCount = 0;
    float tmpUsageCount = 0.0;
    callback2.getStatisticsForWarpedResidual<double>(&tmpGoodCount,
                                                     &tmpBadCount,
                                                     &tmpUsageCount);
    //this->debugInfo_last_residual_ = (float)last_residual_;
    
    sum_lastGoodCount = sum_lastGoodCount + tmpGoodCount;
    sum_lastBadCount = sum_lastBadCount + tmpBadCount;
    sum_lastUsageCount = sum_lastUsageCount + tmpUsageCount;
        
    printf("WRITING.....\n");
    callback2.writeDebugImage();

    // Debug finish
    calcResidualAndBuffers_debugFinish(frameList->back()->width(level));

#if 0 // DEBUG
    //--------------------------------------------------------------------------
    // Print out several pose warping result for DEBUGGING
    {

        int offset = 0;
#if 1//TEMP testing
        referenceAt = referenceList->front().get();
#endif
        double *qycoord = qycoordList.front().get();

        // Initialize variables for semi-dense points
        refPoint = referenceAt->posData[level];
        refColourVariance = referenceAt->colorAndVarData[level];
        refPoint_max = refPoint + referenceAt->numData[level];

        printf("{\n");
        for (int idx = 0;
             refPoint < refPoint_max;
             idx++, refPoint++, refColourVariance++)
        {

            const int indexQyCoord = idx;

            // Scale the reference 3D point by the scale factor
            // from tracking reference, which might be updated by
            // keyframe changes
            Eigen::Vector3f scaled_refPoint = (*refPoint)*scale;

#if 0 // refPoint is a reference point wrt world


            // RefPoint should be transformed by its reference pose with respect
            // to the world coordinate system
            Eigen::Matrix<double,3,3> rotRef =
                    se3FromSim3(referenceAt->keyframe
                               ->getScaledCamToWorld()).rotationMatrix();
            Eigen::Matrix<double,3,1> transRef =
                    se3FromSim3(referenceAt->keyframe
                                ->getScaledCamToWorld()).translation();

            Eigen::Vector3f new_refPoint =
                    rotRef.cast<float>()*scaled_refPoint
                    + transRef.cast<float>();


#else // refPoint is a reference point wrt reference

            Eigen::Vector3f new_refPoint = scaled_refPoint;
#endif

#if 0//DEBUG
            printf("REF_Frame%08d_Level%d: %f %f %f\n",
                   frameList->front()->id(), level,
                   (*refPoint)[0], (*refPoint)[1], (*refPoint)[2]);

            printf("NEWREF_Frame%08d_Level%d: %f %f %f\n",
                   frameList->front()->id(), level,
                   (new_refPoint)[0], (new_refPoint)[1], (new_refPoint)[2]);
#endif


            // Print out $0
            printf("SOL_Frame%08d_Level%d: ", frameList->front()->id(), level);
            double x = ((new_refPoint)[0]/(new_refPoint)[2]) *
            frameList->front()->K(level)(0,0) +
            frameList->front()->K(level)(0,2);
            double y = ((new_refPoint)[1]/(new_refPoint)[2]) *
            frameList->front()->K(level)(1,1) +
            frameList->front()->K(level)(1,2);

            // $1 $2
            printf("%f %f ", x, y);

            // x y z utime for $3 $4 $5 $6
            printf("%f %f %f ", (new_refPoint)[0],
                    (new_refPoint)[1], (new_refPoint)[2]);
            printf("%f ", qycoord[indexQyCoord]);


            //------------------------------------------------------------------
            // Print out all poses on the B-spline

            // For coordinate system conversion
            Sophus::SE3d worldToRef =
                    se3FromSim3(referenceAt->keyframe->pose->getCamToWorld()
                                .cast<double>()).inverse();

            // Copy all control poses
            std::vector<Sophus::SE3d> controlPoses_T;
            for (unsigned int i=0; i<controlPoses.size(); i++)
            {

                controlPoses_T.push_back((controlPoses)[i]);// World to Control

            }

            // Then set the last k control poses with the pose
            for (int i=0; i<SPLINE_K; i++)
            {

#if 0 // USE WORLD_TO_CONTROL pose (ORIGINAL IMPLEMENTATION)
                Eigen::Matrix<double,6,1> pose_6vec;
                for (int j=0; j<6; j++)
                    pose_6vec[j] = pose[6*i + j]; // World to Control
#else // USE REF_TO_CONTROL pose
                Eigen::Matrix<double,6,1> pose_6vec;
                for (int j=0; j<6; j++)
                    pose_6vec[j] = pose_of_kCtrlPts_to_frontKeyframe[6*i + j]; // Ctrl to FK
#endif
                Eigen::Matrix<double,6,1> tmp =
                        Sophus::SE3d::exp(pose_6vec).log();
                std::cout << tmp.transpose() << " ";

                // WorldToControl
                Sophus::SE3d newWorldToControl = Sophus::SE3d::exp(pose_6vec);

                // Replace
                controlPoses_T[idxSplineCtrlPts + i] = newWorldToControl;

            }

            // Create a B-Spline for WorldToControl and a get pose by u param
            lsd_slam::Spline<double> spline(controlPoses_T, SPLINE_K);
            double u_time = qycoord[indexQyCoord];
            int i_spline = idxSplineCtrlPts;
            Eigen::Matrix<double,4,4> pose_on_spline =
                    spline.getPoseOnSpline(u_time, i_spline);
            Sophus::SE3d pose_worldToControl = Sophus::SE3d(pose_on_spline);
            Eigen::Matrix<double,4,4> tmpM = pose_worldToControl.matrix();

            // Print out the pose on B-Spline $31 - $46
            for (int i=0; i<4; i++)
                for (int j=0; j<4; j++)
                    printf("%f ", tmpM(i,j));

            //------------------------------------------------------------------
            // Convert pose_on_spline for RefToControl
            // RefToControl = WorldToControl @ RefToWorld
#if 0 // Use new ref to control
            Sophus::SE3d pose_refToControl = pose_worldToControl
                    * worldToRef;
#else // USE ref to control
            Sophus::SE3d pose_refToControl = pose_worldToControl
                    * worldToRef.inverse();
#endif
            // Warped point
            Eigen::Matrix<double,3,3> rotMat =
                    pose_refToControl.rotationMatrix();
            Eigen::Matrix<double,3,1> transVec =
                    pose_refToControl.translation();
            Eigen::Matrix<double,3,1> pointRef;
            pointRef[0] = (new_refPoint)[0];
            pointRef[1] = (new_refPoint)[1];
            pointRef[2] = (new_refPoint)[2];
            Eigen::Matrix<double,3,1> Wxp = rotMat*pointRef + transVec;
            double fx_l = double(frameList->front()->K(level)(0,0));
            double fy_l = double(frameList->front()->K(level)(1,1));
            double cx_l = double(frameList->front()->K(level)(0,2));
            double cy_l = double(frameList->front()->K(level)(1,2));
            double u_new = (Wxp[0]/Wxp[2])*fx_l + cx_l;
            double v_new = (Wxp[1]/Wxp[2])*fy_l + cy_l;
            double v_new_offset = v_new + offset;

            // $47 $48 $49
            printf("%f %f %f ", u_new, v_new, v_new_offset);

            // Camera info from framelist $50 $51 $52 $53
            printf("%f %f %f %f ", fx_l, fy_l, cx_l, cy_l);

            // Wxp $54 $55 $56
            printf("%f %f %f ", Wxp[0], Wxp[1], Wxp[2]);

            //------------------------------------------------------------------
            // Pose at last row
            Eigen::Matrix<double,4,4> pose_lastRow =
                    spline.getPoseOnSpline(1.0, i_spline);

            // Transform refToControl into WorldToControl
            Sophus::SE3d worldToControl_LastRow =
                    Sophus::SE3d(pose_lastRow) * worldToRef;
            tmpM = worldToControl_LastRow.matrix();

            // Print out the pose on B-Spline $57 - $72
            tmpM = tmpM.inverse().eval();
            for (int i=0; i<4; i++)
                for (int j=0; j<4; j++)
                    printf("%f ", tmpM(i,j));

            //------------------------------------------------------------------
            // Print out the pose on B-Spline $73 - $88
            Eigen::Matrix<double,4,4> pose_firstRow =
                    spline.getPoseOnSpline(0.0, i_spline);

            // Convert
            Sophus::SE3d worldToControl_firstRow =
                    Sophus::SE3d(pose_firstRow) * worldToRef;
            tmpM = worldToControl_firstRow.matrix();

            // Print
            tmpM = tmpM.inverse().eval();
            for (int i=0; i<4; i++)
                for (int j=0; j<4; j++)
                    printf("%f ", tmpM(i,j));

            //------------------------------------------------------------------
            // Print out the pose on B-Spline $89 - $104
            Eigen::Matrix<double,4,4> pose_midRow =
                    spline.getPoseOnSpline(0.5, i_spline);

            // Convert
            Sophus::SE3d worldToControl_midRow =
                    Sophus::SE3d(pose_midRow) * worldToRef;
            tmpM = worldToControl_midRow.matrix();

            // Print
            tmpM = tmpM.inverse().eval();
            for (int i=0; i<4; i++)
                for (int j=0; j<4; j++)
                    printf("%f ", tmpM(i,j));

            // Print out the covariance of the depth $105, $106
            printf("%f %f ", (*refColourVariance)[0], (*refColourVariance)[1]);

            printf("\n");
        }
        printf("}\n");

    }
#endif // END DEBUG

    //--------------------------------------------------------------------------
    // Free
    //delete[] qycoord;
    qycoordList.clear();


    //--------------------------------------------------------------------------
    // Return information
    bool isUsableSolution = summary.IsSolutionUsable();
    if (isUsableSolution == true)
    {

        last_residual_ = final_cost_ / num_residual_blocks_;
        fprintf(stderr, "    last_residual_: %f, "
                        "final_cost_: %f, "
                        "num_res_block_: %d\n",
                        last_residual_,
                        final_cost_,
                        num_residual_blocks_);
        lsd_slam::WarpResidual_BN3::lastGoodCount_ = sum_lastGoodCount;
        lsd_slam::WarpResidual_BN3::lastBadCount_ = sum_lastBadCount;
        lsd_slam::WarpResidual_BN3::usageCount_ = sum_lastUsageCount;
        fprintf(stderr, "    usageCount: %f, "
                        "lastGoodCount: %d, "
                        "lastBadCount: %d\n",
                        lsd_slam::WarpResidual_BN3::usageCount_,
                        lsd_slam::WarpResidual_BN3::lastGoodCount_,
                        lsd_slam::WarpResidual_BN3::lastBadCount_);

        return true;

    }
    else
    {

        return false;

    }

}

} // end of namespace


