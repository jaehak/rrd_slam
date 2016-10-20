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
//------------------------------------------------------------------------------
// SE(3) based tracking of new frame to the reference key frame for a rolling
// shutter camera
//------------------------------------------------------------------------------
// Copyright@2014 University of Adelaide
//
// Author: Jae-Hak Kim (jaehak.kim@adelaide.edu.au)
//------------------------------------------------------------------------------
#include "Tracking/SE3TrackerRollingShutter.h"
#include <opencv2/highgui/highgui.hpp>
#include "DataStructures/Frame.h"
#include "Tracking/TrackingReference.h"
#include "util/globalFuncs.h"
#include "IOWrapper/ImageDisplay.h"
#include "Tracking/least_squares.h"
#include <math.h>

#include "Tracking/WarpResidual.h"
#include "ceres/ceres.h"
#include "Tracking/Spline.h"

#include <unsupported/Eigen/MatrixFunctions>


namespace lsd_slam
{

// Allocation of static debug image
//cv::Mat lsd_slam::SE3TrackerRollingShutter::debugImageResidualsRollingShutter
//    = cv::Mat(480, 640, CV_8UC3);
//cv::Mat lsd_slam::SE3TrackerRollingShutter::debugImageOldImageSourceRollingShutter
//    = cv::Mat(480, 640, CV_8UC3);
//cv::Mat lsd_slam::SE3TrackerRollingShutter::debugImageOldImageWarpedRollingShutter
//    = cv::Mat(480, 640, CV_8UC3);
//cv::Mat lsd_slam::SE3TrackerRollingShutter::debugImageSecondFrameRollingShutter
//    = cv::Mat(480, 640, CV_8UC3);
//cv::Mat lsd_slam::SE3TrackerRollingShutter::debugImageWeightsRollingShutter
//    = cv::Mat(480, 640, CV_8UC3);

#if defined(ENABLE_NEON)
	#define callOptimized(function, arguments) function##NEON arguments
#else
	#if defined(ENABLE_SSE)
        #define callOptimized(function, arguments) (USESSE ? function##SSE arguments : function arguments)
	#else
		#define callOptimized(function, arguments) function arguments
	#endif
#endif


SE3TrackerRollingShutter::SE3TrackerRollingShutter(
        int w, int h, Eigen::Matrix3f K,
        Undistorter *undistorter,
        int width_bsplineSegImg, int height_bsplineSegImg,
        int methodRollingShutter)
: SE3Tracker(w, h, K)
{

    this->undistorter = undistorter;
    this->methodRollingShutter = methodRollingShutter;

    // Set B-Spline segment image size
    this->width_bsplineSegImg = width_bsplineSegImg;
    this->height_bsplineSegImg = height_bsplineSegImg;

    // Initialize the debug image
    //debugImageOriginalInput = cv::Mat(undistorter->getInputHeight(),
    //                              undistorter->getInputWidth(), CV_8UC3);
    debugImageWeightsRollingShutter =
            cv::Mat(h, w, CV_8UC3);
    debugImageResidualsRollingShutter =
            cv::Mat(h, w, CV_8UC3);
    debugImageSecondFrameRollingShutter =
            cv::Mat(h, w, CV_8UC3);
    debugImageOldImageWarpedRollingShutter =
            cv::Mat(h, w, CV_8UC3);
    debugImageOldImageSourceRollingShutter =
            cv::Mat(h, w, CV_8UC3);

    debugGradientXImage0_RollingShutter =
            cv::Mat(h, w, CV_8UC3);
    debugGradientXImage1_RollingShutter =
            cv::Mat(h, w, CV_8UC3);

    debugRefImage0 =
            cv::Mat(h, w, CV_8UC3);
    debugRefImage1 =
            cv::Mat(h, w, CV_8UC3);
    debugRefImage2 =
            cv::Mat(h, w, CV_8UC3);
    
    // Create internal spline structure
    std::vector<Sophus::SE3d> controlPoses;
    for (int i=0;i<SPLINE_K;i++)
    {
        
        controlPoses.push_back(Sophus::SE3d(Eigen::Matrix4d::Identity()));
        
    }
    internalSpline = std::make_shared<Spline<double> >(controlPoses, SPLINE_K);
    
}

SE3TrackerRollingShutter::~SE3TrackerRollingShutter()
{

    debugImageResidualsRollingShutter.release();
    debugImageWeightsRollingShutter.release();
    debugImageSecondFrameRollingShutter.release();
    debugImageOldImageSourceRollingShutter.release();
    debugImageOldImageWarpedRollingShutter.release();

    debugGradientXImage0_RollingShutter.release();
    debugGradientXImage1_RollingShutter.release();

    debugRefImage0.release();
    debugRefImage1.release();
    debugRefImage2.release();

    // Set boolean for indicating windows moved or not
    isWindowMoved = false;
    isWindowMoved2 = false;

}


    
// tracks a frame.
// first_frame has depth, second_frame DOES NOT have depth.
SE3 SE3TrackerRollingShutter::trackFrame(
                           TrackingReference* reference,
                           Frame* frame,
                           const SE3& frameToReference_initialEstimate)
{
    
    boost::shared_lock<boost::shared_mutex> lock = frame->getActiveLock();
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

        const float* frameImage = frame->image();
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
    Sophus::SE3f referenceToFrame =
            frameToReference_initialEstimate.inverse().cast<float>();
    NormalEquationsLeastSquares ls;
    
    int numCalcResidualCalls[PYRAMID_LEVELS];
    int numCalcWarpUpdateCalls[PYRAMID_LEVELS];
    
    float last_residual = 0;
    
    for (int lvl=SE3TRACKING_MAX_LEVEL-1; lvl >= SE3TRACKING_MIN_LEVEL; lvl--)
    {

        numCalcResidualCalls[lvl] = 0;
        numCalcWarpUpdateCalls[lvl] = 0;
        
        reference->makePointCloud(lvl);
        
        callOptimized(calcResidualAndBuffers,
                        (reference->posData[lvl],
                         reference->colorAndVarData[lvl],
                         SE3TRACKING_MIN_LEVEL == lvl ?
                            reference->pointPosInXYGrid[lvl] : 0,
                         reference->numData[lvl],
                         frame,
                         referenceToFrame,
                         lvl,
                         (plotTracking && lvl == SE3TRACKING_MIN_LEVEL)));

        if (buf_warped_size <
                MIN_GOODPERALL_PIXEL_ABSMIN * (width>>lvl)*(height>>lvl)*0.25)
        {

            diverged = true; //THIS DIVERGED SETTING SHOULD BE RECOMPUTED BASED ON A WIDTH SIZE OF ROW
            trackingWasGood = false;
            return SE3();

        }
        
        if (useAffineLightningEstimation)
        {

            affineEstimation_a = affineEstimation_a_lastIt;
            affineEstimation_b = affineEstimation_b_lastIt;

        }
        float lastErr =
                callOptimized(calcWeightsAndResidual,(referenceToFrame));
        
        numCalcResidualCalls[lvl]++;
        
        
        float LM_lambda = settings.lambdaInitial[lvl];
        
        for (int iteration=0; iteration < settings.maxItsPerLvl[lvl];
             iteration++)
        {
            
            callOptimized(calculateWarpUpdate,(ls));
            
            numCalcWarpUpdateCalls[lvl]++;
            
            iterationNumber = iteration;
            
            int incTry=0;
            while (true)
            {

                // solve LS system with current lambda
                Vector6 b = -ls.b;
                Matrix6x6 A = ls.A;
                for(int i=0;i<6;i++) A(i,i) *= 1+LM_lambda;
                Vector6 inc = A.ldlt().solve(b);
                incTry++;
                
                // apply increment. pretty sure this way round is correct,
                // but hard to test.
                Sophus::SE3f new_referenceToFrame =
                        Sophus::SE3f::exp((inc)) * referenceToFrame;
                //Sophus::SE3f new_referenceToFrame =
                //      referenceToFrame * Sophus::SE3f::exp((inc));
                
                // re-evaluate residual
                callOptimized(calcResidualAndBuffers,
                              (reference->posData[lvl],
                               reference->colorAndVarData[lvl],
                               SE3TRACKING_MIN_LEVEL == lvl ?
                                   reference->pointPosInXYGrid[lvl] : 0,
                               reference->numData[lvl],
                               frame,
                               new_referenceToFrame,
                               lvl,
                               (plotTracking && lvl == SE3TRACKING_MIN_LEVEL)));

                if (buf_warped_size < MIN_GOODPERALL_PIXEL_ABSMIN *
                        (width>>lvl)*(height>>lvl)*0.25)
                {

                    diverged = true; //THIS DIVERGED SETTING SHOULD BE RECOMPUTED BASED ON A WIDTH SIZE OF ROW
                    trackingWasGood = false;
                    return SE3();

                }
                
                float error = callOptimized(calcWeightsAndResidual,
                                            (new_referenceToFrame));
                numCalcResidualCalls[lvl]++;
                
                
                // accept inc?
                if( error < lastErr)
                {

                    // accept inc
                    referenceToFrame = new_referenceToFrame;
                    if(useAffineLightningEstimation)
                    {

                        affineEstimation_a = affineEstimation_a_lastIt;
                        affineEstimation_b = affineEstimation_b_lastIt;

                    }
                    
                    
                    if (enablePrintDebugInfo && printTrackingIterationInfo)
                    {

                        // debug output
                        printf("(%d-%d): ACCEPTED increment of "
                               "%f with lambda %.1f, residual: %f -> %f\n",
                               lvl,iteration, sqrt(inc.dot(inc)),
                               LM_lambda, lastErr, error);
                        
                        printf("         p=%.4f %.4f %.4f %.4f %.4f %.4f\n",
                                referenceToFrame.log()[0],
                                referenceToFrame.log()[1],
                                referenceToFrame.log()[2],
                                referenceToFrame.log()[3],
                                referenceToFrame.log()[4],
                                referenceToFrame.log()[5]);

                    }
                    
                    // converged?
                    if (error / lastErr > settings.convergenceEps[lvl])
                    {

                        if (enablePrintDebugInfo && printTrackingIterationInfo)
                        {

                            printf("(%d-%d): FINISHED pyramid level "
                                   "(last residual reduction too small).\n",
                                   lvl,iteration);

                        }
                        iteration = settings.maxItsPerLvl[lvl];

                    }
                    
                    last_residual = lastErr = error;
                    
                    
                    if (LM_lambda <= 0.2)
                        LM_lambda = 0;
                    else
                        LM_lambda *= settings.lambdaSuccessFac;
                    
                    break;
                }
                else
                {

                    if (enablePrintDebugInfo && printTrackingIterationInfo)
                    {

                        printf("(%d-%d): REJECTED increment of %f with "
                               "lambda %.1f, (residual: %f -> %f)\n",
                               lvl,iteration, sqrt(inc.dot(inc)),
                               LM_lambda, lastErr, error);

                    }
                    
                    if (!(inc.dot(inc) > settings.stepSizeMin[lvl]))
                    {

                        if (enablePrintDebugInfo && printTrackingIterationInfo)
                        {

                            printf("(%d-%d): FINISHED pyramid level "
                                   "(stepsize too small).\n",
                                   lvl,iteration);

                        }
                        iteration = settings.maxItsPerLvl[lvl];
                        break;

                    }
                    
                    if(LM_lambda == 0)
                        LM_lambda = 0.2;
                    else
                        LM_lambda *= std::pow(settings.lambdaFailFac, incTry);

                }

            }

        }

    }
    
    
    if(plotTracking)
        Util::displayImage("TrackingResidual", debugImageResidualsRollingShutter, false);
    
    
    if(enablePrintDebugInfo && printTrackingIterationInfo)
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
    
    saveAllTrackingStagesInternal = false;
    
    lastResidual = last_residual;
    
    trackingWasGood = !diverged &&
            lastGoodCount / (frame->width(SE3TRACKING_MIN_LEVEL) *
                             frame->height(SE3TRACKING_MIN_LEVEL) * 0.25)
            > MIN_GOODPERALL_PIXEL &&
            lastGoodCount / (lastGoodCount + lastBadCount) >
            MIN_GOODPERGOODBAD_PIXEL;

    if (!trackingWasGood)
    {


//        Tracking was NOT good
//            lastGoodCount = 3022.000000
//            (frame->width(SE3TRACKING_MIN_LEVEL) *                frame->height(SE3TRACKING_MIN_LEVEL)) = 76800
//            MIN_GOODPERALL_PIXEL = 0.040000
//            (lastGoodCount + lastBadCount) = 3026.000000
//            MIN_GOODPERGOODBAD_PIXEL = 0.500000
//        TRACKING LOST for frame 292 (3.93% good Points, which is                99.87% of available points, NOT DIVERGED)!




        printf("Tracking was NOT good\n");
        printf("    lastGoodCount = %f\n", lastGoodCount);
        printf("    (frame->width(SE3TRACKING_MIN_LEVEL) * \
               frame->height(SE3TRACKING_MIN_LEVEL)) = %d\n",
                frame->width(SE3TRACKING_MIN_LEVEL) *
                frame->height(SE3TRACKING_MIN_LEVEL));
        printf("    MIN_GOODPERALL_PIXEL = %f\n", MIN_GOODPERALL_PIXEL);
        printf("    (lastGoodCount + lastBadCount) = %f\n",
               (lastGoodCount + lastBadCount));
        printf("    MIN_GOODPERGOODBAD_PIXEL = %f\n",
               MIN_GOODPERGOODBAD_PIXEL);

    }
    
    if (trackingWasGood)
    {

        reference->keyframe->numFramesTrackedOnThis++;

    }
    
    frame->initialTrackedResidual = lastResidual / pointUsage;
    frame->pose->thisToParent_raw =
            sim3FromSE3(toSophus(referenceToFrame.inverse()),1);

    frame->pose->trackingParent = reference->keyframe->pose;
    return toSophus(referenceToFrame.inverse());

}

float
SE3TrackerRollingShutter::calcResidualAndBuffers(
        const Eigen::Vector3f* refPoint,
        const Eigen::Vector2f* refColVar,
        int* idxBuf,
        int refNum,
        Frame* frame,
        const Sophus::SE3f& referenceToFrame,
        int level,
        bool plotResidual)
{

    calcResidualAndBuffers_debugStart();
    
    if (plotResidual)
    {

        debugImageResidualsRollingShutter.setTo(0);

    }
    
    int w = frame->width(level);
    int h = frame->height(level);
    Eigen::Matrix3f KLvl = frame->K(level);
    float fx_l = KLvl(0,0);
    float fy_l = KLvl(1,1);
    float cx_l = KLvl(0,2);
    float cy_l = KLvl(1,2);
    
    Eigen::Matrix3f rotMat = referenceToFrame.rotationMatrix();
    Eigen::Vector3f transVec = referenceToFrame.translation();
    
    const Eigen::Vector3f* refPoint_max = refPoint + refNum;
    
    const Eigen::Vector4f* frame_gradients = frame->gradients(level);
    
    int idx=0;
    
    float sumResUnweighted = 0;
    
    bool* isGoodOutBuffer = idxBuf != 0 ? frame->refPixelWasGood() : 0;
    
    int goodCount = 0;
    int badCount = 0;
    
    float sumSignedRes = 0;
    
    float sxx=0,syy=0,sx=0,sy=0,sw=0;
    
    float usageCount = 0;
    
    for (;refPoint<refPoint_max; refPoint++, refColVar++, idxBuf++)
    {
        
        Eigen::Vector3f Wxp = rotMat * (*refPoint) + transVec;
        float u_new = (Wxp[0]/Wxp[2])*fx_l + cx_l;
        float v_new = (Wxp[1]/Wxp[2])*fy_l + cy_l;
        
        // step 1a: coordinates have to be in image:
        // (inverse test to exclude NANs)
        if (!(u_new > 1 && v_new > 1 && u_new < w-2 && v_new < h-2))
        {

            if(isGoodOutBuffer != 0)
            {

                isGoodOutBuffer[*idxBuf] = false;

            }
            continue;

        }
        
        // resInterp[0] : interpolated gradient x
        // resInterp[1] : interpolated gradient y
        // resInterp[2] : intensity/colour value
        Eigen::Vector3f resInterp =
                getInterpolatedElement43(frame_gradients, u_new, v_new, w);

        
        float c1 = affineEstimation_a * (*refColVar)[0] + affineEstimation_b;
        float c2 = resInterp[2]; // Interpolated colour value
        

        //----------------------------------------------------------------------
        // Compute residual depending on the row of the warped pixel points

        // Get distorted image point coordinates to check within a scanline
        float u_new_tmp = u_new * pow(2, level);
        float v_new_tmp = v_new * pow(2, level);
        float dist_u_new;
        float dist_v_new;
        undistorter->distortPoint(u_new_tmp, v_new_tmp,
                                  &dist_u_new, &dist_v_new);

        // Check if dist_v_new is in the scanline determined by the time
        // For testing, here a certain size of strip is chosen
        float residual = 0.0;
        if ((dist_v_new >= undistorter->getInputHeight()*0.5) &&
                (dist_v_new <= undistorter->getInputHeight()*0.7))
        {

            residual = c1 - c2;

        }
        else
        {
            // Warped pixel is not in the region of rows, so do not consider
            continue;
        
        }

        //----------------------------------------------------------------------
        
        float weight = fabsf(residual) < 5.0f ? 1 : 5.0f / fabsf(residual);
        sxx += c1*c1*weight;
        syy += c2*c2*weight;
        sx += c1*weight;
        sy += c2*weight;
        sw += weight;
        
        bool isGood = residual*residual /
                (MAX_DIFF_CONSTANT +
                 MAX_DIFF_GRAD_MULT*(resInterp[0]*resInterp[0] +
                resInterp[1]*resInterp[1])) < 1;
        
        if (isGoodOutBuffer != 0)
        {

            isGoodOutBuffer[*idxBuf] = isGood;

        }
        
        *(buf_warped_x+idx) = Wxp(0);
        *(buf_warped_y+idx) = Wxp(1);
        *(buf_warped_z+idx) = Wxp(2);
        
        *(buf_warped_dx+idx) = fx_l * resInterp[0];
        *(buf_warped_dy+idx) = fy_l * resInterp[1];
        *(buf_warped_residual+idx) = residual;
        
        *(buf_d+idx) = 1.0f / (*refPoint)[2];
        *(buf_idepthVar+idx) = (*refColVar)[1];
        idx++;
        
        
        if (isGood)
        {

            sumResUnweighted += residual*residual;
            sumSignedRes += residual;
            goodCount++;

        }
        else
        {

            badCount++;

        }
        
        // If depth becomes larger: pixel becomes "smaller",
        // hence count it less.
        float depthChange = (*refPoint)[2] / Wxp[2];
        usageCount += depthChange < 1 ? depthChange : 1;
        
        
        // DEBUG STUFF
        if (plotTrackingIterationInfo || plotResidual)
        {

            // for debug plot only: find x,y again.
            // horribly inefficient, but who cares at this point...
            Eigen::Vector3f point = KLvl * (*refPoint);
            int x = point[0] / point[2] + 0.5f;
            int y = point[1] / point[2] + 0.5f;
            
            if(plotTrackingIterationInfo || consoleMode)
            {

                setPixelInCvMat(&debugImageOldImageSourceRollingShutter,
                        getGrayCvPixel((float)resInterp[2]),
                        u_new+0.5,v_new+0.5,(width/w));
                setPixelInCvMat(&debugImageOldImageWarpedRollingShutter,
                        getGrayCvPixel((float)resInterp[2]),
                        x,y,(width/w));

            }

            if (isGood)
            {

                setPixelInCvMat(&debugImageResidualsRollingShutter,
                                getGrayCvPixel(residual+128),x,y,(width/w));

            }
            else
            {

                setPixelInCvMat(&debugImageResidualsRollingShutter,
                                cv::Vec3b(0,0,255),x,y,(width/w));
            }

        }

    }
    
    buf_warped_size = idx;
    
    pointUsage = usageCount / (float)refNum;
    lastGoodCount = goodCount;
    lastBadCount = badCount;
    lastMeanRes = sumSignedRes / goodCount;

    affineEstimation_a_lastIt = sqrtf((syy - sy*sy/sw) / (sxx - sx*sx/sw));
    affineEstimation_b_lastIt = (sy - affineEstimation_a_lastIt*sx)/sw;

    calcResidualAndBuffers_debugFinish(w);

    return sumResUnweighted / goodCount;
}

void SE3TrackerRollingShutter::calcResidualAndBuffers_debugStart()
{
    if(plotTrackingIterationInfo || saveAllTrackingStagesInternal || consoleMode)
    {

        int other = saveAllTrackingStagesInternal ? 255 : 0;
        fillCvMat(&debugImageResidualsRollingShutter,cv::Vec3b(other,other,255));
        fillCvMat(&debugImageWeightsRollingShutter,cv::Vec3b(other,other,255));
        fillCvMat(&debugImageOldImageSourceRollingShutter,cv::Vec3b(other,other,255));
        fillCvMat(&debugImageOldImageWarpedRollingShutter,cv::Vec3b(other,other,255));

        fillCvMat(&debugGradientXImage0_RollingShutter,cv::Vec3b(other,other,255));
        fillCvMat(&debugGradientXImage1_RollingShutter,cv::Vec3b(other,other,255));

        fillCvMat(&debugRefImage0,cv::Vec3b(other,other,255));
        fillCvMat(&debugRefImage1,cv::Vec3b(other,other,255));
        fillCvMat(&debugRefImage2,cv::Vec3b(other,other,255));

    }
}

void SE3TrackerRollingShutter::calcResidualAndBuffers_debugFinish(int w, Frame *frame)
{

    if(plotTrackingIterationInfo)
    {
        Util::displayImage( "[RS]Weights", debugImageWeightsRollingShutter );
        Util::displayImage( "[RS]second_frame", debugImageSecondFrameRollingShutter );
        Util::displayImage( "[RS]Intensities of second_frame at transformed positions",
                            debugImageOldImageSourceRollingShutter );
        Util::displayImage( "[RS]Intensities of second_frame at pointcloud in first_frame",
                            debugImageOldImageWarpedRollingShutter );
        Util::displayImage( "[RS]Residuals", debugImageResidualsRollingShutter );

        // Text
        char buf1[200];
        char buf2[200];
        snprintf(buf1,200,"FrameID %d", frame->id());
        snprintf(buf2,200,"ReferenceID %d", frame->referenceID);
        printMessageOnCVImage(debugImageSecondFrameRollingShutter, buf1, buf2);

        // wait for key and handle it
        bool looping = true;
        while(looping)
        {
            int k = Util::waitKey(1);
            if(k == -1)
            {
                if(autoRunWithinFrame)
                    break;
                else
                    continue;
            }

            char key = k;
            if(key == ' ')
                looping = false;
            else
                handleKey(k);
        }
    }

    if(saveAllTrackingStagesInternal)
    {
        char charbuf[500];

        snprintf(charbuf,500,"save/%sresidual-%d-%d.png",packagePath.c_str(),w,iterationNumber);
        cv::imwrite(charbuf,debugImageResidualsRollingShutter);

        snprintf(charbuf,500,"save/%swarped-%d-%d.png",packagePath.c_str(),w,iterationNumber);
        cv::imwrite(charbuf,debugImageOldImageWarpedRollingShutter);

        snprintf(charbuf,500,"save/%sweights-%d-%d.png",packagePath.c_str(),w,iterationNumber);
        cv::imwrite(charbuf,debugImageWeightsRollingShutter);

        printf("saved three images for lvl %d, iteration %d\n",w,iterationNumber);
    }

}

void SE3TrackerRollingShutter::calcResidualAndBuffers_debugFinish(int w)
{
    if(plotTrackingIterationInfo)
    {
        Util::displayImage( "[RS]Weights", debugImageWeightsRollingShutter );
        Util::displayImage( "[RS]second_frame", debugImageSecondFrameRollingShutter );
        Util::displayImage( "[RS]Intensities of second_frame at transformed positions",
                            debugImageOldImageSourceRollingShutter );
        Util::displayImage( "[RS]Intensities of second_frame at pointcloud in first_frame",
                            debugImageOldImageWarpedRollingShutter );
        Util::displayImage( "[RS]Residuals", debugImageResidualsRollingShutter );

        Util::displayImage( "[RS]GradientX Image0",
                            debugGradientXImage0_RollingShutter );
        Util::displayImage( "[RS]GradientX Image1",
                            debugGradientXImage1_RollingShutter );

        Util::displayImage("[RS]TrackingReference 0", debugRefImage0);
        Util::displayImage("[RS]TrackingReference 1", debugRefImage1);
        Util::displayImage("[RS]TrackingReference 2", debugRefImage2);

        // Move windows

        // Top-left corner = [70, 30]
        if (isWindowMoved == false)
        {

            isWindowMoved = true;
            cv::moveWindow("[RS]second_frame",
                           650*0 + 70, 520*2 + 30);
            cv::moveWindow("[RS]Intensities of second_frame "
                           "at pointcloud in first_frame",
                           650*1 + 70, 520*1 + 30);
            cv::moveWindow("[RS]Intensities of second_frame "
                           "at transformed positions",
                           650*1 + 70, 520*2 + 30);
            cv::moveWindow("[RS]Residuals", 650*2 + 70, 520*1 + 30);
            cv::moveWindow("[RS]GradientX Image0", 650*3 + 70, 520*0 + 30);
            cv::moveWindow("[RS]GradientX Image1", 650*3 + 70, 520*1 + 30);

        }

        // wait for key and handle it
        bool looping = true;
        while(looping)
        {
            int k = Util::waitKey(1);
            if(k == -1)
            {
                if(autoRunWithinFrame)
                    break;
                else
                    continue;
            }

            char key = k;
            if(key == ' ')
                looping = false;
            else
                handleKey(k);
        }
    }

    if(saveAllTrackingStagesInternal)
    {
        char charbuf[500];

        snprintf(charbuf,500,"save/%sresidual-%d-%d.png",packagePath.c_str(),w,iterationNumber);
        cv::imwrite(charbuf,debugImageResidualsRollingShutter);

        snprintf(charbuf,500,"save/%swarped-%d-%d.png",packagePath.c_str(),w,iterationNumber);
        cv::imwrite(charbuf,debugImageOldImageWarpedRollingShutter);

        snprintf(charbuf,500,"save/%sweights-%d-%d.png",packagePath.c_str(),w,iterationNumber);
        cv::imwrite(charbuf,debugImageWeightsRollingShutter);

        printf("saved three images for lvl %d, iteration %d\n",w,iterationNumber);
    }
}


SE3 SE3TrackerRollingShutter::trackFrameByRows(
        TrackingReference* reference,
        Frame* frame,
        const SE3& frameToReference_initialEstimate,
        const float startRow,
        const float endRow)
{

    boost::shared_lock<boost::shared_mutex> lock = frame->getActiveLock();
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

        const float* frameImage = frame->image();
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
    Sophus::SE3f referenceToFrame =
            frameToReference_initialEstimate.inverse().cast<float>();
    NormalEquationsLeastSquares ls;

    int numCalcResidualCalls[PYRAMID_LEVELS];
    int numCalcWarpUpdateCalls[PYRAMID_LEVELS];

    float last_residual = 0;

    for (int lvl=SE3TRACKING_MAX_LEVEL-1; lvl >= SE3TRACKING_MIN_LEVEL; lvl--)
    {

        numCalcResidualCalls[lvl] = 0;
        numCalcWarpUpdateCalls[lvl] = 0;

        reference->makePointCloud(lvl);

        callOptimized(calcResidualAndBuffersByRows,
                        (reference->posData[lvl],
                         reference->colorAndVarData[lvl],
                         SE3TRACKING_MIN_LEVEL == lvl ?
                            reference->pointPosInXYGrid[lvl] : 0,
                         reference->numData[lvl],
                         frame,
                         referenceToFrame,
                         lvl,
                         startRow,
                         endRow,
                         (plotTracking && lvl == SE3TRACKING_MIN_LEVEL)));

        if (buf_warped_size <
                MIN_GOODPERALL_PIXEL_ABSMIN * (width>>lvl)*(height>>lvl)*0.25)
        {

            diverged = true;
            trackingWasGood = false;
            return SE3();

        }

        if (useAffineLightningEstimation)
        {

            affineEstimation_a = affineEstimation_a_lastIt;
            affineEstimation_b = affineEstimation_b_lastIt;

        }
        float lastErr =
                callOptimized(calcWeightsAndResidual,(referenceToFrame));

        numCalcResidualCalls[lvl]++;


        float LM_lambda = settings.lambdaInitial[lvl];

        for (int iteration=0; iteration < settings.maxItsPerLvl[lvl];
             iteration++)
        {

            callOptimized(calculateWarpUpdate,(ls));

            numCalcWarpUpdateCalls[lvl]++;

            iterationNumber = iteration;

            int incTry=0;
            while (true)
            {

                // solve LS system with current lambda
                Vector6 b = -ls.b;
                Matrix6x6 A = ls.A;
                for(int i=0;i<6;i++) A(i,i) *= 1+LM_lambda;
                Vector6 inc = A.ldlt().solve(b);
                incTry++;

                // apply increment. pretty sure this way round is correct,
                // but hard to test.
                Sophus::SE3f new_referenceToFrame =
                        Sophus::SE3f::exp((inc)) * referenceToFrame;
                //Sophus::SE3f new_referenceToFrame =
                //      referenceToFrame * Sophus::SE3f::exp((inc));

                // re-evaluate residual
                callOptimized(calcResidualAndBuffersByRows,
                              (reference->posData[lvl],
                               reference->colorAndVarData[lvl],
                               SE3TRACKING_MIN_LEVEL == lvl ?
                                   reference->pointPosInXYGrid[lvl] : 0,
                               reference->numData[lvl],
                               frame,
                               new_referenceToFrame,
                               lvl,
                               startRow,
                               endRow,
                               (plotTracking && lvl == SE3TRACKING_MIN_LEVEL)));

                if (buf_warped_size < MIN_GOODPERALL_PIXEL_ABSMIN *
                        (width>>lvl)*(height>>lvl)*(endRow - startRow))
                {

                    diverged = true;
                    trackingWasGood = false;
                    return SE3();

                }

                float error = callOptimized(calcWeightsAndResidual,
                                            (new_referenceToFrame));
                numCalcResidualCalls[lvl]++;


                // accept inc?
                if( error < lastErr)
                {

                    // accept inc
                    referenceToFrame = new_referenceToFrame;
                    if(useAffineLightningEstimation)
                    {

                        affineEstimation_a = affineEstimation_a_lastIt;
                        affineEstimation_b = affineEstimation_b_lastIt;

                    }


                    if (enablePrintDebugInfo && printTrackingIterationInfo)
                    {

                        // debug output
                        printf("(%d-%d): ACCEPTED increment of "
                               "%f with lambda %.1f, residual: %f -> %f\n",
                               lvl,iteration, sqrt(inc.dot(inc)),
                               LM_lambda, lastErr, error);

                        printf("         p=%.4f %.4f %.4f %.4f %.4f %.4f\n",
                                referenceToFrame.log()[0],
                                referenceToFrame.log()[1],
                                referenceToFrame.log()[2],
                                referenceToFrame.log()[3],
                                referenceToFrame.log()[4],
                                referenceToFrame.log()[5]);

                    }

                    // converged?
                    if (error / lastErr > settings.convergenceEps[lvl])
                    {

                        if (enablePrintDebugInfo && printTrackingIterationInfo)
                        {

                            printf("(%d-%d): FINISHED pyramid level "
                                   "(last residual reduction too small).\n",
                                   lvl,iteration);

                        }
                        iteration = settings.maxItsPerLvl[lvl];

                    }

                    last_residual = lastErr = error;


                    if (LM_lambda <= 0.2)
                    {

                        LM_lambda = 0;

                    }
                    else
                    {

                        LM_lambda *= settings.lambdaSuccessFac;

                    }

                    break;

                }
                else
                {

                    if (enablePrintDebugInfo && printTrackingIterationInfo)
                    {

                        printf("(%d-%d): REJECTED increment of %f with "
                               "lambda %.1f, (residual: %f -> %f)\n",
                               lvl,iteration, sqrt(inc.dot(inc)),
                               LM_lambda, lastErr, error);

                    }

                    if (!(inc.dot(inc) > settings.stepSizeMin[lvl]))
                    {

                        if (enablePrintDebugInfo && printTrackingIterationInfo)
                        {

                            printf("(%d-%d): FINISHED pyramid level "
                                   "(stepsize too small).\n",
                                   lvl,iteration);

                        }
                        iteration = settings.maxItsPerLvl[lvl];
                        break;

                    }

                    if (LM_lambda == 0)
                    {

                        LM_lambda = 0.2;

                    }
                    else
                    {

                        LM_lambda *= std::pow(settings.lambdaFailFac, incTry);

                    }

                }

            } // end of while

        } // end of for (iterations)

    } // end of for (pyramid levels)


    if (plotTracking)
        Util::displayImage("TrackingResidual", debugImageResidualsRollingShutter, false);


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

    saveAllTrackingStagesInternal = false;

    lastResidual = last_residual;

    trackingWasGood = !diverged &&
            lastGoodCount / (frame->width(SE3TRACKING_MIN_LEVEL) *
                             frame->height(SE3TRACKING_MIN_LEVEL) *
                             (endRow - startRow)*0.4)
            > MIN_GOODPERALL_PIXEL &&
            lastGoodCount / (lastGoodCount + lastBadCount) >
            MIN_GOODPERGOODBAD_PIXEL;

    if (!trackingWasGood)
    {

        printf("Tracking was NOT good\n");
        printf("    lastGoodCount = %f\n", lastGoodCount);
        printf("    (frame->width(SE3TRACKING_MIN_LEVEL) * "
               "frame->height(SE3TRACKING_MIN_LEVEL)) = %d\n",
                frame->width(SE3TRACKING_MIN_LEVEL) *
                frame->height(SE3TRACKING_MIN_LEVEL));
        printf("    (endRow - startRow) = %f\n", endRow - startRow);
        printf("    MIN_GOODPERALL_PIXEL = %f\n", MIN_GOODPERALL_PIXEL);
        printf("    (lastGoodCount + lastBadCount) = %f\n",
               (lastGoodCount + lastBadCount));
        printf("    MIN_GOODPERGOODBAD_PIXEL = %f\n",
               MIN_GOODPERGOODBAD_PIXEL);

    }

    if (trackingWasGood)
    {

        reference->keyframe->numFramesTrackedOnThis++;

    }

    frame->initialTrackedResidual = lastResidual / pointUsage;
    frame->pose->thisToParent_raw =
            sim3FromSE3(toSophus(referenceToFrame.inverse()),1);

    frame->pose->trackingParent = reference->keyframe->pose;
    return toSophus(referenceToFrame.inverse());

}

float
SE3TrackerRollingShutter::calcResidualAndBuffersByRows(const Eigen::Vector3f* refPoint,
        const Eigen::Vector2f* refColVar,
        int* idxBuf,
        int refNum,
        Frame* frame,
        const Sophus::SE3f& referenceToFrame,
        int level,
        const float startRow,
        const float endRow,
        bool plotResidual)
{

    calcResidualAndBuffers_debugStart();

    if (plotResidual)
    {

        debugImageResidualsRollingShutter.setTo(0);

    }

    int w = frame->width(level);
    int h = frame->height(level);
    Eigen::Matrix3f KLvl = frame->K(level);
    float fx_l = KLvl(0,0);
    float fy_l = KLvl(1,1);
    float cx_l = KLvl(0,2);
    float cy_l = KLvl(1,2);

    Eigen::Matrix3f rotMat = referenceToFrame.rotationMatrix();
    Eigen::Vector3f transVec = referenceToFrame.translation();

    const Eigen::Vector3f* refPoint_max = refPoint + refNum;

    const Eigen::Vector4f* frame_gradients = frame->gradients(level);

    int idx=0;

    float sumResUnweighted = 0;

    bool* isGoodOutBuffer = idxBuf != 0 ? frame->refPixelWasGood() : 0;

    int goodCount = 0;
    int badCount = 0;

    float sumSignedRes = 0;

    float sxx=0,syy=0,sx=0,sy=0,sw=0;

    float usageCount = 0;

    for (;refPoint<refPoint_max; refPoint++, refColVar++, idxBuf++)
    {

        Eigen::Vector3f Wxp = rotMat * (*refPoint) + transVec;
        float u_new = (Wxp[0]/Wxp[2])*fx_l + cx_l;
        float v_new = (Wxp[1]/Wxp[2])*fy_l + cy_l;

        // step 1a: coordinates have to be in image:
        // (inverse test to exclude NANs)
        if (!(u_new > 1 && v_new > 1 && u_new < w-2 && v_new < h-2))
        {

            if(isGoodOutBuffer != 0)
            {

                isGoodOutBuffer[*idxBuf] = false;

            }
            continue;

        }

        // resInterp[0] : interpolated gradient x
        // resInterp[1] : interpolated gradient y
        // resInterp[2] : intensity/colour value
        Eigen::Vector3f resInterp =
                getInterpolatedElement43(frame_gradients, u_new, v_new, w);


        float c1 = affineEstimation_a * (*refColVar)[0] + affineEstimation_b;
        float c2 = resInterp[2]; // Interpolated colour value


        //----------------------------------------------------------------------
        // Compute residual depending on the row of the warped pixel points

        // Get distorted image point coordinates to check within a scanline
        float u_new_tmp = u_new * pow(2, level);
        float v_new_tmp = v_new * pow(2, level);
        float dist_u_new;
        float dist_v_new;
        undistorter->distortPoint(u_new_tmp, v_new_tmp,
                                  &dist_u_new, &dist_v_new);


//        // Plot distorted points on a debug window
//        if (plotTrackingIterationInfo || plotResidual)
//        {

//            setPixelInCvMat(&debugImageOriginalInput,
//                getGrayCvPixel((float)c2),
//                dist_u_new, dist_v_new, 1);

//        }

        // Check if dist_v_new is in the scanline determined by the time
        // For testing, here a certain size of strip is chosen
        float residual = 0.0;
        float startRowInOriginalImage = undistorter->getInputHeight()*startRow;
        float endRowInOriginalImage = undistorter->getInputHeight()*endRow;
        if ((dist_v_new >= startRowInOriginalImage) &&
                (dist_v_new <= endRowInOriginalImage))
        {

            residual = c1 - c2;

        }
        else
        {
            // Warped pixel is not in the region of rows, so do not consider
            continue;

        }

        //----------------------------------------------------------------------

        float weight = fabsf(residual) < 5.0f ? 1 : 5.0f / fabsf(residual);
        sxx += c1*c1*weight;
        syy += c2*c2*weight;
        sx += c1*weight;
        sy += c2*weight;
        sw += weight;

        bool isGood = residual*residual /
                (MAX_DIFF_CONSTANT +
                 MAX_DIFF_GRAD_MULT*(resInterp[0]*resInterp[0] +
                resInterp[1]*resInterp[1])) < 1;

        if (isGoodOutBuffer != 0)
        {

            isGoodOutBuffer[*idxBuf] = isGood;

        }

        *(buf_warped_x+idx) = Wxp(0);
        *(buf_warped_y+idx) = Wxp(1);
        *(buf_warped_z+idx) = Wxp(2);

        *(buf_warped_dx+idx) = fx_l * resInterp[0];
        *(buf_warped_dy+idx) = fy_l * resInterp[1];
        *(buf_warped_residual+idx) = residual;

        *(buf_d+idx) = 1.0f / (*refPoint)[2];
        *(buf_idepthVar+idx) = (*refColVar)[1];
        idx++;


        if (isGood)
        {

            sumResUnweighted += residual*residual;
            sumSignedRes += residual;
            goodCount++;

        }
        else
        {

            badCount++;

        }

        // If depth becomes larger: pixel becomes "smaller",
        // hence count it less.
        float depthChange = (*refPoint)[2] / Wxp[2];
        usageCount += depthChange < 1 ? depthChange : 1;


        // DEBUG STUFF
        if (plotTrackingIterationInfo || plotResidual)
        {

            // for debug plot only: find x,y again.
            // horribly inefficient, but who cares at this point...
            Eigen::Vector3f point = KLvl * (*refPoint);
            int x = point[0] / point[2] + 0.5f;
            int y = point[1] / point[2] + 0.5f;

            if(plotTrackingIterationInfo || consoleMode)
            {

                setPixelInCvMat(&debugImageOldImageSourceRollingShutter,
                        getGrayCvPixel((float)resInterp[2]),
                        u_new+0.5,v_new+0.5,(width/w));
                setPixelInCvMat(&debugImageOldImageWarpedRollingShutter,
                        getGrayCvPixel((float)resInterp[2]),
                        x,y,(width/w));

            }

            if (isGood)
            {

                setPixelInCvMat(&debugImageResidualsRollingShutter,
                                getGrayCvPixel(residual+128),x,y,(width/w));

            }
            else
            {

                setPixelInCvMat(&debugImageResidualsRollingShutter,
                                cv::Vec3b(0,0,255),x,y,(width/w));
            }

        }

    }

    buf_warped_size = idx;

    pointUsage = usageCount / (float)refNum;
    lastGoodCount = goodCount;
    lastBadCount = badCount;
    lastMeanRes = sumSignedRes / goodCount;

    affineEstimation_a_lastIt = sqrtf((syy - sy*sy/sw) / (sxx - sx*sx/sw));
    affineEstimation_b_lastIt = (sy - affineEstimation_a_lastIt*sx)/sw;

    calcResidualAndBuffers_debugFinish(w);

    return sumResUnweighted / goodCount;
}

#if defined(ENABLE_SSE)
float SE3TrackerRollingShutter::calcResidualAndBuffersByRowsSSE(
        const Eigen::Vector3f* refPoint,
        const Eigen::Vector2f* refColVar,
        int* idxBuf,
        int refNum,
        Frame* frame,
        const Sophus::SE3f& referenceToFrame,
        int level,
        const float startRow,
        const float endRow,
        bool plotResidual)
{

    return calcResidualAndBuffersByRows(refPoint, refColVar, idxBuf, refNum,
                                        frame, referenceToFrame, level,
                                        startRow, endRow, plotResidual);

}
#endif

#if 1//WarpResidual related
//------------------------------------------------------------------------------
// Track frame as original SE3 tracker that uses the whole size of input image
// However, this version uses Ceres for optimization part
// This implementation is purpose of testing Ceres optimization for the
// integration with LSD-SLAM
SE3 SE3TrackerRollingShutter::trackFrame_Ceres(
        TrackingReference* reference,
        Frame* frame,
        const SE3& frameToReference_initialEstimate,
        const float startRow,
        const float endRow)
{

    boost::shared_lock<boost::shared_mutex> lock = frame->getActiveLock();
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

        const float* frameImage = frame->image();
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
    Sophus::SE3f referenceToFrame =
            frameToReference_initialEstimate.inverse().cast<float>();

    int numCalcResidualCalls[PYRAMID_LEVELS];
    int numCalcWarpUpdateCalls[PYRAMID_LEVELS];
    //float last_residual = 0;

#define USE_CERES_SPLINE 1
#if 1//USE_CERES

    bool isOptimized = optimizeUsingCeres(reference,
                                          frame,
                                          &referenceToFrame,
                                          startRow,
                                          endRow);


#elif USE_CERES_SPLINE

//    // Compute a transformation of reference with respect to the world
//    Sophus::SE3d refToWorld =
//            se3FromSim3(reference->keyframe->pose->getCamToWorld());
//    Sophus::SE3d frameToWorld =
//            se3FromSim3(frame->pose->getCamToWorld());

    // Set all initial 4 control poses with the reference to frame
    std::vector<Sophus::SE3d> controlPoses;
    for (int i=0; i<4; i++)
    {

        controlPoses.push_back(referenceToFrame.cast<double>());

    }

    printf("To process keyframe %d and frame %d...\n",
           reference->frameID, frame->id());

    // Run optimization to refine the 4 control poses
    bool isOptimized = optimizeUsingCeresBSpline(reference,
                                                frame,
                                                &referenceToFrame,
                                                controlPoses);

#else

    bool isOptimized =
            optimzeByNormalEquationLeastSquares(numCalcResidualCalls,
                                        numCalcWarpUpdateCalls,
                                        &lastResidual,
                                        reference,
                                        frame,
                                        &referenceToFrame,
                                        startRow,
                                        endRow);

#endif


    if (isOptimized == false)
    {
        diverged = true;
        trackingWasGood = false;

        printf("Optimization failed.. DIVERGED.\n");

        return SE3();

    }



    if (plotTracking)
        Util::displayImage("TrackingResidual", debugImageResidualsRollingShutter, false);


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

    trackingWasGood = !diverged &&
            lastGoodCount / (frame->width(SE3TRACKING_MIN_LEVEL) *
                             frame->height(SE3TRACKING_MIN_LEVEL) *
                             (endRow - startRow))
            > MIN_GOODPERALL_PIXEL &&
            lastGoodCount / (lastGoodCount + lastBadCount) >
            MIN_GOODPERGOODBAD_PIXEL;

    if (!trackingWasGood)
    {

        printf("Tracking was NOT good\n");
        printf("    lastGoodCount = %f\n", lastGoodCount);
        printf("    (frame->width(SE3TRACKING_MIN_LEVEL) * "
               "frame->height(SE3TRACKING_MIN_LEVEL)) = %d\n",
                frame->width(SE3TRACKING_MIN_LEVEL) *
                frame->height(SE3TRACKING_MIN_LEVEL));
        printf("    (endRow - startRow) = %f\n", endRow - startRow);
        printf("    MIN_GOODPERALL_PIXEL = %f\n", MIN_GOODPERALL_PIXEL);
        printf("    (lastGoodCount + lastBadCount) = %f\n",
               (lastGoodCount + lastBadCount));
        printf("    MIN_GOODPERGOODBAD_PIXEL = %f\n",
               MIN_GOODPERGOODBAD_PIXEL);

    }

    if (trackingWasGood)
    {

        reference->keyframe->numFramesTrackedOnThis++;

    }

    frame->initialTrackedResidual = lastResidual / pointUsage;
    frame->pose->thisToParent_raw =
            sim3FromSE3(toSophus(referenceToFrame.inverse()),1);

    frame->pose->trackingParent = reference->keyframe->pose;
    return toSophus(referenceToFrame.inverse());

}
#endif
bool
SE3TrackerRollingShutter::optimzeByNormalEquationLeastSquares(
        int *numCalcResidualCalls,
        int *numCalcWarpUpdateCalls,
        float *last_residual,
        TrackingReference* reference,
        Frame* frame,
        Sophus::SE3f *referenceToFrame,
        float startRow,
        float endRow)
{

    NormalEquationsLeastSquares ls;

//    float last_residual = 0;

    for (int lvl=SE3TRACKING_MAX_LEVEL-1; lvl >= SE3TRACKING_MIN_LEVEL; lvl--)
    {

        numCalcResidualCalls[lvl] = 0;
        numCalcWarpUpdateCalls[lvl] = 0;

        reference->makePointCloud(lvl);

        callOptimized(calcResidualAndBuffersByRows,
                        (reference->posData[lvl],
                         reference->colorAndVarData[lvl],
                         SE3TRACKING_MIN_LEVEL == lvl ?
                            reference->pointPosInXYGrid[lvl] : 0,
                         reference->numData[lvl],
                         frame,
                         *referenceToFrame,
                         lvl,
                         startRow,
                         endRow,
                         (plotTracking && lvl == SE3TRACKING_MIN_LEVEL)));

        if (buf_warped_size <
                MIN_GOODPERALL_PIXEL_ABSMIN * (width>>lvl)*(height>>lvl)*0.25)
        {

            diverged = true;
            trackingWasGood = false;
            return false;

        }

        if (useAffineLightningEstimation)
        {

            affineEstimation_a = affineEstimation_a_lastIt;
            affineEstimation_b = affineEstimation_b_lastIt;

        }
        float lastErr =
                callOptimized(calcWeightsAndResidual,(*referenceToFrame));

        numCalcResidualCalls[lvl]++;


        float LM_lambda = settings.lambdaInitial[lvl];

        for (int iteration=0; iteration < settings.maxItsPerLvl[lvl];
             iteration++)
        {

            callOptimized(calculateWarpUpdate,(ls));

            numCalcWarpUpdateCalls[lvl]++;

            iterationNumber = iteration;

            int incTry=0;
            while (true)
            {

                // solve LS system with current lambda
                Vector6 b = -ls.b;
                Matrix6x6 A = ls.A;
                for(int i=0;i<6;i++) A(i,i) *= 1+LM_lambda;
                Vector6 inc = A.ldlt().solve(b);
                incTry++;

                // apply increment. pretty sure this way round is correct,
                // but hard to test.
                Sophus::SE3f new_referenceToFrame =
                        Sophus::SE3f::exp((inc))*(*referenceToFrame);
                //Sophus::SE3f new_referenceToFrame =
                //      referenceToFrame * Sophus::SE3f::exp((inc));

                // re-evaluate residual
                callOptimized(calcResidualAndBuffersByRows,
                              (reference->posData[lvl],
                               reference->colorAndVarData[lvl],
                               SE3TRACKING_MIN_LEVEL == lvl ?
                                   reference->pointPosInXYGrid[lvl] : 0,
                               reference->numData[lvl],
                               frame,
                               new_referenceToFrame,
                               lvl,
                               startRow,
                               endRow,
                               (plotTracking && lvl == SE3TRACKING_MIN_LEVEL)));

                if (buf_warped_size < MIN_GOODPERALL_PIXEL_ABSMIN *
                        (width>>lvl)*(height>>lvl)*(endRow - startRow))
                {

                    diverged = true;
                    trackingWasGood = false;
                    return false;

                }

                float error = callOptimized(calcWeightsAndResidual,
                                            (new_referenceToFrame));
                numCalcResidualCalls[lvl]++;


                // accept inc?
                if( error < lastErr)
                {

                    // accept inc
                    *referenceToFrame = new_referenceToFrame;
                    if(useAffineLightningEstimation)
                    {

                        affineEstimation_a = affineEstimation_a_lastIt;
                        affineEstimation_b = affineEstimation_b_lastIt;

                    }


                    if (enablePrintDebugInfo && printTrackingIterationInfo)
                    {

                        // debug output
                        printf("(%d-%d): ACCEPTED increment of "
                               "%f with lambda %.1f, residual: %f -> %f\n",
                               lvl,iteration, sqrt(inc.dot(inc)),
                               LM_lambda, lastErr, error);

                        printf("         p=%.4f %.4f %.4f %.4f %.4f %.4f\n",
                                referenceToFrame->log()[0],
                                referenceToFrame->log()[1],
                                referenceToFrame->log()[2],
                                referenceToFrame->log()[3],
                                referenceToFrame->log()[4],
                                referenceToFrame->log()[5]);

                    }

                    // converged?
                    if (error / lastErr > settings.convergenceEps[lvl])
                    {

                        if (enablePrintDebugInfo && printTrackingIterationInfo)
                        {

                            printf("(%d-%d): FINISHED pyramid level "
                                   "(last residual reduction too small).\n",
                                   lvl,iteration);

                        }
                        iteration = settings.maxItsPerLvl[lvl];

                    }

                    *last_residual = lastErr = error;


                    if (LM_lambda <= 0.2)
                    {

                        LM_lambda = 0;

                    }
                    else
                    {

                        LM_lambda *= settings.lambdaSuccessFac;

                    }

                    break;

                }
                else
                {

                    if (enablePrintDebugInfo && printTrackingIterationInfo)
                    {

                        printf("(%d-%d): REJECTED increment of %f with "
                               "lambda %.1f, (residual: %f -> %f)\n",
                               lvl,iteration, sqrt(inc.dot(inc)),
                               LM_lambda, lastErr, error);

                    }

                    if (!(inc.dot(inc) > settings.stepSizeMin[lvl]))
                    {

                        if (enablePrintDebugInfo && printTrackingIterationInfo)
                        {

                            printf("(%d-%d): FINISHED pyramid level "
                                   "(stepsize too small).\n",
                                   lvl,iteration);

                        }
                        iteration = settings.maxItsPerLvl[lvl];
                        break;

                    }

                    if (LM_lambda == 0)
                    {

                        LM_lambda = 0.2;

                    }
                    else
                    {

                        LM_lambda *= std::pow(settings.lambdaFailFac, incTry);

                    }

                }

            } // end of while

        } // end of for (iterations)

    } // end of for (pyramid levels)


    double poseOut[6];
    for (int i=0; i<6; i++)
    {

        poseOut[i] = referenceToFrame->log()[i];

    }
    printf("[Frame %d] Finale pose = %f %f %f %f %f %f\n",
           frame->id(),
           poseOut[0], poseOut[1], poseOut[2], poseOut[3], poseOut[4], poseOut[5]);

    return true;

}

#if 1//WarpResidual related
//------------------------------------------------------------------------------
// Optimize SE(3) tracker using Ceres-solver. Yet another implementation
// of SE(3) tracker in LSD-SLAM using Ceres-solver's automatic differentiation
//
bool
SE3TrackerRollingShutter::optimizeUsingCeres(
        TrackingReference* reference,
        Frame* frame,
        Sophus::SE3f *referenceToFrame,
        float startRow,
        float endRow)
{

    // DEBUG

    //--------------------------------------------------------------------------
    // Set the initial pose 6-vector in Lie algebra elements of SE(3)
    double pose[6];
    for (int i=0; i<6; i++)
    {

        pose[i] = referenceToFrame->log()[i];

    }

#if DEBUG_MODE
    printf("--------------------------------------\n");
    printf("\n\nInitial pose = %f %f %f %f %f %f\n",
           pose[0], pose[1], pose[2], pose[3], pose[4], pose[5]);
#endif

    //--------------------------------------------------------------------------
    // For each pyramid level, image warping is computed the the result pose
    // is estimated and reused for the next pyramid level
    bool isSolved = false;
    for (int pyramidLevel = SE3TRACKING_MAX_LEVEL - 1;
         pyramidLevel >= SE3TRACKING_MIN_LEVEL;
         pyramidLevel--)
    {

        // Compute pose for image warping for the given pyramid level
        // NOTE: this function calls Ceres-solver, so calls
        //       the operator() in WarpResidual up to maximum iteration times.
        isSolved = optimizeUsingCeresPyramidLevelAt(pyramidLevel, pose,
                                         reference, frame);

        //----------------------------------------------------------------------
        // Store statistics and check the tracking was good

        // Total number of warpped pixel points from source to target image
        int totalWarpedPixels =
                WarpResidual::lastGoodCount_ +
                WarpResidual::lastBadCount_;

#if DEBUG_MODE
        printf("\n");
        printf("  TotalWarpedPixels = %d\n", totalWarpedPixels);
        printf("  Good %d, Bad %d\n",
               WarpResidual::lastGoodCount_,
               WarpResidual::lastBadCount_);

#endif

        // Store point usage
        pointUsage = WarpResidual::usageCount_ /
                (float) reference->numData[pyramidLevel];

#if DEBUG_MODE
        printf("\n\n");
        printf("  Reference frame ID = %d\n", reference->frameID);
        printf("  new image id = %d\n", frame->id());
        printf("  pointUsage = %f\n", pointUsage);
        printf("  WarpResidual::usageCount_ = %f\n", WarpResidual::usageCount_);

        printf("  reference->numData[pyramidLevel] = %f\n",
               (float)reference->numData[pyramidLevel]);
        printf("  at pyramidLevel = %d\n", pyramidLevel);
#endif

#if DEBUG_TEST_KEY_SELECTION_CRITERIA

        if ((pointUsage == 0 || pointUsage > 1.0) && reference->frameID > 0)
        {

        exit(1);

        }
#endif

        // Check tracking valid
        if (totalWarpedPixels <
                MIN_GOODPERALL_PIXEL_ABSMIN *
                (frame->width(pyramidLevel)*(frame->height(pyramidLevel))))
        {

            // Too less number of pixels is warped
            diverged = true;
            trackingWasGood = false;

            printf("\nDIVERGED -- Set divered = TRUE\n\n");

            isSolved = false;
            return isSolved;

        }

#if DEBUG_MODE
        printf("[Frame %d/Level %d] Estimated pose = %f %f %f %f %f %f\n",
               frame->id(),
               pyramidLevel,
               pose[0], pose[1], pose[2], pose[3], pose[4], pose[5]);
#endif

        //----------------------------------------------------------------------

    } // end of pyramid

    //--------------------------------------------------------------------------
    // Store counters

    lastGoodCount = WarpResidual::lastGoodCount_;
    lastBadCount = WarpResidual::lastBadCount_;

#if DEBUG_MODE
    printf("lastGood = %f <== %d\n",
           lastGoodCount, WarpResidual::lastGoodCount_);
    printf("lasBad = %f <=== %d\n",
           lastBadCount, WarpResidual::lastBadCount_);
#endif

    //--------------------------------------------------------------------------
    // Copy the final estimated pose to referenceToFrame
    Vector6 poseOut;
    for (int i=0; i<6; i++)
    {

        poseOut[i] = pose[i];

    }
    *referenceToFrame = Sophus::SE3f::exp(poseOut);

#if DEBUG_MODE
    printf("[Frame %d] Finale pose = %f %f %f %f %f %f\n",
           frame->id(),
           poseOut[0], poseOut[1], poseOut[2], poseOut[3], poseOut[4], poseOut[5]);
#endif

    // Return if the optimization is solved or not
    return isSolved;

}
#endif
#if 1 // WarpResidual related
bool
SE3TrackerRollingShutter::optimizeUsingCeresPyramidLevelAt(
        int level,
        double *pose,
        TrackingReference* reference,
        Frame *frame)
{

    // Debug start
    calcResidualAndBuffers_debugStart();

    if (lsd_slam::plotTracking)
    {

        debugImageResidualsRollingShutter.setTo(0);

    }

    reference->makePointCloud(level);

    //--------------------------------------------------------------------------
    // Set the optimization using Ceres-solver
    //--------------------------------------------------------------------------
    Problem problem;

    // For each pixel in the reference frame (keyframe), a residual block
    // computing the error of image warping is computed
    // NOTE: the image is from pyramid, so their size varies on each level

    // Pointer to first semi-dense point for the given pyramid level
    // It consists of (x, y, z) for x-coordinate, y-coordinate and
    // z inverse depth.
    Eigen::Vector3f* refPoint = reference->posData[level];

    // Pointer to last semi-dense point in the reference image
    Eigen::Vector3f* refPoint_max = refPoint + reference->numData[level];

    // Pointer to first semi-dense point containing colour and variance
    // information in the given pyramid level
    Eigen::Vector2f* refColourVariance = reference->colorAndVarData[level];

    // For each available semi-dense point, a residual block is built
    for (; refPoint<refPoint_max; refPoint++, refColourVariance++)
    {

        // Values 0 ... 1 indicating start and end of rows
        // E.g. For entire image, use startRow = 0 and endRow = 1.0;
        float startRow = 0.0;
        float endRow = 1.0;

        // Check if this semi-dense point has a chance to be warped at
        // the out of target image boundary
        //bool isOutOfImage =
        //        checkSemiDensePointWarpBeOutOfImage(refPoint,
        //                                            reference, frame, pose);
        bool isOutOfImage = false;

        // If they are mapped within the target image
        if (isOutOfImage == false)
        {

            // Add a residual block for each semi-dense pixel
            // computing warping error and returning the estimate of pose
            problem.AddResidualBlock(
                new AutoDiffCostFunction<WarpResidual, 1, 6>
                (new WarpResidual(level, refPoint, refColourVariance,
                reference, frame,
                undistorter, startRow, endRow)),
                new ceres::HuberLoss(0.001),
                pose);

        }

    }

    //--------------------------------------------------------------------------
    // Solve using Ceres-solver

    // Set options for solver
    ceres::Solver::Options options;

    // Set a maximum number of iterations for each pyramid level
    //const int maxIterations[6] = {2, 4, 10, 20, 20, 20};
    options.max_num_iterations = settings.maxItsPerLvl[level];
    //options.max_num_iterations = 25;

    options.linear_solver_type = ceres::DENSE_SCHUR;
    //options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    //options.linear_solver_type = ceres::CGNR;

    options.minimizer_progress_to_stdout = true;

    // Number of threads for Jacobian computation
    options.num_threads = 16;

    // Add user callback examining variables at every iteration end
    UserCallback callback(level, reference, frame, pose,
                          &debugImageWeightsRollingShutter,
                          &debugImageResidualsRollingShutter,
                          &debugImageSecondFrameRollingShutter,
                          &debugImageOldImageWarpedRollingShutter,
                          &debugImageOldImageSourceRollingShutter);
    options.callbacks.push_back(&callback);
    options.update_state_every_iteration = true;

    // Summary
    ceres::Solver::Summary summary;

#if DEBUG_MODE
    printf("    cesres:Solving..\n");
#endif

    // Solve
    ceres::Solve(options, &problem, &summary);
    //std::cout << summary.FullReport() << "\n";

    // Update counters using callback
    //UserCallback callback2(level, reference, frame, pose);
    //callback2.getStatisticsForWarpedResidual<double>(
    callback.getStatisticsForWarpedResidual<double>(
                &(lsd_slam::WarpResidual::lastGoodCount_),
                &(lsd_slam::WarpResidual::lastBadCount_),
                &(lsd_slam::WarpResidual::usageCount_));

#if DEBUG_MODE
    printf("\n\n");
    printf("    ceres::Solve done\n");
    printf("    goodCount = %d\n",
           lsd_slam::WarpResidual::goodWarpPixelCountInPyramidImage_);
    printf("    badCount = %d\n",
           lsd_slam::WarpResidual::badWarpPixelCountInPyramidImage_);
    printf("    lastGood %d, lastBad %d\n",
           lsd_slam::WarpResidual::lastGoodCount_,
           lsd_slam::WarpResidual::lastBadCount_);
    printf("    divered = %d\n", lsd_slam::SE3TrackerRollingShutter::diverged);
    printf("    trackingWasGood = %d\n\n",
           lsd_slam::SE3TrackerRollingShutter::trackingWasGood);
#endif

    // Debug finish
    calcResidualAndBuffers_debugFinish(frame->width(level), frame);

    //--------------------------------------------------------------------------
    // Store the final value of cost function
    lastResidual = summary.final_cost / (float) summary.num_residual_blocks;

#if DEBUG_MODE
    if (summary.termination_type == ceres::CONVERGENCE)
    {

        printf("    ceres::CONVERGENCE\n");

    }
    else
    {

        printf("    ceres::NON-CONVERGENCE\n");

    }
#endif

    if (summary.IsSolutionUsable() == true)
    {

        return true;

    }
    else
    {
        return false;

    }

}
#endif


} // end of namespace
