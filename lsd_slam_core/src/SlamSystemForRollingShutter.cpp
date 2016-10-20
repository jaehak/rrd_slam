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

#include "SlamSystemForRollingShutter.h"

#include "DataStructures/Frame.h"
//#include "Tracking/SE3Tracker.h"
#include "Tracking/SE3TrackerRollingShutter.h"
#include "Tracking/Sim3Tracker.h"
#include "DepthEstimation/DepthMap.h"
#include "Tracking/TrackingReference.h"
#include "util/globalFuncs.h"
#include "GlobalMapping/KeyFrameGraph.h"
#include "GlobalMapping/TrackableKeyFrameSearch.h"
#include "GlobalMapping/g2oTypeSim3Sophus.h"
#include "IOWrapper/ImageDisplay.h"
#include "IOWrapper/Output3DWrapper.h"
#include <g2o/core/robust_kernel_impl.h>
#include "DataStructures/FrameMemory.h"
#include "deque"

// for mkdir
#include <sys/types.h>
#include <sys/stat.h>

#ifdef ANDROID
#include <android/log.h>
#endif

#include "opencv2/opencv.hpp"

using namespace lsd_slam;

SlamSystemForRollingShutter::
SlamSystemForRollingShutter(int width_output,
                            int height_output,
                            Eigen::Matrix3f K,
                            Undistorter *undistorter,
                            int width_bsplineSegImg,
                            int height_bsplineSegImg,
                            int rowPerSegment,
                            int methodRollingShutter,
                            bool enableSLAM,
                            bool useGTmotion)
: SlamSystem(width_output, height_output, K, enableSLAM)
{
    
    this->rowPerSegment = rowPerSegment;

    // Set B-spline segment image size
    this->width_bsplineSegImg = width_bsplineSegImg;
    this->height_bsplineSegImg = height_bsplineSegImg;

    // Use a rolling shutter version of SE3 tracker
    tracker = new SE3TrackerRollingShutter(width_output,
                                           height_output,
                                           K,
                                           undistorter,
                                           width_bsplineSegImg,
                                           height_bsplineSegImg,
                                           methodRollingShutter);

    this->undistorter = undistorter;

#if 1// TURN_OFF SIM3 TRACKER
    if(SLAMEnabled)
    {
        constraintSE3Tracker =
                new SE3TrackerRollingShutter(width_output,
                                             height_output,
                                             K,
                                             undistorter,
                                             width_bsplineSegImg,
                                             height_bsplineSegImg,
                                             methodRollingShutter);
    }
#endif
    
    // Reallocate a map since the SlamSystem constructor already allocated
    // a map with DepthMap class
    //thread_mapping.join();
    

    // Depth-map
    switch (methodRollingShutter)
    {
            
        case 0: // Rolling-shutter
            map = new DepthMapForRollingShutter(width_output, height_output, K, useGTmotion);
            break;
            
        default: // (Default 1) Radial-Rolling-shutter
            map = new DepthMapForRadialRolling(width_output, height_output, K,
                                               undistorter->getOutputDistort(),
                                               useGTmotion);
            break;
            
    }

    //thread_mapping = boost::thread(&SlamSystemForRollingShutter::mappingThreadLoop, this);

}

SlamSystemForRollingShutter::~SlamSystemForRollingShutter()
{

    delete tracker;

}

void
SlamSystemForRollingShutter::trackFrame_GS(
                            uchar* image,
                            cv::Mat imageMask,
                            unsigned int frameID,
                            bool blockUntilMapped,
                            double timestamp,
                            int height_rollingShutter_inputImage,
                            std::string imageFilename,
                            std::string depthFilename)
{
    
    // Create new frame
    std::shared_ptr<Frame> trackingNewFrame(
                new Frame(frameID, width, height, K,
                          timestamp,
                          image, imageMask.data,
                          height_rollingShutter_inputImage));

    // Set info for this new frame
    trackingNewFrame.get()->imageFilename = imageFilename;
    trackingNewFrame.get()->depthFilename = depthFilename;
    
    // Set internal spline for every tracking new frame
    trackingNewFrame.get()->setSpline(tracker->internalSpline);
    trackingNewFrame.get()->setUndistorter(this->undistorter);
    
    if(!trackingIsGood)
    {
        relocalizer.updateCurrentFrame(trackingNewFrame);
        
        unmappedTrackedFramesMutex.lock();
        unmappedTrackedFramesSignal.notify_one();
        unmappedTrackedFramesMutex.unlock();
        return;
    }
    
    currentKeyFrameMutex.lock();
    bool my_createNewKeyframe = createNewKeyFrame;	// pre-save here, to make decision afterwards.
    if(trackingReference->keyframe != currentKeyFrame.get() || currentKeyFrame->depthHasBeenUpdatedFlag)
    {
        trackingReference->importFrame(currentKeyFrame.get());
        currentKeyFrame->depthHasBeenUpdatedFlag = false;
        trackingReferenceFrameSharedPT = currentKeyFrame;
    }
    
    FramePoseStruct* trackingReferencePose = trackingReference->keyframe->pose;
    currentKeyFrameMutex.unlock();
    
    // DO TRACKING & Show tracking result.
    if(enablePrintDebugInfo && printThreadingInfo)
        printf("TRACKING %d on %d\n", trackingNewFrame->id(), trackingReferencePose->frameID);
    
    
    poseConsistencyMutex.lock_shared();
    SE3 frameToReference_initialEstimate = se3FromSim3(
                                                       trackingReferencePose->getCamToWorld().inverse() * keyFrameGraph->allFramePoses.back()->getCamToWorld());
    poseConsistencyMutex.unlock_shared();
    
#if 0//DEBUG_frameToReference_Info
    
    printf("frameToReference_initialEstimate=\n");
    std::cout << frameToReference_initialEstimate.matrix() << std::endl;
    printf("RefToWorld=\n");
    std::cout << se3FromSim3(trackingReferencePose->getCamToWorld()).matrix()
    << std::endl;
    printf("FrameToWorld=\n");
    std::cout << keyFrameGraph->allFramePoses.back()->getCamToWorld().matrix()
    << std::endl;
    
#endif
    struct timeval tv_start, tv_end;
    gettimeofday(&tv_start, NULL);
    
    SE3 newRefToFrame_poseUpdate = tracker->trackFrame(
                                                       trackingReference,
                                                       trackingNewFrame.get(),
                                                       frameToReference_initialEstimate);
    
    
    gettimeofday(&tv_end, NULL);
    msTrackFrame = 0.9*msTrackFrame + 0.1*((tv_end.tv_sec-tv_start.tv_sec)*1000.0f + (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
    nTrackFrame++;
    
    tracking_lastResidual = tracker->lastResidual;
    tracking_lastUsage = tracker->pointUsage;
    tracking_lastGoodPerBad = tracker->lastGoodCount / (tracker->lastGoodCount + tracker->lastBadCount);
    tracking_lastGoodPerTotal = tracker->lastGoodCount / (trackingNewFrame->width(SE3TRACKING_MIN_LEVEL)*trackingNewFrame->height(SE3TRACKING_MIN_LEVEL));
    
    
    if(manualTrackingLossIndicated || tracker->diverged || (keyFrameGraph->keyframesAll.size() > INITIALIZATION_PHASE_COUNT && !tracker->trackingWasGood))
    {
        printf("TRACKING LOST for frame %d (%1.2f%% good Points, which is %1.2f%% of available points, %s)!\n",
               trackingNewFrame->id(),
               100*tracking_lastGoodPerTotal,
               100*tracking_lastGoodPerBad,
               tracker->diverged ? "DIVERGED" : "NOT DIVERGED");
        
        trackingReference->invalidate();
        
        trackingIsGood = false;
        nextRelocIdx = -1;
        
        unmappedTrackedFramesMutex.lock();
        unmappedTrackedFramesSignal.notify_one();
        unmappedTrackedFramesMutex.unlock();
        
        manualTrackingLossIndicated = false;
        return;
    }
    
    
    
    if(plotTracking)
    {
        Eigen::Matrix<float, 20, 1> data;
        data.setZero();
        data[0] = tracker->lastResidual;
        
        data[3] = tracker->lastGoodCount / (tracker->lastGoodCount + tracker->lastBadCount);
        data[4] = 4*tracker->lastGoodCount / (width*height);
        data[5] = tracker->pointUsage;
        
        data[6] = tracker->affineEstimation_a;
        data[7] = tracker->affineEstimation_b;
        outputWrapper->publishDebugInfo(data);
    }
    
    keyFrameGraph->addFrame(trackingNewFrame.get());
    
    
    //Sim3 lastTrackedCamToWorld = mostCurrentTrackedFrame->getScaledCamToWorld();//  mostCurrentTrackedFrame->TrackingParent->getScaledCamToWorld() * sim3FromSE3(mostCurrentTrackedFrame->thisToParent_SE3TrackingResult, 1.0);
    if (outputWrapper != 0)
    {
        outputWrapper->publishTrackedFrame(trackingNewFrame.get());
    }
    
    
    // Keyframe selection
    latestTrackedFrame = trackingNewFrame;
    if (!my_createNewKeyframe && currentKeyFrame->numMappedOnThisTotal > MIN_NUM_MAPPED)
    {
        Sophus::Vector3d dist = newRefToFrame_poseUpdate.translation() * currentKeyFrame->meanIdepth;
        float minVal = fmin(0.2f + keyFrameGraph->keyframesAll.size() * 0.8f / INITIALIZATION_PHASE_COUNT,1.0f);
        
        if(keyFrameGraph->keyframesAll.size() < INITIALIZATION_PHASE_COUNT)	minVal *= 0.7;
        
        lastTrackingClosenessScore = trackableKeyFrameSearch->getRefFrameScore(dist.dot(dist), tracker->pointUsage);
        
        if (lastTrackingClosenessScore > minVal)
        {
            createNewKeyFrame = true;
            
            if(enablePrintDebugInfo && printKeyframeSelectionInfo)
                printf("SELECT %d on %d! dist %.3f + usage %.3f = %.3f > 1\n",trackingNewFrame->id(),trackingNewFrame->getTrackingParent()->id(), dist.dot(dist), tracker->pointUsage, trackableKeyFrameSearch->getRefFrameScore(dist.dot(dist), tracker->pointUsage));
        }
        else
        {
            if(enablePrintDebugInfo && printKeyframeSelectionInfo)
                printf("SKIPPD %d on %d! dist %.3f + usage %.3f = %.3f > 1\n",trackingNewFrame->id(),trackingNewFrame->getTrackingParent()->id(), dist.dot(dist), tracker->pointUsage, trackableKeyFrameSearch->getRefFrameScore(dist.dot(dist), tracker->pointUsage));
            
        }
    }
    
    
    unmappedTrackedFramesMutex.lock();
    if(unmappedTrackedFrames.size() < 50 || (unmappedTrackedFrames.size() < 100 && trackingNewFrame->getTrackingParent()->numMappedOnThisTotal < 10))
        unmappedTrackedFrames.push_back(trackingNewFrame);
    unmappedTrackedFramesSignal.notify_one();
    unmappedTrackedFramesMutex.unlock();
    
    // implement blocking
    if(blockUntilMapped && trackingIsGood)
    {
        boost::unique_lock<boost::mutex> lock(newFrameMappedMutex);
        while(unmappedTrackedFrames.size() > 0)
        {
            //printf("TRACKING IS BLOCKING, waiting for %d frames to finish mapping.\n", (int)unmappedTrackedFrames.size());
            newFrameMappedSignal.wait(lock);
        }
        lock.unlock();
    }
}

bool SlamSystemForRollingShutter::doMappingIteration()
{

    if (currentKeyFrame == 0)
    {

        return false;

    }

    if (!doMapping && currentKeyFrame->idxInKeyframes < 0)
    {

        if (currentKeyFrame->numMappedOnThisTotal >= MIN_NUM_MAPPED)
        {

            finishCurrentKeyframe();

        }
        else
        {

            discardCurrentKeyframe();

        }

        map->invalidate();
        printf("Finished KF %d as Mapping got disabled!\n",
               currentKeyFrame->id());

        // No create and forced
        changeKeyframe(true, true, 1.0f);

    }

#if 0 // TURN OFF MERGE OPTIMIZATION OFFSET (Jae-Hak Kim)
#else // TURN ON
    mergeOptimizationOffset();
#endif
    addTimingSamples();

    if (dumpMap)
    {

        keyFrameGraph->dumpMap(packagePath+"/save");
        dumpMap = false;

    }


    // set mappingFrame
    if (trackingIsGood)
    {

        if (!doMapping)
        {

            //printf("tryToChange refframe, lastScore %f!\n",
            //       lastTrackingClosenessScore);
            if (lastTrackingClosenessScore > 1)
            {

                // No creation but not forced
                changeKeyframe(true, false, lastTrackingClosenessScore * 0.75);

            }

            if (displayDepthMap || depthMapScreenshotFlag)
            {

                debugDisplayDepthMap();

            }

            return false;

        }


        if (createNewKeyFrame)
        {

            finishCurrentKeyframe();
            // Creation and forced
            changeKeyframe(false, true, 1.0f);


            if (displayDepthMap || depthMapScreenshotFlag)
            {

                debugDisplayDepthMap();

            }

        }
        else
        {

            bool didSomething = updateKeyframe();

            if (displayDepthMap || depthMapScreenshotFlag)
            {

                debugDisplayDepthMap();

            }

            if (!didSomething)
            {

                return false;

            }

        }

        return true;

    }
    else
    {
        // invalidate map if it was valid.
        if (map->isValid())
        {

            if(currentKeyFrame->numMappedOnThisTotal >= MIN_NUM_MAPPED)
            {

                finishCurrentKeyframe();

            }
            else
            {

                discardCurrentKeyframe();

            }

            map->invalidate();

        }

#if 0// TURN OFF RELOCALIZER (Jae-Hak Kim, 2015-04-17)
#else// TURN ON RELOCALIZER
        // start relocalizer if it isnt running already
        if (!relocalizer.isRunning)
        {

            relocalizer.start(keyFrameGraph->keyframesAll);

        }

        // did we find a frame to relocalize with?
        if (relocalizer.waitResult(50))
        {

            takeRelocalizeResult();

        }
#endif
        return true;

    }

}

void SlamSystemForRollingShutter::gtDepthInit(uchar* image,
                                              uchar* imageMask,
                                              float* depth,
                                              double timeStamp,
                                              int id)
{
    printf("Doing GT initialization!\n");

    currentKeyFrameMutex.lock();

    // Update for rolling shutter
    Frame *newFrame = new Frame(id, width, height, K, timeStamp, image, imageMask);
    newFrame->startRowImage = 0;
    newFrame->startRowSeg = 0;
    newFrame->rowPerSegment = this->rowPerSegment;
    newFrame->numExtraCtrlPts_ = SPLINE_K - 1;
    newFrame->setUndistorter(undistorter);

    currentKeyFrame.reset(newFrame);
    currentKeyFrame->setDepthFromGroundTruth(depth);

    map->initializeFromGTDepth(currentKeyFrame.get());
    keyFrameGraph->addFrame(currentKeyFrame.get());
    
    //--------------------------------------------------------------------------
    // Make mean of inverse depth be one.
    map->makeMeanInvDepthOne();

    currentKeyFrameMutex.unlock();

    if(doSlam)
    {
        keyFrameGraph->idToKeyFrameMutex.lock();
        keyFrameGraph->idToKeyFrame.insert(std::make_pair(currentKeyFrame->id(), currentKeyFrame));
        keyFrameGraph->idToKeyFrameMutex.unlock();
    }
    if(continuousPCOutput && outputWrapper != 0) outputWrapper->publishKeyframe(currentKeyFrame.get());

    printf("Done GT initialization!\n");
}

// Groud truth depth with an initial noise
void SlamSystemForRollingShutter::gtDepthInitWithNoise(uchar* image, float* depth,
                                              double timeStamp, int id)
{
    printf("Doing GT initialization!\n");
    
    currentKeyFrameMutex.lock();
    
    // Update for rolling shutter
    Frame *newFrame = new Frame(id, width, height, K, timeStamp, image);
    newFrame->startRowImage = 0;
    newFrame->startRowSeg = 0;
    newFrame->rowPerSegment = this->rowPerSegment;
    newFrame->numExtraCtrlPts_ = SPLINE_K - 1;
    newFrame->setUndistorter(undistorter);
    
    currentKeyFrame.reset(newFrame);
    currentKeyFrame->setDepthFromGroundTruth(depth);
    
    map->initializeFromGTDepthWithNoise(currentKeyFrame.get());
    keyFrameGraph->addFrame(currentKeyFrame.get());
    
    currentKeyFrameMutex.unlock();
    
    if(doSlam)
    {
        keyFrameGraph->idToKeyFrameMutex.lock();
        keyFrameGraph->idToKeyFrame.insert(std::make_pair(currentKeyFrame->id(), currentKeyFrame));
        keyFrameGraph->idToKeyFrameMutex.unlock();
    }
    if(continuousPCOutput && outputWrapper != 0) outputWrapper->publishKeyframe(currentKeyFrame.get());
    
    printf("Done GT initialization!\n");
}


void SlamSystemForRollingShutter::gtDepthInitWithName(
          uchar* image, float* depth,
          std::string depthFilename,
          std::string motionFilename,
          double timeStamp, int id)
{
    printf("Doing GT initialization!\n");
    
    currentKeyFrameMutex.lock();
    
    // Update for rolling shutter
    Frame *newFrame = new Frame(id, width, height, K, timeStamp, image);
    newFrame->startRowImage = 0;
    newFrame->startRowSeg = 0;
    newFrame->rowPerSegment = this->rowPerSegment;
    newFrame->depthFilename = depthFilename;
    newFrame->motionFilename = motionFilename;
    newFrame->setGTMotion();
    newFrame->numExtraCtrlPts_ = SPLINE_K - 1;
    newFrame->setUndistorter(undistorter);
    newFrame->setSpline(tracker->internalSpline);
    
    currentKeyFrame.reset(newFrame);
    currentKeyFrame->setDepthFromGroundTruth(depth);
    
    map->initializeFromGTDepth(currentKeyFrame.get());
    keyFrameGraph->addFrame(currentKeyFrame.get());
    
    currentKeyFrameMutex.unlock();
    
    if(doSlam)
    {
        keyFrameGraph->idToKeyFrameMutex.lock();
        keyFrameGraph->idToKeyFrame.insert(std::make_pair(currentKeyFrame->id(), currentKeyFrame));
        keyFrameGraph->idToKeyFrameMutex.unlock();
    }
    if(continuousPCOutput && outputWrapper != 0) outputWrapper->publishKeyframe(currentKeyFrame.get());
    
    printf("Done GT initialization!\n");
}


bool SlamSystemForRollingShutter::updateKeyframe()
{

    std::shared_ptr<Frame> reference = nullptr;
    std::deque< std::shared_ptr<Frame> > references;

    unmappedTrackedFramesMutex.lock();

#if 0 // THIS IS NOT NECESSARY FOR MULTIPLE KEYFRAMES IN DEQUE
    
    // remove frames that have a different tracking parent.
    while(unmappedTrackedFrames.size() > 0 &&
          (!unmappedTrackedFrames.front()->hasTrackingParent() ||
           unmappedTrackedFrames.front()->getTrackingParent() !=
           currentKeyFrame.get()))
    {

        unmappedTrackedFrames.front()->clear_refPixelWasGood();
        unmappedTrackedFrames.pop_front();

    }
    
#else
    
    // Remove frames that has no tracking parent.
    // Remove a frame that is identical to the current keyframe
    while(unmappedTrackedFrames.size() > 0 &&
          (!unmappedTrackedFrames.front()->hasTrackingParent() ||
           unmappedTrackedFrames.front()->id() ==
           currentKeyFrame.get()->id()))
    {
        
        unmappedTrackedFrames.front()->clear_refPixelWasGood();
        unmappedTrackedFrames.pop_front();
        
    }
    
#endif

    // clone list
    if (unmappedTrackedFrames.size() > 0)
    {

        for(unsigned int i=0;i<unmappedTrackedFrames.size(); i++)
        {

            references.push_back(unmappedTrackedFrames[i]);

        }

        std::shared_ptr<Frame> popped = unmappedTrackedFrames.front();
        unmappedTrackedFrames.pop_front();
        unmappedTrackedFramesMutex.unlock();

        if (enablePrintDebugInfo && printThreadingInfo)
        {

            printf("MAPPING %d on %d to %d (%d frames)\n",
                   currentKeyFrame->id(), references.front()->id(),
                   references.back()->id(), (int)references.size());

        }
#if 0 // TURN OFF UPDATING DEPTH IN KEYFRAME
#else
        // TURN ON UPDATING DEPTH IN KEYFRAME (Default)
        map->updateKeyframe(references);
#endif
        popped->clear_refPixelWasGood();
        references.clear();

    }
    else
    {

        unmappedTrackedFramesMutex.unlock();
        return false;

    }

    if (enablePrintDebugInfo && printRegularizeStatistics)
    {

        Eigen::Matrix<float, 20, 1> data;
        data.setZero();
        data[0] = runningStats.num_reg_created;
        data[2] = runningStats.num_reg_smeared;
        data[3] = runningStats.num_reg_deleted_secondary;
        data[4] = runningStats.num_reg_deleted_occluded;
        data[5] = runningStats.num_reg_blacklisted;

        data[6] = runningStats.num_observe_created;
        data[7] = runningStats.num_observe_create_attempted;
        data[8] = runningStats.num_observe_updated;
        data[9] = runningStats.num_observe_update_attempted;


        data[10] = runningStats.num_observe_good;
        data[11] = runningStats.num_observe_inconsistent;
        data[12] = runningStats.num_observe_notfound;
        data[13] = runningStats.num_observe_skip_oob;
        data[14] = runningStats.num_observe_skip_fail;

        outputWrapper->publishDebugInfo(data);

    }

    if (outputWrapper != 0 && continuousPCOutput && currentKeyFrame != 0)
    {

        outputWrapper->publishKeyframe(currentKeyFrame.get());

    }

    return true;

}

void SlamSystemForRollingShutter::randomInitDepthRange(cv::Mat image,
                                             cv::Mat imageMask,
                                             double timeStamp, int id,
                                             float minInitDepth,
                                             float maxInitDepth)
{
    printf("Doing Random initialization!\n");

    if (!doMapping)
    {

        printf("WARNING: mapping is disabled, but we just initialized... "
               "THIS WILL NOT WORK! Set doMapping to true.\n");

    }

    currentKeyFrameMutex.lock();
    Frame *newFrame = new Frame(id, image.cols, image.rows, K,
                                timeStamp, image.data, imageMask.data);
    newFrame->startRowImage = 0;
    newFrame->startRowSeg = 0;
    newFrame->rowPerSegment = this->rowPerSegment;
    newFrame->numExtraCtrlPts_ = SPLINE_K - 1;
    newFrame->setUndistorter(undistorter);
    
    currentKeyFrame.reset(newFrame);
    
    map->initializeRandomly(currentKeyFrame.get(),
                            minInitDepth,
                            maxInitDepth);
    keyFrameGraph->addFrame(currentKeyFrame.get());

    //--------------------------------------------------------------------------
    // Make mean of inverse depth be one.
    map->makeMeanInvDepthOne();
    
    currentKeyFrameMutex.unlock();

    if (doSlam)
    {

        keyFrameGraph->idToKeyFrameMutex.lock();
        keyFrameGraph->idToKeyFrame.insert(
                    std::make_pair(currentKeyFrame->id(), currentKeyFrame));
        keyFrameGraph->idToKeyFrameMutex.unlock();

    }

    if (continuousPCOutput && outputWrapper != 0)
    {

        outputWrapper->publishKeyframe(currentKeyFrame.get());

    }

    if (displayDepthMap || depthMapScreenshotFlag)
    {

        debugDisplayDepthMap();

    }

    printf("Done Random initialization!\n");

}


void
SlamSystemForRollingShutter::
depthTransferFrom(Frame *firstKeyframe)
{
    printf("Doing depth transfer initialization from GS SLAM!\n");
    
    if (!doMapping)
    {
        
        printf("WARNING: mapping is disabled, but we just initialized... "
               "THIS WILL NOT WORK! Set doMapping to true.\n");
        
    }
    
    currentKeyFrameMutex.lock();
    
    // Copy first key frame
    Frame *newFrame = new Frame(*firstKeyframe);
    
    // Re-init newFrame except depth
    newFrame->initializeExceptDepth(firstKeyframe->id(),
                                    firstKeyframe->width(),
                                    firstKeyframe->height(),
                                    firstKeyframe->K(),
                                    firstKeyframe->timestamp(),
                                    firstKeyframe->height_rollingShutter_inputImage);
    
    // Reset
    newFrame->startRowImage = 0;
    newFrame->startRowSeg = 0;
    newFrame->rowPerSegment = this->rowPerSegment;
    newFrame->numExtraCtrlPts_ = SPLINE_K - 1;
    newFrame->setUndistorter(undistorter);
    
    currentKeyFrame.reset(newFrame);
    
    map->initializeFromKeyframe(currentKeyFrame.get()); // Transfer map
    keyFrameGraph->addFrame(currentKeyFrame.get());
    
    //--------------------------------------------------------------------------
    // Make mean of inverse depth be one.
    map->makeMeanInvDepthOne();
    
    currentKeyFrameMutex.unlock();
    
    if (doSlam)
    {
        
        keyFrameGraph->idToKeyFrameMutex.lock();
        keyFrameGraph->idToKeyFrame.insert(
            std::make_pair(currentKeyFrame->id(), currentKeyFrame));
        keyFrameGraph->idToKeyFrameMutex.unlock();
        
    }
    
    if (continuousPCOutput && outputWrapper != 0)
    {
        
        outputWrapper->publishKeyframe(currentKeyFrame.get());
        
    }
    
    if (displayDepthMap || depthMapScreenshotFlag)
    {
        
        debugDisplayDepthMap();
        
    }
    
    printf("Done depth transfer initialization!\n");
    
}



void SlamSystemForRollingShutter::trackFrameRS(cv::Mat image,
                                             unsigned int frameID,
                                             bool blockUntilMapped,
                                             double timestamp)
{
    // Create new frame
    std::shared_ptr<Frame> trackingNewFrame(
                new Frame(frameID, image.cols, image.rows, K,
                          timestamp, image.data));

    if (!trackingIsGood)
    {

        relocalizer.updateCurrentFrame(trackingNewFrame);

        unmappedTrackedFramesMutex.lock();
        unmappedTrackedFramesSignal.notify_one();
        unmappedTrackedFramesMutex.unlock();
        return;

    }

    currentKeyFrameMutex.lock();
    // pre-save here, to make decision afterwards.
    bool my_createNewKeyframe = createNewKeyFrame;
    if (trackingReference->keyframe != currentKeyFrame.get() ||
            currentKeyFrame->depthHasBeenUpdatedFlag)
    {

        trackingReference->importFrame(currentKeyFrame.get());
        currentKeyFrame->depthHasBeenUpdatedFlag = false;
        trackingReferenceFrameSharedPT = currentKeyFrame;

    }

    FramePoseStruct* trackingReferencePose = trackingReference->keyframe->pose;
    currentKeyFrameMutex.unlock();

    // DO TRACKING & Show tracking result.
    if (enablePrintDebugInfo && printThreadingInfo)
    {

        printf("TRACKING %d on %d\n", trackingNewFrame->id(),
               trackingReferencePose->frameID);

    }


    poseConsistencyMutex.lock_shared();
    SE3 frameToReference_initialEstimate =
            se3FromSim3(trackingReferencePose->getCamToWorld()).inverse() *
            se3FromSim3(keyFrameGraph->allFramePoses.back()->getCamToWorld());
    poseConsistencyMutex.unlock_shared();



    struct timeval tv_start, tv_end;
    gettimeofday(&tv_start, NULL);

    // Testing the original using Ceres
    SE3 newRefToFrame_poseUpdate =
            tracker->trackFrame_Ceres(
                trackingReference,
                trackingNewFrame.get(),
                frameToReference_initialEstimate,
                0.0,
                1.0);

//    //--------------------------------------------------------------------------
//    // Tracking rows in new frame with respect to the keyframe

//    // Globally
//    std::shared_ptr<Frame> trackingNewFrame_copy(
//                new Frame(frameID, width, height, K, timestamp, image));
//    SE3 global_frameToReference_initialEstimate =
//            tracker->trackFrameByRows(
//                trackingReference,
//                trackingNewFrame_copy.get(),
//                frameToReference_initialEstimate,
//                0.0,
//                1.0);

//    // Update the initial estimate
//    frameToReference_initialEstimate =
//            global_frameToReference_initialEstimate;

//    // Locallly in terms of rows
//    SE3 newRefToFrame_poseUpdate;

//    // Thickness of the row band between start and end of row
//    // e.g. 0.5 means the half of the image height
//    float thicknessRowBand = 0.2;

//    // Increase rate of rows between previous and current row band
//    // If rowIncRate is less than thicknessRowBand, adjacent row band
//    // will overlap. If they are equal, there will be no overlap
//    float rowIncRate = 0.2;

//    // For each band of rows, process tracking frame by the rows
////    for (float startRow = 0.0, endRow = startRow + thicknessRowBand;
////         startRow <= 1.0 - thicknessRowBand && endRow <= 1.0;
////         startRow += rowIncRate, endRow += rowIncRate)
////    {
//        float startRow = 0.5;
//        float endRow = 0.6;
//        newRefToFrame_poseUpdate =
//                tracker->trackFrameByRows(
//                    trackingReference,
//                    trackingNewFrame.get(),
//                    frameToReference_initialEstimate,
//                    startRow,
//                    endRow);

////    }

    //--------------------------------------------------------------------------


    gettimeofday(&tv_end, NULL);
    msTrackFrame = 0.9*msTrackFrame +
            0.1*((tv_end.tv_sec-tv_start.tv_sec)*1000.0f +
                 (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
    nTrackFrame++;

    tracking_lastResidual = tracker->lastResidual;
    tracking_lastUsage = tracker->pointUsage;
    tracking_lastGoodPerBad = tracker->lastGoodCount /
            (tracker->lastGoodCount + tracker->lastBadCount);
    tracking_lastGoodPerTotal = tracker->lastGoodCount /
            (trackingNewFrame->width(SE3TRACKING_MIN_LEVEL)*
             trackingNewFrame->height(SE3TRACKING_MIN_LEVEL));


    if (manualTrackingLossIndicated || tracker->diverged ||
            (keyFrameGraph->keyframesAll.size() >
             INITIALIZATION_PHASE_COUNT && !tracker->trackingWasGood))
    {

        printf("TRACKING LOST for frame %d (%1.2f%% good Points, which is "
               "%1.2f%% of available points, %s)!\n",
                trackingNewFrame->id(),
                100*tracking_lastGoodPerTotal,
                100*tracking_lastGoodPerBad,
                tracker->diverged ? "DIVERGED" : "NOT DIVERGED");

        trackingReference->invalidate();

        trackingIsGood = false;
        nextRelocIdx = -1;

        unmappedTrackedFramesMutex.lock();
        unmappedTrackedFramesSignal.notify_one();
        unmappedTrackedFramesMutex.unlock();

        manualTrackingLossIndicated = false;
        return;

    }



    if (plotTracking)
    {

        Eigen::Matrix<float, 20, 1> data;
        data.setZero();
        data[0] = tracker->lastResidual;

        data[3] = tracker->lastGoodCount / (tracker->lastGoodCount +
                                            tracker->lastBadCount);
        data[4] = 4*tracker->lastGoodCount / (width*height);
        data[5] = tracker->pointUsage;

        data[6] = tracker->affineEstimation_a;
        data[7] = tracker->affineEstimation_b;
        outputWrapper->publishDebugInfo(data);

    }

    keyFrameGraph->addFrame(trackingNewFrame.get());


    if (outputWrapper != 0)
    {

        outputWrapper->publishTrackedFrame(trackingNewFrame.get());

    }


    // Keyframe selection
    latestTrackedFrame = trackingNewFrame;
    if (!my_createNewKeyframe && currentKeyFrame->numMappedOnThisTotal
            > MIN_NUM_MAPPED)
    {

        Sophus::Vector3d dist = newRefToFrame_poseUpdate.translation() *
                currentKeyFrame->meanIdepth;
        float minVal = fmin(0.2f + keyFrameGraph->keyframesAll.size() * 0.8f /
                            INITIALIZATION_PHASE_COUNT,1.0f);

        if(keyFrameGraph->keyframesAll.size() < INITIALIZATION_PHASE_COUNT)
            minVal *= 0.7;

        lastTrackingClosenessScore =
                trackableKeyFrameSearch->getRefFrameScore(dist.dot(dist),
                                                          tracker->pointUsage);

        if (lastTrackingClosenessScore > minVal)
        {

#if DEBUG_MODE
            printf("lastTrackingClosenessScore %f > minVal %f\n",
                   lastTrackingClosenessScore, minVal);
            printf("tracker->pointUsage = %f\n", tracker->pointUsage);
#endif
            createNewKeyFrame = true;

            if (enablePrintDebugInfo && printKeyframeSelectionInfo)
            {

                printf("SELECT %d on %d! dist %.3f + usage %.3f = %.3f > 1\n",
                       trackingNewFrame->id(),
                       trackingNewFrame->getTrackingParent()->id(),
                       dist.dot(dist),
                       tracker->pointUsage,
                       trackableKeyFrameSearch->getRefFrameScore(
                           dist.dot(dist), tracker->pointUsage));

            }

        }
        else
        {

            if (enablePrintDebugInfo && printKeyframeSelectionInfo)
            {

                printf("SKIPPED %d on %d! dist %.3f + usage %.3f = %.3f > 1\n",
                       trackingNewFrame->id(),
                       trackingNewFrame->getTrackingParent()->id(),
                       dist.dot(dist),
                       tracker->pointUsage,
                       trackableKeyFrameSearch->getRefFrameScore(
                           dist.dot(dist), tracker->pointUsage));

            }


        }

    }

    unmappedTrackedFramesMutex.lock();
    if(unmappedTrackedFrames.size() < 50 ||
            (unmappedTrackedFrames.size() < 100 &&
             trackingNewFrame->getTrackingParent()->numMappedOnThisTotal < 10))
        unmappedTrackedFrames.push_back(trackingNewFrame);
    unmappedTrackedFramesSignal.notify_one();
    unmappedTrackedFramesMutex.unlock();

    // implement blocking
    if (blockUntilMapped && trackingIsGood)
    {

        boost::unique_lock<boost::mutex> lock(newFrameMappedMutex);
        while(unmappedTrackedFrames.size() > 0)
        {

//            printf("TRACKING IS BLOCKING, waiting for %d frames to "
//                   "finish mapping.\n", (int)unmappedTrackedFrames.size());
            newFrameMappedSignal.wait(lock);

        }
        lock.unlock();

    }
}

void SlamSystemForRollingShutter::debugDisplayDepthMap()
{

    map->debugPlotDepthMap();
    double scale = 1;
    if(currentKeyFrame != 0 && currentKeyFrame != 0)
        scale = currentKeyFrame->getScaledCamToWorld().scale();
    // debug plot depthmap
    char buf1[200];
    char buf2[200];


    snprintf(buf1,200,"Map: Upd %3.0fms (%2.4fHz); Trk %3.0fms (%2.4fHz); %d / %d / %d",
            map->msUpdate, map->nAvgUpdate,
            msTrackFrame, nAvgTrackFrame,
            currentKeyFrame->numFramesTrackedOnThis, currentKeyFrame->numMappedOnThis, (int)unmappedTrackedFrames.size());

    snprintf(buf2,200,"dens %2.0f%%; good %2.0f%%; scale %2.2f; res %2.1f/; usg %2.0f%%; Map: %d F, %d KF, %d E, %.1fm Pts",
            100*currentKeyFrame->numPoints/(float)(width*height),
            100*tracking_lastGoodPerBad,
            scale,
            tracking_lastResidual,
            100*tracking_lastUsage,
            (int)keyFrameGraph->allFramePoses.size(),
            keyFrameGraph->totalVertices,
            (int)keyFrameGraph->edgesAll.size(),
            1e-6 * (float)keyFrameGraph->totalPoints);
    
//    if (keyFrameGraph->totalVertices > 0)
//    {
//        
//        fprintf(stderr, "Paused. Press any key to continue");
//        char n;
//        std::cin >> n;
//        
//    }


    if(onSceenInfoDisplay)
        printMessageOnCVImage(map->debugImageDepth, buf1, buf2);
    if ((displayDepthMap) && (consoleMode == false))
        Util::displayImage( "DebugWindow DEPTH", map->debugImageDepth, false );
    
//    // Display generalised stereo pair
//    bool displayGenStereoPair = true;
//    if (displayGenStereoPair)
//    {
//        Util::displayImage("GenStereo Src", map->debugImageSource, false);
//        Util::displayImage("GenStereo Trg", map->debugImageTarget, false);
//    }


    // Save a depth map image as a file
    if ((displayDepthMap) || (consoleMode == true))
    {

        char buf[255];
        sprintf(buf, "depth_key%04d_cur%04d.jpg",
                currentKeyFrame->id(),
                tracker->getFrameNumberToDisplay());
        if (!isFileExist(buf))
            cv::imwrite(buf, map->debugImageDepth);
  
#if 0 // SAVE GEN STEREO PAIR
        // Save general stereo pair as file
        sprintf(buf, "genstereo_src_%06d.png", tracker->getIterationNumber());
        if (!isFileExist(buf))
            cv::imwrite(buf, map->debugImageSource);
        sprintf(buf, "genstereo_trg_%06d.png", tracker->getIterationNumber());
        if (!isFileExist(buf))
            cv::imwrite(buf, map->debugImageTarget);
#endif
        
#if 1 // SAVE depth as file
        sprintf(buf, "depth_key%04d_cur%04d.txt",
                currentKeyFrame->id(),
                tracker->getFrameNumberToDisplay());
        if (!isFileExist(buf))
            map->saveDepthAsFile(buf);
#endif

    }

    int pressedKey = Util::waitKey(1);
    handleKey(pressedKey);

}

//------------------------------------------------------------------------------
// Track frame for rolling shutter
// Version BN3
void SlamSystemForRollingShutter::trackFrameRS_BN3(
                                     cv::Mat image,
                                     cv::Mat imageMask,
                                     unsigned int frameID,
                                     bool blockUntilMapped,
                                     double timestamp,
                                     int startRowImage,
                                     int startRowSeg,
                                     double rowPerSegment,
                                     double framePerSegment,
                                     int height_rollingShutter_inputImage,
                                     std::vector<std::string> imageFiles,
                                     std::vector<std::string> depthFiles,
                                     std::vector<std::string> motionFiles)
{

    // You will need to use timestamp intelligently
    // since it could have a frame drop, so do not rely on row per segment only
    // Well, or just we let assume no frame drop

    // Maximum capacity of the list above
    //double framePerSegment = rowPerSegment / image.rows;
    int numToPop = (int)framePerSegment; // Number of pop from the list
    this->maxCapList_ = (int)(framePerSegment*SPLINE_K); // For past window
    //this->maxCapList_ = (int)(framePerSegment*(SPLINE_K + 1)); // For past + 1 window

#if 0    //FORCE for DEBUG - Use only one image per segment
    this->maxCapList_ = 1;
#endif
    // Create new frame
    std::shared_ptr<Frame> trackingNewFrame(
                new Frame(frameID, image.cols, image.rows, K,
                          timestamp, image.data, imageMask.data,
                          height_rollingShutter_inputImage));
    // Set info for this new frame
    trackingNewFrame.get()->startRowImage = startRowImage;
    trackingNewFrame.get()->startRowSeg = startRowSeg;
    trackingNewFrame.get()->rowPerSegment = rowPerSegment;
    trackingNewFrame.get()->numExtraCtrlPts_ = SPLINE_K - 1;

    if (imageFiles.size() >= frameID + 1)
        trackingNewFrame.get()->imageFilename = imageFiles[frameID];
    if (depthFiles.size() >= frameID + 1)
        trackingNewFrame.get()->depthFilename = depthFiles[frameID];
    if (motionFiles.size() >= frameID + 1)
    {
        trackingNewFrame.get()->motionFilename = motionFiles[frameID];
        trackingNewFrame.get()->setGTMotion();
    }
    
#if 1 // INTERNAL_SPLINE for every frame
    
    // Set internal spline for every tracking new frame
    trackingNewFrame.get()->setSpline(tracker->internalSpline);
    trackingNewFrame.get()->setUndistorter(this->undistorter);
    
#endif

    //--------------------------------------------------------------------------
    // Push to the queue and pop if exceeds the maximum capacity of the list
    trackingNewFrameList_.push_back(trackingNewFrame);
    if ((int)trackingNewFrameList_.size() > maxCapList_)
    {

        // Pop the number of elements in the list where
        // the number of elements to be popped determined by
        // the frame number per segment (numToPop)
        for (int i=0; i<numToPop; i++)
        {
        
            trackingNewFrameList_.pop_front();
            
        }

    }
    startRowImageList_.push_back(startRowImage);
    if ((int)startRowImageList_.size() > maxCapList_)
    {

        // Pop the number of elements in the list where
        // the number of elements to be popped determined by
        // the frame number per segment (numToPop)
        for (int i=0; i<numToPop; i++)
        {
            
            startRowImageList_.pop_front();
        
        }

    }

    startRowSegList_.push_back(startRowSeg);
    if ((int)startRowSegList_.size() > maxCapList_)
    {

        // Pop the number of elements in the list where
        // the number of elements to be popped determined by
        // the frame number per segment (numToPop)
        for (int i=0; i<numToPop; i++)
        {
            
            startRowSegList_.pop_front();
            
        }

    }

    printf("trackingNewFrameList size = %d\n",
           (int)trackingNewFrameList_.size());



#if 0
    if (referenceList_.size() != 0)
    {

        referenceList_.pop_back();
        // Copy new one
        TrackingReference *refNew = new TrackingReference();
        *refNew = *trackingReference;
        printf("Copy dest one Dim_pyd:%d %d %d %d %d\n",
               refNew->numData[0],
               refNew->numData[1],
               refNew->numData[2],
               refNew->numData[3],
               refNew->numData[4]);
        printf("Copy source one Dim_pyd:%d %d %d %d %d\n",
               trackingReference->numData[0],
               trackingReference->numData[1],
               trackingReference->numData[2],
               trackingReference->numData[3],
               trackingReference->numData[4]);
        referenceList_.push_back(refNew);

    }
#endif


    //--------------------------------------------------------------------------

    if (!trackingIsGood)
    {

        relocalizer.updateCurrentFrame(trackingNewFrame);

        unmappedTrackedFramesMutex.lock();
        unmappedTrackedFramesSignal.notify_one();
        unmappedTrackedFramesMutex.unlock();
        return;

    }

    currentKeyFrameMutex.lock();
    // pre-save here, to make decision afterwards.
    bool my_createNewKeyframe = createNewKeyFrame;

    if (trackingReference->keyframe != currentKeyFrame.get() ||
            currentKeyFrame->depthHasBeenUpdatedFlag)
    {

        trackingReference->importFrame(currentKeyFrame.get());
        currentKeyFrame->depthHasBeenUpdatedFlag = false;
        trackingReferenceFrameSharedPT = currentKeyFrame;

    }


    FramePoseStruct* trackingReferencePose = trackingReference->keyframe->pose;

#if 1 // Within currentKeyFrame Mutex
    //--------------------------------------------------------------------------
    // Push to the dequeue and pop if exceeds the maximum capacity of the list
    //std::shared_ptr<TrackingReference> trackingReferencePtr(trackingReference);

    printf("Push new one Dim_pyd:%d %d %d %d %d\n",
           trackingReference->numData[0],
           trackingReference->numData[1],
           trackingReference->numData[2],
           trackingReference->numData[3],
           trackingReference->numData[4]);
    std::shared_ptr<TrackingReference> refNew(new TrackingReference());
    *(refNew.get()) = *trackingReference;
    referenceList_.push_back(refNew);
    if ((int)referenceList_.size() > maxCapList_)
    {

        for (int i=0; i<numToPop; i++)
        {
            
            // Pop the number of elements in the list where
            // the number of elements to be popped determined by
            // the frame number per segment (numToPop)
            referenceList_.pop_front();
            
        }

    }
#endif

    currentKeyFrameMutex.unlock();


    //----------------------------------------------------------------------
    // Check if all lists are in the same size
    //----------------------------------------------------------------------
    if (referenceList_.size() != trackingNewFrameList_.size())
    {
        fprintf(stderr,"CHECKSIZE: referenceList_.size: %d != "
               "trackingNewFrameList_.size: %d\n",
               (int)referenceList_.size(),
               (int)trackingNewFrameList_.size());
        exit(1);
    }
    if (referenceList_.size() != startRowImageList_.size())
    {
        fprintf(stderr,"CHECKSIZE: referenceList_.size: %d != "
               "startRowImageList_.size: %d\n",
               (int)referenceList_.size(),
               (int)startRowImageList_.size());
        exit(1);
    }
    if (referenceList_.size() != startRowSegList_.size())
    {
        fprintf(stderr,"CHECKSIZE: referenceList_.size: %d != "
               "startRowSegList_.size: %d\n",
               (int)referenceList_.size(),
               (int)startRowSegList_.size());
        exit(1);
    }

    //----------------------------------------------------------------------

#if 0// CHECK REF
    //----------------------------------------------------------------------
    // Check ref and frame images
    // Am I really using correct images?
    //----------------------------------------------------------------------
    if (referenceList_.size() > 1)
    for (int i=0; i<(int)referenceList_.size(); i++)
    {

        printf("Checking reference %d\n", referenceList_[i].get()->frameID);
        // Load an image of the reference frame ID
        cv::Mat refImg1 = cv::imread(imageFiles[referenceList_[i].get()->frameID],
                                    CV_LOAD_IMAGE_GRAYSCALE);
        // Compare with the reference's keyframe image
        float sum_val = 0.0;
        for (int j=0;j<refImg1.rows;j++)
        {
            for (int k=0;k<refImg1.cols;k++)
            {

                float val1 = (float)refImg1.at<uchar>(j,k);
                float val2 = (referenceList_[i].get()
                            ->keyframe->image(0)[j*refImg1.cols + k]);
                float val = val1 - val2;

                sum_val = sum_val + val*val;
            }
        }
        float mean_val = sum_val / (float) (refImg1.rows*refImg1.cols);
        if (sum_val != 0.0)
        {
            fprintf(stderr,"CHECKREF ref and keyframe difference: "
                   "Frame %d - track %d, Sum_val %f; Mean_val %f\n",
                   referenceList_[i].get()->frameID,
                   trackingNewFrameList_[i].get()->id(),
                   sum_val, mean_val);
            exit(1);
        }

    }
#endif

    // DO TRACKING & Show tracking result.
    if (enablePrintDebugInfo && printThreadingInfo)
    {

        printf("TRACKING %d on %d\n", trackingNewFrame->id(),
               trackingReferencePose->frameID);

    }

    //--------------------------------------------------------------------------
    // Push to the dequeue and pop if exceeds the maximum capacity of the list
#if 0
    //std::shared_ptr<TrackingReference> trackingReferencePtr(trackingReference);

    printf("Push new one Dim_pyd:%d %d %d %d %d\n",
           trackingReference->numData[0],
           trackingReference->numData[1],
           trackingReference->numData[2],
           trackingReference->numData[3],
           trackingReference->numData[4]);
    referenceList_.push_back(trackingReference);
#endif
#if 0
    if ((int)referenceList_.size() > maxCapList_)
    {

        referenceList_.pop_front();

    }
#endif

#if 0
    // If key frame changes, replace front with an old reference to front
    if (referenceList_.size() > 1)
    if (trackingReference->keyframe->id() !=
            referenceList_[referenceList_.size() - 2]->keyframe->id())
    {

        referenceList_.pop_front();
        trackingNewFrameList_.pop_front();
        startRowImageList_.pop_front();
        startRowSegList_.pop_front();

        TrackingReference *refNew = new TrackingReference();
        *refNew = *(referenceList_[referenceList_.size() - 2]);
        referenceList_.push_front(refNew);
        trackingNewFrameList_.push_front(trackingNewFrame);
        int startRowImage = startRowImageList_[startRowImageList_.size() - 1];
        startRowImageList_.push_front(startRowImage);
        int startRowSeg = startRowSegList_[startRowSegList_.size() - 1];
        startRowImageList_.push_front(startRowSeg);

    }
#endif

#if 0// POP_FRONT_ALL_BUT_LAST
    // If key frame changes, pop front all but last
    if (referenceList_.size() > 1)
    if (trackingReference->keyframe->id() !=
            referenceList_[referenceList_.size()-2]->keyframe->id())
    {
        while (trackingNewFrameList_.size() > 1)
        {
            trackingNewFrameList_.pop_front();
            startRowImageList_.pop_front();
            startRowSegList_.pop_front();
            referenceList_.pop_front();
        }
    }
#endif

    //--------------------------------------------------------------------------

    // Transformation from frame coordiante system to Reference coordinate
    // system
    // frameToRef = WorldToRef @ frameToWorld
    poseConsistencyMutex.lock_shared();
#if 0//USE FAKE WORLDTOREF only for the first time
    SE3 frameToReference_initialEstimate;
    if (frameID == 0)
    {

        trackingReference->keyframe->pose->thisToParent_raw
                = sim3FromSE3(
                    SE3::exp(Eigen::Matrix<double,6,1>(1,2,3,-1,-3,-5)),
                    1.0);
        trackingReference->keyframe->pose->trackingParent
        frameToReference_initialEstimate
            = se3FromSim3(trackingReferencePose->getCamToWorld()).inverse() *
              se3FromSim3(keyFrameGraph->allFramePoses.back()->getCamToWorld());

    }
    frameToReference_initialEstimate =
              se3FromSim3(trackingReferencePose->getCamToWorld()).inverse() *
              se3FromSim3(keyFrameGraph->allFramePoses.back()->getCamToWorld());
#else
    SE3 frameToReference_initialEstimate =
            se3FromSim3(trackingReferencePose->getCamToWorld()).inverse() *
            se3FromSim3(keyFrameGraph->allFramePoses.back()->getCamToWorld());
#endif
    poseConsistencyMutex.unlock_shared();

#if 1// DEBUG_FRAME_TO_REFERENCE
    // Print out frameToReference
    std::cout << "[REF] frameToReference_initialEstimate = ";
    Eigen::Matrix<double,6,1> tmp = frameToReference_initialEstimate.log();
    std::cout << tmp.transpose() << std::endl;


    std::cout << "[REF]     RefToWorld_SE3 = ";
    SE3 M = se3FromSim3(trackingReferencePose->getCamToWorld());
    tmp = M.log();
    std::cout << tmp.transpose() << std::endl;
    std::cout << M.matrix() << std::endl;

    std::cout << "[REF]     RefToWorld_Sim3 = ";
    std::cout << trackingReferencePose->getCamToWorld().log().transpose() << std::endl;
    std::cout << trackingReferencePose->getCamToWorld().matrix() << std::endl;

    double scale = trackingReferencePose->getCamToWorld().scale();
    std::cout << "scale = " << scale << std::endl;

    std::cout << "[REF]     LastFrameToWorld_SE3 = ";
    M = se3FromSim3(keyFrameGraph->allFramePoses.back()->getCamToWorld());
    tmp = M.log();
    std::cout << tmp.transpose() << std::endl;
    std::cout << M.matrix() << std::endl;

    std::cout << "[REF]     LastFrameToWorld_Sim3 = ";
    std::cout << keyFrameGraph->allFramePoses.back()->getCamToWorld().log().transpose()
              << std::endl;
    std::cout << keyFrameGraph->allFramePoses.back()->getCamToWorld().matrix()
              << std::endl;



#endif

    struct timeval tv_start, tv_end;
    gettimeofday(&tv_start, NULL);


#if 0 // TEST
    SE3 newRefToFrame_poseUpdate = SE3();
#else

    if ((int)trackingNewFrameList_.size() < this->maxCapList_)
    {
        
        // Not enough frames in the deque
        // Do nothing.
        return;
        
    }
    
    // Testing the original using Ceres
    struct timeval tic, toc;
    gettimeofday(&tic, NULL);
    SE3 newRefToFrame_poseUpdate =
            tracker->trackFrame_Ceres_BN3(
                trackingReference,
                &referenceList_,
                &trackingNewFrameList_,
                frameToReference_initialEstimate,
                startRowImageList_,
                startRowSegList_,
                rowPerSegment,
                framePerSegment,
                1.0,
                imageFiles);
    gettimeofday(&toc, NULL);

    fprintf(stderr, "Tracking frame %d completed (%5.4f sec)\n",
            trackingNewFrameList_.front().get()->id(),
            (double)(toc.tv_sec - tic.tv_sec) +
            (double)(toc.tv_usec - tic.tv_usec)*1e-06);
#endif

//    //--------------------------------------------------------------------------
//    // Tracking rows in new frame with respect to the keyframe

//    // Globally
//    std::shared_ptr<Frame> trackingNewFrame_copy(
//                new Frame(frameID, width, height, K, timestamp, image));
//    SE3 global_frameToReference_initialEstimate =
//            tracker->trackFrameByRows(
//                trackingReference,
//                trackingNewFrame_copy.get(),
//                frameToReference_initialEstimate,
//                0.0,
//                1.0);

//    // Update the initial estimate
//    frameToReference_initialEstimate =
//            global_frameToReference_initialEstimate;

//    // Locallly in terms of rows
//    SE3 newRefToFrame_poseUpdate;

//    // Thickness of the row band between start and end of row
//    // e.g. 0.5 means the half of the image height
//    float thicknessRowBand = 0.2;

//    // Increase rate of rows between previous and current row band
//    // If rowIncRate is less than thicknessRowBand, adjacent row band
//    // will overlap. If they are equal, there will be no overlap
//    float rowIncRate = 0.2;

//    // For each band of rows, process tracking frame by the rows
////    for (float startRow = 0.0, endRow = startRow + thicknessRowBand;
////         startRow <= 1.0 - thicknessRowBand && endRow <= 1.0;
////         startRow += rowIncRate, endRow += rowIncRate)
////    {
//        float startRow = 0.5;
//        float endRow = 0.6;
//        newRefToFrame_poseUpdate =
//                tracker->trackFrameByRows(
//                    trackingReference,
//                    trackingNewFrame.get(),
//                    frameToReference_initialEstimate,
//                    startRow,
//                    endRow);

////    }

    //--------------------------------------------------------------------------

    gettimeofday(&tv_end, NULL);
    msTrackFrame = 0.9*msTrackFrame +
            0.1*((tv_end.tv_sec-tv_start.tv_sec)*1000.0f +
                 (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
    nTrackFrame++;

    tracking_lastResidual = tracker->lastResidual;
    tracking_lastUsage = tracker->pointUsage;
    tracking_lastGoodPerBad = tracker->lastGoodCount /
            (tracker->lastGoodCount + tracker->lastBadCount);
    tracking_lastGoodPerTotal = tracker->lastGoodCount /
            (trackingNewFrame->width(SE3TRACKING_MIN_LEVEL)*
             trackingNewFrame->height(SE3TRACKING_MIN_LEVEL)*
             lsd_slam::pixelDensityForTracking);

    if (manualTrackingLossIndicated == true)
    {
        
        printf("TRACKING LOST for frame %d by MANUALLY (Key L pressed)\n",
               trackingNewFrame->id());
        
        
    }

    if (manualTrackingLossIndicated || tracker->diverged ||
            (keyFrameGraph->keyframesAll.size() >
             INITIALIZATION_PHASE_COUNT && !tracker->trackingWasGood))
    {

        printf("TRACKING LOST for frame %d (%1.2f%% good Points, which is "
               "%1.2f%% of available points, %s)!\n",
                trackingNewFrame->id(),
                100*tracking_lastGoodPerTotal,
                100*tracking_lastGoodPerBad,
                tracker->diverged ? "DIVERGED" : "NOT DIVERGED");

        trackingReference->invalidate();

        trackingIsGood = false;
        nextRelocIdx = -1;

        unmappedTrackedFramesMutex.lock();
        unmappedTrackedFramesSignal.notify_one();
        unmappedTrackedFramesMutex.unlock();

        manualTrackingLossIndicated = false;
        return;

    }



    if (plotTracking)
    {

        Eigen::Matrix<float, 20, 1> data;
        data.setZero();
        data[0] = tracker->lastResidual;

        data[3] = tracker->lastGoodCount / (tracker->lastGoodCount +
                                            tracker->lastBadCount);
        data[4] = 4*tracker->lastGoodCount / (width*height);
        data[5] = tracker->pointUsage;

        data[6] = tracker->affineEstimation_a;
        data[7] = tracker->affineEstimation_b;
        outputWrapper->publishDebugInfo(data);

    }

    //keyFrameGraph->addFrame(trackingNewFrame.get());
#if 0 // USE BACK

    keyFrameGraph->addFrame(trackingNewFrameList_.back().get());
    
#else // USE FRONT
    
//    for (int i=0;i<numToPop;i++)
//    {
//        // Add numToPop front frames into keyframe graph
//        keyFrameGraph->addFrame(trackingNewFrameList_[i].get());
//        
//    }
    
    
#endif
    printf("Adding a pose into keyFrameGraph..\n");
    std::cout << " keyFrameGraph.pose = "
              << trackingNewFrame.get()->pose->getCamToWorld().log().transpose()
              << std::endl;


    if (outputWrapper != 0)
    {

        //outputWrapper->publishTrackedFrame(trackingNewFrame.get());
#if 0 // USE BACK
        outputWrapper->publishTrackedFrame(trackingNewFrameList_.back().get());
#else // USE FRONT
        
        for (int i=0;i<numToPop;i++)
        {
            
            // Add numToPop front frames into outputWrapper
            outputWrapper->publishTrackedFrame(trackingNewFrameList_[i].get());
        
        }
                
#endif

    }


    // Keyframe selection
    //latestTrackedFrame = trackingNewFrame;
#if 0 // USE BACK FOR LAST_TRACKED_FRAME
    
    latestTrackedFrame = trackingNewFrameList_.back();
    
#else // USE FRONT FOR LAST_TRACKED_FRAME
    
#if 1 // SET INTERNAL SPLINE
    trackingNewFrameList_.front().get()->setSpline(tracker->internalSpline);
    trackingNewFrameList_.front().get()->setUndistorter(this->undistorter);
#endif
    latestTrackedFrame = trackingNewFrameList_.front();
    
#endif
#if 0 //1 for FORCE CREATE NEW KEY FRAME
    if (!my_createNewKeyframe && currentKeyFrame->numMappedOnThisTotal
            > MIN_NUM_MAPPED)
    {
        createNewKeyFrame = true;
    }
#else
#if 0 // FORCE DO NOT CREATE KEY FRAME
    my_createNewKeyframe = true; // NO CREATION OF KEYFRAME WHATEVER YOU DO
#endif
    if (!my_createNewKeyframe && currentKeyFrame->numMappedOnThisTotal
            > MIN_NUM_MAPPED)
    {

        Sophus::Vector3d dist = newRefToFrame_poseUpdate.translation() *
                currentKeyFrame->meanIdepth;
        float minVal = fmin(0.2f + keyFrameGraph->keyframesAll.size() * 0.8f /
                            INITIALIZATION_PHASE_COUNT,1.0f);

        if(keyFrameGraph->keyframesAll.size() < INITIALIZATION_PHASE_COUNT)
            minVal *= 0.7;

        lastTrackingClosenessScore =
                trackableKeyFrameSearch->getRefFrameScore(dist.dot(dist),
                                     tracker->pointUsage*
                                     lsd_slam::pixelDensityForTracking);

        if (lastTrackingClosenessScore > minVal)
        {

#if 0//DEBUG_MODE
            printf("lastTrackingClosenessScore %f > minVal %f\n",
                   lastTrackingClosenessScore, minVal);
            printf("tracker->pointUsage = %f\n", tracker->pointUsage);
#endif
            createNewKeyFrame = true;

            if (enablePrintDebugInfo && printKeyframeSelectionInfo)
            {

                printf("SELECT %d on %d! dist %.3f + usage %.3f = %.3f > 1\n",
                       trackingNewFrame->id(),
                       trackingNewFrame->getTrackingParent()->id(),
                       dist.dot(dist),
                       tracker->pointUsage,
                       trackableKeyFrameSearch->getRefFrameScore(
                           dist.dot(dist),
                           tracker->pointUsage*
                           lsd_slam::pixelDensityForTracking));

            }

        }
        else
        {
#if 0//DEBUG
            printf("lastTrackingClosenessScore %f <= minVal %f\n",
                   lastTrackingClosenessScore, minVal);
            printf("tracker->pointUsage = %f\n", tracker->pointUsage);
#endif

            if (enablePrintDebugInfo && printKeyframeSelectionInfo)
            {

                printf("SKIPPED %d on %d! dist %.3f + usage %.3f = %.3f > 1\n",
                       trackingNewFrame->id(),
                       trackingNewFrame->getTrackingParent()->id(),
                       dist.dot(dist),
                       tracker->pointUsage,
                       trackableKeyFrameSearch->getRefFrameScore(
                           dist.dot(dist),
                           tracker->pointUsage*
                           lsd_slam::pixelDensityForTracking));

            }


        }

    }
#endif

    unmappedTrackedFramesMutex.lock();
    if(unmappedTrackedFrames.size() < 50 ||
            (unmappedTrackedFrames.size() < 100 &&
             trackingNewFrame->getTrackingParent()->numMappedOnThisTotal < 10))
#if 1 // PUSH THE LAST FRAME INTO UNMAPPED_TRACKED_FRAMES
    {
        trackingNewFrameList_.back()->setSpline(tracker->internalSpline);
        trackingNewFrameList_.back()->setUndistorter(this->undistorter);
        unmappedTrackedFrames.push_back(trackingNewFrameList_.back());
    }
#else // PUSH THE FIRST FRAME INTO UNMAPPED_TRACKED_FRAMES
    
    {
        
        for (int i=0;i<numToPop;i++)
        {
        
#if 1 // INTERNAL SPLINE FOR UNMAPPED TRACKED FRAMES
            // Make a trackingNewFrameList element have the internal spline
            trackingNewFrameList_[i].get()->setSpline(tracker->internalSpline);
            trackingNewFrameList_[i].get()->setUndistorter(this->undistorter);
#endif
            // First numToPop front into the unmapped tracked frames
            unmappedTrackedFrames.push_back(trackingNewFrameList_[i]);
            
            printf("UNMAPPED FRAME %d\n", trackingNewFrameList_[i].get()->id());
        
        }
    
    }
    
#endif
    unmappedTrackedFramesSignal.notify_one();
    unmappedTrackedFramesMutex.unlock();

    // implement blocking
    if (blockUntilMapped && trackingIsGood)
    {

        boost::unique_lock<boost::mutex> lock(newFrameMappedMutex);
        while(unmappedTrackedFrames.size() > 0)
        {

//            printf("TRACKING IS BLOCKING, waiting for %d frames to "
//                   "finish mapping.\n", (int)unmappedTrackedFrames.size());
            newFrameMappedSignal.wait(lock);

        }
        lock.unlock();

    }

}

#if 0
void SlamSystemForRollingShutter::mappingThreadLoop()
{
    printf("Started mapping thread!\n");
    while(keepRunning)
    {
        if (!doMappingIteration())
        {
            boost::unique_lock<boost::mutex> lock(unmappedTrackedFramesMutex);
            unmappedTrackedFramesSignal.timed_wait(lock,boost::posix_time::milliseconds(200));	// slight chance of deadlock otherwise
            lock.unlock();
        }
        
        newFrameMappedMutex.lock();
        newFrameMappedSignal.notify_all();
        newFrameMappedMutex.unlock();
    }
    printf("Exited mapping thread \n");
}
#endif

void SlamSystemForRollingShutter::finishCurrentKeyframe()
{
    if(enablePrintDebugInfo && printThreadingInfo)
        printf("FINALIZING KF %d\n", currentKeyFrame->id());
    
    map->finalizeKeyFrame(); // Use depth map for rolling shutter
    
    if(SLAMEnabled)
    {
        mappingTrackingReference->importFrame(currentKeyFrame.get());
        currentKeyFrame->setPermaRef(mappingTrackingReference);
        mappingTrackingReference->invalidate();
        
        if(currentKeyFrame->idxInKeyframes < 0)
        {
            keyFrameGraph->keyframesAllMutex.lock();
            currentKeyFrame->idxInKeyframes = keyFrameGraph->keyframesAll.size();
            keyFrameGraph->keyframesAll.push_back(currentKeyFrame.get());
            keyFrameGraph->totalPoints += currentKeyFrame->numPoints;
            keyFrameGraph->totalVertices ++;
            keyFrameGraph->keyframesAllMutex.unlock();
            
            newKeyFrameMutex.lock();
            newKeyFrames.push_back(currentKeyFrame.get());
            newKeyFrameCreatedSignal.notify_all();
            newKeyFrameMutex.unlock();
        }
    }
    
    if(outputWrapper!= 0)
        outputWrapper->publishKeyframe(currentKeyFrame.get());
}

void SlamSystemForRollingShutter::createNewCurrentKeyframe(std::shared_ptr<Frame> newKeyframeCandidate)
{
    if(enablePrintDebugInfo && printThreadingInfo)
        printf("CREATE NEW KF %d from %d\n", newKeyframeCandidate->id(), currentKeyFrame->id());
    
    
    if(SLAMEnabled)
    {
        // add NEW keyframe to id-lookup
        keyFrameGraph->idToKeyFrameMutex.lock();
        keyFrameGraph->idToKeyFrame.insert(std::make_pair(newKeyframeCandidate->id(), newKeyframeCandidate));
        keyFrameGraph->idToKeyFrameMutex.unlock();
    }
    
    // propagate & make new.
    map->createKeyFrame(newKeyframeCandidate.get());
    
    if(printPropagationStatistics)
    {
        
        Eigen::Matrix<float, 20, 1> data;
        data.setZero();
        data[0] = runningStats.num_prop_attempts / ((float)width*height);
        data[1] = (runningStats.num_prop_created + runningStats.num_prop_merged) / (float)runningStats.num_prop_attempts;
        data[2] = runningStats.num_prop_removed_colorDiff / (float)runningStats.num_prop_attempts;
        
        outputWrapper->publishDebugInfo(data);
    }
    
    currentKeyFrameMutex.lock();
    currentKeyFrame = newKeyframeCandidate;
    currentKeyFrameMutex.unlock();
}

void SlamSystemForRollingShutter::addTimingSamples()
{
    map->addTimingSample();
    struct timeval now;
    gettimeofday(&now, NULL);
    float sPassed = ((now.tv_sec-lastHzUpdate.tv_sec) + (now.tv_usec-lastHzUpdate.tv_usec)/1000000.0f);
    if(sPassed > 1.0f)
    {
        nAvgTrackFrame = 0.8*nAvgTrackFrame + 0.2*(nTrackFrame / sPassed); nTrackFrame = 0;
        nAvgOptimizationIteration = 0.8*nAvgOptimizationIteration + 0.2*(nOptimizationIteration / sPassed); nOptimizationIteration = 0;
        nAvgFindReferences = 0.8*nAvgFindReferences + 0.2*(nFindReferences / sPassed); nFindReferences = 0;
        
        if(trackableKeyFrameSearch != 0)
        {
            trackableKeyFrameSearch->nAvgTrackPermaRef = 0.8*trackableKeyFrameSearch->nAvgTrackPermaRef + 0.2*(trackableKeyFrameSearch->nTrackPermaRef / sPassed); trackableKeyFrameSearch->nTrackPermaRef = 0;
        }
        nAvgFindConstraintsItaration = 0.8*nAvgFindConstraintsItaration + 0.2*(nFindConstraintsItaration / sPassed); nFindConstraintsItaration = 0;
        nAvgOptimizationIteration = 0.8*nAvgOptimizationIteration + 0.2*(nOptimizationIteration / sPassed); nOptimizationIteration = 0;
        
        lastHzUpdate = now;
        
        
        if(enablePrintDebugInfo && printOverallTiming)
        {
            printf("MapIt: %3.1fms (%.1fHz); Track: %3.1fms (%.1fHz); Create: %3.1fms (%.1fHz); FindRef: %3.1fms (%.1fHz); PermaTrk: %3.1fms (%.1fHz); Opt: %3.1fms (%.1fHz); FindConst: %3.1fms (%.1fHz);\n",
                   map->msUpdate, map->nAvgUpdate,
                   msTrackFrame, nAvgTrackFrame,
                   map->msCreate+map->msFinalize, map->nAvgCreate,
                   msFindReferences, nAvgFindReferences,
                   trackableKeyFrameSearch != 0 ? trackableKeyFrameSearch->msTrackPermaRef : 0, trackableKeyFrameSearch != 0 ? trackableKeyFrameSearch->nAvgTrackPermaRef : 0,
                   msOptimizationIteration, nAvgOptimizationIteration,
                   msFindConstraintsItaration, nAvgFindConstraintsItaration);
        }
    }
    
}

void SlamSystemForRollingShutter::discardCurrentKeyframe()
{
    if(enablePrintDebugInfo && printThreadingInfo)
        printf("DISCARDING KF %d\n", currentKeyFrame->id());
    
    if(currentKeyFrame->idxInKeyframes >= 0)
    {
        printf("WARNING: trying to discard a KF that has already been added to the graph... finalizing instead.\n");
        finishCurrentKeyframe();
        return;
    }
    
    
    map->invalidate();
    
    keyFrameGraph->allFramePosesMutex.lock();
    for(FramePoseStruct* p : keyFrameGraph->allFramePoses)
    {
        if(p->trackingParent != 0 && p->trackingParent->frameID == currentKeyFrame->id())
            p->trackingParent = 0;
    }
    keyFrameGraph->allFramePosesMutex.unlock();
    
    
    keyFrameGraph->idToKeyFrameMutex.lock();
    keyFrameGraph->idToKeyFrame.erase(currentKeyFrame->id());
    keyFrameGraph->idToKeyFrameMutex.unlock();
    
}

void SlamSystemForRollingShutter::loadNewCurrentKeyframe(Frame* keyframeToLoad)
{
    if(enablePrintDebugInfo && printThreadingInfo)
        printf("RE-ACTIVATE KF %d\n", keyframeToLoad->id());
    
    map->setFromExistingKF(keyframeToLoad);
    
    if(enablePrintDebugInfo && printRegularizeStatistics)
        printf("re-activate frame %d!\n", keyframeToLoad->id());
    
    currentKeyFrameMutex.lock();
    currentKeyFrame = keyFrameGraph->idToKeyFrame.find(keyframeToLoad->id())->second;
    currentKeyFrame->depthHasBeenUpdatedFlag = false;
    currentKeyFrameMutex.unlock();
}

void SlamSystemForRollingShutter::changeKeyframe(bool noCreate, bool force, float maxScore)
{
    Frame* newReferenceKF=0;
    std::shared_ptr<Frame> newKeyframeCandidate = latestTrackedFrame;
    if(doKFReActivation && SLAMEnabled)
    {
        struct timeval tv_start, tv_end;
        gettimeofday(&tv_start, NULL);
        newReferenceKF = trackableKeyFrameSearch->findRePositionCandidate(newKeyframeCandidate.get(), maxScore);
        gettimeofday(&tv_end, NULL);
        msFindReferences = 0.9*msFindReferences + 0.1*((tv_end.tv_sec-tv_start.tv_sec)*1000.0f + (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
        nFindReferences++;
    }
    
    if(newReferenceKF != 0)
        loadNewCurrentKeyframe(newReferenceKF);
    else
    {
        if(force)
        {
            if(noCreate)
            {
                trackingIsGood = false;
                nextRelocIdx = -1;
                printf("mapping is disabled & moved outside of known map. Starting Relocalizer!\n");
            }
            else
                createNewCurrentKeyframe(newKeyframeCandidate);
        }
    }
    
    
    createNewKeyFrame = false;


#if 0    // Clear deques related to tracking
    trackingNewFrameList_.clear();
    referenceList_.clear();
    startRowImageList_.clear();
    startRowSegList_.clear();
#endif
    
}