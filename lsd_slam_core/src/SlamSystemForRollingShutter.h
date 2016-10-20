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

#pragma once
#include "SlamSystem.h"

#include <vector>
#include <boost/thread.hpp>
#include <boost/thread/shared_mutex.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/thread/locks.hpp>
#include "util/settings.h"
#include "IOWrapper/Timestamp.h"
#include "opencv2/core/core.hpp"

#include "util/SophusUtil.h"

#include "Tracking/Relocalizer.h"
#include "util/Undistorter.h"
#include "DepthEstimation/DepthMapForRollingShutter.h"
#include "DepthEstimation/DepthMapForRadialRolling.h"

namespace lsd_slam
{

class TrackingReference;
class KeyFrameGraph;
class SE3TrackerRollingShutter;
class Sim3Tracker;
class DepthMap;
class Frame;
class DataSet;
class LiveSLAMWrapper;
class Output3DWrapper;
class TrackableKeyFrameSearch;
class FramePoseStruct;
struct KFConstraintStruct;

typedef Eigen::Matrix<float, 7, 7> Matrix7x7;

struct TrackingReferenceData
{
    // Structure keeping a list of tracking references, which contain 3D points,
    // image gradient, image and variance, with the starting row in image
    // the starting row in B-Spline segument for corresponding the tracking
    // reference

    std::shared_ptr<TrackingReference> reference;
    int startRowImage;
    int startRowSeg;

};

struct TrackingFrameData
{

    std::shared_ptr<Frame> frame;
    int startRowImage;
    int startRowSeg;

};

class SlamSystemForRollingShutter : public SlamSystem
{

    friend class IntegrationTest;

private:
    
    int rowPerSegment; // Number of rows per B-Spline segment

    // Image size of B-Spline segment image
    int width_bsplineSegImg;
    int height_bsplineSegImg;

    // tracking new frame list
    std::deque<std::shared_ptr<Frame>> trackingNewFrameList_;
    std::deque<int> startRowImageList_;
    int maxCapList_;
    std::deque<int> startRowSegList_;

    // reference list
    std::deque<std::shared_ptr<TrackingReference> > referenceList_;

public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    SlamSystemForRollingShutter(int w, int h, Eigen::Matrix3f K,
                                Undistorter *undistorter,
                                int width_bsplineSegImg,
                                int height_bsplineSegImg,
                                int rowPerSegment,
                                int methodRollingShutter,
                                bool enableSLAM = true,
                                bool useGTmotion = false);

	~SlamSystemForRollingShutter();

    // Set undistortion class for rolling shutter
    void setUndistorter(Undistorter *undistorter)
    {

        this->undistorter = undistorter;

    }

    // Get undistortion
    Undistorter *getUndistorter()
    {

        return undistorter;

    }
    
    // Get tracker
    SE3TrackerRollingShutter* getTracker()
    {
        return tracker;
    }

    // Mapping for rolling shutter camera
    bool doMappingIteration();
    bool updateKeyframe();

    // Random depth initialization
    using SlamSystem::randomInit;
    void randomInitDepthRange(cv::Mat image,
                    cv::Mat imageMask,
                    double timeStamp, int id,
                    float minInitDepth = 0.6666f,
                    float maxInitDepth = 2.0f);
    
    // Depth transfer from first keyframe of GS SLAM (e.g. LSD-SLAM) system
    void depthTransferFrom(Frame* firstKeyFrame);
    
    // GT depth init
    void gtDepthInit(uchar* image, uchar* imageMask,
                     float* depth, double timeStamp, int id);
    void gtDepthInitWithName(uchar* image, float* depth,
                             std::string depthFilename,
                             std::string motionFilename,
                             double timeStamp, int id);
    
    // GT depth with an initial noise
    void gtDepthInitWithNoise(uchar* image, float* depth,
                              double timeStamp, int id);

    // Track frame
    void trackFrame_GS(uchar* image,
                    cv::Mat imageMask,
                    unsigned int frameID,
                    bool blockUntilMapped,
                    double timestamp,
                    int height_rollingShutter_inputImage,
                    std::string imageFilename = "",
                    std::string depthFilename = "");
    void trackFrameRS(cv::Mat image, unsigned int frameID,
                      bool blockUntilMapped, double timestamp);

    // Track frame with B-spline neigbhour segment version 3
    void trackFrameRS_BN3(cv::Mat image,
                          cv::Mat imageMask,
                          unsigned int frameID,
                          bool blockUntilMapped,
                          double timestamp,
                          int startRowImage, int startRowSeg,
                          double rowPerSegment,
                          double framePerSegment,
                          int height_rollingShutter_inputImage,
                          std::vector<std::string> imageFiles,
                          std::vector<std::string> depthFiles,
                          std::vector<std::string> motionFiles);

protected:

    // Image undistortion pointer (to find a row in the original image)
    Undistorter *undistorter;

    // ============= EXCLUSIVELY TRACKING THREAD (+ init) ===============
	SE3TrackerRollingShutter* tracker;
    
#if 1// WORKING for Depth map ROLLING SHUTTER
    
    // ============= EXCLUSIVELY MAPPING THREAD (+ init) =============
    DepthMapForRollingShutter* map;
    
#endif

    // ============= EXCLUSIVELY FIND-CONSTRAINT THREAD (+ init) ========
    SE3TrackerRollingShutter* constraintSE3Tracker;


    // ============= ETC
    void debugDisplayDepthMap();
#if 0
	void mappingThreadLoop();    
#endif
    void finishCurrentKeyframe();
	void changeKeyframe(bool noCreate, bool force, float maxScore);    
	void createNewCurrentKeyframe(std::shared_ptr<Frame> newKeyframeCandidate);    
	void addTimingSamples();
    void discardCurrentKeyframe();
    void loadNewCurrentKeyframe(Frame* keyframeToLoad);
    
    
};

}
