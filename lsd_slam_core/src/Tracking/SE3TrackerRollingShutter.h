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
// HISTORY
//------------------------------------------------------------------------------
// 2014-11-11 B-Spline cost function for rolling shutter
// 2014-11-09 Ceres version for SE3 tracking
// 2014-10-28 Writing started
//------------------------------------------------------------------------------
#pragma once
#include <opencv2/core/core.hpp>
#include "util/settings.h"
#include "util/EigenCoreInclude.h"
#include "util/SophusUtil.h"
#include "Tracking/least_squares.h"

#include "Tracking/SE3Tracker.h"
#include "util/Undistorter.h"
#include "DataStructures/Frame.h"

#include "Tracking/Spline.h"

#include <memory>

#include "ceres/ceres.h"

#include <deque>

namespace lsd_slam
{

#define SPLINE_K 2 // Linear B-spline (degree 1)
//#define SPLINE_K 3 // Quadratic B-spline (degree 2)
//#define SPLINE_K 4 // Cubic B-spline (degree 3)
    
//------------------------------------------------------------------------------
#define INIT_NUM_CTRLPTS (3*SPLINE_K - 2)
#define START_INDEX_CTRLPTS (SPLINE_K - 1)
#define NUM_EXTRA_CTRLPTS (SPLINE_K - 1)

// Given the order of B-spline k, the initial number of control points
// being estimated for a first tracking frame is 3*k - 2. This number
// was computed from two assumptions that (1) the starting point
// on the B-spline curve coincides at the location of
// the first control point of B-spline; (2) all control points being estimated
// are determined by B-spline segments affected by the first k control points.
// With the first assumption, in B-spline curves, these extra and static
// control points are required to make B-spline pass through a specific point.
// The number of these extra control points are k - 1.
// For the second assumption, the number of segments affected by k control
// points is k, and each segment shares k - 1 neighbour control points. Thus,
// In total, the number of control points being estimated at first becomes
// 2*k - 1. Now, with these numbers, we may proceed adding a new control point
// for the next B-spline segment.
// A condition adding a new control point is when i + k - 1 >= n, where
// n is the number of all control points having been estimated,
// k is the order of B-spline and
// i is the index of a starting control point for a new set of k control points
//   to be estimated.
//
// Example)
//
// For B-spline with the order k = 2,
// the number of extra/static control points is 1,
// the initial number of control points is 4 and
// the starting index of control points is 1.
//
//    0 [1  2]
//         [2  3]
//
// Thus, [1] is the starting index from all control points being estimated,
// and [0] is the extra control point. The condition adding a new set of
// control points [3 4] happens when the index i = 3 since 3 + 2 - 1 >= 4.
//
// For B-spline with the order k = 3,
// the number of extra/static control points is 2,
// the initial number of control points is 7 and
// the starting index of control points is 2.
//
//    0  1 [2  3  4]
//            [3  4  5]
//               [4  5  6]
//
// Thus, [2] is the starting index from all control points being estimated,
// and [0] and [1] are the extra control points. The condition adding a new set
// of control points [5 6 7] happens when the index i = 5 since 5 + 3 - 1 >= 7.
//
// For B-spline with the order k = 4,
// the number of extra/static control points is 3,
// the initial number of control points is 10 and
// the starting index of control points is 3.
//
//    0  1  2 [3  4  5  6]
//               [4  5  6  7]
//                  [5  6  7  8]
//                     [6  7  8  9]
//
// Thus, [3] is the starting index from all control points being estimated,
// and [0], [1] and [2] are the extra control points. The condition adding
// a new set of control points [7 8 9 10] happens when the index i = 7
// since 7 + 4 - 1 >= 10.
//
//------------------------------------------------------------------------------

class TrackingReference;
class Frame;

class SE3TrackerRollingShutter : public SE3Tracker
{
public:
    
    int spline_k = SPLINE_K;
    
    // Debug info
    float debugInfo_last_residual_;

    // A vector of control poses
    std::vector<Sophus::SE3d> controlPoses; // (WorldToControl)

    // Ceres problems holding for each pyramid level
    // problems[0] : ceres::Problem for Level 0
    // problems[1] : ceres::Problem for Level 1
    std::vector<ceres::Problem*> problemList;


    // file list
    std::vector<std::string> files_;
    
    
    // Internal spline structure
    std::shared_ptr<Spline<double> > internalSpline;
    
    // method used for rolling-shutter
    int methodRollingShutter; // 0: RSD-SLAM, 1: RRD-SLAM

private:

    // Size of B-Spline segment image
    int width_bsplineSegImg;
    int height_bsplineSegImg;

    // Info optimiztaion result
    double final_cost_;   // final cost from Ceres::Solver::Summary
    double last_residual_; // Last residual (sumResSquare/num)
    int num_residual_blocks_; // Number of residual blocks

public:

    // debug image (original input image)
    //cv::Mat debugImageOriginalInput;

    // Declare them as static
    cv::Mat debugImageResidualsRollingShutter;
    cv::Mat debugImageWeightsRollingShutter;
    cv::Mat debugImageSecondFrameRollingShutter;
    cv::Mat debugImageOldImageSourceRollingShutter;
    cv::Mat debugImageOldImageWarpedRollingShutter;

    // Gradient images
    cv::Mat debugGradientXImage0_RollingShutter;
    cv::Mat debugGradientXImage1_RollingShutter;

    // Tracking reference image
    cv::Mat debugRefImage0;
    cv::Mat debugRefImage1;
    cv::Mat debugRefImage2;

    bool isWindowMoved;
    bool isWindowMoved2;

    // Frame number to be used in displaying or saving images
    int frameNumberToDisplay;
    int getFrameNumberToDisplay() { return frameNumberToDisplay; }

public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    SE3TrackerRollingShutter(int w, int h, Eigen::Matrix3f K,
                             Undistorter *undistorter,
                             int width_segImg, int height_bsplineSegImg,
                             int methodRollingShutter);
    ~SE3TrackerRollingShutter();

    //--------------------------------------------------------------------------
    // Track the whole size frame as SE3Tracker using Ceres solver
    SE3 trackFrame_Ceres(
            TrackingReference* reference,
            Frame* frame,
            const SE3& frameToReference_initialEstimate,
            const float startRow,
            const float endRow);

    //--------------------------------------------------------------------------
    // Track frame between reference image and current frame by a row
    //--------------------------------------------------------------------------
    //
    //    startRow : the starting row value for warping between 0 and 1
    //               0 indicates the first row, 1 indiates the last
    //               row
    //      endRow : the end row value for warping between 0 and 1
    //               0 means the first, and 1 means the last row
    //               E.g. startRow = 0.4 and endRow = 0.5
    SE3 trackFrameByRows(TrackingReference* reference, Frame* frame,
                   const SE3& frameToReference_initialEstimate,
                   const float startRow,
                   const float endRow);

    //--------------------------------------------------------------------------
    // Track frame between reference image and current frame
    // Experimental method testing a set of predefined rows
    SE3 trackFrame(TrackingReference* reference, Frame* frame,
                   const SE3& frameToReference_initialEstimate);

    //--------------------------------------------------------------------------
    // Debug related methods
    void calcResidualAndBuffers_debugStart();
    void calcResidualAndBuffers_debugFinish(int w);
    void calcResidualAndBuffers_debugFinish(int w, Frame *frame);

    //--------------------------------------------------------------------------
    // B-spline with neighour segments version 3 with full window
    SE3 trackFrame_Ceres_BN3(TrackingReference* reference,
            std::deque<std::shared_ptr<TrackingReference> > *referenceList,
            std::deque<std::shared_ptr<Frame> > *frameList, // for tracking new frame
            const SE3& frameToReference_initialEstimate,
            std::deque<int> startRowImageList,
            std::deque<int> startRowSegList,
            double rowPerSegment,
            double framePerSegment,
            double scale,
            std::vector<std::string> files);

protected:
    
    //--------------------------------------------------------------------------
    // Compute a residual error of warping between keyframe and current frame
    // by a set of rows for rolling shutter
    float calcResidualAndBuffersByRows(const Eigen::Vector3f* refPoint,
                                 const Eigen::Vector2f* refColVar,
                                 int* idxBuf,
                                 int refNum,
                                 Frame* frame,
                                 const Sophus::SE3f& referenceToFrame,
                                 int level,
                                 const float startRow,
                                 const float endRow,
                                 bool plotResidual = false);

#if defined(ENABLE_SSE)

    float calcResidualAndBuffersByRowsSSE(
            const Eigen::Vector3f* refPoint,
            const Eigen::Vector2f* refColVar,
            int* idxBuf,
            int refNum,
            Frame* frame,
            const Sophus::SE3f& referenceToFrame,
            int level,
            const float startRow,
            const float endRow,
            bool plotResidual = false);

#endif

    //--------------------------------------------------------------------------
    // Compute a residual error of warping between keyframe and current frame
    float calcResidualAndBuffers(const Eigen::Vector3f* refPoint,
                                 const Eigen::Vector2f* refColVar,
                                 int* idxBuf,
                                 int refNum,
                                 Frame* frame,
                                 const Sophus::SE3f& referenceToFrame,
                                 int level,
                                 bool plotResidual = false);

private:

    Undistorter *undistorter;

    // Optimization by normal equation (LSD-SLAM way for optimization)
    bool optimzeByNormalEquationLeastSquares(int *numCalcResidualCalls,
            int *numCalcWarpUpdateCalls,
            float *last_residual,
            TrackingReference* reference,
            Frame* frame,
            Sophus::SE3f *referenceToFrame, float startRow, float endRow);

    //--------------------------------------------------------------------------
    // Optimization by Ceres-solver
    //--------------------------------------------------------------------------
    // Optimize SE(3) tracker using Ceres-solver. Yet another implementation
    // of SE(3) tracker in LSD-SLAM using automatic differentiation
    bool optimizeUsingCeres(TrackingReference* reference,
                            Frame* frame,
                            Sophus::SE3f *referenceToFrame,
                            float startRow,
                            float endRow);

    // Optimization by Ceres-solver for the given level of image
    bool
    optimizeUsingCeresPyramidLevelAt(
            int level,
            double *pose,
            TrackingReference* reference,
            Frame *frame);

    //--------------------------------------------------------------------------
    // Using Ceres for B-Spline

    // Track a new frame and return control poses as B-Spline basis
    bool optimizeUsingCeresBSpline(TrackingReference* reference,
                                   Frame* frame,
                                   Sophus::SE3f *referenceToFrame,
                                   std::vector<Sophus::SE3d> controlPoses);

    // Optimize at pyramid level
    bool optimizeUsingCeresBSplinePyramidLevelAt(int level,
            double *pose,
            TrackingReference* reference,
            Frame *frame, const int numControlPoses);

    //--------------------------------------------------------------------------
    // BN3
    //--------------------------------------------------------------------------
    // Track a new frame and return control poses as B-Spline basis
    bool optimizeUsingCeres_BN3(TrackingReference *reference,
            std::deque<std::shared_ptr<TrackingReference> > *referenceList,
            std::deque<std::shared_ptr<Frame> > *frameList,
            Sophus::SE3d *referenceToFrame,
            std::vector<Sophus::SE3d> *controlPoses,
            std::deque<int> startRowImageList,
            std::deque<int> startRowSegList,
            double rowPerSegment,
            double scale);

    // Optimize at pyramid level
    bool optimizeUsingCeres_BN3_PyramidLevelAt(
            TrackingReference *reference,
            int level,
            double *pose,
            std::deque<std::shared_ptr<TrackingReference> > *referenceList,
            std::deque<std::shared_ptr<Frame> > *frameList,
            const int numControlPoses,
            std::deque<int> startRowImageList,
            std::deque<int> startRowSegList,
            double rowPerSegment,
            double lowerBound_uValue,
            double upperBound_uValue,
            const std::vector<Sophus::SE3d> controlPoses,
            int idxSplineCtrlPts,
            double scale);


}; // end of class

} // end of namespace
