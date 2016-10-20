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
// SE(3) based tracking of new frame to the reference key frame for a rolling
// shutter camera with B-Spline control pose tuning
//
// Copyright@2015 University of Adelaide
// Author: Jae-Hak Kim (jaehak.kim@adelaide.edu.au)
//==============================================================================
#ifndef WARPRESIDUAL_BN3_H
#define WARPRESIDUAL_BN3_H

#include "sophus/se3.hpp"
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "ceres/rotation.h"
#include "Tracking/JetOps.h"

#include "Tracking/SE3TrackerRollingShutter.h"
#include "Tracking/TrackingReference.h"
#include "DataStructures/Frame.h"
#include "util/Undistorter.h"
#include "util/globalFuncs.h"

#include "Tracking/Spline.h"

#include <opencv2/opencv.hpp>

#define OUTLIER_RESIDUAL_INTENSITY_BN3 NULL // (10*5)/1.0 Outlier residual for intensity difference
#define OUTLIER_RESIDUAL_ROWDIFF_BN3 NULL // (# rows)/1.0 Outlier residual for row difference
#define ALPHA_BN3 500 // 500 // Weight for the second residual

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;

namespace lsd_slam
{

class TrackingReference;
class Frame;

struct WarpResidual_BN3
{

public:

    // Resulf of warping statistics
    static float usageCount_;
    static int lastGoodCount_;
    static int lastBadCount_;

    static int goodWarpPixelCountInPyramidImage_;
    static int badWarpPixelCountInPyramidImage_;

private:

    // Pyramid level
    const int level_;

    // X,Y,Z point in the reference (source) image
    const Eigen::Vector3f refPoint_;

    // Intensity/Colour and Variance for inverse depth of the point
    // in the reference (source) image
    const Eigen::Vector2f refColourVariance_;

    // tracking reference containing the reference (source)image,
    // gradients and inverse depth
    //TrackingReference *reference_;

    // A row pose of an image point in the reference image
    Sophus::SE3d rowPose_RefPoint_;

    // Current image (target) as a vector of Frame pointers
    //std::deque<std::shared_ptr<Frame> > *frameList_;

    // Current frame
    Frame *frame_;
    //std::shared_ptr<Frame> frame_;

    // Image undistortion/distortion object
    Undistorter *undistorter;

    // Index for rows
    const int startRowImage_;
    const int startRowSeg_;
    const double rowPerSegment_;

    // All control points
    const std::vector<Sophus::SE3d> controlPoses_; // WorldToControl
    int idxSC_frontKF_; // Index of control point for front keyframe
    int idxSC_RefPoint_; // Index of control point for each keyframe
    int idxSplineCtrlPts_; // Index of control point to be estimated
    int numExtraCtrlPts_; // Number of extra contrl points for beginning

    // Scale for point
    double scale_;
    int idx_;

    const DenseDepthTrackerSettings settings_;

    // file list
    const std::vector<std::string> files_;
    
    // rolling-shutter method used
    const int methodRollingShutter_;

public:

    //--------------------------------------------------------------------------
    // WarpResidual_BN3 CostFunctor for Ceres-solver
    //--------------------------------------------------------------------------
    // NOTE: It warps a point from "reference" image to each image in
    //       the "frameList"
    //
    //         offsetRow : [int] Row offset of source segment image
    //     idxSrcSegment : [int] frame index for source of B-Spline segment
    //             level : [int] pyramid level such that
    //                     SE3TRACKING_MIN_LEVEL <= level &&
    //                     level <= SE3TRACKER_MAX_LEVEL - 1
    //          refPoint : [3x1 float] semi-dense point of (x, y, z) such that
    //                     x coords, y coords and inverse depth z in the
    //                     reference image
    // refColourVariance : [2x1 float] semi-dense point of (I, V) such that
    //                     image intensity and variance
    //         reference : [TrackingReference] Keyframe image (source image)
    //         frameList : [Frame *] Current new image (target image) in list
    //  referenceToFrame : [SE3f] Key frame pose wrt new image pose
    //       undistorter : [Undistorter] Undistortion object
    //                     used for finding warpd image coordinates in original
    //                     input image
    //      controlPoses : [vector<SE3d> *] pointer to all control poses
    //                         respective to the world
    //     idxSC_frontKF : [int] index of control point for the front keyframe
    //    idxSC_RefPoint : [int] index of control point for each keyframe
    //  idxSplineCtrlPts : [int] index of starting B-Spline control points
    //             scale : [double] Scale for reference point
    //
    WarpResidual_BN3 (
            const int level,
            const Eigen::Vector3f refPoint,
            const Eigen::Vector2f refColourVariance,
            TrackingReference *reference,
            Sophus::SE3d rowPose_RefPoint,
            Frame *frame,
            Undistorter *undistorter,
            const int startRowImage,
            const int startRowSeg,
            const double rowPerSegment,
            const std::vector<Sophus::SE3d> controlPoses,
            int idxSC_frontKF,
            int idxSC_RefPoint,
            int idxSplineCtrlPts,
            const double scale,
            const int idx,
            const DenseDepthTrackerSettings settings,
            const std::vector<std::string> files,
            const int methodRollingShutter)
        : level_(level),
          refPoint_(refPoint),
          refColourVariance_(refColourVariance),
          //reference_(reference),
          rowPose_RefPoint_(rowPose_RefPoint),
          frame_(frame),
          startRowImage_(startRowImage),
          startRowSeg_(startRowSeg),
          rowPerSegment_(rowPerSegment),
          controlPoses_(controlPoses),
          idxSC_frontKF_(idxSC_frontKF),
          idxSC_RefPoint_(idxSC_RefPoint),
          idxSplineCtrlPts_(idxSplineCtrlPts),
          scale_(scale),
          idx_(idx),
          settings_(settings),
          files_(files),
          methodRollingShutter_(methodRollingShutter)
    {

        this->undistorter = undistorter;
        this->numExtraCtrlPts_ = SPLINE_K - 1;

#if 0// CHECK REF
        //----------------------------------------------------------------------
        // Check ref and frame images
        // Am I really using correct images?
        //----------------------------------------------------------------------

        // Load an image of the reference frame ID
        cv::Mat refImg1 = cv::imread(files_[reference_->frameID],
                                    CV_LOAD_IMAGE_GRAYSCALE);
        // Compare with the reference's keyframe image
        float sum_val = 0.0;
        for (int i=0;i<refImg1.rows;i++)
        {
            for (int j=0;j<refImg1.cols;j++)
            {

                float val = (float)refImg1.at<uchar>(i,j) -
                (reference_->keyframe->image(0)[i*refImg1.cols + j]);

                sum_val = sum_val + val*val;
            }
        }
        float mean_val = sum_val / (float) (refImg1.rows*refImg1.cols);
        if (sum_val != 0.0)
        {
            fprintf(stderr, "CHECKREF ref and keyframe difference: "
               "Frame %d, Sum_val %f; Mean_val %f\n",
               reference_->frameID,
               sum_val, mean_val);
            exit(1);
        }
#endif


    }
    
    //--------------------------------------------------------------------------
    // Operator computing the residual error for warping a pixel from the
    // reference image to the new image
    //--------------------------------------------------------------------------
    //      pose : [6*k double] 6*k parameters for k control poses
    //                          Ctrl pts to Front Keyframe
    //    u_time : [1x1 double] u value of the warped point q
    //  residual : [2x1 double] residual error
    //
    template <typename T> bool
    operator() (const T* const pose, const T* const uTimes, T* residual) const
    {

        // Define SE3Group based on template typename T
        typedef Sophus::SE3Group<T> SE3T;

        //----------------------------------------------------------------------
        // Declare Spline in template T


#if 1// USE ROW POSE for keyframe

        // Get refToWorld for template T
        SE3 refToWorld_se3 = rowPose_RefPoint_.inverse(); // row ref to world

#else // USE global pose for keyframe

        // Get refToWorld for template Tf
        SE3 refToWorld_se3 = se3FromSim3(
                    reference_->keyframe->pose->getCamToWorld());

#endif
        Eigen::Matrix<T,6,1> refToWorld_se3T_log;
#if 0//USE LOG AND EXP
        for (int i=0; i<6; i++)
        {

            refToWorld_se3T_log[i] = T(refToWorld_se3.log()[i]);

        }
        SE3T refToWorldT = SE3T::exp(refToWorld_se3T_log);
#else // DO NOT USE LOG OR EXP

        Eigen::Matrix<double,4,4> tmpM = refToWorld_se3.matrix();
        Eigen::Matrix<T,4,4> tmpMT;
        tmpMT << T(tmpM(0,0)), T(tmpM(0,1)), T(tmpM(0,2)), T(tmpM(0,3)),
                 T(tmpM(1,0)), T(tmpM(1,1)), T(tmpM(1,2)), T(tmpM(1,3)),
                 T(tmpM(2,0)), T(tmpM(2,1)), T(tmpM(2,2)), T(tmpM(2,3)),
                 T(tmpM(3,0)), T(tmpM(3,1)), T(tmpM(3,2)), T(tmpM(3,3));
        SE3T refToWorldT = SE3T(tmpMT);

#endif

        // Copy all control poses
        std::vector<SE3T> controlPoses_T;
        for (int i=0; i<(int)controlPoses_.size(); i++)
        {
#if 0//USE LOG AND EXP
            // World To Control in template T
            Eigen::Matrix<T,6,1> controlPoses_se3T_log;
            for (int j=0; j<6; j++)
            {

                controlPoses_se3T_log[j] = T((controlPoses_)[i].log()[j]);

            }
            SE3T cpT = SE3T::exp(controlPoses_se3T_log);
#else //DO NOT USE LOG OR EXP

            Eigen::Matrix<double,4,4> M = (controlPoses_)[i].matrix();
            Eigen::Matrix<T,4,4> MT;
            MT << T(M(0,0)), T(M(0,1)), T(M(0,2)), T(M(0,3)),
                  T(M(1,0)), T(M(1,1)), T(M(1,2)), T(M(1,3)),
                  T(M(2,0)), T(M(2,1)), T(M(2,2)), T(M(2,3)),
                  T(M(3,0)), T(M(3,1)), T(M(3,2)), T(M(3,3));
            SE3T cpT = SE3T(MT); // World to Control

#endif

#if 0 // TEST new ref to control from the assumption about getScaledCamToWorld is actually world to cam
            SE3T cpT_refToControl = cpT * refToWorld.inverse();
            controlPoses_T.push_back(cpT_refToControl);
#else // USE ref to control
            // Convert to refToControl
            // RefToControl = WorldToControl * RefToWorld
            SE3T cpT_refToControl = cpT * refToWorldT;


            // Store (RefToControl)
            controlPoses_T.push_back(cpT_refToControl);
#endif

        }

        // Copy 4 control poses from array to SE3T vector
#if 0//BACKUP
//        std::vector<SE3T> controlPoses_T;
//        for (int i=0; i<4; i++)
//        {

//            Eigen::Matrix<T,6,1> pose_6vec;
//            for (int j=0; j<6; j++)
//                pose_6vec[j] = pose[6*i + j]; // ReferenceToControl

//            controlPoses_T.push_back(SE3T::exp(pose_6vec));

//        }
#endif
        // Determine the number of active control points to be estimated
        int idxForFrameI = (int)floor(startRowImage_ / rowPerSegment_)
                           + NUM_EXTRA_CTRLPTS;
        int numActiveCtrlPts = SPLINE_K - (idxSplineCtrlPts_ - idxForFrameI);
        
#if 0 // USE WORLD To CONTROL pose
#else // USE CONTROl TO FRONT KEYFRAME
        
        // World to Front keyframe
        Eigen::Matrix<T,6,1> worldToFrontKeyframe_se3T_log;
        for (int j=0; j<6; j++)
        {
            
            // ControlPoses_[idxSC_frontKF_] is a control point related to
            // the front keyframe
            worldToFrontKeyframe_se3T_log[j]
                = T((controlPoses_)[idxSC_frontKF_].log()[j]);
        }
        SE3T worldToFrontKeyframe_T = SE3T::exp(worldToFrontKeyframe_se3T_log);
#endif

#if 0// USE active number of control points
        for (int i=0; i<numActiveCtrlPts; i++)
#else // Use all control points
        for (int i=0; i<SPLINE_K; i++)
#endif
        {

#if 0 // USE WORLD_TO_CONTROL pose (ORIGINAL IMPLEMENTATION)

            Eigen::Matrix<T,6,1> pose_6vec;
            for (int j=0; j<6; j++)
            {
                
                pose_6vec[j] = pose[6*i + j]; // WorldToControl
                
            }
            
#else // USE CONTROL TO FRONT KEYFRAME
                
            Eigen::Matrix<T,6,1> pose_6vec;
            for (int j=0; j<6; j++)
            {
                
                pose_6vec[j] = pose[6*i + j]; // Control to Front Keyframe
                
            }
            
#endif


#if 0// USE new ref to control
            controlPoses_T[idxSplineCtrlPts_ + i] // new RefToControl
                    = SE3T::exp(pose_6vec) * refToWorld_T.inverse();
#else// USE ref to control
            
#if 0 // USE WORLD_TO_CONTROL for pose (ORIGINAL IMPLEMENTATION)
            the
            controlPoses_T[idxSplineCtrlPts_ + i] // RefToControl
                    = SE3T::exp(pose_6vec) * refToWorldT; // WorldToCntl*RefToWorld
            
#else // USE CONTROL TO FRONT KEYFRAME
            
            // ControlToFrontKeyframe^{-1} is front_keyframe to Control
            // RefToControl
            //     = FrontKeyframeToControl * WorldToFrontKeyFrame * RefToWorld
            SE3T fkToControl_T = SE3T::exp(pose_6vec).inverse();
            controlPoses_T[idxSplineCtrlPts_ + i] // RefToControl
                    = fkToControl_T * worldToFrontKeyframe_T * refToWorldT;
            
#endif
            
#endif

        }


        // Instantiate Spline object with typename T (Ref To Cntl)
        lsd_slam::Spline<T> spline(controlPoses_T, SPLINE_K);

        //----------------------------------------------------------------------
        // Determine the rows for each control pose
        // "i" is the interval segment in B-Spline, and predefined as
        // i - 1, i, i + 1 and i + 2 become 0, height/3, height/3*2 and height,
        // where "height" is the image height at each pyramid level

        // Determine the interval time u and i parameter for Bspline
        // from y-coordinate
//        T ycoord = qycoord[0];
//        T u_time;
//        int i_spline;
//        determineTimeIntervalFromYcoord(ycoord, &u_time, &i_spline);
        T u_time = uTimes[0];

        //----------------------------------------------------------------------

#if 0
//        double u_time_scalar = ceres::JetOps<T>::GetScalar(u_time);

//        if ((offsetRow_ != 0) &&
//                ((u_time_scalar < 0) || (u_time_scalar > 1.0)))
        if ((offsetRow_ != 0) &&
            (fabs((*frameList_)[idxSrcSegment_]->timestamp() - 0.06) < 0.001))
        {

            printf("timestamp %f\n", (*frameList_)[idxSrcSegment_]->timestamp());
            printf("idxSrcSegment_ %d\n", idxSrcSegment_);
            printf("level_ %d\n", level_);
            printf("offsetRow %d\n", offsetRow_);
            printf("ycoord %f; u_time %f; i_spline %d\n",
                   ceres::JetOps<T>::GetScalar(ycoord),
                   ceres::JetOps<T>::GetScalar(u_time), i_spline);
            printf("refPoint_ %f %f %f\n",
                   (*refPoint_)[0], (*refPoint_)[1], (*refPoint_)[2]);
            printf("refPoint_ in image %f %f\n",
                   ((*refPoint_)[0]/(*refPoint_)[2]) *
                   (*frameList_)[idxSrcSegment_]->K(level_)(0,0) +
                   (*frameList_)[idxSrcSegment_]->K(level_)(0,2),
                   ((*refPoint_)[1]/(*refPoint_)[2]) *
                   (*frameList_)[idxSrcSegment_]->K(level_)(1,1) +
                   (*frameList_)[idxSrcSegment_]->K(level_)(1,2));
            printf("y-image point + offset  = %f\n",
                   ((*refPoint_)[1]/(*refPoint_)[2]) *
                   (*frameList_)[idxSrcSegment_]->K(level_)(1,1) +
                   (*frameList_)[idxSrcSegment_]->K(level_)(1,2) +
                   offsetRow_);


        }
#endif

        // u_time for Spline is out of range, therefore return maximum residual
        if (u_time < T(0.0) || u_time > T(1.0))
        {
#if 1// MAX RESIDUAL
            
#if 1 // Squared distance
            residual[0] = T(OUTLIER_RESIDUAL_INTENSITY_BN3);
            residual[1] = T(OUTLIER_RESIDUAL_ROWDIFF_BN3);

#else // L1-norm
            residual[0] = T(sqrt(OUTLIER_RESIDUAL_INTENSITY_BN3));
            residual[1] = T(sqrt(OUTLIER_RESIDUAL_ROWDIFF_BN3));
#endif
            //residual[1] = T(OUTLIER_RESIDUAL_ROWDIFF_BN3);
#else // MIN RESIDUAL
            residual[0] = T(0);
            residual[1] = T(0);
#endif

#if DEBUG
            printf("residual[0] = %f [Out of range u_time]\n",
                ceres::JetOps<T>::GetScalar(residual[0]));
#endif
            return true;

        }
#if 1 //DEBUG - USE SPLINE
        //----------------------------------------------------------------------
        // Get a pose on the spline

        // Determine i_spline by startRow
        // NOTE: idxSplineCtrlPts contains index of control points for the
        //       currently optimizing new frame
        //       idxForFrame contains index of control point for
        //       stored frame which depends on the starting row
        int idxForFrame = (int)floor(startRowImage_ / rowPerSegment_)
                          + NUM_EXTRA_CTRLPTS;
#if 1 // COMPUTE FROM STARTROW AND ROWPERSEG
        int i_spline = idxForFrame;
#else // USE DIRECTLY FRON PARAMETER
        int i_spline = idxSplineCtrlPts_;
#endif

#ifdef DEBUG
        if ((frame_->id() >= 7) && (frame_->id() <= 11) &&
            (idx_ == 0))
        {

            printf("RefID %d; frameID %d; Level %d; idxForFrame = %d, "
                   "idxSplineCtrlPts_ = %d\n",
                   reference_->frameID,
                   frame_->id(),
                   level_,
                   idxForFrame, idxSplineCtrlPts_);

            for (int i=0; i<(int)controlPoses_T.size(); i++)
            {

                Eigen::Matrix<T,6,1> val = controlPoses_T[i].log();
                printf(" c%03d: %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f ",
                        i,
                        ceres::JetOps<T>::GetScalar(val[0]),
                        ceres::JetOps<T>::GetScalar(val[1]),
                        ceres::JetOps<T>::GetScalar(val[2]),
                        ceres::JetOps<T>::GetScalar(val[3]),
                        ceres::JetOps<T>::GetScalar(val[4]),
                        ceres::JetOps<T>::GetScalar(val[5]));

                // Print points (world wrt control)
                Eigen::Matrix<T,6,1> val2 =
                        (controlPoses_T[i]*refToWorldT.inverse()).log();
                printf(" |[w] %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
                       ceres::JetOps<T>::GetScalar(val2[0]),
                       ceres::JetOps<T>::GetScalar(val2[1]),
                       ceres::JetOps<T>::GetScalar(val2[2]),
                       ceres::JetOps<T>::GetScalar(val2[3]),
                       ceres::JetOps<T>::GetScalar(val2[4]),
                       ceres::JetOps<T>::GetScalar(val2[5]));

            }

        }
#endif
       //----------------------------------------------------------------------

        //int i_spline = idxSplineCtrlPts_;
        Eigen::Matrix<T,4,4> pose_on_spline =
                spline.getPoseOnSpline(u_time, i_spline); // (Ref to control)

        //----------------------------------------------------------------------
#else
        // FOR DEBUGGING - NO SPLINE USED
        Eigen::Matrix<T,6,1> pose_6vec;
        for (int j=0; j<6; j++)
            pose_6vec[j] = pose[j];
        Eigen::Matrix<T,4,4> pose_on_spline = SE3T::exp(pose_6vec).matrix();
#endif

        // Compute a residual by warping
        computeResidualByWarping(pose_on_spline, u_time, residual);

#if DEBUG
        printf("residual[0] = %f [Return]\n",
               ceres::JetOps<T>::GetScalar(residual[0]));
#endif
        return true;

    }

    template<class T>
    bool distortPoint_Template(T x_undist, T y_undist,
                               T *x_dist_out, T *y_dist_out) const
    {
        
        int out_width = undistorter->getOutputWidth();
        int out_height = undistorter->getOutputHeight();
        if (ceres::JetOps<T>::GetScalar(x_undist) < 0 ||
            ceres::JetOps<T>::GetScalar(x_undist) > out_width - 1 ||
            ceres::JetOps<T>::GetScalar(y_undist) < 0 ||
            ceres::JetOps<T>::GetScalar(y_undist) > out_height - 1)
        {
            
            // points are not in the image
            //printf("Errror: points are not inside the image\n");
            *x_dist_out = T(-1);
            *y_dist_out = T(-1);
            return false;
            
        }
        
        cv::Mat K_in = undistorter->getOriginalK();
        cv::Mat K_out = undistorter->getK();
        
        // Camera parameters
        // Distorted
        T fx_dist = T(K_in.at<double>(0, 0));
        T fy_dist = T(K_in.at<double>(1, 1));
        T cx_dist = T(K_in.at<double>(0, 2));
        T cy_dist = T(K_in.at<double>(1, 2));
        T fx_undist = T(K_out.at<double>(0, 0)); // Undistorted
        T fy_undist = T(K_out.at<double>(1, 1));
        T cx_undist = T(K_out.at<double>(0, 2));
        T cy_undist = T(K_out.at<double>(1, 2));
        T w = T(undistorter->getInputDistort());
        T d2t = T(2.0)*T(tan(w/T(2.0)));
        
        // Normalized distorted points
        T x_undist_norm = (x_undist - cx_undist)/fx_undist;
        T y_undist_norm = (y_undist - cy_undist)/fy_undist;
        
        // Radius of the normalized undistort pixels
        T ru = sqrt(x_undist_norm*x_undist_norm + y_undist_norm*y_undist_norm);
        
        // Computer factor ratio rd/ru to convert from undistort to distort
        T fac = T(1.0);
        if ((ru == T(0)) || (w == T(0)))
        {
            
            fac = T(1.0);
            
        }
        else
        {
            
            fac = T(atan(ru*d2t)/(w*ru)); // Ratio rd/ru
            
        }
        
        // Distorted pixels (input)
        T x_dist = fx_dist*fac*x_undist_norm + cx_dist;
        T y_dist = fy_dist*fac*y_undist_norm + cy_dist;
        
        
        // Output
        *x_dist_out = x_dist;
        *y_dist_out = y_dist;
        return true;
        
    }
    
    template<typename T>
    bool distortPointLevel_Template(T x_undist, T y_undist,
                                    T *x_dist, T *y_dist,
                                    int level) const
    {
        
        // Image scale by the level
        T scale = T(pow(2.0, level));
        
        // Input to original scale
#if 0
        T x_undist_org = (x_undist + T(0.5))*scale - T(0.5);
        T y_undist_org = (y_undist + T(0.5))*scale - T(0.5);
#else
        T x_undist_org = (x_undist)*scale;
        T y_undist_org = (y_undist)*scale;
#endif
        
        T x_dist_org, y_dist_org;
        bool found = distortPoint_Template<T>(x_undist_org, y_undist_org,
                              &x_dist_org, &y_dist_org);
        if (found == true)
        {
            
            // Output to level size
#if 0
            *x_dist = (x_dist_org + T(0.5))/scale - T(0.5);
            *y_dist = (y_dist_org + T(0.5))/scale - T(0.5);
#else
            *x_dist = (x_dist_org)/scale;
            *y_dist = (y_dist_org)/scale;
#endif
            return true;
            
        }
        else
        {
            return false;
        }
        
    }
    
    //--------------------------------------------------------------------------
    template <typename T>
    void computeResidualByWarping(Eigen::Matrix<T,4,4> pose_on_spline,
                                T u_time, T* residual) const
    {
        double alpha = ALPHA_BN3;

        residual[0] = T(0.0);
        residual[1] = T(0.0);
        
        // Define SE3Group based on template typename T
        typedef Sophus::SE3Group<T> SE3T;

        //----------------------------------------------------------------------
        // A 3D point warped by the pose vector

        SE3T pose_SE3T = SE3T(pose_on_spline);

        // Warp by rotation and translation
        Eigen::Matrix<T,3,3> rotMat = pose_SE3T.rotationMatrix();
        Eigen::Matrix<T,3,1> transVec = pose_SE3T.translation();
        Eigen::Matrix<T,3,1> pointRef;
        pointRef[0] = T((refPoint_)[0])*T(scale_);
        pointRef[1] = T((refPoint_)[1])*T(scale_);
        pointRef[2] = T((refPoint_)[2])*T(scale_);
#if 0 // NO_ROTATION_ESTIMATE
        // NO ROTATION - only translation
        Eigen::Matrix<T,3,1> Wxp = pointRef + transVec;
#else
        Eigen::Matrix<T,3,1> Wxp = rotMat*pointRef + transVec;
#endif

#if  DEBUG
        // Print out point
        printf("pointRef[0] = %f\n", ceres::JetOps<T>::GetScalar(pointRef[0]));
        printf("pointRef[1] = %f\n", ceres::JetOps<T>::GetScalar(pointRef[1]));
        printf("pointRef[2] = %f\n", ceres::JetOps<T>::GetScalar(pointRef[2]));

        // Print out rot and trans
        printf("rotMat = \n");
        for (int i=0; i<3; i++)
        {
            for (int j=0; j<3; j++)
            {
                printf("%f ", ceres::JetOps<T>::GetScalar(rotMat(i,j)));
            }
            printf("\n");
        }
        printf("trans = %f %f %f\n",
               ceres::JetOps<T>::GetScalar(transVec(0)),
               ceres::JetOps<T>::GetScalar(transVec(1)),
               ceres::JetOps<T>::GetScalar(transVec(2)));
#endif

        //----------------------------------------------------------------------
        // Projection of the 3D point in the (target) image

        // Camera parameters in the given pyramid level
        T fx_l = T(frame_->K(level_)(0,0));
        T fy_l = T(frame_->K(level_)(1,1));
        T cx_l = T(frame_->K(level_)(0,2));
        T cy_l = T(frame_->K(level_)(1,2));

#if 0
        printf("camera parameters: \n");
        printf(" fx_l = %f\n", ceres::JetOps<T>::GetScalar(fx_l));
        printf(" fy_l = %f\n", ceres::JetOps<T>::GetScalar(fy_l));
        printf(" cx_l = %f\n", ceres::JetOps<T>::GetScalar(cx_l));
        printf(" cy_l = %f\n", ceres::JetOps<T>::GetScalar(cy_l));
#endif

        // New image coordinates warped and projected in the target image
        T u_new = (Wxp[0]/Wxp[2])*fx_l + cx_l;
        T v_new = (Wxp[1]/Wxp[2])*fy_l + cy_l;

        // Translate the new image coordiantes according to the source
        // segment image
        double offsetRow_ = 0.0;
        T v_new_offset = v_new + T((double)offsetRow_);

#if DEBUG
        printf("u_new = %f\n", ceres::JetOps<T>::GetScalar(u_new));
        printf("v_new = %f\n", ceres::JetOps<T>::GetScalar(v_new));
        printf("offsetRow_ = %f\n", (double)offsetRow_);
        printf("v_new_offset = %f\n", ceres::JetOps<T>::GetScalar(v_new_offset));
#endif

        //----------------------------------------------------------------------
        // Check if the warped image point is wihtin the image frame

        // Image width and height for the given pyramid level
        int width_l = frame_->width(level_);
        int height_l = frame_->height(level_);

        // Boundary check (NOTE: the boundary is defined as
        // [1:width_l-1, 1:height_l-1] excluding the one-pixel boundary
        // in the original level image size.
        double u = ceres::JetOps<T>::GetScalar(u_new);
        double v_offset = ceres::JetOps<T>::GetScalar(v_new_offset);
        if ((u <= 1) ||
            (v_offset <= 1) ||
            (u >= (frame_->width(level_)) - 2) ||
            (v_offset >= (frame_->height(level_) - 2)))
        {

            // Out of target image boundary
            // We return zero (the lowest intensity) in this case.
            //WarpResidual::numWarpedPixelOutOfImage++;
#if 0//DEBUG
            printf("Invalid warping - out of boundary"
                   "(u:%f v_offset:%f) for u:[%f %f] v_offset[%f %f]\n",
                   u, v_offset,
                   1.0, (double)(imageSeg_.cols >> level_) - 2,
                   1.0, (double)(imageSeg_.rows >> level_) - 2);
#endif

#if 0//DEBUG
            printf("    pointRef[0] = %f\n", ceres::JetOps<T>::GetScalar(pointRef[0]));
            printf("    pointRef[1] = %f\n", ceres::JetOps<T>::GetScalar(pointRef[1]));
            printf("    pointRef[2] = %f\n", ceres::JetOps<T>::GetScalar(pointRef[2]));
            printf("    u_new = %f\n", ceres::JetOps<T>::GetScalar(u_new));
            printf("    v_new = %f\n", ceres::JetOps<T>::GetScalar(v_new));
            printf("    offsetRow_ = %f\n", (double)offsetRow_);
            printf("    v_new_offset = %f\n", ceres::JetOps<T>::GetScalar(v_new_offset));
#endif
            //residual[0] = T(0);
            //residual[1] = T(0); //v_new_offset - T((imageSeg_.rows >> level_) - 1)*u_time;
#if 1// MAX RESIDUAL
            
#if 1 // Squared distance
            residual[0] = T(OUTLIER_RESIDUAL_INTENSITY_BN3);
            residual[1] = T(OUTLIER_RESIDUAL_ROWDIFF_BN3);

#else // L1-norm
            residual[0] = T(sqrt(OUTLIER_RESIDUAL_INTENSITY_BN3));
            residual[1] = T(sqrt(OUTLIER_RESIDUAL_ROWDIFF_BN3));
#endif

#else // MIN RESIDUAL
            residual[0] = T(0);
            residual[1] = T(0);
#endif

            return;

        }

        // Check if any NaN value in u and v
        if ((isnan(u) == true)|| (isnan(v_offset) == true))
        {

            //printf("NaN value detected\n");
#if 1 // MAX RESIDUAL
            
#if 1 // Squared-distance
            residual[0] = T(OUTLIER_RESIDUAL_INTENSITY_BN3);
            residual[1] = T(OUTLIER_RESIDUAL_ROWDIFF_BN3);

#else // L1-norm
            residual[0] = T(sqrt(OUTLIER_RESIDUAL_INTENSITY_BN3));
            residual[1] = T(sqrt(OUTLIER_RESIDUAL_ROWDIFF_BN3));
#endif
            //residual[1] = T(OUTLIER_RESIDUAL_ROWDIFF_BN3);
#else // MIN RESIDUAL
            residual[0] = T(0);
            residual[1] = T(0);
#endif
            return;

        }

        //----------------------------------------------------------------------
        // Get intensity values in the source and target image

        // Gradient for the given pyramid level in the reference (source) image
        const Eigen::Vector4f* frame_gradients = frame_->gradients(level_);

        //----------------------------------------------------------------------
        // Interpolated gradients and intensity of the warped point in
        // the new (target) image.
        // Returns:
        //    resInterp[0] : interpolated gradient x
        //    resInterp[1] : interpolated gradient y
        //    resInterp[2] : intensity/colour value

        // v coordinate for source image is determined
        // If offset is positive, v_offset needs a correction
        // E.g. v = 160, offset = 320, v_offset = 480
        //      but v_source should be 160.
        T v_new_source = (offsetRow_ > 0) ?
                    v_new_offset - T(offsetRow_) : v_new_offset;
        double v_source = ceres::JetOps<T>::GetScalar(v_new_source);

        // Check if v_source is in the boundary
        if ((v_source < 1) || (v_source > height_l - 2))
        {

            //printf("V_source is out of boundary\n");
#if 1 // MAX RESIDUAL
            
#if 1 // Squared distance
            residual[0] = T(OUTLIER_RESIDUAL_INTENSITY_BN3);
            residual[1] = T(OUTLIER_RESIDUAL_ROWDIFF_BN3);
            
#else // L1-norm
            residual[0] = T(sqrt(OUTLIER_RESIDUAL_INTENSITY_BN3));
            residual[1] = T(sqrt(OUTLIER_RESIDUAL_ROWDIFF_BN3));
#endif
            
#else // MIN RESIDUAL
            residual[0] = T(0);
            residual[1] = T(0);
#endif
            return;

        }

        Eigen::Vector3f resInterp =
                getInterpolatedElement43(frame_gradients,
                                         (float)u, (float)v_source, width_l);

        // Intensity value of the pixel in the source image
        float c1 = (refColourVariance_)[0];

        // Intensity value of the pixel in the target image
        // Take a chain rule
        double sample[3];
        sample[0] = double(resInterp[2]); // Intensity
        sample[1] = double(resInterp[0]); // Gradient x
        sample[2] = double(resInterp[1]); // Gradient y
        T uv[2] = { u_new, v_new_source };
        // f = sample[0]
        // dfdx = sample + 1
        // z = uv
        // ==> dfdx * dx*d(uv) : Derivative of intensity with respect to
        //                       the warped coordinates u and v
        T c2 = ceres::Chain<double, 2, T>::Rule(sample[0], sample + 1, uv);

#if 0
        // Check if c2 is valid, otherwise return a maximum residual
        if (ceres::JetOps<T>::GetScalar(c2) > 255 ||
            ceres::JetOps<T>::GetScalar(c2) < 0)
        {

            residual = T(OUTLIER_RESIDUAL_INTENSITY_BN3);
            return residual;

        }
#endif

        //----------------------------------------------------------------------
        // Inverse depth variance
        T c3 = T(refColourVariance_[1]);

        //----------------------------------------------------------------------
        // Partial derivative of intensity error wrt depth
        T tx = transVec[0];
        T ty = transVec[1];
        T tz = transVec[2];
        T px = Wxp[0];
        T py = Wxp[1];
        T pz = Wxp[2];
        T d = T(1.0/refPoint_[2]);
        T g0 = (tx * pz - tz * px) / (pz*pz*d);
        T g1 = (ty * pz - tz * py) / (pz*pz*d);
        T gx = fx_l*T(resInterp[0]);
        T gy = fy_l*T(resInterp[1]);
        T drpdd = gx * g0 + gy * g1;

        // Partial derivative of row-difference error wrt depth
        T Ns = T(rowPerSegment_);
        T drrdd = T(alpha)/(Ns*d)*g1;
        
        //----------------------------------------------------------------------
        // This is an error for various segHz and camHz
//        double v_new_offset_scalar = ceres::JetOps<T>::GetScalar(v_new_offset);
//        double uaccForRow = (v_new_offset_scalar + 0
//                             + (((int)startRowImage_) >> level_)
//                             - (((int)startRowSeg_) >> level_))
//        / (((int)rowPerSegment_) >> level_);
//        double uValue = uaccForRow;
        
        // Update v_new if radial distortion exists
        double u_new_scalar = ceres::JetOps<T>::GetScalar(u_new);
        double v_new_scalar = ceres::JetOps<T>::GetScalar(v_new);
        T v_new_update = v_new;
#if 0
        if (undistorter->isNoRectification() == false &&
            methodRollingShutter_ == 1)
        {
            
            // Update yRef via undistortion
            float x_dist, y_dist;
            undistorter->distortPointLevel(u_new_scalar, v_new_scalar,
                                      &x_dist, &y_dist, level_);
            v_new_update = T(y_dist);
        }
#else // USE Template
        bool isRowWithinImage = true;
        if (undistorter->isNoRectification() == false &&
            methodRollingShutter_ == 1)
        {
            
            // Update yRef via undistortion
            T x_dist, y_dist;
            isRowWithinImage = distortPointLevel_Template(u_new, v_new,
                                       &x_dist, &y_dist, level_);
            v_new_update = y_dist;
            
        }
#endif
        
        T uaccForRow_T = (v_new_update + T(0)
                          + T(((int)startRowImage_) >> level_)
                          - T(((int)startRowSeg_) >> level_))
        / T(((int)rowPerSegment_) >> level_);
        T uValue_T = uaccForRow_T;
        double uValue = ceres::JetOps<T>::GetScalar(uValue_T);

        // Partial derivative of row-difference error wrt time
        T drrdu = T(-alpha*uValue);
        
        //----------------------------------------------------------------------
        // A Huber weight to residual
        // Sigma_residual^2 = cameraPixelNoise2 + s1*drpdd*drpdd
        //                    + s1*c3*drrdd
        //                    + s2*drrdu
        // where s1 = var_weight_depth * variance_inverse_depth
        //       s2 = var_weight_row_time * variance_row_time
        
        // Normalized variance of photometric error wrt depth
        T s1 = T(settings_.var_weight) * c3;
        float VAR_ROW_TIME = 1.0f / ((int)rowPerSegment_ >> level_);
        T s2 = T(settings_.var_weight * VAR_ROW_TIME);
        
#if 0 // NORMALIZATION BY INTENSITY VARIANCE
        
        T w_p = T(1.0)/(T(cameraPixelNoise2) + s1*drpdd*drpdd; // 1/Sigma_res^2
        
#else // NORMALIZATION BY INTENSITY AND TIME VARIANCE
        
        T w_p = T(1.0)/(T(cameraPixelNoise2) + s1*drpdd*drpdd
                        + s1*drrdd*drrdd
                        + s2*drrdu*drrdu); // 1/Sigma_res^2
        
#endif
        
        T rp = (T(c1) - c2);
        T weighted_rp2 = rp*rp*w_p; // rp^2/Sigma_res^2
        T weighted_rp = sqrt(weighted_rp2); // rp/Sigma_res
        T wh;
        if (ceres::JetOps<T>::GetScalar(weighted_rp) < settings_.huber_d/2.0)
        {

            // wh is one if (rp/Sigma_res < hd/2),
            wh = T(1.0);

        }
        else
        {

            // wh is (hd/2)*(Sigma_res/rp), otherwise
            wh = (T(settings_.huber_d/2.0) / weighted_rp);

        }
        
        //----------------------------------------------------------------------
        // Residual error of intensity between the source and target pixel
        // normalized by the inverse depth variance
#if 0 // USE intensity residual with normalized variance

#if 0// LSD-SLAM::HUBER + NULL loss function
        
        // intensityResidual is rp/Sigma_res if (rp/Sigma_res < hd/2)
        //                   is (hd/2)*(Sigma_res/rp) otherwise
        T intensityResidual = sqrt(wh*w_p*rp*rp);
        // Check: sqrt(wh*wp*rp*rp) = sqrt((1.0)*(1/Sigma_res^2)*rp*rp)
        //                          = rp/Sigma_res, if (rp/Sigma_res < hd/2)
        //        sqrt(wh*wp*rp*rp) = sqrt((hd/2)*(Sigma_res/rp)
        //                                  *(1/Sigma_res^2)*rp*rp)
        //                          = sqrt((hd/2)*(rp/Sigma_res), otherwise.
        //
        // Note: NULL loss required, since err = intensityResidual^2
        
#else // Ceres::Huber
        
        // intensityResidual is rp*sqrt(w_p) = rp*sqrt(1/Sigma_res^2)
        //                                   = rp/Sigma_res
        //
        // Note: ceres::HuberLoss((hd/2)^2) is required in AddResidualBlock()
        //       since err = (rp/Sigma_res)^2 if (rp/Sigma_res)^2 < (hd/2)^2
        //                 or
        //                 = 2*sqrt((rp/Sigma_res)^2) - 1
        //                 = 2*(rp/Sigma_res) - 1 if (rp/Sigma_res)^2 > (hd/2)^2
        //       NOTE: this is different from the Huber in the original LSDSLAM
        //             implementation
        T intensityResidual = rp*sqrt(w_p); // rp/Sigma_res
        
#endif

#else // USE Huber intensity residual without normalized variance
        // Note: use ceres::HuberLoss() for loss function in AddResidualBlock

        T intensityResidual = rp;

#endif
        //residual[0] = intensityResidual;

#if 0// PRINT WEIGHTVAR
        //----------------------------------------------------------------------
        // Print out variance for pixels
        printf("WEIGHTVAR_Frame%08d_Level%d_Idx%d %f %f %f\n",
               frame_->id(),
               level_,
               idx_,
               u, v_source, ceres::JetOps<T>::GetScalar(weightVar));
#endif
        //----------------------------------------------------------------------
        // Residual error for difference of y-coordinate and utime

#if 0 // USE Huber row difference residual with normalized variance

        T rr = T(alpha)*(T(uValue) - u_time);
        T weighted_rr2 = rr*rr*w_p; // rr^2/Sigma_res^2
        T weighted_rr = sqrt(weighted_rr2); // rr/Sigma_res
        float HUBER_ROWDIFF = alpha*VAR_ROW_TIME/2.0f; // Huber paramter for row difference
        T wh_rr;
        if (ceres::JetOps<T>::GetScalar(weighted_rr) < HUBER_ROWDIFF/2.0)
        {
            
            // wh_r is one if (rp/Sigma_res < hd/2),
            wh_rr = T(1.0);
            
        }
        else
        {
            
            // wh is (hd/2)*(Sigma_res/rp), otherwise
            wh_rr = (T(HUBER_ROWDIFF/2.0) / weighted_rr);
            
        }
        
#if 0   // LSD-SLAM HUBER + Ceres::NULL
        T rowDiffResidual = sqrt(wh_rr*w_p*rr*rr);
        // NOTE: rowDiffResidual is sqrt(1.0*w_p*rr*rr)
        //                          = sqrt((1/Sigma_res^2)*(rr)*(rr))
        //                          = rr/Sigma_res
        //                       If (rr/Sigma_res < hd/2).
        //
        //                       is sqrt(wh_rr*w_p*rr*rr)
        //                          = sqrt((hd/2)/(rr/Sigma_res)
        //                            *(1/Sigma_res^2)*(rr)*(rr))
        //                          = sqrt((hd/2)*(Sigma_res/rr)
        //                            *(1/Sigma_res^2)*(rr)*(rr))
        //                          = sqrt((hd/2)*(rr/Sigma_res))
        //                       Otherwise.
        
#else   // Use Ceres::HuberLoss
        T rowDiffResidual = rr*sqrt(w_p);
        // rowDiffResidual is rr*sqrt(w_p) = rr*sqrt(1/Sigma_res^2)
        //                                 = rr/Sigma_res
        //
        // Note: ceres::HuberLoss((hd/2)^2) is required in AddResidualBlock()
        //       since err = (rr/Sigma_res)^2 if (rr/Sigma_res)^2 < (hd/2)^2
        //                 or
        //                 = 2*sqrt((rr/Sigma_res)^2) - 1
        //                 = 2*(rr/Sigma_res) - 1 if (rr/Sigma_res)^2 > (hd/2)^2
        //       NOTE: this is different from the Huber in the original LSDSLAM
        //             implementation
#endif

        
#else // USE row difference residual without normalized variance
                        
        T rowDiffResidual;
        if (isRowWithinImage == true)
        {
        
            rowDiffResidual = T(alpha)*(uValue_T - u_time);
            
        }
        else
        {
            
            intensityResidual = T(OUTLIER_RESIDUAL_INTENSITY_BN3);
            rowDiffResidual = T(OUTLIER_RESIDUAL_ROWDIFF_BN3);
            
        }

#endif
        //residual[1] = rowDiffResidual;
        
#if 1  // Squared distance

        residual[0] = intensityResidual;
        residual[1] = rowDiffResidual;

#else   // L1-norm
        // NOTE: Ceres cannot handle sqrt(0.0) therefore we need a case
        //       handling this by a condition checking
        float res0 = ceres::JetOps<T>::GetScalar(intensityResidual);
        float res1 = ceres::JetOps<T>::GetScalar(rowDiffResidual);
        if (res0 < 0.0)
        {
            res0 = -res0;
            intensityResidual = -intensityResidual;
        }
        if (res1 < 0.0)
        {
            res1 = -res1;
            rowDiffResidual = -rowDiffResidual;
        }
        residual[0] = (res0 < 1e-12) ?
            T(0.0) : sqrt(intensityResidual);
        residual[1] = (res1 < 1e-12) ?
            T(0.0) : sqrt(rowDiffResidual);
        
        residual[0] = sqrt(abs(intensityResidual));
        residual[1] = sqrt(abs(rowDiffResidual));

#endif

#if 0//DEBUG
        printf("v_new_offset_scalar = %f\n", v_new_offset_scalar);
        printf("uaccForRow = %f\n", uaccForRow);
        printf("           = (%f + (%d >> %d))/(%f >> %d)\n",
               v_new_offset_scalar, startRowImage_, level_,
               rowPerSegment_, level_);
        printf("  rowPerSegment = %f\n", rowPerSegment_);
        printf("  startRowImage = %d\n", startRowImage_);
        printf("uValue = %f\n", uValue);
        printf(" u_time: %f\n",
               (double)ceres::JetOps<T>::GetScalar(u_time));

#endif


#if 0//DEBUG

        double residual_0 = ceres::JetOps<T>::GetScalar(residual[0]);
        //double residual_1 = ceres::JetOps<T>::GetScalar(residual[1]);

        // Print out
        printf("residual = %f (c1: %f - c2: %f) at u:%f, v_source:%f",
               (double)residual_0,
               (double)c1,
               (double)ceres::JetOps<T>::GetScalar(c2),
               (double)u,
               (double)v_source);
        printf(" - v: %f, v_offset: %f, offset: %f",
               (double)ceres::JetOps<T>::GetScalar(v_new),
               (double)v_offset,
               (double)offsetRow_);
        printf(" height(level):%d, width(level):%d", width_l, level_);
        printf(" uValue: %f, u_time: %f",
               uValue, ceres::JetOps<T>::GetScalar(u_time));
        printf(" utime_diff: %.2f\n",
               uValue - ceres::JetOps<T>::GetScalar(u_time));
        //printf(" residual[1]: %.2f\n", residual_1);
        
#endif

        return;

    }

}; // end of class

struct se3PlusBN3 {
    template<typename T>
    bool operator()(const T* x, const T* delta, T* x_plus_delta) const
    {
        // x_plus_delta = plus(x, delta)
        typedef Sophus::SE3Group<T> SE3T;

        for (int i=0; i<SPLINE_K; i++)
        {
            Eigen::Matrix<T,6,1> x_vec;
            x_vec << T(x[6*i + 0]),
                     T(x[6*i + 1]),
                     T(x[6*i + 2]),
                     T(x[6*i + 3]),
                     T(x[6*i + 4]),
                     T(x[6*i + 5]);
            SE3T lieX = SE3T::exp(x_vec);

            Eigen::Matrix<T,6,1> delta_vec;
            delta_vec << T(delta[6*i + 0]),
                         T(delta[6*i + 1]),
                         T(delta[6*i + 2]),
                         T(delta[6*i + 3]),
                         T(delta[6*i + 4]),
                         T(delta[6*i + 5]);
            SE3T lieDelta = SE3T::exp(delta_vec);

            SE3T lie_x_plus_delta = lieX*lieDelta; // (Control to world)[Right]
            //SE3T lie_x_plus_delta = lieDelta*lieX; // (World to control) [Left]
            //SE3T lie_x_plus_delta = lieX; // DO NOTHING

            Eigen::Matrix<T,6,1> Xi_x_plus_delta = lie_x_plus_delta.log();
            x_plus_delta[6*i + 0] = Xi_x_plus_delta(0);
            x_plus_delta[6*i + 1] = Xi_x_plus_delta(1);
            x_plus_delta[6*i + 2] = Xi_x_plus_delta(2);
            x_plus_delta[6*i + 3] = Xi_x_plus_delta(3);
            x_plus_delta[6*i + 4] = Xi_x_plus_delta(4);
            x_plus_delta[6*i + 5] = Xi_x_plus_delta(5);
        }

        return true;
    }
};


// This UserCallback class runs at the end of iteration.
// This is for examining variables and status at every iteration end.
class BN3_UserCallback : public ceres::IterationCallback
{

public:

    const int level_;
    const TrackingReference *reference_;
    Frame *frame_;
    const std::vector<Sophus::SE3d> controlPoses_; // WorldToControl
    const double *pose_;
    const int idxSC_frontKF_;
    const int idxSplineCtrlPts_;
    const double *qycoord_;
    const double rowPerSegment_;
    const double startRowImage_;
    const double startRowSeg_;
    const double scale_;
    const DenseDepthTrackerSettings settings_;
    cv::Mat *debugImageWeightRollingShutter_;
    cv::Mat *debugImageResidualsRollingShutter_;
    cv::Mat *debugImageSecondFrameRollingShutter_;
    cv::Mat *debugImageOldImageWarpedRollingShutter_;
    cv::Mat *debugImageOldImageSourceRollingShutter_;
    // Image undistortion/distortion object
    Undistorter *undistorter_;
    const int methodRollingShutter_;
    const std::vector<int> bin_; // Bin for sparse pixels
    bool saveToFile_;

    
    // Debug info
    double debugInfo_sumSquareResOut;
    double debugInfo_lastResidual;
    double debugInfo_resOut;
    double debugInfo_weightVar;
    double debugInfo_rp;
    double debugInfo_source_x;
    double debugInfo_source_y;
    double debugInfo_target_u;
    double debugInfo_target_v;
    
    int numExtraCtrlPts_;
    
    double debug_xRef;
    double debug_yRef;
    double debug_yRef_update;
    double debug_uTime_refPoint;
    int debug_idxSC_frontKF;

public:

    BN3_UserCallback(const int pyramidLevel,
                 const TrackingReference *reference,
                 Frame *frame,
                 const std::vector<Sophus::SE3d> controlPoses, // WorldToControl
                 const double *pose,
                 const int idxSC_frontKF,
                 const int idxSplineCtrlPts,
                 const double *qycoord,
                 const double rowPerSegment,
                 const double startRowImage,
                 const double startRowSeg,
                 const double scale,
                 const DenseDepthTrackerSettings settings,
                 cv::Mat *debugImageWeightRollingShutter,
                 cv::Mat *debugImageResidualsRollingShutter,
                 cv::Mat *debugImageSecondFrameRollingShutter,
                 cv::Mat *debugImageOldImageWarpedRollingShutter,
                 cv::Mat *debugImageOldImageSourceRollingShutter,
                 Undistorter *undistorter,
                 const int methodRollingShutter,
                 const std::vector<int> bin, // Bin for indices of sparse pixels
                 bool saveToFile = false)
        : level_(pyramidLevel),
          reference_(reference),
          frame_(frame),
          controlPoses_(controlPoses),
          pose_(pose),
          idxSC_frontKF_(idxSC_frontKF),
          idxSplineCtrlPts_(idxSplineCtrlPts),
          qycoord_(qycoord),
          rowPerSegment_(rowPerSegment),
          startRowImage_(startRowImage),
          startRowSeg_(startRowSeg),
          scale_(scale),
          settings_(settings),
          debugImageWeightRollingShutter_(debugImageWeightRollingShutter),
          debugImageResidualsRollingShutter_(debugImageResidualsRollingShutter),
          debugImageSecondFrameRollingShutter_(debugImageSecondFrameRollingShutter),
          debugImageOldImageWarpedRollingShutter_(debugImageOldImageWarpedRollingShutter),
          debugImageOldImageSourceRollingShutter_(debugImageOldImageSourceRollingShutter),
          undistorter_(undistorter),
          methodRollingShutter_(methodRollingShutter),
          bin_(bin),
          saveToFile_(saveToFile)
    {
        
        this->numExtraCtrlPts_ = SPLINE_K - 1;

    }

    ceres::CallbackReturnType operator()
        (const ceres::IterationSummary& summary)

    {

        if (summary.step_is_successful == false)
        {

            // Nothing to do
            return ceres::SOLVER_CONTINUE;

        }
        else if (summary.iteration == 0)
        {
            lsd_slam::WarpResidual_BN3::lastGoodCount_ = 0;
            lsd_slam::WarpResidual_BN3::lastBadCount_ = 0;
            lsd_slam::WarpResidual_BN3::usageCount_ = 0.0;

            return ceres::SOLVER_CONTINUE;
        }
        else if (summary.step_is_successful == true)
        {

            // Evaluate the warp error and returns the statistics
            // related to number of good and badly warped pixels, and
            // the number of depth changed pixels

            // Counters for each source segment
            int goodCount = 0;
            int badCount = 0;
            float usageFractus = 0.0;
            getStatisticsForWarpedResidual<double>(&goodCount,
                                                   &badCount,
                                                   &usageFractus);
            lsd_slam::WarpResidual_BN3::lastGoodCount_ = goodCount;
            lsd_slam::WarpResidual_BN3::lastBadCount_ = badCount;
            lsd_slam::WarpResidual_BN3::usageCount_ = usageFractus;

#if DEBUG_MODE
            printf("      %d: cost:%.15g, WarpResidual::usageCount_ = %f\n",
                   summary.iteration, summary.cost, WarpResidual::usageCount_);
#endif

            return ceres::SOLVER_CONTINUE;
        }

        return ceres::SOLVER_CONTINUE;

    } // end of operator

    // Get statistics for warped residual given reference frame,
    // current frame and updated pose
    template <typename T>
    T getStatisticsForWarpedResidual(int *goodCount,
                                        int *badCount,
                                        float *usageFractus)
    {
        // Define SE3Group based on template typename T
        typedef Sophus::SE3Group<T> SE3T;

        // Set to zeros
        *goodCount = 0;
        *badCount = 0;
        *usageFractus = 0;

        Eigen::Vector3f* refPoint = reference_->posData[level_];
        Eigen::Vector3f* refPoint_max = refPoint +
                reference_->numData[level_];
        Eigen::Vector2f* refColourVariance =
                reference_->colorAndVarData[level_];

        const int numSemiPtsInLevel = reference_->numData[level_];
        const int numSourceImage = 1;

        T residual;

        double sumSquareResOut = 0.0;
        this->debugInfo_sumSquareResOut = sumSquareResOut;
        printf("{\n");
        
#if 0 // USE_SPARSE_STEP

        // Determine sparse skip pixels
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
            frame_->width(level_)*
            frame_->height(level_)*
            lsd_slam::pixelDensityForTracking)
        {
            
            SPARSE_STEP = 1;
            numSparsePts = numSemiPtsInLevel;
            
        }
        
        // For each available pixels, 
        for (int idx = 0;
             refPoint < refPoint_max - SPARSE_STEP;
             idx+=SPARSE_STEP,
             refPoint+=SPARSE_STEP,
             refColourVariance+=SPARSE_STEP)
        {
            
#else // USE_BIN_SORTED_MERGED
            
        int idx = -1;
        for (auto it: bin_)
        {

            // Init pointer to point and color with variance
            refPoint = reference_->posData[level_];
            refColourVariance = reference_->colorAndVarData[level_];
            
            // Pointing the sample out of the bin
            idx++;
            refPoint = refPoint + it;
            refColourVariance = refColourVariance + it;
            
#endif

            // Define SE3Group based on template typename T
            typedef Sophus::SE3Group<T> SE3T;

            // Scale the reference point by a scalar (e.g. 1.0)
            Eigen::Vector3f scaled_refPoint = (*refPoint)*(float)scale_;
            
            // Check if the refPoint contains nan
            if (isnan(scaled_refPoint[0]) == true ||
                isnan(scaled_refPoint[1]) == true ||
                isnan(scaled_refPoint[2]) == true)
            {
                
                // 3D point contains a NaN value therefore do not consider
                // in optimization for image warping
                continue;
                
            }

#if 0// Point wrt world

            // RefPoint should be transformed by its reference pose with respect
            // to the world coordinate system
            Eigen::Matrix<double,3,3> rotRef =
                    se3FromSim3(reference_->keyframe
                               ->getScaledCamToWorld()).rotationMatrix();
            Eigen::Matrix<double,3,1> transRef =
                    se3FromSim3(reference_->keyframe
                                ->getScaledCamToWorld()).translation();

            Eigen::Vector3f new_refPoint =
                    rotRef.cast<float>()*scaled_refPoint
                    + transRef.cast<float>();
#else // Point wrt reference
            Eigen::Vector3f new_refPoint = scaled_refPoint;
#endif
            
            //------------------------------------------------------------------
            // Get refToWorld for template T
            Spline<double> splineRef(controlPoses_, SPLINE_K); // World to control
            double xRef = ((new_refPoint)[0]/(new_refPoint)[2]) *
                            frame_->K(level_)(0,0) +
                            frame_->K(level_)(0,2);
            double yRef = ((new_refPoint)[1]/(new_refPoint)[2]) *
                            frame_->K(level_)(1,1) +
                            frame_->K(level_)(1,2);
            double yRef_update = yRef;
            if (this->undistorter_->isNoRectification() == false &&
                methodRollingShutter_ == 1)
            {
                
                // Update yRef via undistortion
                float x_dist, y_dist;
                this->undistorter_->distortPointLevel((float)xRef, (float)yRef,
                                                      &x_dist, &y_dist,
                                                      level_);
                yRef_update = (double)y_dist;
                
            }
            
            double uTime_RefPoint =
                    (yRef_update + 0
                     + ((int) reference_->keyframe->startRowImage >> level_)
                     - ((int) reference_->keyframe->startRowSeg >> level_))
                    / ((int) reference_->keyframe->rowPerSegment >> level_);
            int idxSC_frontKF =
                    (int)floor(reference_->keyframe->startRowImage
                       / reference_->keyframe->rowPerSegment)
                       + NUM_EXTRA_CTRLPTS;
            
            this->debug_xRef = xRef;
            this->debug_yRef = yRef;
            this->debug_yRef_update = yRef_update;
            this->debug_uTime_refPoint = uTime_RefPoint;
            this->debug_idxSC_frontKF = idxSC_frontKF;
            
            Sophus::SE3d rowPose_RefPoint = // World to RowRef
                    Sophus::SE3d(splineRef.getPoseOnSpline(uTime_RefPoint,
                                                idxSC_frontKF));
            SE3 refToWorld_se3 = rowPose_RefPoint.inverse(); // row ref to world
            Eigen::Matrix<T,6,1> refToWorld_se3T_log;
            for (int i=0; i<6; i++)
            {
    
                refToWorld_se3T_log[i] = T(refToWorld_se3.log()[i]);
    
            }
            SE3T refToWorld = SE3T::exp(refToWorld_se3T_log);
    
            //------------------------------------------------------------------
            // Declare Spline in template T

            // Copy all control poses from array to SE3T vector
            std::vector<SE3T> controlPoses_T;
            for (int i=0; i<(int)controlPoses_.size(); i++)
            {

                // World To Control in template T
                Eigen::Matrix<T,6,1> controlPoses_se3T_log;
                for (int j=0; j<6; j++)
                {

                    controlPoses_se3T_log[j] = T((controlPoses_)[i].log()[j]);

                }
                SE3T cpT = SE3T::exp(controlPoses_se3T_log);

                // Convert to refToControl
                SE3T cpT_refToControl = cpT * refToWorld;

                // Store
                controlPoses_T.push_back(cpT_refToControl);

            }
            
#if 0 // USE WORLD To CONTROL pose
#else // USE CONTROL TO FRONT KEYFRAME
            
            // World to Front keyframe
            Eigen::Matrix<T,6,1> worldToFrontKeyframe_se3T_log;
            for (int j=0; j<6; j++)
            {
                
                // ControlPoses_[idxSC_frontKF_] is a control point related to
                // the front keyframe
                worldToFrontKeyframe_se3T_log[j]
                = T((controlPoses_)[idxSC_frontKF_].log()[j]);
                
            }
            SE3T worldToFrontKeyframe_T
                = SE3T::exp(worldToFrontKeyframe_se3T_log);
#endif

            // Copy k control poses from array to SE3T vector
            for (int i=0; i<SPLINE_K; i++)
            {

#if 0 // USE WORLD_TO_CONTROL for pose (ORIGINAL)
                
                Eigen::Matrix<T,6,1> pose_6vec;
                for (int j=0; j<6; j++)
                    pose_6vec[j] = pose_[6*i + j]; // WorldToControl
                
                // Set control poses as RefToControl
                controlPoses_T[idxSplineCtrlPts_ + i] // RefToControl
                        = SE3T::exp(pose_6vec) * refToWorld;
                
#else // USE CONTRL TO FRONT KEYFRAME
                
                Eigen::Matrix<T,6,1> pose_6vec;
                for (int j=0; j<6; j++)
                    pose_6vec[j] = pose_[6*i + j]; // Control to Front keyframe
                SE3T fkToControl_T = SE3T::exp(pose_6vec).inverse();
                
                // Set control poses as RefToControl
                controlPoses_T[idxSplineCtrlPts_ + i] // RefToControl
                    = fkToControl_T * worldToFrontKeyframe_T * refToWorld;
                
#endif

            }

            // Instantiate Spline object with typename T
            lsd_slam::Spline<T> spline(controlPoses_T, SPLINE_K);

//            // Determine the interval time u and i parameter for Bspline
//            // from y-coordinate
//            T ycoord = T(qycoord_[idxSrcSegment*numSemiPtsInLevel + idx]);
//            T u_time;
//            int i_spline  = 0;
//            determineTimeIntervalFromYcoord<T>(idxSrcSegment,
//                                              ycoord, &u_time, &i_spline);
            T u_time = T(qycoord_[idx]);

            // u_time for Spline is out of range, therefore return maximum residual
            if (u_time < T(0.0) || u_time > T(1.0))
            {

                residual = T(0);
                continue;

            }

#if 1 // USE SPLINE
            // Get a pose on the spline
            int i_spline = (int)floor(startRowImage_ / rowPerSegment_)
                           + NUM_EXTRA_CTRLPTS;
            
#if 0 // USE IDX SPLINE
            Eigen::Matrix<T,4,4> pose_on_spline =
            spline.getPoseOnSpline(u_time, idxSplineCtrlPts_);
#else
            Eigen::Matrix<T,4,4> pose_on_spline =
                    spline.getPoseOnSpline(u_time, i_spline);
#endif
#else
            Eigen::Matrix<T,6,1> pose_6vec;
            for (int j=0; j<6; j++)
                pose_6vec[j] = pose_[j];
            Eigen::Matrix<T,4,4> pose_on_spline = SE3T::exp(pose_6vec).matrix();
#endif
            //----------------------------------------------------------------------
            // A 3D point warped by the pose vector

            SE3T pose_SE3T = SE3T(pose_on_spline);

            // Warp by rotation and translation
            Eigen::Matrix<T,3,3> rotMat = pose_SE3T.rotationMatrix();
            Eigen::Matrix<T,3,1> transVec = pose_SE3T.translation();
            Eigen::Matrix<T,3,1> pointRef;
            pointRef[0] = T((new_refPoint)[0]);
            pointRef[1] = T((new_refPoint)[1]);
            pointRef[2] = T((new_refPoint)[2]);
#if 0 //NO_ROTATION_EST
            Eigen::Matrix<T,3,1> Wxp = pointRef + transVec;
#else
            Eigen::Matrix<T,3,1> Wxp = rotMat*pointRef + transVec;
#endif
            //------------------------------------------------------------------
            // Projection of the 3D point in the (target) image

            // Camera parameters in the given pyramid level
            T fx_l = T(frame_->K(level_)(0,0));
            T fy_l = T(frame_->K(level_)(1,1));
            T cx_l = T(frame_->K(level_)(0,2));
            T cy_l = T(frame_->K(level_)(1,2));

            // New image coordinates warped and projected in the target image
            T u_new = (Wxp[0]/Wxp[2])*fx_l + cx_l;
            T v_new = (Wxp[1]/Wxp[2])*fy_l + cy_l;

            // Translate the new image coordiantes according to the source
            // segment image
//            int offset = imageSeg_.offsetRows[idxSrcSegment] >> level_;
            int offset = 0;
#ifdef DEBUG
            printf("v_new = %f\n", v_new);
            printf("offset = %d\n", offset);
#endif
            T v_new_offset = v_new + T(offset);
#ifdef DEBUG
            printf("v_new + offset = %f\n", v_new);
            exit(1);
#endif
            //------------------------------------------------------------------
            // Check if the warped image point is wihtin the image frame

            // Image width and height for the given pyramid level
            int width_l = frame_->width(level_);
            int height_l = frame_->height(level_);

            // Boundary check (NOTE: the boundary is defined as
            // [1:width_l-1, 1:height_l-1] excluding the one-pixel boundary
            // in the original level image size.
            double u = T(u_new);
            double v = T(v_new);
            double v_offset = T(v_new_offset);
            if ((u <= 1) ||
                (v_offset <= 1) ||
                (u >= (frame_->width(level_) - 2)) ||
                (v_offset >= (frame_->height(level_) - 2)))
            {

                // Out of target image boundary
                // We return zero (the lowest intensity) in this case.
                residual = T(0);

                continue;

            }

            // Check if any NaN value in u and v
            if ((isnan(u) == true)|| (isnan(v_offset) == true))
            {
#if DEBUG
                printf("NaN value detected\n");
                printf(" pointRef: %f %f %f, ",
                       (double)pointRef[0],
                       (double)pointRef[1],
                       (double)pointRef[2]);
                printf("Dim of qycoord = %d",
                       (int)(numSemiPtsInLevel*numSourceImage));
                printf(" qycoord[%d] = %f\n",
                       (int)(idx),
                       qycoord_[idx]);
                printf("u_time: %f, i_spline: %f\n",
                       (double)u_time,
                       (double)i_spline);
                printf("pose_on_spline = \n");
                std::cout << pose_on_spline << std::endl;
                printf("rotMat = ");
                std::cout << rotMat << std::endl;
                printf("transVec = ");
                std::cout << transVec << std::endl;
                printf(" Wxp: %f %f %f, ",
                       (double)Wxp[0],
                       (double)Wxp[1],
                       (double)Wxp[2]);
                printf("  u: %f, v:%f, v_offset:%f\n",
                       (double)u,
                       (double)v,
                       (double)v_offset);

                for (int i=0; i<SPLINE_K; i++)
                printf("pose %d: %.2f %.2f %.2f %.2f %.2f %.2f\n", i,
                   pose_[6*i+0], pose_[6*i+1], pose_[6*i+2],
                   pose_[6*i+3], pose_[6*i+4], pose_[6*i+5]);
#endif
                residual = T(0);

                continue;

            }

            //------------------------------------------------------------------
            // Get intensity values in the source and target image

            // Gradient for the given pyramid level in the reference image
            const Eigen::Vector4f* frame_gradients =
                    frame_->gradients(level_);

            //------------------------------------------------------------------
            // Interpolated gradients and intensity of the warped point in
            // the new (target) image.
            // Returns:
            //    resInterp[0] : interpolated gradient x
            //    resInterp[1] : interpolated gradient y
            //    resInterp[2] : intensity/colour value

            // v coordinate for source image is determined
            // If offset is positive, v_offset needs a correction
            // E.g. v = 160, offset = 320, v_offset = 480
            //      but v_source should be 160.
            double v_source = (offset > 0) ? v_offset - offset : v_offset;

            // Check if v_source is in the boundary
            if ((v_source < 1) || (v_source > height_l - 2))
            {

                residual = T(0);
                continue;

            }

            Eigen::Vector3f resInterp =
                    getInterpolatedElement43(frame_gradients,
                                             (float)u, (float)v_source, width_l);

            //------------------------------------------------------------------
            // Residual between intensities from source and target
            T c1 = T((*refColourVariance)[0]);
            T c2 = T(resInterp[2]);
            
#if 1 // RESIDUAL WITHOUT WEIGHT OR HUBER
            
            residual = c1 - c2;
            
#else // RESIDUAL WITH WEIGHT OR HUBER
            
            //------------------------------------------------------------------
            // Inverse depth variance
            T c3 = T((*refColourVariance)[1]);
            
            //------------------------------------------------------------------
            // Partial derivative of intensity error wrt depth
            T tx = transVec[0];
            T ty = transVec[1];
            T tz = transVec[2];
            T px = Wxp[0];
            T py = Wxp[1];
            T pz = Wxp[2];
            T d = T(1.0/pointRef[2]);
            T g0 = (tx * pz - tz * px) / (pz*pz*d);
            T g1 = (ty * pz - tz * py) / (pz*pz*d);
            T gx = fx_l*T(resInterp[0]);
            T gy = fy_l*T(resInterp[1]);
            T drpdd = gx * g0 + gy * g1;
            
            //------------------------------------------------------------------
            // A Huber weight to residual
            // Normalized variance
            T w_p = T(1.0)/(T(cameraPixelNoise2) + c3*drpdd*drpdd);
            T rp = (T(c1) - c2);
            T weighted_rp2 = rp*rp*w_p;
            T weighted_rp = sqrt(weighted_rp2);
            T wh;
            if (weighted_rp < settings_.huber_d/2)
            {
                
                wh = T(1);
                
            }
            else
            {
                
                wh = (T(settings_.huber_d/2) / weighted_rp);
                
            }
            
            //T weightVar = 1.0/sqrt(w_p);
            if (wh < 0) wh = -wh;
            
            T weightVar = wh * w_p;
            
            residual = weightVar*rp;
            
#endif

#if 0
            // Check if c2 is valid, otherwise return a maximum residual
            if (ceres::JetOps<T>::GetScalar(c2) > 255 ||
                ceres::JetOps<T>::GetScalar(c2) < 0)
            {

                residual = T(OUTLIER_RESIDUAL_INTENSITY_BN3);
                continue;

            }
#endif

#if 0//DEBUG
            // Print out
            printf("residual = %.4f (c1: %.4f - c2: %.4f) at u:%.4f, v_source:%.4f",
                   (double)residual,
                   (double)c1,
                   (double)c2,
                   (double)u,
                   (double)v_source);
            printf(" - v: %.4f, v_offset: %.4f, offset: %.4f",
                   (double)v,
                   (double)v_offset,
                   (double)offset);
            printf(" height(level):%d, width(level):%d",
                   (frame_->height(level_)),
                   (frame_->width(level_)));
            printf(" u_time: %.4f, i_spline: %.4f",
                   (double)u_time,
                   (double)i_spline);

            printf(" utime*(rowPerSeg_level - 1) = %.4f ",
                   u_time*(((int)rowPerSegment_ >> level_) - 1));

            printf(" (rowPerSeg_level = %d ",
                   ((int)rowPerSegment_ >> level_));
/*
            printf(" utime_diff: %.4f",
                   v_offset + 1 - u_time*((int)rowPerSegment_ >> level_));
*/
#if 0//SAVE
            double uaccForRow = (v_new_offset +
                                 (((int)startRowImage_) >> level_))
                                /((((int)rowPerSegment_) >> level_) - 1);
            printf(" v_min_uaccForRow*rowPerSeg_level = %f\n",
                   v_new_offset + ((int)startRowImage_ >> level_)
                   - uaccForRow*(((int)rowPerSegment_ >> level_) - 1));
#endif
            double uaccForRow = (v_new_offset + 1
                                 + (((int)startRowImage_) >> level_)
                                 - (((int)startRowSeg_) >>level_))
                                / (((int)rowPerSegment_) >> level_);

            printf(" utime_diff: %.4f",
                   u_time - uaccForRow);

            printf(" v_min_uaccForRow*rowPerSeg_level = %f\n",
                   v_new_offset + 1 + ((int)startRowImage_ >> level_)
                   - (((int)startRowSeg_) >>level_)
                   - uaccForRow*(((int)rowPerSegment_ >> level_)));


#endif
            //------------------------------------------------------------------
            // Determine if the pixel is good or badly warped
            T resOut = residual;
            T squareResOut = resOut*resOut;
            sumSquareResOut = sumSquareResOut + squareResOut;
            
            //this->debugInfo_weightVar = weightVar;
            //this->debugInfo_rp = rp;
            this->debugInfo_resOut = resOut;
            this->debugInfo_sumSquareResOut = sumSquareResOut;
            this->debugInfo_source_x = (new_refPoint[0]/new_refPoint[2])*fx_l + cx_l;
            this->debugInfo_source_y = (new_refPoint[1]/new_refPoint[2])*fy_l + cy_l;
            this->debugInfo_target_u = u_new;
            this->debugInfo_target_v = v_new;
            
            bool isGood = resOut*resOut / (MAX_DIFF_CONSTANT +
                          MAX_DIFF_GRAD_MULT*(resInterp[0]*resInterp[0] +
                          resInterp[1]*resInterp[1])) < 1.0;
            if (isGood == true)
            {

                WarpResidual_BN3::goodWarpPixelCountInPyramidImage_++;
                (*goodCount)++;

            }
            else
            {

                WarpResidual_BN3::badWarpPixelCountInPyramidImage_++;
                (*badCount)++;

            }

            //------------------------------------------------------------------
            // If depth becomes larger, pixels becomes smaller
            // hence it count less
            float depthChange = (float)((new_refPoint)[2] / Wxp[2]);
            *usageFractus += depthChange < 1 ? depthChange : 1;

            //------------------------------------------------------------------
            // Plot debug images
            // Display debug images
            if (lsd_slam::plotTracking || consoleMode)
            {

                // Force to plot only for idxSrcSegment == 1
//                if (idxSrcSegment == 1)
//                {

                    // Force resInterp same as the input
                    //resInterp[2] = (*refColourVariance)[0];
                plotPixelToDebugImages<T>(u_new,
                                          v_new_offset,
                                          residual,
                                          &new_refPoint,
                                          resInterp, isGood);
//                }

                //WORK BACK resInterp may not be correct for v_new + offset!!


            }

        } // end of for
        
        this->debugInfo_sumSquareResOut = sumSquareResOut;
        this->debugInfo_lastResidual = sumSquareResOut / (double) (*goodCount);

#if DEBUG
        printf("sum of squared res out = %f (Mean %f)\n", sumSquareResOut,
                sumSquareResOut / (double) (*goodCount + *badCount));
        printf("}\n");
#endif

        return sumSquareResOut / (double) (*goodCount + *badCount);

    } // end of function

    // Display debug images
    template <typename T>
    void plotPixelToDebugImages(T u_new, T v_new,
                                T residual,
                         const Eigen::Vector3f *refPoint,
                         const Eigen::Vector3f resInterp,
                         bool isGood)
    {

        //----------------------------------------------------------------------
        // Debug


        // for debug plot only: find x,y again.
        // horribly inefficient, but who cares at this point...
        int width = frame_->width(0);
        Eigen::Matrix3f KLvl = frame_->K(level_);
        Eigen::Vector3f point = KLvl * (*refPoint);
        int x = (int)(point[0] / point[2] + 0.5f);
        int y = (int)(point[1] / point[2] + 0.5f);
        int width_l = frame_->width(level_);

        setPixelInCvMat(
            debugImageOldImageSourceRollingShutter_,
                        getGrayCvPixel((float)resInterp[2]),
                (int)(T(u_new)+0.5),
                (int)(T(v_new)+0.5),
                (width/width_l));
        setPixelInCvMat(
            debugImageOldImageWarpedRollingShutter_,
                        getGrayCvPixel((float)resInterp[2]),
                x,y,(width/width_l));

        if (isGood == true)
        {

            // Good warp pixel in gray with residual brightness
            // Larger residual : brighter gray
            // Less residual : darker gray
            setPixelInCvMat(
                debugImageResidualsRollingShutter_,
                getGrayCvPixel(
                    (float)T(residual)+128),
                            x,y,(width/width_l));
        }
        else
        {

            // Bad warp pixel in red colour
            setPixelInCvMat(
                 debugImageResidualsRollingShutter_,
                            cv::Vec3b(0,0,255),x,y,(width/width_l));

        }

    } // end of function

    void writeDebugImage()
    {

        if (saveToFile_ == true)
        {

            char buf[255];
            sprintf(buf, "output_source_%06d.jpg", frame_->id());
            cv::imwrite(buf, *debugImageOldImageSourceRollingShutter_);
            
#if 1 // DEBUG
            sprintf(buf, "output_warp_%06d.jpg", frame_->id());
            cv::imwrite(buf, *debugImageOldImageWarpedRollingShutter_);

            sprintf(buf, "output_residual_%06d.jpg", frame_->id());
            cv::imwrite(buf, *debugImageResidualsRollingShutter_);
#endif

        }

    }


}; // end of class UserCallback
    


} // end of namespace

#endif // WARPRESIDUAL_BN3_H

