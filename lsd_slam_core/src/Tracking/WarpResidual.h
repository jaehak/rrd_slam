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
#ifndef WARPRESIDUAL_H
#define WARPRESIDUAL_H

#include "sophus/se3.hpp"
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "ceres/rotation.h"

#include "Tracking/SE3TrackerRollingShutter.h"
#include "Tracking/TrackingReference.h"
#include "DataStructures/Frame.h"
#include "util/Undistorter.h"
#include "util/globalFuncs.h"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

#include "Tracking/JetOps.h"

namespace lsd_slam
{

class TrackingReference;
class Frame;

// Automatic differentiation for warping residual computation
struct WarpResidual
{

    //--------------------------------------------------------------------------
    // WarpResidual CostFunctor for Ceres-solver
    //--------------------------------------------------------------------------
    // NOTE: It warps a point from "reference" image to "frame" image
    //
    //             level : [int] pyramid level such that
    //                     SE3TRACKING_MIN_LEVEL <= level &&
    //                     level <= SE3TRACKER_MAX_LEVEL - 1
    //          refPoint : [3x1 float] semi-dense point of (x, y, z) such that
    //                     x coords, y coords and inverse depth z in the
    //                     reference image
    // refColourVariance : [2x1 float] semi-dense point of (I, V) such that
    //                     image intensity and variance
    //         reference : [TrackingReference] Keyframe image (source image)
    //             frame : [Frame *] Current new image (target image)
    //  referenceToFrame : [SE3f] Key frame pose wrt new image pose
    //       undistorter : [Undistorter] Undistortion object
    //                     used for finding warpd image coordinates in original
    //                     input image
    //          startRow : [float] starting row in 0..1
    //            endRow : [float] ending row in 0..1
    WarpResidual (const int level,
                  const Eigen::Vector3f *refPoint,
                  const Eigen::Vector2f *refColourVariance,
                  TrackingReference *reference,
                  Frame *frame,
                  Undistorter *undistorter,
                  float startRow,
                  float endRow)
        : level_(level),
          refPoint_(refPoint),
          refColourVariance_(refColourVariance),
          reference_(reference),
          frame_(frame),
          startRow_(startRow),
          endRow_(endRow)
    {

        this->undistorter = undistorter;

    }

    //--------------------------------------------------------------------------
    // Operator computing the residual error for warping a pixel from the
    // reference image to the new image
    //
    //--------------------------------------------------------------------------
    // Notes below are no longer correct. They were not a good idea at all.
    //       Since Ceres-solver somehow calls this operator() more than once
    //       therefore any counters may be invalid.
    //
    // NOTE: In Ceres-solver, operator() should be "const" and const-correctness
    //       to all sub-functions, but we can access and change static or
    //       mutable variables - we take an advantage of this to report
    //       statistics and debug image.
    //
    // NOTE2: This operator() is being called for every residual block and
    //        for every iteration. Therefore, any counter storing
    //        goodPixel, badPixel and usageCount should be set ZERO
    //        inside the iteration callback routine provided by
    //        UserCallback class.
    template <typename T> bool
    operator() (const T* const pose, T* residual) const
    {

        // Define SE3Group based on template typename T
        typedef Sophus::SE3Group<T> SE3T;

        //----------------------------------------------------------------------
        // A 3D point warped by the pose vector

        // Convert input pose variable to SE3T type
        Eigen::Matrix<T,6,1> pose_6vec;
        for (int i=0; i<6; i++)
            pose_6vec[i] = pose[i];

        SE3T pose_SE3T = SE3T::exp(pose_6vec);


        // Warp by rotation and translation
        Eigen::Matrix<T,3,3> rotMat = pose_SE3T.rotationMatrix();
        Eigen::Matrix<T,3,1> transVec = pose_SE3T.translation();
        Eigen::Matrix<T,3,1> pointRef;
        pointRef[0] = T((*refPoint_)[0]);
        pointRef[1] = T((*refPoint_)[1]);
        pointRef[2] = T((*refPoint_)[2]);
        Eigen::Matrix<T,3,1> Wxp = rotMat*pointRef + transVec;

        //----------------------------------------------------------------------
        // Projection of the 3D point in the (target) image

        // Camera parameters in the given pyramid level
        T fx_l = T(frame_->K(level_)(0,0));
        T fy_l = T(frame_->K(level_)(1,1));
        T cx_l = T(frame_->K(level_)(0,2));
        T cy_l = T(frame_->K(level_)(1,2));

        // New image coordinates warped and projected in the target image
        T u_new = (Wxp[0]/Wxp[2])*fx_l + cx_l;
        T v_new = (Wxp[1]/Wxp[2])*fy_l + cy_l;

        //----------------------------------------------------------------------
        // Check if the warped image point is wihtin the image frame

        // Image width and height for the given pyramid level
        int width_l = frame_->width(level_);
        int height_l = frame_->height(level_);

        // Boundary check (NOTE: the boundary is defined as
        // [1:width_l-1, 1:height_l-1] excluding the one-pixel boundary
        // in the original level image size.
        double u = ceres::JetOps<T>::GetScalar(u_new);
        double v = ceres::JetOps<T>::GetScalar(v_new);
        if ((u <= 1) ||
            (v <= 1) ||
            (u >= width_l - 2) ||
            (v >= height_l - 2))
        {

            // Out of target image boundary
            // We return zero (the lowest intensity) in this case.
            //WarpResidual::numWarpedPixelOutOfImage++;
            residual[0] = T(255);

            return true;

        }

        if ((isnan(u) == true) || (isnan(v) == true))
        {

            residual[0] = T(255);
            return true;

        }

        //----------------------------------------------------------------------
        // Get intensity values in the source and target image

        // Gradient for the given pyramid level in the reference (source) image
        const Eigen::Vector4f* frame_gradients = frame_->gradients(level_);

        // Interpolated gradients and intensity of the warped point in
        // the new (target) image.
        // Returns:
        //    resInterp[0] : interpolated gradient x
        //    resInterp[1] : interpolated gradient y
        //    resInterp[2] : intensity/colour value
        Eigen::Vector3f resInterp =
                getInterpolatedElement43(frame_gradients,
                                         u, v, width_l);

        // Intensity value of the pixel in the source image with applying
        // affine estimation
//        float c1 = affineEstimation_a * (*refColourVariance_)[0]
//                + affineEstimation_a;
        float c1 = (*refColourVariance_)[0];

        // Intensity value of the pixel in the target image
        // Take a chain rule
        double sample[3];
        sample[0] = double(resInterp[2]); // Intensity
        sample[1] = double(resInterp[0]); // Gradient x
        sample[2] = double(resInterp[1]); // Gradient y
        T uv[2] = { u_new, v_new };
        // f = sample[0]
        // dfdx = sample + 1
        // z = uv
        // ==> dfdx * dx*d(uv) : Derivative of intensity with respect to
        //                       the warped coordinates u and v
        T c2 = ceres::Chain<double, 2, T>::Rule(sample[0], sample + 1, uv);

        //----------------------------------------------------------------------
        // Compute residual depending on the row of the warped pixel points

//        // Get distorted image point coordinates to check within a scanline
//        T u_new_tmp = u_new * pow(2, level_);
//        T v_new_tmp = v_new * pow(2, level_);
//        float dist_u_new;
//        float dist_v_new;
//        undistorter->distortPoint(u_new_tmp, v_new_tmp,
//                                  &dist_u_new, &dist_v_new);

//        // Check if dist_v_new is in the scanline determined by the time
//        // For testing, here a certain size of strip is chosen
//        float startRowInOriginalImage = undistorter->getInputHeight()*startRow_;
//        float endRowInOriginalImage = undistorter->getInputHeight()*endRow_;
//        if ((dist_v_new >= startRowInOriginalImage) &&
//                (dist_v_new <= endRowInOriginalImage))
//        {

            // Residual error of intensity between the source and target pixel
            residual[0] = T(c1) - c2;

//        }
//        else
//        {
//            // The warped pixel is not in the region of rows, so do not consider
//            return false;

//        }

#if 0
        //----------------------------------------------------------------------
        // Was this point tracked good?
        double resOut = ceres::JetOps<T>::GetScalar(residual[0]);
        bool isGood = resOut*resOut / (MAX_DIFF_CONSTANT +
                      MAX_DIFF_GRAD_MULT*(resInterp[0]*resInterp[0] +
                      resInterp[1]*resInterp[1])) < 1.0;
#endif

        //----------------------------------------------------------------------
        // Compute weight on the residual
#if 0
        T tx = transVec[0];
        T ty = transVec[1];
        T tz = transVec[2];
        T px = Wxp[0];
        T py = Wxp[1];
        T pz = Wxp[2];
        T d = T(1.0)/pointRef[2];
        T rp = residual[0];
        T gx = fx_l*T(resInterp[0]);
        T gy = fy_l*T(resInterp[1]);
        T s = T(1.0)*T((*refColourVariance_)[1]);

        T g0 = (tx*pz - tz*px)/(pz*pz*d);
        T g1 = (ty*pz - tz*py)/(pz*pz*d);

        T drpdd = gx*g0 + gy*g1;
        T w_p = T(1.0) / (T((lsd_slam::cameraPixelNoise2)) + s*drpdd*drpdd);

        // Multiply the weight to the residual
        T weighted_rp = ceres::abs(rp*ceres::sqrt(w_p));
        T wh = ceres::abs(weighted_rp < T(3.0/2) ?
                               T(1.0) : T(3.0/2) / weighted_rp);
        residual[0] = ceres::sqrt(wh * w_p) * rp;
#else
        // My computation of weight
        // Derivative of the residual with respect to inverse depth
        // f = residual[0].a
        // dfdx = residual[0].v
        // z = depth
        //T a = ceres::Chain<double, 1, T>::Rule(sample[0], sample + 1, Wxp[2]);

        //T s = T(1.0)*T((*refColourVariance_)[1]);
        //T w_p = T(1.0)/(T(lsd_slam::cameraPixelNoise2) + s*drpdd*drpdd);
#endif

#if 0
        printf("depthChange = %f, refPoint2 = %f, Wxp[2] =%f, usageCount_=%f\n",
               depthChange,
               (*refPoint_)[2],
                (float) ceres::JetOps<T>::GetScalar(Wxp[2]),
                usageCount_);
#endif

#if 0
        //----------------------------------------------------------------------
        // Debug
        if (lsd_slam::plotTrackingIterationInfo == true)
        {

            // for debug plot only: find x,y again.
            // horribly inefficient, but who cares at this point...
            int width = frame_->width(0);
            Eigen::Matrix3f KLvl = frame_->K(level_);
            Eigen::Vector3f point = KLvl * (*refPoint_);
            int x = point[0] / point[2] + 0.5f;
            int y = point[1] / point[2] + 0.5f;


            setPixelInCvMat(
                &lsd_slam::SE3TrackerRollingShutter::debugImageOldImageSourceRollingShutter,
                            getGrayCvPixel((float)resInterp[2]),
                    ceres::JetOps<T>::GetScalar(u_new)+0.5,
                    ceres::JetOps<T>::GetScalar(v_new)+0.5,
                    (width/width_l));
            setPixelInCvMat(
                &lsd_slam::SE3TrackerRollingShutter::debugImageOldImageWarpedRollingShutter,
                            getGrayCvPixel((float)resInterp[2]),
                    x,y,(width/width_l));

            if (isGood == true)
            {

                // Good warp pixel in gray with residual brightness
                // Larger residual : brighter gray
                // Less residual : darker gray
                setPixelInCvMat(
                    &lsd_slam::SE3TrackerRollingShutter::debugImageResidualsRollingShutter,
                    getGrayCvPixel(
                        (float)ceres::JetOps<T>::GetScalar(residual[0])+128),
                                x,y,(width/width_l));
            }
            else
            {

                // Bad warp pixel in red colour
                setPixelInCvMat(
                     &lsd_slam::SE3TrackerRollingShutter::debugImageResidualsRollingShutter,
                                cv::Vec3b(0,0,255),x,y,(width/width_l));

            }

        }
#endif

        return true;

    } // end of operator()

    void debugStart() const;
    void debugFinish() const;

public:

    // Result of warping
    static int goodWarpPixelCountInPyramidImage_;
    static int badWarpPixelCountInPyramidImage_;
    static float usageCount_;
    static int lastGoodCount_;
    static int lastBadCount_;
    static int numWarpedPixelOutOfImage;

    // Debug images
    static cv::Mat debugImageResiduals;
    static cv::Mat debugImageWeights;
    static cv::Mat debugImageOldImageSource;
    static cv::Mat debugImageOldImageWarped;

private:

    // Pyramid level
    const int level_;

    // X,Y,Z point in the reference (source) image
    const Eigen::Vector3f *refPoint_;

    // Intensity/Colour and Variance for inverse depth of the point
    // in the reference (source) image
    const Eigen::Vector2f *refColourVariance_;

    // tracking reference containing the reference (source)image,
    // gradients and inverse depth
    TrackingReference *reference_;

    // Current image (target)
    Frame *frame_;

    // Values between 0 and 1 indicating the start and end of rows in
    // the target image for warping
    // NOTE: designed for specifying the range of rows
    float startRow_;
    float endRow_;

    // Variables used in warping
    float affineEstimation_a = 1.0;
    float affineEstimation_b = 0.0;

    // Image undistortion/distortion object
    Undistorter *undistorter;

    // Warped image buffer
    float* buf_warped_residual;
    float* buf_warped_dx;
    float* buf_warped_dy;
    float* buf_warped_x;
    float* buf_warped_y;
    float* buf_warped_z;
    float* buf_d;
    float* buf_idepthVar;
    int buf_warped_size;

    // Result of warping
    float lastMeanRes;
    float affineEstimation_a_lastIt;
    float affineEstimation_b_lastIt;

};

// This UserCallback class runs at the end of iteration of Ceres-solver.
// This is for examining variables and status at every iteration end.
class UserCallback : public ceres::IterationCallback
{

public:

    const int level_;
    const TrackingReference *reference_;
    Frame *frame_;
    const double *pose_;
    cv::Mat *debugImageWeightRollingShutter_;
    cv::Mat *debugImageResidualsRollingShutter_;
    cv::Mat *debugImageSecondFrameRollingShutter_;
    cv::Mat *debugImageOldImageWarpedRollingShutter_;
    cv::Mat *debugImageOldImageSourceRollingShutter_;

public:

    UserCallback(const int pyramidLevel,
                 const TrackingReference *reference,
                 Frame *frame,
                 const double *pose,
                 cv::Mat *debugImageWeightRollingShutter,
                 cv::Mat *debugImageResidualsRollingShutter,
                 cv::Mat *debugImageSecondFrameRollingShutter,
                 cv::Mat *debugImageOldImageWarpedRollingShutter,
                 cv::Mat *debugImageOldImageSourceRollingShutter)
        : level_(pyramidLevel),
          reference_(reference),
          frame_(frame),
          pose_(pose),
          debugImageWeightRollingShutter_(debugImageWeightRollingShutter),
          debugImageResidualsRollingShutter_(debugImageResidualsRollingShutter),
          debugImageSecondFrameRollingShutter_(debugImageSecondFrameRollingShutter),
          debugImageOldImageWarpedRollingShutter_(debugImageOldImageWarpedRollingShutter),
          debugImageOldImageSourceRollingShutter_(debugImageOldImageSourceRollingShutter)
    {

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
            lsd_slam::WarpResidual::lastGoodCount_ = 0;
            lsd_slam::WarpResidual::lastBadCount_ = 0;
            lsd_slam::WarpResidual::usageCount_ = 0.0;

            return ceres::SOLVER_CONTINUE;
        }
        else if (summary.step_is_successful == true)
        {

            // Evaluate the warp error and returns the statistics
            // related to number of good and badly warped pixels, and
            // the number of depth changed pixels
            int goodCount = 0;
            int badCount = 0;
            float usageFractus = 0.0;
            getStatisticsForWarpedResidual<double>(&goodCount,
                                           &badCount,
                                           &usageFractus);
            lsd_slam::WarpResidual::lastGoodCount_ = goodCount;
            lsd_slam::WarpResidual::lastBadCount_ = badCount;
            lsd_slam::WarpResidual::usageCount_ = usageFractus;

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
    void getStatisticsForWarpedResidual(int *goodCount,
                                        int *badCount,
                                        float *usageFractus)
    {

        // Set to zeros
        *goodCount = 0;
        *badCount = 0;
        *usageFractus = 0;


        Eigen::Vector3f* refPoint = reference_->posData[level_];
        Eigen::Vector3f* refPoint_max = refPoint +
                reference_->numData[level_];
        Eigen::Vector2f* refColourVariance =
                reference_->colorAndVarData[level_];

        T residual;

        for (; refPoint<refPoint_max; refPoint++, refColourVariance++)
        {

            // Define SE3Group based on template typename T
            typedef Sophus::SE3Group<T> SE3T;

            //------------------------------------------------------------------
            // A 3D point warped by the pose vector

            // Convert input pose variable to SE3T type
            Eigen::Matrix<T,6,1> pose_6vec;
            for (int i=0; i<6; i++)
                pose_6vec[i] = pose_[i];

            SE3T pose_SE3T = SE3T::exp(pose_6vec);

            // Warp by rotation and translation
            Eigen::Matrix<T,3,3> rotMat = pose_SE3T.rotationMatrix();
            Eigen::Matrix<T,3,1> transVec = pose_SE3T.translation();
            Eigen::Matrix<T,3,1> pointRef;
            pointRef[0] = T((*refPoint)[0]);
            pointRef[1] = T((*refPoint)[1]);
            pointRef[2] = T((*refPoint)[2]);
            Eigen::Matrix<T,3,1> Wxp = rotMat*pointRef + transVec;

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
            if ((u <= 1) ||
                (v <= 1) ||
                (u >= width_l - 2) ||
                (v >= height_l - 2))
            {

                // Out of target image boundary
                // We return zero (the lowest intensity) in this case.
                WarpResidual::numWarpedPixelOutOfImage++;
                residual = T(0);

                continue;

            }

            //------------------------------------------------------------------
            // Get intensity values in the source and target image

            // Gradient for the given pyramid level in the reference image
            const Eigen::Vector4f* frame_gradients = frame_->gradients(level_);

            // Interpolated gradients and intensity of the warped point in
            // the new (target) image.
            // Returns:
            //    resInterp[0] : interpolated gradient x
            //    resInterp[1] : interpolated gradient y
            //    resInterp[2] : intensity/colour value
            Eigen::Vector3f resInterp =
                    getInterpolatedElement43(frame_gradients,
                                             u, v, width_l);

            //------------------------------------------------------------------
            // Residual between intensities from source and target
            T c1 = T((*refColourVariance)[0]);
            T c2 = T(resInterp[2]);
            residual = c1 - c2;

            //------------------------------------------------------------------
            // Determine if the pixel is good or badly warped
            T resOut = residual;
            bool isGood = resOut*resOut / (MAX_DIFF_CONSTANT +
                          MAX_DIFF_GRAD_MULT*(resInterp[0]*resInterp[0] +
                          resInterp[1]*resInterp[1])) < 1.0;
            if (isGood == true)
            {

                WarpResidual::goodWarpPixelCountInPyramidImage_++;
                (*goodCount)++;

            }
            else
            {

                WarpResidual::badWarpPixelCountInPyramidImage_++;
                (*badCount)++;

            }

            //------------------------------------------------------------------
            // If depth becomes larger, pixels becomes smaller
            // hence it count less
            float depthChange = (*refPoint)[2] / Wxp[2];
            *usageFractus += depthChange < 1 ? depthChange : 1;

            //------------------------------------------------------------------
            // Plot debug images
            // Display debug images
            if (lsd_slam::plotTrackingIterationInfo == true)
            {

                plotPixelToDebugImages<double>(u_new, v_new, residual,
                                               refPoint, resInterp, isGood);

            }

        } // end of for

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
        int x = point[0] / point[2] + 0.5f;
        int y = point[1] / point[2] + 0.5f;
        int width_l = frame_->width(level_);


        setPixelInCvMat(
            debugImageOldImageSourceRollingShutter_,
                        getGrayCvPixel((float)resInterp[2]),
                T(u_new)+0.5,
                T(v_new)+0.5,
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

}; // end of class UserCallback


} // end of namespace

#endif // WARPRESIDUAL_H
