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
//  DepthMapForRadialRolling.h
//
//  Created by Jae-Hak Kim on 8/06/2016.
//
//

#ifndef __lsd_slam_core__DepthMapForRadialRolling__
#define __lsd_slam_core__DepthMapForRadialRolling__

#include <stdio.h>
#include "DepthEstimation/DepthMapForRollingShutter.h"
#include "util/settings.h"

namespace lsd_slam
{
    
class DepthMapForRadialRolling : public DepthMapForRollingShutter
{
public:
    
    float dist_param_; // Lens distortion parameter
    
    double debug_x_sol;
    double debug_y_sol;
    double debug_a;
    double debug_b;
    double debug_c;
    double debug_y_dist;
    double debug_x_init;
    double debug_y_init;
    
public:
    
    Eigen::Matrix<double,3,3> skewsym(Eigen::Matrix<double,3,1> v)
    {
        
        Eigen::Matrix<double,3,3> M;
        
        M << 0, -v(2), v(1),
             v(2), 0, -v(0),
             -v(1), v(0), 0;
        
        return M;
        
    }
    
    DepthMapForRadialRolling(int w, int h, const Eigen::Matrix3f& K,
                             float dist_param,
                             bool useGTmotion = false)
    : DepthMapForRollingShutter(w, h, K)
    {
        
        debugImageSource = cv::Mat(h, w, CV_8UC3);
        debugImageTarget = cv::Mat(h, w, CV_8UC3);
        useGTmotion_ = useGTmotion;
        dist_param_ = dist_param;
    }
    
protected:
    
    // Extract image points on generalised epipolar curve with degenerate
    // handing from radial-rolling stereo images
    bool extractPointsOnGeneralEpipolarLineBand(const int x,
                                                const int y,
                                                Frame* const ref,
                                                std::list<PointEpiCurve>
                                                *ptsEpiCurve,
                                                RunningStats* const stats,
                                                bool useGT = false);
    
    // Add points to Epipolar Curve when NO radial distortion exists
    std::list<PointEpiCurve>
    addPtsToEpiCurveForNoRadial(double src_x_undist,
                                double src_y_undist,
                                Sophus::SE3d Xi_ySource_to_rowTarget,
                                Eigen::Matrix<double,3,1> eline,
                                double epx,
                                double epy,
                                double aRow_dist,
                                RunningStats* const stats);
public:
    // Add points to Epipolar Curve when radial and rolling distortion exists
    std::list<PointEpiCurve>
    addPtsToEpiCurveForRadialRolling(double src_x_undist,
                                     double src_y_undist,
                                     Sophus::SE3d Xi_ySource_to_rowTarget,
                                     Eigen::Matrix<double,3,1> eline,
                                     double epx,
                                     double epy,
                                     double tgt_y_dist,
                                     RunningStats* const stats);
public:    
    void build_UniformSamplePoints(std::list<PointEpiCurve> *ptsEpiCurve,
											  Frame *sourceImage,
											  Frame *targetImage);

protected:
    // Find intersections for epipolar curve
    bool findInterEpipolar_Optimization(Eigen::Matrix<double,3,1> epline,
                                        double y_dist,
                                        double x_init,
                                        double y_init,
                                        double *pInter_x,
                                        double *pInter_y);
    
    // Check rescale factor
    bool checkRescale(double u,
                      double v,
                      double nepx,
                      double nepy,
                      Eigen::Matrix<double,3,3> R,
                      Eigen::Matrix<double,3,1> t,
                      double prior_depth,
                      double u2,
                      double v2);
    
    void propagateDepth(Frame* new_keyframe);
    
    void processPropagationForPixel(int x, int y,
                                    const Frame *old_keyframe,
                                    const Frame *new_keyframe,
                                    const DepthMapPixelHypothesis *source,
                                    const float* newKFMaxGrad,
                                    const bool *trackingWasGood,
                                    const float* activeKFImageData,
                                    const float* newKFImageData);
    
    void
    processPropagationForPixelRow(
                                  int x, int y,
                                  const Frame *old_keyframe,
                                  Frame *new_keyframe,
                                  const DepthMapPixelHypothesis *source,
                                  int yMin,
                                  int yMax,
                                  RunningStats* stats);
    
    void
    processPropagationForPixelRow_Distort(
                                          int x_undist, int y_undist,
                                          const Frame *old_keyframe,
                                          Frame *new_keyframe,
                                          const DepthMapPixelHypothesis *source,
                                          int yMin_dist,
                                          int yMax_dist,
                                          RunningStats* stats);

    public: void initializeFromGTDepth(Frame* new_frame);

    
}; // end of class
    
}; // end of namespace lsd_slam

#endif /* defined(__lsd_slam_core__DepthMapForRadialRolling__) */
