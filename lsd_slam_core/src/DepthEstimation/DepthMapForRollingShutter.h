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
//  DepthMapForRollingShutter.h
//
//  Created by Jae-Hak Kim on 14/05/2015.
//
//

#ifndef __DepthMapForRollingShutter__
#define __DepthMapForRollingShutter__

#include <stdio.h>
#include "DepthEstimation/DepthMap.h"

namespace lsd_slam
{
    
struct StereoDebugInfo
{
    
    float alpha; // Result alpha
    float sampleDist; // Result sample distance
    float geoDispError; // Result geometric disparity error
    float photoDispError; // Result photo disparity error
    float rescaleFactor; // rescale factor return
    float gradAlongLine; // Return of gradient along line
    float trackingErrorFac; // Result of tracking error factor
    float geoDispError_denom; // geoDispError denominator
    float gradsInterpX; // Gradient interpolation x
    float gradsInterpY; // Gradient interpolation y
    float initialTrackedResidual;
    
};
    
struct PointEpiCurve
{
    
    double xSource; // x-coord of the point in the source
    double ySource; // y-coord of the point in the source
    double xTarget; // x-coord of the point on epipolar curve in target
    double yTarget; // y-coord of the point on epipolar curve in target
    double nepx; // x-coord of the normal on the point
    double nepy; // y-coord of the normal on the point
    Sophus::SE3d Xi; // motion, row y source to row target
    Eigen::Matrix<float,5,1> vals; // Intensities on (x,y) in the target
    bool isValsAvailable = false; // Check the intensity values available
    float rescaleFactor; // Rescale factor
    float invDepth; // Inverse depth of the point in the source
    float alpha; // Alpha for variance
    Eigen::Vector3f pClose; // Close point
    Eigen::Vector3f pFar; // Far point
    float incx; // Inc x
    float incy; // Inc y
    float eplLength; // Epipolar line length
    float incx_epl; // Inc x used for epipolar lengh computation
    float incy_epl; // Inc y used for epipolar lengh computation
    bool isMotionAvailable = false;
    bool isIncxIncyAvailable = false;
    
    PointEpiCurve()
    {
        // Init variables
        xSource = 0;
        ySource = 0;
        xTarget = 0;
        yTarget = 0;
        nepx = 0;
        nepy = 0;
        rescaleFactor = 0;
        invDepth = 0;
        alpha = 0;
        incx = 0;
        incy = 0;
        eplLength = 0;
        incx_epl = 0;
        incy_epl = 0;
        isMotionAvailable = false;
        isIncxIncyAvailable = false;
    
    };
    
    PointEpiCurve(double xSource, double ySource,
                  double xTarget, double yTarget,
                  double nepx, double nepy,
                  Sophus::SE3d motion)
    {
        
        // Init from parameters
        this->xSource = xSource;
        this->ySource = ySource;
        this->xTarget = xTarget;
        this->yTarget = yTarget;
        this->nepx = nepx;
        this->nepy = nepy;
        this->Xi = motion;
        this->isMotionAvailable = true;
        
        // Init the rest
        this->rescaleFactor = 0;
        this->invDepth = 0;
        this->alpha = 0;
        this->incx = 0;
        this->incy = 0;
        this->isIncxIncyAvailable = true;
        this->eplLength = 0;
        this->incx_epl = 0;
        this->incy_epl = 0;
        
    }
    
    void setSE3(Sophus::SE3d Xi)
    {
        
        this->Xi = Xi;
        this->isMotionAvailable = true;
        
    }
    
    void setVals(Eigen::Matrix<float,5,1> vals)
    {
        
        this->vals = vals;
        isValsAvailable = true;
        
    }
    
    void setRescaleFactor(float rescaleFactor)
    {
        
        this->rescaleFactor = rescaleFactor;
        
    }
    
    void setPointCloseFar(Eigen::Vector3f pClose, Eigen::Vector3f pFar)
    {
     
        this->pClose = pClose;
        this->pFar = pFar;
        
    }
    
    void setIncXandY(float incx, float incy)
    {
        
        this->incx = incx;
        this->incy = incy;
        this->isIncxIncyAvailable = true;
        
    }
    
    void setEplLength(float eplLength)
    {
        
        this->eplLength = eplLength;
        
    }
    
    void setIncXandYForEpl(float incx_epl, float incy_epl)
    {
        
        this->incx_epl = incx_epl;
        this->incy_epl = incy_epl;
        
    }
    
    // Compute inverse depth and set its member variable
    void computeInvDepth(float fxi, float fyi, float cxi, float cyi)
    {
        
        if (this != NULL)
        {
            // Check
            if ((this->incx == 0) && (this->incy == 0))
            {
                
                fprintf(stderr,
                        "Error computeInvDepth: "
                        "incx and incy are not available.\n");
                exit(1);
                return;
                
            }
            
            // Check
            if (this->isMotionAvailable == false)
            {
                
                fprintf(stderr,
                        "Error computeInvDepth: "
                        "Xi is not available.\n");
                exit(1);
                return;
                
            }
            
            // Check
            if (this->isIncxIncyAvailable == false)
            {

                fprintf(stderr,
                        "Error computeInvDepth: "
                        "incX and incY are not available.\n");
                exit(1);
                return;
                
            }
            
#if 0 // Incx and incy from source and target point coordinates
            
            float incx = xSource - xTarget;
            float incy = ySource - yTarget;

#else // Incx and Incy from stored PointEpiCurve structure

            float incxST = (float)(xSource - xTarget);
            float incyST = (float)(ySource - yTarget);
            float incx = this->incx;
            float incy = this->incy;

#endif
            float iDepth;
            float alphaValue;
            
            if (incxST*incxST > incyST*incyST)
            {

                // Other is Source; This is Target.
                Sophus::SE3d otherToThis = this->Xi;
                Sophus::SE3d sourceToTarget = this->Xi;
                Eigen::Vector3f KinvP = Eigen::Vector3f(fxi*(float)xSource + cxi,
                                                        fyi*(float)ySource + cyi,
                                                        1.0f);
                Eigen::Matrix<float,3,1> t =
                    otherToThis.translation().cast<float>();
                Eigen::Matrix<float,3,3> R =
                    otherToThis.rotationMatrix().cast<float>();
                float oldX = fxi*(float)xTarget + cxi;
                float nominator = (oldX*t[2] - t[0]);
                float dot0 = KinvP.dot(R.row(0));
                float dot2 = KinvP.dot(R.row(2));
                iDepth = (dot0 - oldX*dot2) / nominator;
                alphaValue = incx*fxi*(dot0*t[2] - dot2*t[0])
                    / (nominator*nominator);
                
            }
            else
            {
                
                Sophus::SE3d otherToThis = this->Xi;
                Eigen::Vector3f KinvP = Eigen::Vector3f(fxi*(float)xSource + cxi,
                                                        fyi*(float)ySource + cyi,
                                                        1.0f);
                Eigen::Matrix<float,3,1> t =
                otherToThis.translation().cast<float>();
                Eigen::Matrix<float,3,3> R =
                otherToThis.rotationMatrix().cast<float>();
                float oldY = fyi*(float)yTarget + cyi;
                float nominator = (oldY*t[2] - t[1]);
                float dot1 = KinvP.dot(R.row(1));
                float dot2 = KinvP.dot(R.row(2));
                iDepth = (dot1 - oldY*dot2) / nominator;
                alphaValue = incy*fyi*(dot1*t[2] - dot2*t[1])
                    / (nominator*nominator);
                
            }
            
//            // Make inverse depth positive always
//            if (iDepth < 0)
//            {
//                iDepth = -iDepth;
//            }
            
            // Set invDepth member variable
            this->invDepth = iDepth;
            
            // Set alpha
            this->alpha = alphaValue;
            
        }
        else
        {
            
            fprintf(stderr, "Error: PointEpiCurve object is not available\n");
            
        }
        
    }
    
    // The less-than operator for sorting a vector of PointEpiCurve structures
    // by their inverse depth in the ascending order.
    bool operator < (const PointEpiCurve& pt) const
    {
        
        return (invDepth < pt.invDepth);
    
    }
    
    
    
};

class DepthMapForRollingShutter : public DepthMap
{
    
public:

    DepthMapForRollingShutter(int w, int h, const Eigen::Matrix3f& K,
                              bool useGTmotion = false)
    : DepthMap(w, h, K)
    {
    
        debugImageSource = cv::Mat(h, w, CV_8UC3);
        debugImageTarget = cv::Mat(h, w, CV_8UC3);
        useGTmotion_ = useGTmotion;
    }
    
    DepthMapForRollingShutter(const DepthMapForRollingShutter&) = delete;
    
    DepthMapForRollingShutter& operator=(const DepthMapForRollingShutter&) = delete;
    
    ~DepthMapForRollingShutter()
    {
        
        debugImageSource.release();
        debugImageTarget.release();
        
    }
    
    /**
     * does obervation and regularization only.
     **/
    void updateKeyframe(std::deque< std::shared_ptr<Frame> > referenceFrames);
    
    /**
     * does propagation and whole-filling-regularization 
     * (no observation, for that need to call updateKeyframe()!)
     **/
    void createKeyFrame(Frame* new_keyframe);
    
    // Save the depth as a text file
    void saveDepthAsFile(char *filename);
    
    // Debug images
    cv::Mat debugImageSource;
    cv::Mat debugImageTarget;
    float debugInfo_rescaleFactor;
    float debugInfo_epxn;
    float debugInfo_epyn;
    float debugInfo_alpha;
    float debugInfo_sampleDist;
    float debugInfo_geoDispError;
    float debugInfo_photoDispError;
    float debugInfo_result_var;
    float debugInfo_trackingErrorFac;
    float debugInfo_gradAlongLine;
    float debugInfo_gradsInterp_0;
    float debugInfo_gradsInterp_1;
    float debugInfo_referenceFrame_initialTrackedResidual;
    float debugInfo_geoDispError_first;
    float debugInfo_xSource;
    float debugInfo_ySource;
    float debugInfo_xTarget;
    float debugInfo_yTarget;
    float debugInfo_best_match_err;
        
    Eigen::Vector3d debugInfo_old2DPoint;
    Eigen::Vector3d debugInfo_new2DPoint;
    Sophus::SE3d debugInfo_worldToNew;
    Sophus::SE3d debugInfo_worldToOld;
    Sophus::SE3d debugInfo_oldToNew;
    int debugInfo_rowNew;

    // Flag using GTmotion or not
    bool useGTmotion_;
    
    // Initialize map from ground truth with noise
    void initializeFromGTDepthWithNoise(Frame* new_frame);

protected:
    
    // ============ internal functions ==================================================
    // does the line-stereo seeking.
    // takes a lot of parameters, because they all have been pre-computed before.
    inline float
    doGeneralEpipolarLineStereo(
        const std::list<PointEpiCurve> *ptsEpiCurve, // Pts on Epipolar curve
        const float min_idepth, // Minimum inverse depth
        const float prior_idepth, // Prior inverse depth
        float max_idepth, // Maximum inverse depth
        lsd_slam::Frame *const referenceFrame, // refFrame (target)
        const float *referenceFrameImage, // tracking image
        float &result_idepth, // Result inverse depth
        float &result_var, // Result variance
        float &result_eplLength, // Result epipolar length
        lsd_slam::StereoDebugInfo &stereoDebugInfo,
        PointEpiCurve &best_ptEC, // Best point on epipolar curve
        lsd_slam::RunningStats *const stats // Return status
    );
    
    inline float doLineStereo(
                              const float u, const float v, const float epxn, const float epyn,
                              const float min_idepth, const float prior_idepth, float max_idepth,
                              Frame* const referenceFrame, const float* referenceFrameImage,
                              float &result_idepth, float &result_var, float &result_eplLength,
                              RunningStats* const stats);
    inline float doLineStereo_ver2(
                              const float u, const float v, const float epxn, const float epyn,
                              const float min_idepth, const float prior_idepth, float max_idepth,
                              Frame* const referenceFrame, const float* referenceFrameImage,
                              Sophus::Sim3d Xi_worldToKF,
                              float &result_idepth, float &result_var, float &result_eplLength,
                              RunningStats* const stats);
    
    virtual void propagateDepth(Frame* new_keyframe);
    virtual void processPropagationForPixel(int x, int y,
                                    const Frame *old_keyframe,
                                    const Frame *new_keyframe,
                                    const DepthMapPixelHypothesis *source,
                                    const float* newKFMaxGrad,
                                    const bool *trackingWasGood,
                                    const float* activeKFImageData,
                                    const float* newKFImageData);
    
    bool observeDepthCreate(const int &x, const int &y, const int &idx,
                            RunningStats* const &stats,
                            bool debugPlot);
    bool observeDepthUpdate(const int &x, const int &y, const int &idx, const float* keyFrameMaxGradBuf, RunningStats* const &stats, bool debugPlot);
    bool makeAndCheckEPL(const int x, const int y, const Frame* const ref,
                         float* pepx, float* pepy, RunningStats* const stats);
    bool makeAndCheckEPL_ver2(const int x, const int y, const Frame* const ref,
                         Sophus::Sim3d Xi_worldToKF,
                         float* pepx, float* pepy, RunningStats* const stats);
    void observeDepthRow(int yMin, int yMax, RunningStats* stats);
    void observeDepth();
    
    void regularizeDepthMapFillHoles();
    void buildRegIntegralBuffer();
    void regularizeDepthMapFillHolesRow(int yMin, int yMax, RunningStats* stats);
    void buildRegIntegralBufferRow1(int yMin, int yMax, RunningStats* stats);

    
    void regularizeDepthMap(bool removeOcclusion, int validityTH);
    template<bool removeOcclusions> void regularizeDepthMapRow(int validityTH, int yMin, int yMax, RunningStats* stats);
    
    
    void setDepthFromGT();
    void setDepthFromGTRow(int yMin, int yMax, RunningStats* stats);
    std::shared_ptr<float> readDepthFile(std::string filename);
    
    bool extractPointsOnGeneralEpipolarLine(const int x,
                                const int y,
                                Frame* const ref,
                                std::vector<PointEpiCurve> *ptsEpiCurve,
                                RunningStats* const stats,
                                bool useGT = false);

    // Extract image points on generalised epipolar curve with degenerate
    // handling
    virtual bool extractPointsOnGeneralEpipolarLineBand(const int x,
                                            const int y,
                                            Frame* const ref,
                                            std::list<PointEpiCurve>
                                                *ptsEpiCurve,
                                            RunningStats* const stats,
                                            bool useGT = false);
    
    int copyImageBandsToBuffer(const float min_idepth,
                               const float prior_idepth,
                               const float max_idepth,
                                const float *targetImage,
                                std::list<PointEpiCurve> *ptsEpiCurve,
                                RunningStats* const stats);
    
    bool makeAndCheckEPL_forGeneralEpipolarLine(
           const int x, // x-coord of semi-dense points in the source/key image
           const int y, // y-coord of semi-dense points in the source/key image
           Eigen::Matrix<double,3,1> thisToOther_t, // translation this to other
           double* pepx,
           double* pepy,
           RunningStats* const stats);
    
    Eigen::Matrix<float,5,1>
    getSourceValsFromEpiPtCurve(const PointEpiCurve* ptEC);
    
    void sortPointsOnGeneralEpipolarLineByIDepth(std::vector<PointEpiCurve>
                                                *ptsEpiCurve);
    
    cv::Mat plotGeneralStereoImages(Frame* const referenceFrame,
                                    Frame* const targetFrame,
                                    const std::list<PointEpiCurve>
                                        *ptsEpiCurve,
                                    PointEpiCurve *best_ptEC);

    // Plot general stereo images with epipolar lines of estimate
    // and the ground truth.
    // Estimate: ptsEpiCurve
    // Ground Truth estimate: GT_ptsEpiCurve which can be obtained
    // from extractPointsOnGeneralEpipolarLine() with useGT = true
    cv::Mat plotGeneralStereoImagesWithGT(Frame* const referenceFrame,
                                    Frame* const targetFrame,
                                    const std::list<PointEpiCurve>
                                    *ptsEpiCurve,
                                    const std::list<PointEpiCurve>
                                    *GT_ptsEpiCurve,
                                    PointEpiCurve *best_ptEC);
    
public:
    
    // Push the point into Epipolar Curve Structure with depth and
    // normal information
    void pushPointToEpiCurve(double xs, double ys,
                             double xt, double yt,
                             double nepx, double nepy,
                             Sophus::SE3d motion,
                             std::list<PointEpiCurve> *ptsEpiCurve);

    // Convert a vector into homogeneous coordinates
    template <typename Derived>
    static
    Eigen::MatrixBase<Derived>& hcoord(const Eigen::MatrixBase<Derived>& x)
    {
        
        int m = x.rows();
        int n = x.cols();
        
        Eigen::MatrixBase<Derived>& y =
        const_cast< Eigen::MatrixBase<Derived>& >(x);
        for(int j=0;j<n;j++)
        {
            
            if (x(m - 1, j) <= (1.0e-06) &&
                x(m - 1, j) >= (-1.0e-06))
            {
                // Skip for very small last coordinate
                continue;
                
            }
            else
            {
                
                for (int i=0;i<m;i++)
                {
                    
                    // Divide by the last coordinate
                    y(i,j) = x(i,j)/x(m - 1, j);
                    
                }
                
            }
            
        }
        
        return y;
    }
    
    std::list<PointEpiCurve>
    buildCurveForPointAtInfinity(double x, double y,
                                 double nepx, double nepy,
                                 Sophus::SE3d motion,
                                 int rowIndex,
                                 Eigen::Matrix<double,3,1> upper_pInter,
                                 Eigen::Matrix<double,3,1> lower_pInter)
    {
        
        std::list<PointEpiCurve> ptsEpiCurve;
    
        // Points
        double x1 = upper_pInter[0];
        double z1 = upper_pInter[2];
        double x2 = lower_pInter[0];
        double z2 = lower_pInter[2];
        double x1z1 = x1/z1;
        double x2z2 = x2/z2;
    
        // Check if there are within the band
        if (x1z1 > x2z2)
        {
            
            if ((x2z2 > width) || (x1z1 < 0))
            {
                
                // Discard
                // continue;
                
            }
            else
            {
                // When 0 < x2z2 < x1z1 < width
                
                // Choose pixels max{0, x2z2} up to min{x1z1, width}
                double max_0_x2z2 = x2z2 > 0 ? x2z2 : 0;
                double min_x1z1_width = x1z1 < width ? x1z1 : width;
                
                // Push all pixels in the row
                for (int i=(int)round(max_0_x2z2);i<=min_x1z1_width;i++)
                {
                    
                    pushPointToEpiCurve(x, y,
                                        i, rowIndex,
                                        nepx, nepy,
                                        motion,
                                        &ptsEpiCurve);
                    
                }
                
            }
            
        }
        else
        {
                
            // When x2z2 is on the right to x1z1
            
            
            if ((x1z1 > width) || (x2z2 < 0))
            {
                
                // Discard
                // continue;
                
            }
            else
            {
                // When 0 < x1z1 < x2z2 < width
                
                // Choose pixels max{0, x1z1} up to min{x2z2, width}
                double max_0_x1z1 = x1z1 > 0 ? x1z1 : 0;
                double min_x2z2_width = x2z2 < width ? x2z2 : width;
                
                // Push all pixels in the row
                for (int i=(int)round(max_0_x1z1);i<=min_x2z2_width;i++)
                {
                    
                    pushPointToEpiCurve(x, y,
                                        i, rowIndex,
                                        nepx, nepy,
                                        motion,
                                        &ptsEpiCurve);
                    
                }
                
            }
            
        } // end of if
        
        return ptsEpiCurve;
        
    }
    
}; // end of class DepthMapForRollingShutter
    
} // end of namespace lsd_slam

#endif /* defined(__DepthMapForRollingShutter__) */
