/**
* This file is part of LSD-SLAM.
*
* Copyright 2013 Jakob Engel <engelj at in dot tum dot de> (Technical University of Munich)
* For more information see <http://vision.in.tum.de/lsdslam> 
*
* LSD-SLAM is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* LSD-SLAM is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with LSD-SLAM. If not, see <http://www.gnu.org/licenses/>.
*/

// Modified by Jae-Hak Kim <jaehak.kim@adelaide.edu.au.>
// Copyrignt@2015, University of Adelaide
//------------------------------------------------------------------------------
// HISTORY
//------------------------------------------------------------------------------
// 2015-05-06 New member variables for rolling shutter frame
//------------------------------------------------------------------------------

#pragma once
#include "util/SophusUtil.h"
#include "util/settings.h"
#include <boost/thread/recursive_mutex.hpp>
#include <boost/thread/shared_mutex.hpp>
#include "DataStructures/FramePoseStruct.h"
#include "DataStructures/FrameMemory.h"
#include "unordered_set"
#include "util/settings.h"

#include "Tracking/Spline.h"
#include <stdio.h>
#include "util/MyMats.h"

#include "util/Undistorter.h"
#define EIGEN_USE_NEW_STDVECTOR
#include<Eigen/StdVector>

namespace lsd_slam
{


class DepthMapPixelHypothesis;
class TrackingReference;
/**
 */

class Frame
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	friend class FrameMemory;


    //==========================================================================
    // Change for rolling shutter by Jae-Hak Kim
    //--------------------------------------------------------------------------
    // Members
    //--------------------------------------------------------------------------
    int startRowImage = 0; // Starting row for image (Rolling shutter)
    int startRowSeg = 0; // Starting row for segment (Rolling shutter)
    int rowPerSegment = 0; // Number of row per segment (for rolling shutter)
    bool isSplineValid = false;
    std::string depthFilename;
    std::string imageFilename;
    std::string motionFilename;
    std::vector<Eigen::Matrix<double,4,4>,
    Eigen::aligned_allocator<Eigen::Matrix<double,4,4>> > RT_motion; // World to Cam
    bool is_RT_motion_Available = false;
    
    // Spline pointer
    std::shared_ptr<lsd_slam::Spline<double> > spline_;
    
    // Undistorter
    lsd_slam::Undistorter *undistorter_;
    
    
    // Number of extra control points to make B-spline pass through
    // a first starting point of the curve, which is determined by
    // the order of B-spline k.
    int numExtraCtrlPts_ = 0;
    
    // Set GT motion
    void setGTMotion();
    
    //--------------------------------------------------------------------------
    // Methods
    //--------------------------------------------------------------------------
    // A time parameter associated with the row of a pixel in the rolling
    // shutter image
    double uTimeAtRow(double y, int level = 0) const
    {

        return ((double)y + 0
                + ((int) startRowImage >> level)
                - ((int) startRowSeg >> level))
        / ((int) rowPerSegment >> level);
        
    }
    
    //--------------------------------------------------------------------------
    // A time parameter associate with the row of a pixel in the radially
    // disotrted rolling-shutter image. It uses an undistorter object
    // to distort undistorted pixel locations to be used in finding the time
    // parameter. It requires both x and y coordinates.
    // INPUT:   x, y : [float, float] undistorted pixels
    // OUTPUT:   out : [float] u-time
    double uTimeAtRow_Distort(double x, double y, int level = 0) const
    {
        // Update y according to the undistorter object
        float y_new = (float) y;
        if ((undistorter_->isValid() == true) &&
            (undistorter_->isNoRectification() == false))
        {
            
            // Update via undistortion
            float x_dist, y_dist;
            undistorter_->distortPointLevel((float)x, (float)y,
                                            &x_dist, &y_dist,
                                            level);
            //printf("DISTORT_POINT_LEVEL: x, y = %f %f\n", x, y);
            //printf("DISTORT_POINT_LEVEL: x_dist, y_dist = %f %f\n", x_dist, y_dist);
            //printf("DISTORT_POINT_LEVEL: level = %d\n", level);
            y_new = y_dist;
            
        }
        
        double u = ((double)y_new + 0
                + ((int) startRowImage >> level)
                - ((int) startRowSeg >> level))
                / ((int) rowPerSegment >> level);
//        printf("DISTORT_POINT_LEVEL: x, y, y_new, startRowImg, startRowSeg "
//               "rowPerSeg level u = %f %f %f %d %d %d %d %f\n",
//               x, y, y_new, startRowImage, startRowSeg, rowPerSegment, level, u);
        
        return u;
        
    }
    
    //--------------------------------------------------------------------------
    // An index of spline control points associated with a rolling shutter image
    int indexSplineCtrlPoint() const
    {
        
        int index = 1;
        if (isSplineValid == false)
        {
            
            // Error: no spline available
            ///fprintf(stderr, "ERROR: No spline available in Frame\n");
            //fprintf(stderr, "       imageFilename = %s\n",
            //        imageFilename.c_str());
            //fprintf(stderr, "       startRowImage = %d\n", startRowImage);
            //fprintf(stderr, "       startRowSeg = %d\n", startRowSeg);
            //fprintf(stderr, "       rowPerSegment = %d\n", rowPerSegment);
            //exit(1);
            
        }
        else
        {
        
            index = (int)floor(startRowImage/rowPerSegment) + numExtraCtrlPts_;
    
        }
        return index;
    }
    
    //--------------------------------------------------------------------------
    // Prepare stereo with pixel point for rolling shutter image
    void prepareForRollingShutterStereo(const Frame* other,
                                        Sim3 thisToOther,
                                        const int y,
                                        const Eigen::Matrix3f& K,
                                        const int level);
    
    //--------------------------------------------------------------------------
    // Set a B-Spline pointer
    void setSpline(const std::shared_ptr<lsd_slam::Spline<double> > spline)
    {
        spline_ = spline;
        numExtraCtrlPts_ = spline_->order() - 1;
        isSplineValid = true;
    }
    
    //--------------------------------------------------------------------------
    // Set an undistorter object
    void setUndistorter(lsd_slam::Undistorter *undist)
    {
        undistorter_ = undist;
    }
    
    //--------------------------------------------------------------------------
    // Set frame id
    void setID(int id)
    {
        data.id = id;
    }
    
    //--------------------------------------------------------------------------
    // Set timestamp
    void setTimestamp(double timestamp)
    {
        data.timestamp = timestamp;
    }
    
    //--------------------------------------------------------------------------
    // Set image mask
    void setImageMask(const unsigned char* imageMask)
    {

        // Allocate memory of image mask of the data, if empty
        if (data.imageMaskValid[0] == false)
        {
        
            data.imageMask[0] = FrameMemory::getInstance().getFloatBuffer(
                                                data.width[0]*data.height[0]);
            
        }
        float *maxPt = data.imageMask[0] + data.width[0]*data.height[0];
        const unsigned char* ptr_imageMask = imageMask;
        
        // Copy image mask to data
        for(float* pt = data.imageMask[0]; pt < maxPt; pt++)
        {
            
            if (ptr_imageMask == 0)
            {
                
                *pt = 255; // Default 255, when no imageMask is given
                
            }
            else
            {
                
                *pt = *imageMask; // Copy imageMask to data array
                imageMask++;
                
            }
            
        }
        
        data.imageMaskValid[0] = true;
    
    }
    
    //--------------------------------------------------------------------------
    // Set as if it is the first frame but do not init depth
    void initializeExceptDepth(int id, int width, int height,
                               const Eigen::Matrix3f& K, double timestamp,
                               int height_rollingShutter_inputImage)
    {
        
        this->height_rollingShutter_inputImage =
            height_rollingShutter_inputImage;
        data.id = id;
        
        // Init pose
        pose = new FramePoseStruct(this);
        
        data.K[0] = K;
        data.fx[0] = K(0,0);
        data.fy[0] = K(1,1);
        data.cx[0] = K(0,2);
        data.cy[0] = K(1,2);
        
        data.KInv[0] = K.inverse();
        data.fxInv[0] = data.KInv[0](0,0);
        data.fyInv[0] = data.KInv[0](1,1);
        data.cxInv[0] = data.KInv[0](0,2);
        data.cyInv[0] = data.KInv[0](1,2);

        data.timestamp = timestamp;
        
        // Depths are available.
        data.hasIDepthBeenSet = true;
        depthHasBeenUpdatedFlag = true;
        
        referenceID = -1;
        referenceLevel = -1;
        
        numMappablePixels = -1;
        
        // Do not touch depths (they are already available)
//        data.idepthValid[0] = true;
//        data.idepthVarValid[0] = true;
//        release(IDEPTH | IDEPTH_VAR, true, true);
//        data.hasIDepthBeenSet = true;
//        depthHasBeenUpdatedFlag = true;
        
        data.validity_reAct = 0;
        data.idepthVar_reAct = 0;
        data.idepth_reAct = 0;
        
        data.refPixelWasGood = 0;
        
        permaRefNumPts = 0;
        permaRef_colorAndVarData = 0;
        permaRef_posData = 0;
        
        meanIdepth = 1;
        numPoints = 0;
        
        numFramesTrackedOnThis = numMappedOnThis = numMappedOnThisTotal = 0;
        
        idxInKeyframes = -1;
        
        edgeErrorSum = edgesNum = 1;
        
        lastConstraintTrackedCamToWorld = Sim3();
        
        isActive = false;

    }
    
    //--------------------------------------------------------------------------
    // Get a pose of the row
    // returns
    //       pose : [4x4 double] world to cam
    Eigen::Matrix<double,4,4> getPoseAtRow(double y, int level = 0) const
    {
        
        // Lock for accessing
        //boost::shared_lock<boost::shared_mutex> lock = getActiveLock();
        
        double u = uTimeAtRow(y, level);
        int i = indexSplineCtrlPoint();
        
        if (isSplineValid == true &&
            (int)spline_.get()->getControlPoses()->size() == 0)
        {
            
            fprintf(stderr, "How this happens?!!\n");
            fprintf(stderr, "u = %f\n", u);
            fprintf(stderr, "i = %d\n", i);
            fprintf(stderr, "spline_.get()->getControlPoses()->size() = %d\n",
                    (int)spline_.get()->getControlPoses()->size());
            fprintf(stderr, "isSplieValid = %d\n", isSplineValid);
            fprintf(stderr, "y = %f\n", y);
            fprintf(stderr, "level = %d\n", level);
            fprintf(stderr, "spline_.use_counter() = %ld\n", spline_.use_count());
            fprintf(stderr, "spline_.unique() = %d\n", spline_.unique());
 
            fprintf(stdout, "How this happens?!!\n");
            fprintf(stdout, "u = %f\n", u);
            fprintf(stdout, "i = %d\n", i);
            fprintf(stdout, "spline_.get()->getControlPoses()->size() = %d\n",
                    (int)spline_.get()->getControlPoses()->size());
            fprintf(stdout, "isSplieValid = %d\n", isSplineValid);
            fprintf(stdout, "y = %f\n", y);
            fprintf(stdout, "level = %d\n", level);
            fprintf(stdout, "spline_.use_counter() = %ld\n", spline_.use_count());
            fprintf(stdout, "spline_.unique() = %d\n", spline_.unique());
            
        }
        
        if (isSplineValid == true &&
            spline_.get()->numControlPoses() > 0 &&
            (int)spline_.get()->getControlPoses()->size() > 0 &&
            spline_.get()->order() > 0) // Check spline is valid
        {
            
            Eigen::Matrix<double,4,4> M;
            M = spline_.get()->getPoseOnSpline(u, i);
            return M;
            
        }
        else
        {
            
            SE3 Xi_pose = se3FromSim3(pose->getCamToWorld());
            SE3 Xi_pose_inv = Xi_pose.inverse();
            Eigen::Matrix<double,4,4> M;
            M = Xi_pose_inv.matrix();
            return M;
            
        }
        
    }
    
    //--------------------------------------------------------------------------
    // Get a pose of the row of a distorted image corresponding to the given
    // pixel coordinates in the undistorted image.
    //
    // Input
    //        x, y : [double] x and y coordiantes in undistorted image
    //       level : [positive int] image pyramid level
    // returns
    //        pose : [4x4 double] world to cam
    Eigen::Matrix<double,4,4> getPoseAtRow_Distort(double x, double y,
                                           int level = 0) const
    {
        
        // Lock for accessing
        //boost::shared_lock<boost::shared_mutex> lock = getActiveLock();
        
        double u = uTimeAtRow_Distort(x, y, level);
        int i = indexSplineCtrlPoint();

         // printf("DISTORT: x, y, level = %f %f %d; "
         //      "DISTORT: u = %f; "
         //      "DISTORT: i = %d\n", x, y, level, u, i);
        if (isSplineValid == true &&
            spline_.get()->numControlPoses() > 0 &&
            spline_.get()->order() > 0) // Check spline is valid
        {
            
            Eigen::Matrix<double,4,4> M;
            M = spline_.get()->getPoseOnSpline(u, i);
            return M;
            
        }
        else
        {
            
            if (i != 1)
            {
                fprintf(stderr, "WARNING: "
                        "getPoseAtRow_Distort from getCamToWorld()\n");
                fprintf(stdout, "WARNING: "
                        "getPoseAtRow_Distort from getCamToWorld()\n");
            }
            SE3 Xi_pose = se3FromSim3(pose->getCamToWorld());
            SE3 Xi_pose_inv = Xi_pose.inverse();
            Eigen::Matrix<double,4,4> M;
            M = Xi_pose_inv.matrix();
            return M;
            
        }
        
    }
    
    //--------------------------------------------------------------------------
    // Get a ground truth pose of the row
    // returns
    //       pose : [4x4 double] world to cam
    Eigen::Matrix<double,4,4> getGTPoseAtRow(int y, int level = 0) const
    {
        
        // Lock for accessing RT_motion
        //boost::shared_lock<boost::shared_mutex> lock = getActiveLock();
        
        Eigen::Matrix<double,4,4> poseMat;

        // Row index
        int rowIndex = y;
        if ((rowIndex < 0) || (rowIndex >= height()))
        {
            
            fprintf(stderr, "Error: rowIndex in getGTposeAtRow()\n");
            exit(1); 
            
        }
        
        poseMat = RT_motion.at(rowIndex);
        
        // Get the pose from the matrix
        // Return world to cam
        return poseMat;
        
    }
    
    //--------------------------------------------------------------------------
    // Get a ground truth pose of the row
    // returnsd
    //       pose : [4x4 double] world to cam
    Eigen::Matrix<double,4,4> getGTPoseAtRow_Distort(int x, int y,
                                                     int level = 0) const
    {
        
        // Lock for accessing RT_motion
        //boost::shared_lock<boost::shared_mutex> lock = getActiveLock();
        
        Eigen::Matrix<double,4,4> poseMat;
        
        // Distort coordinates
        float x_dist, y_dist;
        undistorter_->distortPointLevel((float)x, (float)y,
                                        &x_dist, &y_dist, level);
        
        // Row index
        int rowIndex = (int)round(y_dist);
        if ((rowIndex < 0) || (rowIndex >= height()))
        {
            
            fprintf(stderr, "Error: rowIndex in getGTposeAtRow_Distort()\n");
            exit(1);
            
        }
        
        poseMat = RT_motion.at(rowIndex);
        
        // Get the pose from the matrix
        // Return world to cam
        return poseMat;
        
    }
    
    //--------------------------------------------------------------------------
    // Back-projection
    //
    //      pt : [3x1 double] 2D input image point
    //   level : [int] Pyramid level 0 (original) to 5 (coarse)
    //
    // returns
    //
    //      X : [4x1 double] 3D point wrt the world
    //
    Eigen::Matrix<double,4,1> backProjection(Eigen::Matrix<double,3,1> pt,
                                             int level = 0)
    {
        
        Eigen::Matrix<double,4,1> X;
        double y = pt[1]/pt[2];
        Eigen::Matrix<double,4,4> M = getPoseAtRow(y, level); // World to Cam
        Eigen::Matrix<double,3,3> R = M.block(0,0,3,3);
        Eigen::Matrix<double,3,1> t = M.block(0,3,3,1);
        X << R.inverse()*(KInv(level).cast<double>()*pt - t),
             1.0;
        
        return X;
    }
    
    //--------------------------------------------------------------------------
    // Projection of a 3D point X
    //      X : [4x1 double] a 3D point wrt the world
    // returns
    //      xs : [3xn double] n 2D image points
    //
    Eigen::Matrix<double,3,Eigen::Dynamic>
    projection(Eigen::Matrix<double,4,1> X, int level = 0)
    {
        
        Eigen::Matrix<double,3,Eigen::Dynamic> xs;
        
        int height = this->height(level);
        for (int i=0;i<height;i++)
        {
            
            // Motion at row i
            Eigen::Matrix<double,4,4> M = getPoseAtRow(i, level); // WorldToCam
            Eigen::Matrix<double,3,4> Rt = M.block(0,0,3,4);
            
            // Projection
            const Eigen::Matrix<double,3,1> m = K(level).cast<double>()*Rt*X;
            Eigen::Matrix<double,3,1> x = JHK::MyMats::hcoord(m);
            
            // Check if the projected point x belongs to the row i
            if (fabs(x[1] - i) < 0.5) // A half pixel
            {
                // Append
                xs << xs, x;
            
            }
            
        } // end for
        
        return xs;
        
    }
    
    //==========================================================================


	Frame(int id,
          int width,
          int height,
          const Eigen::Matrix3f& K,
          double timestamp,
          const unsigned char* image,
          const unsigned char* imageMask = 0,
          const int height_rollingShutter_inputImage = 0);

	Frame(int id,
          int width,
          int height,
          const Eigen::Matrix3f& K,
          double timestamp,
          const float* image,
          const float* imageMask = 0,
          const int height_rollingShutter_inputImage = 0);

	~Frame();

    // Copy constructor
    Frame(const Frame &copyOfFrame);
    
    // Copy assignment operator, DISABLED
    Frame &operator= (const Frame & copyOfFrame) = delete;
	
	
	/** Sets or updates idepth and idepthVar on level zero. Invalidates higher levels. */
	void setDepth(const DepthMapPixelHypothesis* newDepth);

	/** Calculates mean information for statistical purposes. */
	void calculateMeanInformation();
	
	/** Sets ground truth depth (real, not inverse!) from a float array on level zero. Invalidates higher levels. */
	void setDepthFromGroundTruth(const float* depth, float cov_scale = 1.0f);
	
	/** Prepares this frame for stereo comparisons with the other frame (computes some intermediate values that will be needed) */
	void prepareForStereoWith(Frame* other, Sim3 thisToOther, const Eigen::Matrix3f& K, const int level);

	

	// Accessors
	/** Returns the unique frame id. */
	inline int id() const;
	
	/** Returns the frame's image width. */
	inline int width(int level = 0) const;
	/** Returns the frame's image height. */
	inline int height(int level = 0) const;
	
	/** Returns the frame's intrinsics matrix. */
	inline const Eigen::Matrix3f& K(int level = 0) const;
	/** Returns the frame's inverse intrinsics matrix. */
	inline const Eigen::Matrix3f& KInv(int level = 0) const;
	/** Returns K(0, 0). */
	inline float fx(int level = 0) const;
	/** Returns K(1, 1). */
	inline float fy(int level = 0) const;
	/** Returns K(0, 2). */
	inline float cx(int level = 0) const;
	/** Returns K(1, 2). */
	inline float cy(int level = 0) const;
	/** Returns KInv(0, 0). */
	inline float fxInv(int level = 0) const;
	/** Returns KInv(1, 1). */
	inline float fyInv(int level = 0) const;
	/** Returns KInv(0, 2). */
	inline float cxInv(int level = 0) const;
	/** Returns KInv(1, 2). */
	inline float cyInv(int level = 0) const;
	
	/** Returns the frame's recording timestamp. */
	inline double timestamp() const;
	
	inline float* image(int level = 0);
    inline float* imageMask(int level = 0);
	inline const Eigen::Vector4f* gradients(int level = 0);
	inline const float* maxGradients(int level = 0);
	inline bool hasIDepthBeenSet() const;
	inline const float* idepth(int level = 0);
	inline const float* idepthVar(int level = 0);
	inline const unsigned char* validity_reAct();
	inline const float* idepth_reAct();
	inline const float* idepthVar_reAct();

	inline bool* refPixelWasGood();
	inline bool* refPixelWasGoodNoCreate();
	inline void clear_refPixelWasGood();

	/** Flags for use with require() and requirePyramid(). See the Frame class
	  * documentation for their exact meaning. */
	enum DataFlags
	{
		IMAGE			= 1<<0,
		GRADIENTS		= 1<<1,
		MAX_GRADIENTS	= 1<<2,
		IDEPTH			= 1<<3,
		IDEPTH_VAR		= 1<<4,
		REF_ID			= 1<<5,
        IMAGE_MASK      = 1<<6,
		
		ALL = IMAGE | GRADIENTS | MAX_GRADIENTS | IDEPTH | IDEPTH_VAR | REF_ID | IMAGE_MASK
	};
	

	void setPermaRef(TrackingReference* reference);
	void takeReActivationData(DepthMapPixelHypothesis* depthMap);


	// shared_lock this as long as any minimizable arrays are being used.
	// the minimizer will only minimize frames after getting
	// an exclusive lock on this.
	inline boost::shared_lock<boost::shared_mutex> getActiveLock()
	{
		return FrameMemory::getInstance().activateFrame(this);
	}


	/*
	 * ==================================================================================
	 * Here are ALL central pose and scale informations.
	 * generally, everything is stored relative to the frame
	 */
	FramePoseStruct* pose;
	Sim3 getScaledCamToWorld(int num=0) { return pose->getCamToWorld();}
	bool hasTrackingParent() { return pose->trackingParent != nullptr;}
	Frame* getTrackingParent() { return pose->trackingParent->frame;}

	Sim3 lastConstraintTrackedCamToWorld;



	/** Pointers to all adjacent Frames in graph. empty for non-keyframes.*/
	std::unordered_set< Frame* > neighbors;

	/** Multi-Map indicating for which other keyframes with which initialization tracking failed.*/
	std::unordered_multimap< Frame*, Sim3 > trackingFailed;


	// flag set when depth is updated.
	bool depthHasBeenUpdatedFlag;


	// Tracking Reference for quick test. Always available, never taken out of memory.
	// this is used for re-localization and re-Keyframe positioning.
	boost::mutex permaRef_mutex;
	Eigen::Vector3f* permaRef_posData;	// (x,y,z)
	Eigen::Vector2f* permaRef_colorAndVarData;	// (I, Var)
	int permaRefNumPts;



	// Temporary values
	int referenceID;
	int referenceLevel;
	float distSquared;
	Eigen::Matrix3f K_otherToThis_R;
	Eigen::Vector3f K_otherToThis_t;
	Eigen::Vector3f otherToThis_t;
	Eigen::Vector3f K_thisToOther_t;
	Eigen::Matrix3f thisToOther_R;
	Eigen::Vector3f otherToThis_R_row0;
	Eigen::Vector3f otherToThis_R_row1;
	Eigen::Vector3f otherToThis_R_row2;
	Eigen::Vector3f thisToOther_t;



	// statistics
	float initialTrackedResidual;
	int numFramesTrackedOnThis;
	int numMappedOnThis;
	int numMappedOnThisTotal;
	float meanIdepth;
	int numPoints;
	int idxInKeyframes;
	float edgeErrorSum, edgesNum;
	int numMappablePixels;
	float meanInformation;
    
    // Rolling shutter height
    int height_rollingShutter_inputImage;

private:

	void require(int dataFlags, int level = 0);
	void release(int dataFlags, bool pyramidsOnly, bool invalidateOnly);

	void initialize(int id, int width, int height,
                    const Eigen::Matrix3f& K, double timestamp,
                    int height_rollingShutter_inputImage = 0);
	void setDepth_Allocate();
	
	void buildImage(int level);
	void releaseImage(int level);

    void buildImageMask(int level);
    void releaseImageMask(int level);

	void buildGradients(int level);
	void releaseGradients(int level);
	
	void buildMaxGradients(int level);
	void releaseMaxGradients(int level);
	
	void buildIDepthAndIDepthVar(int level);
	void releaseIDepth(int level);
	void releaseIDepthVar(int level);
	
	void printfAssert(const char* message) const;
	
	struct Data
	{
		int id;
		
		int width[PYRAMID_LEVELS], height[PYRAMID_LEVELS];

		Eigen::Matrix3f K[PYRAMID_LEVELS], KInv[PYRAMID_LEVELS];
		float fx[PYRAMID_LEVELS], fy[PYRAMID_LEVELS], cx[PYRAMID_LEVELS], cy[PYRAMID_LEVELS];
		float fxInv[PYRAMID_LEVELS], fyInv[PYRAMID_LEVELS], cxInv[PYRAMID_LEVELS], cyInv[PYRAMID_LEVELS];
		
		double timestamp;

		
		float* image[PYRAMID_LEVELS];
		bool imageValid[PYRAMID_LEVELS];

        float* imageMask[PYRAMID_LEVELS]; // Mask from ALPHA channel
        bool imageMaskValid[PYRAMID_LEVELS];
		
		Eigen::Vector4f* gradients[PYRAMID_LEVELS];
		bool gradientsValid[PYRAMID_LEVELS];
		
		float* maxGradients[PYRAMID_LEVELS];
		bool maxGradientsValid[PYRAMID_LEVELS];
		

		bool hasIDepthBeenSet;

		// negative depthvalues are actually allowed, so setting this to -1 does NOT invalidate the pixel's depth.
		// a pixel is valid iff idepthVar[i] > 0.
		float* idepth[PYRAMID_LEVELS];
		bool idepthValid[PYRAMID_LEVELS];
		
		// MUST contain -1 for invalid pixel (that dont have depth)!!
		float* idepthVar[PYRAMID_LEVELS];
		bool idepthVarValid[PYRAMID_LEVELS];

		// data needed for re-activating the frame. theoretically, this is all data the
		// frame contains.
		unsigned char* validity_reAct;
		float* idepth_reAct;
		float* idepthVar_reAct;
		bool reActivationDataValid;


		// data from initial tracking, indicating which pixels in the reference frame ware good or not.
		// deleted as soon as frame is used for mapping.
		bool* refPixelWasGood;
	};
	Data data;


	// used internally. locked while something is being built, such that no
	// two threads build anything simultaneously. not locked on require() if nothing is changed.
	boost::mutex buildMutex;

	boost::shared_mutex activeMutex;
	bool isActive;

	/** Releases everything which can be recalculated, but keeps the minimal
	  * representation in memory. Use release(Frame::ALL, false) to store on disk instead.
	  * ONLY CALL THIS, if an exclusive lock on activeMutex is owned! */
	bool minimizeInMemory();
};



inline int Frame::id() const
{
	return data.id;
}

inline int Frame::width(int level) const
{
	return data.width[level];
}

inline int Frame::height(int level) const
{
	return data.height[level];
}

inline const Eigen::Matrix3f& Frame::K(int level) const
{
	return data.K[level];
}
inline const Eigen::Matrix3f& Frame::KInv(int level) const
{
	return data.KInv[level];
}
inline float Frame::fx(int level) const
{
	return data.fx[level];
}
inline float Frame::fy(int level) const
{
	return data.fy[level];
}
inline float Frame::cx(int level) const
{
	return data.cx[level];
}
inline float Frame::cy(int level) const
{
	return data.cy[level];
}
inline float Frame::fxInv(int level) const
{
	return data.fxInv[level];
}
inline float Frame::fyInv(int level) const
{
	return data.fyInv[level];
}
inline float Frame::cxInv(int level) const
{
	return data.cxInv[level];
}
inline float Frame::cyInv(int level) const
{
	return data.cyInv[level];
}

inline double Frame::timestamp() const
{
	return data.timestamp;
}


inline float* Frame::image(int level)
{
	if (! data.imageValid[level])
		require(IMAGE, level);
	return data.image[level];
}

inline float* Frame::imageMask(int level)
{
    if (! data.imageMaskValid[level])
        require(IMAGE_MASK, level);
    return data.imageMask[level];
}

inline const Eigen::Vector4f* Frame::gradients(int level)
{
	if (! data.gradientsValid[level])
		require(GRADIENTS, level);
	return data.gradients[level];
}
inline const float* Frame::maxGradients(int level)
{
	if (! data.maxGradientsValid[level])
		require(MAX_GRADIENTS, level);
	return data.maxGradients[level];
}
inline bool Frame::hasIDepthBeenSet() const
{
	return data.hasIDepthBeenSet;
}
inline const float* Frame::idepth(int level)
{
	if (! data.hasIDepthBeenSet)
	{
		printfAssert("Frame::idepth(): idepth has not been set yet!");
		return nullptr;
	}
	if (! data.idepthValid[level])
		require(IDEPTH, level);
	return data.idepth[level];
}
inline const unsigned char* Frame::validity_reAct()
{
	if( !data.reActivationDataValid)
		return 0;
	return data.validity_reAct;
}
inline const float* Frame::idepth_reAct()
{
	if( !data.reActivationDataValid)
		return 0;
	return data.idepth_reAct;
}
inline const float* Frame::idepthVar_reAct()
{
	if( !data.reActivationDataValid)
		return 0;
	return data.idepthVar_reAct;
}
inline const float* Frame::idepthVar(int level)
{
	if (! data.hasIDepthBeenSet)
	{
		printfAssert("Frame::idepthVar(): idepth has not been set yet!");
		return nullptr;
	}
	if (! data.idepthVarValid[level])
		require(IDEPTH_VAR, level);
	return data.idepthVar[level];
}


inline bool* Frame::refPixelWasGood()
{
	if( data.refPixelWasGood == 0)
	{
		boost::unique_lock<boost::mutex> lock2(buildMutex);

		if(data.refPixelWasGood == 0)
		{
			int width = data.width[SE3TRACKING_MIN_LEVEL];
			int height = data.height[SE3TRACKING_MIN_LEVEL];
			data.refPixelWasGood = (bool*)FrameMemory::getInstance().getBuffer(sizeof(bool) * width * height);

			memset(data.refPixelWasGood, 0xFFFFFFFF, sizeof(bool) * (width * height));
		}
	}
	return data.refPixelWasGood;
}


inline bool* Frame::refPixelWasGoodNoCreate()
{
	return data.refPixelWasGood;
}

inline void Frame::clear_refPixelWasGood()
{
	FrameMemory::getInstance().returnBuffer(reinterpret_cast<float*>(data.refPixelWasGood));
	data.refPixelWasGood=0;
}


}
