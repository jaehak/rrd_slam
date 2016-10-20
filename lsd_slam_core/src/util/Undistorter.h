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

#ifndef _UNDISTORTER_HPP_
#define _UNDISTORTER_HPP_

#include <opencv2/core/core.hpp>



namespace lsd_slam
{

class Undistorter
{
public:
	virtual ~Undistorter();
	
	/**
	 * Undistorts the given image and returns the result image.
	 */
	virtual void undistort(const cv::Mat& image, cv::OutputArray result) const = 0;
	
	/**
	 * Returns the intrinsic parameter matrix of the undistorted images.
	 */
	virtual const cv::Mat& getK() const = 0;
	
	/**
	 * Returns the intrinsic parameter matrix of the original images,
	 */
	virtual const cv::Mat& getOriginalK() const = 0;
	
	/**
	 * Returns the width of the undistorted images in pixels.
	 */
	virtual int getOutputWidth() const = 0;
	
	/**
	 * Returns the height of the undistorted images in pixels.
	 */
	virtual int getOutputHeight() const = 0;
	
	/**
	 * Returns the width of the input images in pixels.
	 */
	virtual int getInputWidth() const = 0;

	/**
	 * Returns the height of the input images in pixels.
	 */
	virtual int getInputHeight() const = 0;

    // Distortion parameters
    virtual float getInputDistort() const = 0;
    
    virtual float getOutputDistort() const = 0;
    
	/**
	 * Returns if the undistorter was initialized successfully.
	 */
	virtual bool isValid() const = 0;
	
	/**
	 * Creates and returns an Undistorter of the type used by the given
	 * configuration file. If the format is not recognized, returns nullptr.
	 */
	static Undistorter* getUndistorterForFile(const char* configFilename);

    /**
     * Distort a pixel point and returns the distorted pixel point
     */
    virtual void distortPoint(float px, float py,
                              float *px_out, float *py_out) const = 0;
    virtual void distortPointLevel(float px, float py,
                              float *px_out, float *py_out,
                              int level) const = 0;
    virtual void distortGray(const cv::Mat& image,
                             cv::OutputArray result) const = 0;
    
    virtual void distortRGB(const cv::Mat& image,
                            cv::OutputArray result) const = 0;
    
    virtual void undistortRGB(const cv::Mat& image,
                              cv::OutputArray result) const = 0;

    // Undistort point
    virtual void undistortPoint(float x_dist, float y_dist,
                        float *x_undist, float *y_undist) const = 0;
    
    /**
     * Returns if the original input is already rectified
     */
    virtual bool isNoRectification() const = 0;

};

class UndistorterPTAM : public Undistorter
{
public:
	/**
	 * Creates an Undistorter by reading the distortion parameters from a file.
	 * 
	 * The file format is as follows:
	 * d1 d2 d3 d4 d5
	 * inputWidth inputHeight
	 * crop / full / none
	 * outputWidth outputHeight
	 */
	UndistorterPTAM(const char* configFileName);
	
	/**
	 * Destructor.
	 */
	~UndistorterPTAM();
	
	UndistorterPTAM(const UndistorterPTAM&) = delete;
	UndistorterPTAM& operator=(const UndistorterPTAM&) = delete;
	
	/**
	 * Undistorts the given image and returns the result image.
	 */
	void undistort(const cv::Mat& image, cv::OutputArray result) const;
	
	/**
	 * Returns the intrinsic parameter matrix of the undistorted images.
	 */
	const cv::Mat& getK() const;
	
	/**
	 * Returns the intrinsic parameter matrix of the original images,
	 */
	const cv::Mat& getOriginalK() const;
	
	/**
	 * Returns the width of the undistorted images in pixels.
	 */
	int getOutputWidth() const;
	
	/**
	 * Returns the height of the undistorted images in pixels.
	 */
	int getOutputHeight() const;
	
	/**
	 * Returns the width of the input images in pixels.
	 */
	int getInputWidth() const;

	/**
	 * Returns the height of the input images in pixels.
	 */
	int getInputHeight() const;


    // Distortion parameters
    float getInputDistort() const;
    
    float getOutputDistort() const;
    
	/**
	 * Returns if the undistorter was initialized successfully.
	 */
	bool isValid() const;
    

    /**
     * Distort a pixel point and returns the distorted pixel point
     */
    void distortPoint(float px, float py,
                              float *px_out, float *py_out) const;
    void distortPointLevel(float px, float py,
                           float *px_out, float *py_out,
                           int level) const;
    void distortGray(const cv::Mat& image,
                     cv::OutputArray result) const;
    void distortRGB(const cv::Mat& image,
                    cv::OutputArray result) const;
    
    /**
     * Undistort RGB image
     */
    void undistortRGB(const cv::Mat& image, cv::OutputArray result) const;

    // Undistort point
    void undistortPoint(float x_dist, float y_dist,
                        float *x_undist, float *y_undist) const;
    
    /**
     * Returns if the original input is already rectified
     */
    bool isNoRectification() const;
    
private:
    // NOTE: K_ and K_original are
    //       a transpose of the standard form of calibtration matrix
    //       K_ = [fx 0 0; 0 fy 0; cx cy 1]
    //       K_original = [fx 0 0; 0 fy 0; cx cy 1]
    cv::Mat K_;
	cv::Mat originalK_;
	
	float inputCalibration[5];
	float outputCalibration[5];
	int out_width, out_height;
	int in_width, in_height;
	float* remapX;
	float* remapY;

	
	/// Is true if the undistorter object is valid (has been initialized with
	/// a valid configuration)
	bool valid;
    
    bool isNoRect_;
    
};

class UndistorterOpenCV : public Undistorter
{
public:
	/**
	 * Creates an Undistorter by reading the distortion parameters from a file.
	 * 
	 * The file format is as follows:
	 * fx fy cx cy d1 d2 d3 d4 d5 d6
	 * inputWidth inputHeight
	 * crop / full / none
	 * outputWidth outputHeight
	 */
	UndistorterOpenCV(const char* configFileName);
	
	/**
	 * Destructor.
	 */
	~UndistorterOpenCV();
	
	UndistorterOpenCV(const UndistorterOpenCV&) = delete;
	UndistorterOpenCV& operator=(const UndistorterOpenCV&) = delete;
	
	/**
	 * Undistorts the given image and returns the result image.
	 */
	void undistort(const cv::Mat& image, cv::OutputArray result) const;
	
	/**
	 * Returns the intrinsic parameter matrix of the undistorted images.
	 */
	const cv::Mat& getK() const;
	
	/**
	 * Returns the intrinsic parameter matrix of the original images,
	 */
	const cv::Mat& getOriginalK() const;
	
	/**
	 * Returns the width of the undistorted images in pixels.
	 */
	int getOutputWidth() const;
	
	/**
	 * Returns the height of the undistorted images in pixels.
	 */
	int getOutputHeight() const;
	

	/**
	 * Returns the width of the input images in pixels.
	 */
	int getInputWidth() const;

	/**
	 * Returns the height of the input images in pixels.
	 */
	int getInputHeight() const;

    // Distortion parameters
    float getInputDistort() const;
    
    float getOutputDistort() const;
    
    /**
	 * Returns if the undistorter was initialized successfully.
	 */
	bool isValid() const;

    /**
     * Distort a pixel point and returns the distorted pixel point
     */
    void distortPoint(float px, float py,
                              float *px_out, float *py_out) const;
    void distortPointLevel(float x_undist, float y_undist,
                           float *x_dist, float *y_dist,
                           int level = 0) const;
    void distortGray(const cv::Mat& image,
                     cv::OutputArray result) const;
    void distortRGB(const cv::Mat& image,
                    cv::OutputArray result) const;
    
    // Undistort RGB image - Test required
    void undistortRGB(const cv::Mat& image, cv::OutputArray result) const;

    // Undistort point
    void undistortPoint(float x_dist, float y_dist,
                        float *x_undist, float *y_undist) const;
    
    /**
     * Returns if the original input is already rectified
     */
    bool isNoRectification() const;

private:
	cv::Mat K_;
	cv::Mat originalK_;
	
	float inputCalibration[10];
	float outputCalibration;
	int out_width, out_height;
	int in_width, in_height;
	cv::Mat map1, map2;
	
	/// Is true if the undistorter object is valid (has been initialized with
	/// a valid configuration)
	bool valid;
    
    bool isNoRect_;
    
};

} // end of namespace
#endif
