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

#include "Undistorter.h"

#include <sstream>
#include <fstream>

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include "util/settings.h"

namespace lsd_slam
{

Undistorter::~Undistorter()
{
}

Undistorter* Undistorter::getUndistorterForFile(const char* configFilename)
{
	std::string completeFileName = configFilename;

	printf("Reading Calibration from file %s",completeFileName.c_str());


	std::ifstream f(completeFileName.c_str());
	if (!f.good())
	{
		f.close();

		completeFileName = packagePath+"calib/"+configFilename;
		printf(" ... not found!\n Trying %s", completeFileName.c_str());

		f.open(completeFileName.c_str());

		if (!f.good())
		{
			printf(" ... not found. Cannot operate without calibration, shutting down.\n");
			f.close();
			return 0;
		}
	}

	printf(" ... found!\n");

	std::string l1;
	std::getline(f,l1);
	f.close();



	float ic[10];
	if(std::sscanf(l1.c_str(), "%f %f %f %f %f %f %f %f",
			&ic[0], &ic[1], &ic[2], &ic[3], &ic[4],
			&ic[5], &ic[6], &ic[7]) == 8)
	{
		printf("found OpenCV camera model, building rectifier.\n");
		Undistorter* u = new UndistorterOpenCV(completeFileName.c_str());
		if(!u->isValid()) return 0;
        
		return u;
	}
	else
	{
		printf("found ATAN camera model, building rectifier.\n");
		Undistorter* u = new UndistorterPTAM(completeFileName.c_str());
		if(!u->isValid()) return 0;
		return u;
	}
}


UndistorterPTAM::UndistorterPTAM(const char* configFileName)
{
    isNoRect_ = false;
	valid = true;

	remapX = nullptr;
	remapY = nullptr;

	
	
	// read parameters
	std::ifstream infile(configFileName);
	assert(infile.good());


	std::string l1,l2,l3,l4;

	std::getline(infile,l1);
	std::getline(infile,l2);
	std::getline(infile,l3);
	std::getline(infile,l4);




	// l1 & l2
	if(std::sscanf(l1.c_str(), "%f %f %f %f %f", &inputCalibration[0], &inputCalibration[1], &inputCalibration[2], &inputCalibration[3], &inputCalibration[4]) == 5 &&
			std::sscanf(l2.c_str(), "%d %d", &in_width, &in_height) == 2)
	{
		printf("Input resolution: %d %d\n",in_width, in_height);
		printf("In: %f %f %f %f %f\n",
				inputCalibration[0], inputCalibration[1], inputCalibration[2], inputCalibration[3], inputCalibration[4]);
	}
	else
	{
		printf("Failed to read camera calibration (invalid format?)\nCalibration file: %s\n", configFileName);
		valid = false;
	}

	// l3
	if(l3 == "crop")
	{
		outputCalibration[0] = -1;
		printf("Out: Crop\n");
	}
	else if(l3 == "full")
	{
		outputCalibration[0] = -2;
		printf("Out: Full\n");
	}
    else if(l3 == "fit")
    {
        outputCalibration[0] = -3;
        printf("Out: Fit\n");
    }
	else if(l3 == "none")
	{
		printf("NO RECTIFICATION\n");
        valid = true;
        isNoRect_ = true;
	}
	else if(std::sscanf(l3.c_str(), "%f %f %f %f %f", &outputCalibration[0], &outputCalibration[1], &outputCalibration[2], &outputCalibration[3], &outputCalibration[4]) == 5)
	{
		printf("Out: %f %f %f %f %f\n",
				outputCalibration[0], outputCalibration[1], outputCalibration[2], outputCalibration[3], outputCalibration[4]);
	}
	else
	{
		printf("Out: Failed to Read Output pars... not rectifying.\n");
		valid = false;
	}


	// l4
	if(std::sscanf(l4.c_str(), "%d %d", &out_width, &out_height) == 2)
	{
		printf("Output resolution: %d %d\n",out_width, out_height);
	}
	else
	{
		printf("Out: Failed to Read Output resolution... not rectifying.\n");
		valid = false;
	}




	// prep warp matrices
	if(valid)
	{
		float dist = inputCalibration[4];
		float d2t = 2.0f * tanf(dist / 2.0f);

		// current camera parameters
		float fx = inputCalibration[0] * in_width;
		float fy = inputCalibration[1] * in_height;
		float cx = inputCalibration[2] * in_width - 0.5f;
		float cy = inputCalibration[3] * in_height - 0.5f;
		
		// scale calibration parameters to input size
		float xfactor = in_width / (1.0f * in_width);
		float yfactor = in_height / (1.0f * in_height);
		fx = fx * xfactor;
		fy = fy * yfactor;
		cx = (cx + 0.5f) * xfactor - 0.5f;
		cy = (cy + 0.5f) * yfactor - 0.5f;

		// output camera parameters
		float ofx, ofy, ocx, ocy;

		// find new camera matrix for "crop" and "full"
		if (inputCalibration[4] == 0)
		{
			ofx = inputCalibration[0] * out_width;
			ofy = inputCalibration[1] * out_height;
			ocx = (inputCalibration[2] * out_width) - 0.5f;
			ocy = (inputCalibration[3] * out_height) - 0.5f;
		}
		else if(outputCalibration[0] == -1)	// "crop"
		{
			// find left-most and right-most radius
			float left_radius = (cx)/fx;
			float right_radius = (in_width-1 - cx)/fx;
			float top_radius = (cy)/fy;
			float bottom_radius = (in_height-1 - cy)/fy;

			float trans_left_radius = tanf(left_radius * dist)/d2t;
			float trans_right_radius = tanf(right_radius * dist)/d2t;
			float trans_top_radius = tanf(top_radius * dist)/d2t;
			float trans_bottom_radius = tanf(bottom_radius * dist)/d2t;

			//printf("left_radius: %f -> %f\n", left_radius, trans_left_radius);
			//printf("right_radius: %f -> %f\n", right_radius, trans_right_radius);
			//printf("top_radius: %f -> %f\n", top_radius, trans_top_radius);
			//printf("bottom_radius: %f -> %f\n", bottom_radius, trans_bottom_radius);


			ofy = fy * ((top_radius + bottom_radius) / (trans_top_radius + trans_bottom_radius)) * ((float)out_height / (float)in_height);
			ocy = (trans_top_radius/top_radius) * ofy*(cy + 0.5f)/fy - 0.5f;

			ofx = fx * ((left_radius + right_radius) / (trans_left_radius + trans_right_radius)) * ((float)out_width / (float)in_width);
			ocx = (trans_left_radius/left_radius) * ofx*(cx + 0.5f)/fx - 0.5f;

			printf("new K: %f %f %f %f\n",ofx,ofy,ocx,ocy);
			printf("old K: %f %f %f %f\n",fx,fy,cx,cy);
		}
		else if(outputCalibration[0] == -2)	 // "full"
		{
			float left_radius = cx/fx;
			float right_radius = (in_width-1 - cx)/fx;
			float top_radius = cy/fy;
			float bottom_radius = (in_height-1 - cy)/fy;

			// find left-most and right-most radius
			float tl_radius = sqrtf(left_radius*left_radius + top_radius*top_radius);
			float tr_radius = sqrtf(right_radius*right_radius + top_radius*top_radius);
			float bl_radius = sqrtf(left_radius*left_radius + bottom_radius*bottom_radius);
			float br_radius = sqrtf(right_radius*right_radius + bottom_radius*bottom_radius);

			float trans_tl_radius = tanf(tl_radius * dist)/d2t;
			float trans_tr_radius = tanf(tr_radius * dist)/d2t;
			float trans_bl_radius = tanf(bl_radius * dist)/d2t;
			float trans_br_radius = tanf(br_radius * dist)/d2t;

			//printf("trans_tl_radius: %f -> %f\n", tl_radius, trans_tl_radius);
			//printf("trans_tr_radius: %f -> %f\n", tr_radius, trans_tr_radius);
			//printf("trans_bl_radius: %f -> %f\n", bl_radius, trans_bl_radius);
			//printf("trans_br_radius: %f -> %f\n", br_radius, trans_br_radius);


			float hor = std::max(br_radius,tr_radius) + std::max(bl_radius,tl_radius);
			float vert = std::max(tr_radius,tl_radius) + std::max(bl_radius,br_radius);

			float trans_hor = std::max(trans_br_radius,trans_tr_radius) + std::max(trans_bl_radius,trans_tl_radius);
			float trans_vert = std::max(trans_tr_radius,trans_tl_radius) + std::max(trans_bl_radius,trans_br_radius);

			ofy = fy * ((vert) / (trans_vert)) * ((float)out_height / (float)in_height);
			ocy = std::max(trans_tl_radius/tl_radius,trans_tr_radius/tr_radius) * ofy*cy/fy;

			ofx = fx * ((hor) / (trans_hor)) * ((float)out_width / (float)in_width);
			ocx = std::max(trans_bl_radius/bl_radius,trans_tl_radius/tl_radius) * ofx*cx/fx;

			printf("new K: %f %f %f %f\n",ofx,ofy,ocx,ocy);
			printf("old K: %f %f %f %f\n",fx,fy,cx,cy);
		}
        else if(outputCalibration[0] == -3)	 // "fit"
        {
            
            // crop canvas with fixed aspect ratio of output
            float newHeight = in_height*(out_height/out_width);
            float newWidth = in_width*(out_height/out_width);
            float newLeft = in_width/2.0f - newWidth/2.0f;
            float newRight = in_width/2.0f + newWidth/2.0f;
            float newTop = in_height/2.0f - newHeight/2.0f;
            float newBottom = in_height/2.0f + newHeight/2.0f;
            
            // find left-most and right-most radius
            float left_radius = (cx - newLeft)/fx;
            float right_radius = (newRight-1 - cx)/fx;
            float top_radius = (cy - newTop)/fy;
            float bottom_radius = (newBottom-1 - cy)/fy;
            
            float trans_left_radius = tanf(left_radius * dist)/d2t;
            float trans_right_radius = tanf(right_radius * dist)/d2t;
            float trans_top_radius = tanf(top_radius * dist)/d2t;
            float trans_bottom_radius = tanf(bottom_radius * dist)/d2t;
            
            //printf("left_radius: %f -> %f\n", left_radius, trans_left_radius);
            //printf("right_radius: %f -> %f\n", right_radius, trans_right_radius);
            //printf("top_radius: %f -> %f\n", top_radius, trans_top_radius);
            //printf("bottom_radius: %f -> %f\n", bottom_radius, trans_bottom_radius);
            
            
            ofy = fy * ((top_radius + bottom_radius) / (trans_top_radius + trans_bottom_radius)) * ((float)out_height / (float)in_height);
            ocy = (trans_top_radius/top_radius) * ofy*(cy + 0.5f)/fy - 0.5f;
            
            ofx = fx * ((left_radius + right_radius) / (trans_left_radius + trans_right_radius)) * ((float)out_width / (float)in_width);
            ocx = (trans_left_radius/left_radius) * ofx*(cx + 0.5f)/fx - 0.5f;
            
            printf("new K: %f %f %f %f\n",ofx,ofy,ocx,ocy);
            printf("old K: %f %f %f %f\n",fx,fy,cx,cy);
        
        }
		else
		{
			ofx = outputCalibration[0] * out_width;
			ofy = outputCalibration[1] * out_height;
			ocx = outputCalibration[2] * out_width-0.5f;	// TODO: -0.5 here or not?
			ocy = outputCalibration[3] * out_height-0.5f;
		}

		outputCalibration[0] = ofx / out_width;
		outputCalibration[1] = ofy / out_height;
		outputCalibration[2] = (ocx+0.5f) / out_width;
		outputCalibration[3] = (ocy+0.5f) / out_height;
		outputCalibration[4] = 0;

        // The "remapX/Y" is a mapping from undistort coords to distort coords
        // such that remapX/Y[undist_coords] = distort_coords.
        //
		remapX = new float[out_width * out_height];
		remapY = new float[out_width * out_height];

		for(int y=0;y<out_height;y++)
		{
			for(int x=0;x<out_width;x++)
			{
				float ix = (x - ocx) / ofx;
				float iy = (y - ocy) / ofy;

				float r = sqrtf(ix*ix + iy*iy);
				float fac = (r==0 || dist==0) ? 1 : atanf(r * d2t)/(dist*r);

				ix = fx*fac*ix+cx;
				iy = fy*fac*iy+cy;

				// make rounding resistant.
				if(ix == 0) ix = 0.01f;
				if(iy == 0) iy = 0.01f;
				if(ix == in_width-1) ix = in_width-1.01f;
                
                // Is this correct?? Shoudn't be iy => instead of ix?
                // Commented by Jae-Hak Kim, 20 June 2016
				if(iy == in_height-1) iy = in_height-1.01f;
                
				if(ix > 0 && iy > 0 && ix < in_width-1 &&  iy < in_height-1)
				{
					remapX[x+y*out_width] = ix;
					remapY[x+y*out_width] = iy;
                }
				else
				{
					remapX[x+y*out_width] = -1;
					remapY[x+y*out_width] = -1;
				}
			}
		}

		printf("Prepped Warp matrices\n");
	}
	else
	{
		printf("Not Rectifying\n");
		outputCalibration[0] = inputCalibration[0];
		outputCalibration[1] = inputCalibration[1];
		outputCalibration[2] = inputCalibration[2];
		outputCalibration[3] = inputCalibration[3];
		outputCalibration[4] = inputCalibration[4];
		out_width = in_width;
		out_height = in_height;
	}

	
	originalK_ = cv::Mat(3, 3, CV_64F, cv::Scalar(0));
	originalK_.at<double>(0, 0) = inputCalibration[0] * in_width;
	originalK_.at<double>(1, 1) = inputCalibration[1] * in_height;
	originalK_.at<double>(2, 2) = 1;
	originalK_.at<double>(0, 2) = inputCalibration[2] * in_width - 0.5;
	originalK_.at<double>(1, 2) = inputCalibration[3] * in_height - 0.5;

	K_ = cv::Mat(3, 3, CV_64F, cv::Scalar(0));
	K_.at<double>(0, 0) = outputCalibration[0] * out_width;
	K_.at<double>(1, 1) = outputCalibration[1] * out_height;
	K_.at<double>(2, 2) = 1;
	K_.at<double>(0, 2) = outputCalibration[2] * out_width - 0.5;
	K_.at<double>(1, 2) = outputCalibration[3] * out_height - 0.5;
    
    printf("Out: %f %f %f %f %f\n",
           outputCalibration[0], outputCalibration[1],
           outputCalibration[2], outputCalibration[3],
           outputCalibration[4]);
    
}

UndistorterPTAM::~UndistorterPTAM()
{

	delete[] remapX;
	delete[] remapY;

}

void UndistorterPTAM::undistortPoint(float x_dist, float y_dist,
                                     float *x_undist, float *y_undist) const
{

    // Camera parameters
    float fx = (float)originalK_.at<double>(0, 0);
    float fy = (float)originalK_.at<double>(1, 1);
    float cx = (float)originalK_.at<double>(0, 2);
    float cy = (float)originalK_.at<double>(1, 2);
    float ofx = (float)K_.at<double>(0, 0);
    float ofy = (float)K_.at<double>(1, 1);
    float ocx = (float)K_.at<double>(0, 2);
    float ocy = (float)K_.at<double>(1, 2);
    float w_ = getInputDistort();
    
    // Normalized distroted points
    float xp = (x_dist - cx)/fx;
    float yp = (y_dist - cy)/fy;
    
    // Radius
    float rp = sqrtf(xp*xp + yp*yp);
    float fac = (rp == 0 || w_ == 0) ? 1 : tanf(rp*w_)/(2*tanf(w_/2))/rp;
    
    // Normalized undistorted x and y
    float x = fac*xp;
    float y = fac*yp;
    
    // Undistorted x and y in the image
    *x_undist = ofx*x + ocx;
    *y_undist = ofy*y + ocy;
    
}
    
void UndistorterPTAM::undistort(const cv::Mat& image, cv::OutputArray result) const
{
	if (!valid)
	{
		result.getMatRef() = image;
		return;
	}
	
	if (image.rows != in_height || image.cols != in_width)
	{
		printf("UndistorterPTAM: input image size differs from expected input size! Not undistorting.\n");
		result.getMatRef() = image;
		return;
	}
	
	if (in_height == out_height && in_width == out_width && inputCalibration[4] == 0)
	{
		// No transformation if neither distortion nor resize
		result.getMatRef() = image;
		return;
	}
	
	result.create(out_height, out_width, CV_8U);
	cv::Mat resultMat = result.getMatRef();
	assert(result.getMatRef().isContinuous());
	assert(image.isContinuous());
	
	uchar* data = resultMat.data;

	for(int idx = out_width*out_height-1;idx>=0;idx--)
	{
		// get interp. values
		float xx = remapX[idx];
		float yy = remapY[idx];

		if(xx<0)
			data[idx] = 0;
		else
		{
			// get integer and rational parts
			int xxi = (int)xx;
			int yyi = (int)yy;
			xx -= xxi;
			yy -= yyi;
			float xxyy = xx*yy;

			// get array base pointer
			const uchar* src = (uchar*)image.data + xxi + yyi * in_width;

			// interpolate (bilinear)
			data[idx] =  (uchar)(xxyy * src[1+in_width]
			                    + (yy-xxyy) * src[in_width]
			                    + (xx-xxyy) * src[1]
			                    + (1-xx-yy+xxyy) * src[0]);
		}
	}
}

// Undistort RGB image
void UndistorterPTAM::undistortRGB(const cv::Mat& image, cv::OutputArray result) const
{
    if (!valid)
    {
        result.getMatRef() = image;
        return;
    }

    if (image.rows != in_height || image.cols != in_width)
    {
        printf("UndistorterPTAM: input image size differs from expected input size! Not undistorting.\n");
        result.getMatRef() = image;
        return;
    }

    if (in_height == out_height && in_width == out_width && inputCalibration[4] == 0)
    {
        // No transformation if neither distortion nor resize
        result.getMatRef() = image;
        return;
    }

    result.create(out_height, out_width, CV_8UC3);
    cv::Mat resultMat = result.getMatRef();
    assert(result.getMatRef().isContinuous());
    assert(image.isContinuous());

    uchar* data = resultMat.data;
    int cn = image.channels();

    for(int idx = out_width*out_height-1;idx>=0;idx--)
    {
        // get interp. values
        float xx = remapX[idx];
        float yy = remapY[idx];

        if(xx<0)
            data[idx] = 0;
        else
        {
            // get integer and rational parts
            int xxi = (int)xx;
            int yyi = (int)yy;
            xx -= xxi;
            yy -= yyi;
            float xxyy = xx*yy;

            // get array base pointer
            //const uchar* src = (uchar*)image.data + xxi + yyi * in_width;
            uchar* bSrc = (uchar*)image.data + xxi*cn + yyi*cn*in_width + 0;
            uchar* gSrc = (uchar*)image.data + xxi*cn + yyi*cn*in_width + 1;
            uchar* rSrc = (uchar*)image.data + xxi*cn + yyi*cn*in_width + 2;;

            // interpolate (bilinear)
            //data[idx] =  xxyy * src[1+in_width]
            //                    + (yy-xxyy) * src[in_width]
            //                    + (xx-xxyy) * src[1]
            //                    + (1-xx-yy+xxyy) * src[0];
            data[cn*idx + 0] = (uchar)(xxyy*bSrc[cn*(1 + in_width)]
                                + (yy - xxyy)*bSrc[cn*(in_width)]
                                + (xx - xxyy)*bSrc[cn*1]
                                + (1 - xx - yy + xxyy)*bSrc[0]);

            data[cn*idx + 1] = (uchar)(xxyy*gSrc[cn*(1 + in_width)]
                                + (yy - xxyy)*gSrc[cn*(in_width)]
                                + (xx - xxyy)*gSrc[cn*1]
                                + (1 - xx - yy + xxyy)*gSrc[0]);

            data[cn*idx + 2] = (uchar)(xxyy*rSrc[cn*(1 + in_width)]
                                + (yy - xxyy)*rSrc[cn*(in_width)]
                                + (xx - xxyy)*rSrc[cn*1]
                                + (1 - xx - yy + xxyy)*rSrc[0]);

        }
    }
}

//------------------------------------------------------------------------------
// Distort Gray image
// input image size of out_width x out_height
// Output image, result, size of in_width x in_height
void UndistorterPTAM::distortGray(const cv::Mat& image_undist,
                                  cv::OutputArray result_dist) const
{
    
    int height_dist = in_height;
    int width_dist = in_width;
    int height_undist = out_height;
    int width_undist = out_width;
    
    if (!valid)
    {
        result_dist.getMatRef() = image_undist;
        return;
    }
    
    if (image_undist.rows != height_undist ||
        image_undist.cols != width_undist)
    {
        printf("UndistorterPTAM: input image size differs from"
               " expected input size! Not distorting.\n");
        result_dist.getMatRef() = image_undist;
        return;
    }
    
    if (height_dist == height_undist &&
        width_dist == width_undist &&
        inputCalibration[4] == 0)
    {
        // No transformation if neither distortion nor resize
        result_dist.getMatRef() = image_undist;
        return;
    }
    
    result_dist.create(height_dist, width_dist, CV_8UC1); // Gray
    cv::Mat resultMat = result_dist.getMatRef();
    assert(result_dist.getMatRef().isContinuous());
    assert(image_undist.isContinuous());
    
    uchar* data = resultMat.data;
    int cn = image_undist.channels();
    
    // For pixels of output distorted image
    for (int x_dist=0;x_dist<width_dist;x_dist++)
    {
        
        for (int y_dist=0;y_dist<height_dist;y_dist++)
        {
        
            // Get undistorted coordiante
            float x_undist, y_undist;
            undistortPoint(x_dist, y_dist, &x_undist, &y_undist);
            
            // Index for distorted coordiante
            int idx = y_dist*width_dist + x_dist;
            
            // Interpolation
            if (x_undist < 0 || x_undist >= width_undist ||
                y_undist < 0 || y_undist >= height_undist)
            {
                
                // Out of boundary
                data[idx] = 0;
                
            }
            else
            {
                
                // Get interget and rational parts
                int xxi = (int)x_undist;
                int yyi = (int)y_undist;
                x_undist -= xxi;
                y_undist -= yyi;
                float xxyy = x_undist*y_undist;
                
                // Get array base pointer
                const uchar *src = (uchar*)image_undist.data + xxi*cn
                                + yyi*cn*width_undist;
                
                // Interpolate
                data[idx] = (uchar)(xxyy*src[cn*(1 + width_undist)]
                            + (y_undist - xxyy)*src[cn*(width_undist)]
                            + (x_undist - xxyy)*src[cn*1]
                            + (1 - x_undist - y_undist + xxyy)*src[0]);
                
            }
            
        
        } // End for y_dist
        
    } // End for x_dist

}


//------------------------------------------------------------------------------
// Distort RGB image
// input image size of out_width x out_height
// Output image, result, size of in_width x in_height
void UndistorterPTAM::distortRGB(const cv::Mat& image_undist,
                                   cv::OutputArray result_dist) const
{
    
    int height_dist = in_height;
    int width_dist = in_width;
    int height_undist = out_height;
    int width_undist = out_width;
    
    if (!valid)
    {
        result_dist.getMatRef() = image_undist;
        return;
    }
    
    if (image_undist.rows != height_undist ||
        image_undist.cols != width_undist)
    {
        printf("UndistorterPTAM: input image size differs from"
               " expected input size! Not distorting.\n");
        result_dist.getMatRef() = image_undist;
        return;
    }
    
    if (height_dist == height_undist &&
        width_dist == width_undist &&
        inputCalibration[4] == 0)
    {
        // No transformation if neither distortion nor resize
        result_dist.getMatRef() = image_undist;
        return;
    }
    
    result_dist.create(height_dist, width_dist, CV_8UC3); // Gray
    cv::Mat resultMat = result_dist.getMatRef();
    assert(result_dist.getMatRef().isContinuous());
    assert(image_undist.isContinuous());
    
    uchar* data = resultMat.data;
    int cn = image_undist.channels();
    
    // For pixels of output distorted image
    for (int x_dist=0;x_dist<width_dist;x_dist++)
    {
        
        for (int y_dist=0;y_dist<height_dist;y_dist++)
        {
            
            // Get undistorted coordiante
            float x_undist, y_undist;
            undistortPoint(x_dist, y_dist, &x_undist, &y_undist);
            
            // Index for distorted coordiante
            int idx = y_dist*width_dist + x_dist;
            
            // Interpolation
            if (x_undist < 0 || x_undist >= width_undist ||
                y_undist < 0 || y_undist >= height_undist)
            {
                
                // Out of boundary
                data[cn*idx + 0] = 0;
                data[cn*idx + 1] = 0;
                data[cn*idx + 2] = 0;
                
            }
            else
            {
                
                // Get interget and rational parts
                int xxi = (int)x_undist;
                int yyi = (int)y_undist;
                x_undist -= xxi;
                y_undist -= yyi;
                float xxyy = x_undist*y_undist;
                
                // Get array base pointer
                const uchar *bSrc = (uchar*)image_undist.data + xxi*cn
                                    + yyi*cn*width_undist + 0;
                const uchar *gSrc = (uchar*)image_undist.data + xxi*cn
                                    + yyi*cn*width_undist + 1;
                const uchar *rSrc = (uchar*)image_undist.data + xxi*cn
                                    + yyi*cn*width_undist + 2;
                
                // Interpolate
                data[cn*idx + 0] = (uchar)(xxyy*bSrc[cn*(1 + width_undist)]
                                    + (y_undist - xxyy)*bSrc[cn*(width_undist)]
                                    + (x_undist - xxyy)*bSrc[cn*1]
                                    + (1 - x_undist - y_undist + xxyy)*bSrc[0]);
                data[cn*idx + 1] = (uchar)(xxyy*gSrc[cn*(1 + width_undist)]
                                    + (y_undist - xxyy)*gSrc[cn*(width_undist)]
                                    + (x_undist - xxyy)*gSrc[cn*1]
                                    + (1 - x_undist - y_undist + xxyy)*gSrc[0]);
                data[cn*idx + 2] = (uchar)(xxyy*rSrc[cn*(1 + width_undist)]
                                    + (y_undist - xxyy)*rSrc[cn*(width_undist)]
                                    + (x_undist - xxyy)*rSrc[cn*1]
                                    + (1 - x_undist - y_undist + xxyy)*rSrc[0]);
            
            }
            
            
        } // End for y_dist
        
    } // End for x_dist
    
}
//// TODO: not working yet
//{
//    
//    int height_dist = in_height;
//    int width_dist = in_width;
//    int height_undist = out_height;
//    int width_undist = out_width;
//
//    if (!valid)
//    {
//        result_undist.getMatRef() = image_dist;
//        return;
//    }
//    
//    if (image_dist.rows != height_undist || image_dist.cols != width_undist)
//    {
//        printf("UndistorterPTAM: input image size differs from"
//               " expected input size! Not distorting.\n");
//        result_undist.getMatRef() = image_dist;
//        return;
//    }
//    
//    if (height_dist == height_undist &&
//        width_dist == width_undist &&
//        inputCalibration[4] == 0)
//    {
//        // No transformation if neither distortion nor resize
//        result_undist.getMatRef() = image_dist;
//        return;
//    }
//
//    result_undist.create(height_dist, width_dist, CV_8UC3);
//    cv::Mat resultMat = result_undist.getMatRef();
//    assert(result_undist.getMatRef().isContinuous());
//    assert(image_dist.isContinuous());
//    
//    uchar* data = resultMat.data;
//    int cn = image_dist.channels();
//   
//    // For pixels of output undistorted image
//    for (int x_undist=0;x_undist<width_undist;x_undist++)
//    {
//        for (int y_undist=0;y_undist<height_undist;y_undist++)
//        {
//            
//            // TO DO -- not working yet
//            fprintf(stderr,"TO DO, NOT WORKING YET\n");
//            exit(1);
//            
////            float x_dist, y_dist;
////            distortPoint(x_undist, y_undist, &x_dist, &y_dist);
////            
////            // Check the boundary
////            if ((x_dist >= 0) && (x_dist <= width_dist - 1) &&
////                (y_dist >= 0) && (y_dist <= height_dist - 1))
////            {
////
////                // get integer and rational parts
////                int xxi = x_dist;
////                int yyi = y_dist;
////                x_dist -= xxi;
////                y_dist -= yyi;
////                float xxyy = x_undist*y_undist;
////                
////                // Get array base pointer
////                // const uchar* src = (uchar*)image.data + xxi + yyi * in_width;
////                uchar* bSrc = (uchar*)image_dist.data + xxi*cn + yyi*cn*width_dist + 0;
////                uchar* gSrc = (uchar*)image_dist.data + xxi*cn + yyi*cn*width_dist + 1;
////                uchar* rSrc = (uchar*)image_dist.data + xxi*cn + yyi*cn*width_dist + 2;
////                
////                int idx = xxi*yyi;
////                
////                // interpolate (bilinear)
////                //data[idx] =  xxyy * src[1+in_width]
////                //                    + (yy-xxyy) * src[in_width]
////                //                    + (xx-xxyy) * src[1]
////                //                    + (1-xx-yy+xxyy) * src[0];
////                data[cn*idx + 0] = xxyy*bSrc[cn*(1 + width_dist)]
////                + (yy - xxyy)*bSrc[cn*(width_dist)]
////                + (xx - xxyy)*bSrc[cn*1]
////                + (1 - xx - yy + xxyy)*bSrc[0];
////                
////                data[cn*idx + 1] = xxyy*gSrc[cn*(1 + width_dist)]
////                + (yy - xxyy)*gSrc[cn*(width_dist)]
////                + (xx - xxyy)*gSrc[cn*1]
////                + (1 - xx - yy + xxyy)*gSrc[0];
////                
////                data[cn*idx + 2] = xxyy*rSrc[cn*(1 + width_dist)]
////                + (yy - xxyy)*rSrc[cn*(width_dist)]
////                + (xx - xxyy)*rSrc[cn*1]
////                + (1 - xx - yy + xxyy)*rSrc[0];
////            
////            }
//            
//        }
//    }
//    
//    
//}

const cv::Mat& UndistorterPTAM::getK() const
{
	return K_;
}

const cv::Mat& UndistorterPTAM::getOriginalK() const
{
	return originalK_;
}

int UndistorterPTAM::getOutputWidth() const
{
	return out_width;
}

int UndistorterPTAM::getOutputHeight() const
{
	return out_height;
}
int UndistorterPTAM::getInputWidth() const
{
	return in_width;
}

int UndistorterPTAM::getInputHeight() const
{
	return in_height;
}
    
float UndistorterPTAM::getInputDistort() const
{
    return inputCalibration[4];
}

float UndistorterPTAM::getOutputDistort() const
{
    return outputCalibration[4];
}

bool UndistorterPTAM::isValid() const
{
	return valid;
}
    
    
bool UndistorterPTAM::isNoRectification() const
{
    
    return isNoRect_;
    
}

//void UndistorterPTAM::distortPoint(float px, float py,
//                          float *px_out, float *py_out) const
//{
//
//    // Use remapX and remapY, which contain mappings from
//    // undistorted pixel points to distorted pixel points
//
//    if (!valid)
//    {
//
//        // No distortion
//        *px_out = px;
//        *py_out = py;
//        return;
//
//    }
//
//    if (px < 0 || px >= out_width - 1 || py < 0 || py >= out_height - 1)
//    {
//
//        // points are not in the image
//        //printf("Errror: points are not inside the image\n");
//        *px_out = -1;
//        *py_out = -1;
//        return;
//
//    }
//    
//    // Bilinear interpolation for x coordinate
//    int x1 = floor(px);
//    int y1 = floor(py);
//    int x2 = x1 + 1.0;
//    int y2 = y1 + 1.0;
//    
//    float X11 = remapX[x1 + y1*out_width];
//    float X12 = remapX[x1 + y2*out_width];
//    float X21 = remapX[x2 + y1*out_width];
//    float X22 = remapX[x2 + y2*out_width];
//    
//    float Xxy = 1.0f/((x2 - x1)*(y2 - y1))*
//                (X11*(x2 - px)*(y2 - py) +
//                 X21*(px - x1)*(y2 - py) +
//                 X12*(x2 - px)*(py - y1) +
//                 X22*(px - x1)*(py - y1));
//    
//    float Y11 = remapY[x1 + y1*out_width];
//    float Y12 = remapY[x1 + y2*out_width];
//    float Y21 = remapY[x2 + y1*out_width];
//    float Y22 = remapY[x2 + y2*out_width];
//    
//    float Yxy = 1.0f/((x2 - x1)*(y2 - y1))*
//                (Y11*(x2 - px)*(y2 - py) +
//                 Y21*(px - x1)*(y2 - py) +
//                 Y12*(x2 - px)*(py - y1) +
//                 Y22*(px - x1)*(py - y1));
//    
//    // Output
//    *px_out = Xxy;
//    *py_out = Yxy;
//
//}
    
void UndistorterPTAM::distortPoint(float x_undist, float y_undist,
                                   float *x_dist_out, float *y_dist_out) const
{
    
    if (!valid)
    {
        
        // No distortion
        *x_dist_out = x_undist;
        *y_dist_out = y_undist;
        return;
        
    }
    
    if (x_undist < 0 || x_undist > out_width - 1 || y_undist < 0 || y_undist > out_height - 1)
    {
        
        // points are not in the image
        //printf("Errror: points are not inside the image\n");
        *x_dist_out = -1;
        *y_dist_out = -1;
        return;
        
    }
    
    cv::Mat K_in = getOriginalK();
    cv::Mat K_out = getK();
    
    // Camera parameters
    // Distorted
    float fx_dist = (float)K_in.at<double>(0, 0);
    float fy_dist = (float)K_in.at<double>(1, 1);
    float cx_dist = (float)K_in.at<double>(0, 2);
    float cy_dist = (float)K_in.at<double>(1, 2);
    float fx_undist = (float)K_out.at<double>(0, 0); // Undistorted
    float fy_undist = (float)K_out.at<double>(1, 1);
    float cx_undist = (float)K_out.at<double>(0, 2);
    float cy_undist = (float)K_out.at<double>(1, 2);
    float w = (float)getInputDistort();
    float d2t = 2.0f*tanf(w/2.0f);
    
    // Normalized distorted points
    float x_undist_norm = (x_undist - cx_undist)/fx_undist;
    float y_undist_norm = (y_undist - cy_undist)/fy_undist;
    
    // Radius of the normalized undistort pixels
    float ru = sqrtf(x_undist_norm*x_undist_norm + y_undist_norm*y_undist_norm);
    
    // Computer factor ratio rd/ru to convert from undistort to distort
    float fac = 1.0f;
    if ((ru == 0) || (w == 0))
    {
        
        fac = 1.0f;
        
    }
    else
    {
        
        fac = atanf(ru*d2t)/(w*ru); // Ratio rd/ru
        
    }
    
    // Distorted pixels (input)
    float x_dist = fx_dist*fac*x_undist_norm + cx_dist;
    float y_dist = fy_dist*fac*y_undist_norm + cy_dist;
    
    
    // Output
    *x_dist_out = x_dist;
    *y_dist_out = y_dist;
    
}

void UndistorterPTAM::distortPointLevel(float x_undist, float y_undist,
                                   float *x_dist, float *y_dist,
                                   int level = 0) const
{
    
    // Image scale by the level
    float scale = powf(2, level);
    
    // Input to original scale
    float x_undist_org = (x_undist + 0.5f)*scale - 0.5f;
    float y_undist_org = (y_undist + 0.5f)*scale - 0.5f;
    
    float x_dist_org, y_dist_org;
    distortPoint(x_undist_org, y_undist_org,
                              &x_dist_org, &y_dist_org);
    // Output to level size
    *x_dist = (x_dist_org + 0.5f)/scale - 0.5f;
    *y_dist = (y_dist_org + 0.5f)/scale - 0.5f;
    
}

//void UndistorterPTAM::distortPointLevel(float px, float py,
//                                   float *px_out, float *py_out,
//                                   int level = 0) const
//{
//    
//    // Use remapX and remapY, which contain mappings from
//    // undistorted pixel points to distorted pixel points
//    
//    // px and py are pixels in undistorted image at the level
//    // px_out and py_out are output pixels in distorted image at the level
//    
//    // Change px and py to the original coords by image pyramid level
//    px = px * pow(2.0, level);
//    py = py * pow(2.0, level);
//    
//    if (!valid)
//    {
//        
//        // No distortion
//        *px_out = px;
//        *py_out = py;
//        return;
//        
//    }
//    
//    if (px < 0 || px >= out_width - 1 || py < 0 || py >= out_height - 1)
//    {
//        
//        // points are not in the image
//        //printf("Errror: points are not inside the image\n");
//        *px_out = -1;
//        *py_out = -1;
//        return;
//        
//    }
//    
//    // Bilinear interpolation for x coordinate
//    int x1 = floor(px);
//    int y1 = floor(py);
//    int x2 = x1 + 1.0;
//    int y2 = y1 + 1.0;
//    
//    float X11 = remapX[x1 + y1*out_width];
//    float X12 = remapX[x1 + y2*out_width];
//    float X21 = remapX[x2 + y1*out_width];
//    float X22 = remapX[x2 + y2*out_width];
//    
//    float Xxy = 1.0f/((x2 - x1)*(y2 - y1))*
//                (X11*(x2 - px)*(y2 - py) +
//                 X21*(px - x1)*(y2 - py) +
//                 X12*(x2 - px)*(py - y1) +
//                 X22*(px - x1)*(py - y1));
//    
//    float Y11 = remapY[x1 + y1*out_width];
//    float Y12 = remapY[x1 + y2*out_width];
//    float Y21 = remapY[x2 + y1*out_width];
//    float Y22 = remapY[x2 + y2*out_width];
//    
//    float Yxy = 1.0f/((x2 - x1)*(y2 - y1))*
//                (Y11*(x2 - px)*(y2 - py) +
//                 Y21*(px - x1)*(y2 - py) +
//                 Y12*(x2 - px)*(py - y1) +
//                 Y22*(px - x1)*(py - y1));
//    
//    // Change Xxy and Yxy by the pyramid level with the same undistorted size
//    Xxy = (Xxy*out_width)/(in_width*pow(2.0, level));
//    Yxy = (Yxy*out_height)/(in_height*pow(2.0, level));
//    
//    // Output
//    *px_out = Xxy;
//    *py_out = Yxy;
//    
//}

UndistorterOpenCV::UndistorterOpenCV(const char* configFileName)
{
	valid = true;
	
	// read parameters
	std::ifstream infile(configFileName);
	assert(infile.good());

	std::string l1, l2, l3, l4;

	std::getline(infile,l1);
	std::getline(infile,l2);
	std::getline(infile,l3);
	std::getline(infile,l4);

	// l1 & l2
	if(std::sscanf(l1.c_str(), "%f %f %f %f %f %f %f %f",
        &inputCalibration[0], &inputCalibration[1], &inputCalibration[2],
                   &inputCalibration[3], &inputCalibration[4],
		&inputCalibration[5], &inputCalibration[6], &inputCalibration[7]
  				) == 8 &&
			std::sscanf(l2.c_str(), "%d %d", &in_width, &in_height) == 2)
	{
		printf("Input resolution: %d %d\n",in_width, in_height);
		printf("In: %f %f %f %f %f %f %f %f\n",
                inputCalibration[0], inputCalibration[1], inputCalibration[2],
                inputCalibration[3], inputCalibration[4],
				inputCalibration[5], inputCalibration[6], inputCalibration[7]);
	}
	else
	{
        printf("Failed to read camera \
               calibration (invalid format?)\nCalibration file: %s\n",
               configFileName);
		valid = false;
	}

	// l3
	if(l3 == "crop")
	{
		outputCalibration = -1;
		printf("Out: Crop\n");
	}
	else if(l3 == "full")
	{
		outputCalibration = -2;
		printf("Out: Full\n");
	}
	else if(l3 == "none")
	{
		printf("NO RECTIFICATION\n");
		valid = true;
        isNoRect_ = true;
	}
	else
	{
		printf("Out: Failed to Read Output pars... not rectifying.\n");
		valid = false;
	}

	// l4
	if(std::sscanf(l4.c_str(), "%d %d", &out_width, &out_height) == 2)
	{
		printf("Output resolution: %d %d\n", out_width, out_height);
	}
	else
	{
		printf("Out: Failed to Read Output resolution... not rectifying.\n");
		valid = false;
	}
	
	cv::Mat distCoeffs = cv::Mat::zeros(4, 1, CV_32F);
	for (int i = 0; i < 4; ++ i)
		distCoeffs.at<float>(i, 0) = inputCalibration[4 + i];

	originalK_ = cv::Mat(3, 3, CV_64F, cv::Scalar(0));
	originalK_.at<double>(0, 0) = inputCalibration[0] * in_width;
	originalK_.at<double>(1, 1) = inputCalibration[1] * in_height;
	originalK_.at<double>(2, 2) = 1;
	originalK_.at<double>(0, 2) = inputCalibration[2] * in_width;
	originalK_.at<double>(1, 2) = inputCalibration[3] * in_height;

	if (valid)
	{
		K_ = cv::getOptimalNewCameraMatrix(originalK_, distCoeffs, cv::Size(in_width, in_height), (outputCalibration == -2) ? 1 : 0, cv::Size(out_width, out_height), nullptr, false);
		
		cv::initUndistortRectifyMap(originalK_, distCoeffs, cv::Mat(), K_,
				cv::Size(out_width, out_height), CV_16SC2, map1, map2);
		
		originalK_.at<double>(0, 0) /= in_width;
		originalK_.at<double>(0, 2) /= in_width;
		originalK_.at<double>(1, 1) /= in_height;
		originalK_.at<double>(1, 2) /= in_height;
	}
	
	originalK_ = originalK_.t();
	K_ = K_.t();
}

UndistorterOpenCV::~UndistorterOpenCV()
{
}

void UndistorterOpenCV::undistort(const cv::Mat& image, cv::OutputArray result) const
{
	cv::remap(image, result, map1, map2, cv::INTER_LINEAR);
}

const cv::Mat& UndistorterOpenCV::getK() const
{
	return K_;
}

const cv::Mat& UndistorterOpenCV::getOriginalK() const
{
	return originalK_;
}

int UndistorterOpenCV::getOutputWidth() const
{
	return out_width;
}

int UndistorterOpenCV::getOutputHeight() const
{
	return out_height;
}
int UndistorterOpenCV::getInputWidth() const
{
	return in_width;
}

int UndistorterOpenCV::getInputHeight() const
{
	return in_height;
}

float UndistorterOpenCV::getInputDistort() const
{
    return inputCalibration[4];
}

float UndistorterOpenCV::getOutputDistort() const
{
    // This is not working. Do not use this method.
    // This is an experimental feature.
    fprintf(stderr,"Please use PTAM format for this method.\n");
    exit(1);
    return inputCalibration[4];
}
    
    bool UndistorterOpenCV::isValid() const
{
	return valid;
}

    
bool UndistorterOpenCV::isNoRectification() const
{
    
    return isNoRect_;
    
}

void UndistorterOpenCV::distortPoint(float px, float py,
                          float *px_out, float *py_out) const
{

    // Use map1 and map2 which contain mappings from
    // undistorted pixel points to distorted pixel points

    if (!valid)
    {

        // No distortion
        *px_out = px;
        *py_out = py;
        return;

    }

    if (px < 0 || px >= out_width - 1 || py < 0 || py >= out_height - 1)
    {

        // points are not in the image
        //printf("Errror: points are not inside the image\n");
        *px_out = -1;
        *py_out = -1;
        return;

    }

    // Bilinear interpolation for x coordinate
    int x1 = (int)floor(px);
    int y1 = (int)floor(py);
    int x2 = (int)(x1 + 1.0);
    int y2 = (int)(y1 + 1.0);
    
    float X11 = (float)map1.at<double>(x1, y1);
    float X12 = (float)map1.at<double>(x1, y2);
    float X21 = (float)map1.at<double>(x2, y1);
    float X22 = (float)map1.at<double>(x2, y2);
    
    float Xxy = 1.0f/((x2 - x1)*(y2 - y1))*
    (X11*(x2 - px)*(y2 - py) +
     X21*(px - x1)*(y2 - py) +
     X12*(x2 - px)*(py - y1) +
     X22*(px - x1)*(py - y1));
    
    float Y11 = (float)map2.at<double>(x1, y1);
    float Y12 = (float)map2.at<double>(x1, y2);
    float Y21 = (float)map2.at<double>(x2, y1);
    float Y22 = (float)map2.at<double>(x2, y2);
    
    float Yxy = 1.0f/((x2 - x1)*(y2 - y1))*
    (Y11*(x2 - px)*(y2 - py) +
     Y21*(px - x1)*(y2 - py) +
     Y12*(x2 - px)*(py - y1) +
     Y22*(px - x1)*(py - y1));
    
    // Output
    *px_out = Xxy;
    *py_out = Yxy;

}

void UndistorterOpenCV::distortPointLevel(float px, float py,
                                     float *px_out, float *py_out,
                                     int level) const
{
    
    // Use map1 and map2 which contain mappings from
    // undistorted pixel points to distorted pixel points
    
    // px and py are pixels in undistorted image at the level
    // px_out and py_out are output pixels in distorted image at the level
    
    // Change px and py to the original coords by image pyramid level
    px = px * powf(2.0, level);
    py = py * powf(2.0, level);
    
    if (!valid)
    {
        
        // No distortion
        *px_out = px;
        *py_out = py;
        return;
        
    }
    
    if (px < 0 || px > out_width - 1 || py < 0 || py > out_height - 1)
    {
        
        // points are not in the image
        //printf("Errror: points are not inside the image\n");
        *px_out = -1;
        *py_out = -1;
        return;
        
    }
    
    // Bilinear interpolation for x coordinate
    int x1 = (int)floor(px);
    int y1 = (int)floor(py);
    int x2 = (int)(x1 + 1.0);
    int y2 = (int)(y1 + 1.0);
    
    float X11 = (float)map1.at<double>(x1, y1);
    float X12 = (float)map1.at<double>(x1, y2);
    float X21 = (float)map1.at<double>(x2, y1);
    float X22 = (float)map1.at<double>(x2, y2);
    
    float Xxy = 1.0f/((x2 - x1)*(y2 - y1))*
                (X11*(x2 - px)*(y2 - py) +
                 X21*(px - x1)*(y2 - py) +
                 X12*(x2 - px)*(py - y1) +
                 X22*(px - x1)*(py - y1));
    
    float Y11 = (float)map2.at<double>(x1, y1);
    float Y12 = (float)map2.at<double>(x1, y2);
    float Y21 = (float)map2.at<double>(x2, y1);
    float Y22 = (float)map2.at<double>(x2, y2);
    
    float Yxy = 1.0f/((x2 - x1)*(y2 - y1))*
                (Y11*(x2 - px)*(y2 - py) +
                 Y21*(px - x1)*(y2 - py) +
                 Y12*(x2 - px)*(py - y1) +
                 Y22*(px - x1)*(py - y1));
    
    // Change Xxy and Yxy by the pyramid level with the same undistorted size
    Xxy = (Xxy*out_width)/(in_width*powf(2.0, level));
    Yxy = (Yxy*out_height)/(in_height*powf(2.0, level));
    
    // Output
    *px_out = Xxy;
    *py_out = Yxy;
    
}

void UndistorterOpenCV::undistortRGB(const cv::Mat& image,
                                     cv::OutputArray result) const
// TODO: Test required
{
    cv::remap(image, result, map1, map2, cv::INTER_LINEAR);
}
    
//------------------------------------------------------------------------------
// Distort Gray image
// input image size of out_width x out_height
// Output image, result, size of in_width x in_height
void UndistorterOpenCV::distortGray(const cv::Mat& image_undist,
                                    cv::OutputArray result_dist) const
{
    
    // TODO: Implementation
    
}
    
//------------------------------------------------------------------------------
// Distort RGB image
// input image size of out_width x out_height
// Output image, result, size of in_width x in_height
void UndistorterOpenCV::distortRGB(const cv::Mat& image_dist,
                                 cv::OutputArray result_undist) const
{
    
    // TODO: Implementation
    
}
    
void UndistorterOpenCV::undistortPoint(float x_dist, float y_dist,
                                       float *x_undist, float *y_undist) const
{
    
    // TODO: Implementation required
    
}


}

// end of namespace
