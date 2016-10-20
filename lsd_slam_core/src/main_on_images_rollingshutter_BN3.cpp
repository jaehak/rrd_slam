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
// Rolling Shutter Direct SLAM (BN3 version)
//
// B-Spline
// Neighbour window
// Multiple keyframes for the window
// Rolling shutter keyframe (all row poses in a keyframe)
//
// Author: Jae-Hak Kim (jaehak.kim@adelaide.edu.au)
// Copyright@2015, University of Adelaide
//==============================================================================


#include <boost/thread.hpp>
#include "util/settings.h"
#include "util/globalFuncs.h"
#include "SlamSystemForRollingShutter.h"

#include <sstream>
#include <fstream>
#include <dirent.h>
#include <algorithm>

#include "IOWrapper/ROS/ROSOutput3DWrapper.h"
#include "IOWrapper/ROS/rosReconfigure.h"

#include "util/Undistorter.h"
#include <ros/package.h>

#include "opencv2/opencv.hpp"
#include "Tracking/WarpResidual.h"

#include <glog/logging.h>

#include <thread>
#include <X11/Xlib.h>

std::string &ltrim(std::string &s) {
    s.erase(s.begin(),
        std::find_if(s.begin(),
            s.end(),
                std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}
std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(),
        s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(),
            s.end());
    return s;
}
std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}
int getdir (std::string dir, std::vector<std::string> &files)
{
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL)
    {
        return -1;
    }

    while ((dirp = readdir(dp)) != NULL) {
        std::string name = std::string(dirp->d_name);

        if(name != "." && name != "..")
            files.push_back(name);
    }
    closedir(dp);


    std::sort(files.begin(), files.end());

    if(dir.at( dir.length() - 1 ) != '/') dir = dir+"/";
    for(unsigned int i=0;i<files.size();i++)
    {
        if(files[i].at(0) != '/')
            files[i] = dir + files[i];
    }

    return files.size();
}

int getFile (std::string source, std::vector<std::string> &files)
{
    std::ifstream f(source.c_str());

    if(f.good() && f.is_open())
    {
        while(!f.eof())
        {
            std::string l;
            std::getline(f,l);

            l = trim(l);

            if(l == "" || l[0] == '#')
                continue;

            files.push_back(l);
        }

        f.close();

        size_t sp = source.find_last_of('/');
        std::string prefix;
        if(sp == std::string::npos)
            prefix = "";
        else
            prefix = source.substr(0,sp);

        for(unsigned int i=0;i<files.size();i++)
        {
            if(files[i].at(0) != '/')
                files[i] = prefix + "/" + files[i];
        }

        return (int)files.size();
    }
    else
    {
        f.close();
        return -1;
    }

}

void useGroundTruthDepth(float *depthInput,
                         std::vector<std::string> depthFiles,
                         int runningIDX,
                         int width_output,
                         float fx,
                         float fy,
                         float cx,
                         float cy,
                         int depthFormat,
                         float noiseOnInitGT)
{

    // Use Ground truth depth init
    printf("---------------------------------------------\n");
    printf("Reading %s\n", depthFiles[runningIDX].c_str());

    FILE* fid = fopen(depthFiles[runningIDX].c_str(),"r");
    float *ptrDepth = depthInput;
    float depthVal;
    int u = 0;
    int v = 0;
    while (fscanf(fid, "%f", &depthVal)!=EOF)
    {
        
        if (u == width_output) // u for column index
        {
            v++; // v for row index
            u = 0;
        }
        
        // Set depth value from file via calibration information
        float u_u0_by_fx = ((float)u - cx)/fx;
        float v_v0_by_fy = ((float)v - cy)/fy;
        float u_u0_by_fx_sq = u_u0_by_fx*u_u0_by_fx;
        float v_v0_by_fy_sq = v_v0_by_fy*v_v0_by_fy;
        
        // Refer to compute3Dpositions.m from ICL-NUIM dataset
        float zVal;
        if (depthFormat == 1) // (JHK Synth depth format)
        {
            
            if (u==0 && v==0)
            {
                fprintf(stderr, "ICL-NUIM RGBD depth format:"
                        "(Eucliden distance depth)\n");
                printf("ICL-NUIM RGBD depth format:"
                       "(Eucliden distance depth)\n");
                
            }
            
            zVal = depthVal /
            sqrtf(u_u0_by_fx_sq + v_v0_by_fy_sq + 1.0f);
            
        }
        else  // depthFormat == 0 (Freiburg depth format)
        {
            
            if (u==0 && v==0)
            {
                fprintf(stderr,"TUM RGBD depth format:"
                        "(z-direction depth)\n");
                printf("TUM RGBD depth format:"
                       "(z-direction depth)\n");
            }
            
            zVal = depthVal;
            
        }
        
        // Add a uniform noise on ground truth depth
        if (noiseOnInitGT != 0.0)
        {
            
            zVal = zVal + noiseOnInitGT *
            ((rand() % 100001) / 100000.0f);
            
        }
        
    #if 0 // SCALE DOWN DEPTH
        zVal = zVal * 0.1f;
    #endif
        
        *ptrDepth = zVal;
        ptrDepth++;
        
        u++;
        
    }
    
    fclose(fid);
    
}

float bilinear_interpolation(float x,
                             float y,
                             float Q11,
                             float Q12,
                             float Q21,
                             float Q22,
                             int x1,
                             int y1,
                             int x2,
                             int y2)
{
// Suppose that we want to find the value of the unknown function f
// at the point (x, y). It is assumed that we know the value of f
// at the four points Q11 = (x1, y1), Q12 = (x1, y2), Q21 = (x2, y1),
// and Q22 = (x2, y2).
    
    float Qxy = 1.0f/((x2 - x1)*(y2 - y1))*
                (Q11*(x2 - x)*(y2 - y) +
                 Q21*(x - x1)*(y2 - y) +
                 Q12*(x2 - x)*(y - y1) +
                 Q22*(x - x1)*(y - y1));
    
    return Qxy;
    
}

// Use ground truth for radial images.
// As our depth files are at integer pixel coordinates
// it requires bilinear interpolation of depths
void useGroundTruthDepthRadial(lsd_slam::Undistorter* undistorter,
                               float *depthInput,
                         std::vector<std::string> depthFiles,
                         int runningIDX,
                         int width, int height,
                         float fx,
                         float fy,
                         float cx,
                         float cy,
                         int depthFormat,
                         float noiseOnInitGT)
{
    
    // Use Ground truth depth init
    printf("---------------------------------------------\n");
    printf("Reading %s\n", depthFiles[runningIDX].c_str());
    
    FILE* fid = fopen(depthFiles[runningIDX].c_str(),"r");
    float *ptrDepth = depthInput;
    float depthVal;
    int u = 0;
    int v = 0;
    while (fscanf(fid, "%f", &depthVal)!=EOF)
    {
        
        if (u == width) // u for column index
        {
            v++; // v for row index
            u = 0;
        }

        *ptrDepth = depthVal;
        ptrDepth++;

        u++;
        
    }
    
    // Bilinear interpolation on depths
    float *depthInterpolated = (float*)malloc(sizeof(float)*width*height);
    ptrDepth = depthInterpolated;
    for (int v=0;v<height;v++)
    {
    
        for (int u=0;u<width;u++)
        {
            
            // Distorted pixel
            float x_distorted, y_distorted;
            undistorter->distortPoint(u, v, &x_distorted, &y_distorted);
            if ((x_distorted == -1) || (y_distorted == -1))
            {
                continue;
            }
            if ((x_distorted < 0) || (x_distorted > width - 1) ||
                (y_distorted < 0) || (y_distorted > height - 1))
            {
                continue;
            }
            
            // Four neighbour pixels
            int x1 = (int)floor(x_distorted);
            int y1 = (int)floor(y_distorted);
            int x2 = x1 + 1;
            int y2 = y1 + 1;
            
            // Four neighbour depths
            float Q11 = depthInput[x1 + y1*width];
            float Q12 = depthInput[x1 + y2*width];
            float Q21 = depthInput[x2 + y1*width];
            float Q22 = depthInput[x2 + y2*width];
            //printf("i, j = %d %d => x_d, y_d %f %f\n", i, j, x_distorted, y_distorted);
            //printf("  x1, y1, x2, y2: %d %d %d %d\n", x1, y1, x2, y2);
            
            
            // Bilinear interpolation from neigbhour depths
            float Qxy = bilinear_interpolation(x_distorted, y_distorted,
                                         Q11, Q12, Q21, Q22,
                                         x1, y1, x2, y2);
            
            float depthVal = Qxy;

            // Set depth value from file via calibration information
            float u_u0_by_fx = ((float)u - cx)/fx;
            float v_v0_by_fy = ((float)v - cy)/fy;
            float u_u0_by_fx_sq = u_u0_by_fx*u_u0_by_fx;
            float v_v0_by_fy_sq = v_v0_by_fy*v_v0_by_fy;
            
            // Refer to compute3Dpositions.m from ICL-NUIM dataset
            float zVal;
            if (depthFormat == 1) // (JHK Synth depth format)
            {
                
                if (u==0 && v==0)
                {
                    fprintf(stderr, "ICL-NUIM RGBD depth format:"
                            "(Eucliden distance depth)\n");
                    printf("ICL-NUIM RGBD depth format:"
                           "(Eucliden distance depth)\n");
                    
                }
                
                zVal = depthVal /
                sqrtf(u_u0_by_fx_sq + v_v0_by_fy_sq + 1.0f);
                
            }
            else  // depthFormat == 0 (Freiburg depth format)
            {
                
                if (u==0 && v==0)
                {
                    fprintf(stderr,"TUM RGBD depth format:"
                            "(z-direction depth)\n");
                    printf("TUM RGBD depth format:"
                           "(z-direction depth)\n");
                }
                
                zVal = depthVal;
                
            }
            
            // Add a uniform noise on ground truth depth
            if (noiseOnInitGT != 0.0) // noiseOnInitGT (meter)
            {
                
                zVal = zVal + noiseOnInitGT *
                ((rand() % 100001) / 100000.0f); // noiseOnInitGT * unif(0, 1)
                
            }
            
            depthInterpolated[u + v*width] = zVal;
            
        }
        
    }
    
    // Copy back
    for(int i=0;i<width*height;i++)
        depthInput[i] = depthInterpolated[i];
    
    // Delete
    free(depthInterpolated);
    fclose(fid);
    
}

int main( int argc, char** argv )
{

    for (int i=0; i<argc; i++)
    {
        printf("argv[%d] = %s\n", i, argv[i]);
    }

    google::InitGoogleLogging(argv[0]);

    ros::init(argc, argv, "LSD_SLAM");

    dynamic_reconfigure::Server<lsd_slam_core::LSDParamsConfig>
            srv(ros::NodeHandle("~"));
    srv.setCallback(lsd_slam::dynConfCb);

    dynamic_reconfigure::Server<lsd_slam_core::LSDDebugParamsConfig>
            srvDebug(ros::NodeHandle("~Debug"));
    srvDebug.setCallback(lsd_slam::dynConfCbDebug);

    lsd_slam::packagePath = ros::package::getPath("lsd_slam_core")+"/";


    //--------------------------------------------------------------------------
    fprintf(stderr, "Press any key to start...\n");
    std::cin.ignore( std::numeric_limits <std::streamsize> ::max(), '\n' );
    fprintf(stderr,"Started\n");
    

    // get camera calibration in form of an undistorter object.
    // if no undistortion is required,
    // the undistorter will just pass images through.
    std::string calibFile;
    lsd_slam::Undistorter* undistorter = 0;
    if(ros::param::get("~calib", calibFile))
    {
         undistorter =
            lsd_slam::Undistorter::getUndistorterForFile(calibFile.c_str());
         ros::param::del("~calib");
    }

    if(undistorter == 0)
    {
        printf("need camera calibration file! (set using _calib:=FILE)\n");
        exit(0);
    }

    int width_output = undistorter->getOutputWidth();
    int height_output = undistorter->getOutputHeight();

    int width_input = undistorter->getInputWidth();
    int height_input = undistorter->getInputHeight();

    float fx = (float)undistorter->getK().at<double>(0, 0);
    float fy = (float)undistorter->getK().at<double>(1, 1);
    float cx = (float)undistorter->getK().at<double>(0, 2);
    float cy = (float)undistorter->getK().at<double>(1, 2);
    Sophus::Matrix3f K;
    K << fx, 0.0, cx, 0.0, fy, cy, 0.0, 0.0, 1.0;
    printf("K matrix = \n");
    std::cout << K << std::endl;
    
    
    //--------------------------------------------------------------------------
    // Option to use a simluation of distortion
    // This feature is to emulate distorted images input by distorting the
    // input undistorted images as if they are distroted input images.
    // Also, we distort undistorted depth as well by radial distortion parameter
    // The purpose of this feature is mainly for debugging.
    int emulDistort = 0;
    if (!ros::param::get("~emulDistort", emulDistort))
    {
        printf("Emulate the radial distortion on input undistorted images\n");
        emulDistort = 0;
    }
    ros::param::del("~emulDistort");
    printf("emulDistort = %d\n", emulDistort);
    fprintf(stderr, "emulDistort = %d\n", emulDistort);
    
    //--------------------------------------------------------------------------
    // make output wrapper. just set to zero if no output is required.
    lsd_slam::Output3DWrapper* outputWrapper =
            new lsd_slam::ROSOutput3DWrapper(width_output, height_output);

    //--------------------------------------------------------------------------
    // open image files: first try to open as file.
    std::string source;
    std::vector<std::string> imageFiles;
    if(!ros::param::get("~files", source))
    {
        printf("need source files! (set using _files:=FOLDER)\n");
        exit(0);
    }
    ros::param::del("~files");

    if(getdir(source, imageFiles) >= 0)
    {
        printf("found %d image files in folder %s!\n",
               (int)imageFiles.size(), source.c_str());
    }
    else if(getFile(source, imageFiles) >= 0)
    {
        printf("found %d image files in file %s!\n",
               (int)imageFiles.size(), source.c_str());
    }
    else
    {
        printf("could not load file list! wrong path / file?\n");
    }

    //--------------------------------------------------------------------------
    // open depth files
    std::string depthSource;
    std::vector<std::string> depthFiles;
    bool useGTdepth = false;
    if (!ros::param::get("~depthFiles", depthSource))
    {
        printf("No depth input files! (set using _depthFiles:=FOLDER)\n");
        //exit(0);
        useGTdepth = false;
    }
    ros::param::del("~depthFiles");

    if(getdir(depthSource, depthFiles) >= 0)
    {
        printf("found %d depth files in folder %s!\n",
               (int)depthFiles.size(), depthSource.c_str());
        useGTdepth = true;
    }
    else if(getFile(depthSource, depthFiles) >= 0)
    {
        printf("found %d depth files in file %s!\n",
               (int)depthFiles.size(), depthSource.c_str());
        useGTdepth = true;
    }
    else
    {
        printf("could not load file list! wrong path / file?\n");
        useGTdepth = false;
    }
    
    //--------------------------------------------------------------------------
    // open motion files
    std::string motionSource;
    std::vector<std::string> motionFiles;
    bool useGTmotion = false;
    if (!ros::param::get("~motionFiles", motionSource))
    {
        printf("No motion input files! (set using _motionFiles:=FOLDER)\n");
        //exit(0);
        useGTmotion = false;
    }
    ros::param::del("~motionFiles");
    
    if (getdir(motionSource, motionFiles) >= 0)
    {
        printf("Found %d motion files in folder %s!\n",
               (int)motionFiles.size(), motionSource.c_str());
        useGTmotion = true;
    }
    else if (getFile(motionSource, motionFiles) >= 0)
    {
        printf("found %d motion files in file %s!\n",
               (int)motionFiles.size(), motionSource.c_str());
        useGTmotion = true;
    }
    else
    {
        printf("could not load file list! wrong path / file?\n");
        useGTmotion = false;
    }
    
    //--------------------------------------------------------------------------
    // Option for image mask file
    
    // Initialize undistorted image mask with one
    cv::Mat imageMask = cv::Mat(height_output, width_output, CV_8U, 255);
    
    std::string imageMaskFile;
    bool useImageMask = false;
    cv::Mat imageMaskDist; // Distorted image mask (input)
    if (ros::param::get("~imageMaskFile", imageMaskFile))
    {

        // Read an image mask file (distorted)
        std::ifstream f(imageMaskFile.c_str());
        if (f.good() && f.is_open())
        {
            f.close();
            imageMaskDist = cv::imread(imageMaskFile, CV_LOAD_IMAGE_GRAYSCALE);
            ros::param::del("~imageMaskFile");
            useImageMask = true;
        }
        else
        {

            f.close();
            ros::param::del("~imageMaskFile");
            fprintf(stderr, "No file exist: %s\n", imageMaskFile.c_str());
            printf("No file exist: %s\n", imageMaskFile.c_str());
            exit(-1);
        
        }
        
    }
    else
    {
        
        fprintf(stderr, "imageMaskFile was not set => using the whole image\n");
        printf("imageMaskFile was not set => using the whole image\n");
    
    }
    printf("_imageMaskFile:=%s\n", imageMaskFile.c_str());

    //--------------------------------------------------------------------------
    // get HZ
    double hz = 0;
    if(!ros::param::get("~hz", hz))
        hz = 0;
    ros::param::del("~hz");

    cv::Mat imageInput = cv::Mat(height_output,width_output,CV_8U);
    cv::Mat imageInputRGB = cv::Mat(height_output,width_output,CV_8UC3);
    int runningIDX=0;
    float fakeTimeStamp = 0;

    ros::Rate r(hz);

    //--------------------------------------------------------------------------
    // get B-Spline segment Hz
    double segHz = 0;
    if (!ros::param::get("~segHz", segHz))
    {
        fprintf(stderr, "Error: need spline segment rate in Hz\n"
                        "       e.g. _segHz:=30(Hz)\n");
        exit(0);

    }
    ros::param::del("~segHz");

    printf("Hz = %.2f and segHz = %.2fHz\n", hz, segHz);
    
    //--------------------------------------------------------------------------
    // get readout time in millisecond
    double readoutTime;
    if (!ros::param::get("~readoutTime", readoutTime))
    {
        // Default reading out time in ms same as camera fps
        readoutTime = (1.0/hz)*1000.0;
    }
    ros::param::del("~readoutTime");
    printf("readoutTime = %.2f ms\n", readoutTime);
    
    // Number of rows for one readout
    // This value is used in u_time per image row (u_q)
    double totalRow_per_camHz = ((1.0/hz*1000.0)*height_input)/readoutTime;
    printf("totalRow_per_camHz = %f\n", totalRow_per_camHz);
    
    // Currently this is not used,
    // WORK BACK <- we need numRow_readout in optimization
    // Use numRow_readout in SE3Tracker_BN3
    double delayTime = (1.0/hz)*1000.0 - readoutTime;
    int delayRows = std::round(height_output*(delayTime/readoutTime));
    printf("delayTime = %f ms\n", delayTime);
    printf("delayRows = %d\n", delayRows);
    
    //-----------------------------------------------------------------
    // Get depth file format
    int depthFormat = 0; // Default 0: TUM format; 1: JHK(ICL-NUIM) format
    if (!ros::param::get("~depthFormat", depthFormat))
    {
        
        depthFormat = 0; // Default 0, TUM format
        
    }
    ros::param::del("~depthFormat");

    if ((depthFormat != 0) && (depthFormat != 1)) // Invalid format
    {
        fprintf(stderr, "Error: specify depth file format\n"
                        "E.g. _depthFormat:=0 for JHK or ICL-NUIM format\n"
                        "     _depthFormat:=1 for TUM format\n");
    }

    //--------------------------------------------------------------------------
    // Option for Radial-Rolling or Rolling
    int methodRollingShutter = 1; // 0 : Rolling-shutter (RSD-SLAM)
                                  // 1 : Radial-rolling-shutter (RRD-SLAM)
    if (!ros::param::get("~methodRollingShutter", methodRollingShutter))
    {
        
        methodRollingShutter = 1; // Default 1, Radial-Rolling-shutter
        
    }
    ros::param::del("~methodRollingShutter");
 
    if ((methodRollingShutter != 0) &&
        (methodRollingShutter != 1) ) // Invalid format
    {
        
        fprintf(stderr, "Error: specify method for rolling-shutter\n"
                "E.g. _methodRollingShutter:=0 for Rolling-shutter\n"
                "     _methodRollingShutter:=1 for Radial-Rolling-shutter\n");
        
    }
    
    //--------------------------------------------------------------------------
    // Set noise on GT initial depth
    float noiseOnInitGT = 0.0;
    if (!ros::param::get("~noiseOnInitGT", noiseOnInitGT))
    {
        noiseOnInitGT = 0.0;
    }
    ros::param::del("~noiseOnInitGT");
    fprintf(stderr, "_noiseOnInitGT:=%f\n", noiseOnInitGT);
    printf("_noiseOnInitGT:=%f\n", noiseOnInitGT);
    
    //--------------------------------------------------------------------------
    // Set random init depth range
    float minInitDepth = 0.6666;
    float maxInitDepth = 2.0;
    if (!ros::param::get("~minInitDepth", minInitDepth))
    {
        minInitDepth = 0.6666;
    }
    ros::param::del("~minInitDepth");
    fprintf(stderr, "minInitDepth:=%f\n", minInitDepth);
    printf("minInitDepth:=%f\n", minInitDepth);
    
    if (!ros::param::get("~maxInitDepth", maxInitDepth))
    {
        maxInitDepth = 2.0;
    }
    ros::param::del("~maxInitDepth");
    fprintf(stderr, "maxInitDepth:=%f\n", maxInitDepth);
    printf("maxInitDepth:=%f\n", maxInitDepth);
    
    //--------------------------------------------------------------------------
    // Option for depth display color scheme
    // Setting for debug display
    // Default: 0 (inverse-depth smoothed)
    //          1 (inverse-depth)
    //          2 (validity counter)
    //          3 (inverse-depth variance smoothed)
    //          4 (inverse-depth variance)
    if (!ros::param::get("~debugDisplay", lsd_slam::debugDisplay))
    {
        
        lsd_slam::debugDisplay = 1; // Default, inverse-depth
        
    }
    ros::param::del("~debugDisplay");
    fprintf(stderr, "_debugDisplay:=%d\n", lsd_slam::debugDisplay);
    printf("_debugDisplay:=%d\n", lsd_slam::debugDisplay);
    
    //--------------------------------------------------------------------------
    // Option for sparse tracking
    // Default: 0 (Sparse tracking disabled)
    //          1 (Sparse tracking enabled)
    if (!ros::param::get("~sparseTracking", lsd_slam::sparseTracking))
    {
        
        lsd_slam::sparseTracking = false; // Default: sparse tracking disabled
        
    }
    ros::param::del("~sparseTracking");
    fprintf(stderr, "sparseTracking:=%d\n", lsd_slam::sparseTracking);
    printf("sparseTracking:=%d\n", lsd_slam::sparseTracking);
    if (lsd_slam::sparseTracking == true)
    {
        
        lsd_slam::pixelDensityForTracking = 0.16f; // 16% of image size
        
    }
    fprintf(stderr, "pixelDensityForTracking:=%f\n",
            lsd_slam::pixelDensityForTracking);
    printf("pixelDensityForTracking:=%f\n",
           lsd_slam::pixelDensityForTracking);
    
    //--------------------------------------------------------------------------
    // Option for enable undistortion
    // Default: true
    //
    bool enableUndistort = true;
    if (!ros::param::get("~enableUndistort", enableUndistort))
    {
        
        enableUndistort = true; // Default: enable undistort
        
    }
    ros::param::del("~enableUndistort");
    fprintf(stderr, "enableUndistort:=%d\n", enableUndistort);
    printf("enableUndistort:=%d\n", enableUndistort);
    
    
    //--------------------------------------------------------------------------
    // Option for initialization phase count
    // Default: 5
    //
    if (!ros::param::get("~initializationPhaseCount",
                         lsd_slam::initializationPhaseCount))
    {
        
        // Default: 5
        lsd_slam::initializationPhaseCount = 0;
        
    }
    ros::param::del("~initializationPhaseCount");
    fprintf(stderr, "initializationPhaseCount:=%d\n",
            lsd_slam::initializationPhaseCount);
    printf("initializationPhaseCount:=%d\n",
           lsd_slam::initializationPhaseCount);
    

    //--------------------------------------------------------------------------
    // Set other rate parameters
    double framePerSegment = hz/segHz;
    double rowPerSegment = framePerSegment*(height_input + delayRows);
    printf("rowPerSegment = %f\n", rowPerSegment);

#if 0
    if (segHz/hz - floor(segHz/hz) != 0.0)
    {

        fprintf(stderr, "Error: segHz shoud be divisible by hz.\n"
                        "       e.g. __segHz:=60 and _hz:=30\n");
        exit(0);

    }
#endif

    printf("imageSource = %s\n", source.c_str());
    printf("useGTmotion = %d\n", useGTmotion);
    printf("useGTdepth = %d\n", useGTdepth);
    printf("methodRollingShutter = %d\n", methodRollingShutter);
    fprintf(stderr, "imageSource = %s\n", source.c_str());
    fprintf(stderr, "useGTmotion = %d\n", useGTmotion);
    fprintf(stderr, "useGTdepth = %d\n", useGTdepth);
    fprintf(stderr, "methodRollingShutter = %d\n", methodRollingShutter);

    //--------------------------------------------------------------------------
    // Number of image rows per one segment of B-Spline curve
    int numRowsPerSegment = ceil((hz/segHz)*(height_input + delayRows));
    printf("number of image rows per segment = %d\n", numRowsPerSegment);

    //--------------------------------------------------------------------------
    // make slam system
    
    // Global Shutter SLAM System for initial random depth
    lsd_slam::SlamSystem* systemGS =
            new lsd_slam::SlamSystem(width_output,
                                     height_output,
                                     K,
                                     lsd_slam::doSlam);
    systemGS->setVisualization(outputWrapper);
    
    //--------------------------------------------------------------------------
    
    //--------------------------------------------------------------------------
    // Settings
    
    // Settings for depth
    if (!ros::param::get("~minUseGrad", lsd_slam::minUseGrad))
    {
        lsd_slam::minUseGrad = 5; // Default 5
    }
    ros::param::del("~minUseGrad");
    lsd_slam::cameraPixelNoise2 = 4*4; // Default: 4*4
    if (!ros::param::get("~allowNegativeIdepths", lsd_slam::allowNegativeIdepths))
    {
        lsd_slam::allowNegativeIdepths = true; // Default: true
    }
    ros::param::del("~allowNegativeIdepths");
    lsd_slam::useSubpixelStereo = true; // Default: true

    // Settings for keyframe selection
    if (!ros::param::get("~KFUsageWeight", lsd_slam::KFUsageWeight))
    {
        lsd_slam::KFUsageWeight = 3; //1.2; // Default:3, Max:20
    }
    ros::param::del("~KFUsageWeight");
    if (!ros::param::get("~KFDistWeight", lsd_slam::KFDistWeight))
    {
        lsd_slam::KFDistWeight = 4; //1.2; // Default:4, Max: 20
    }
    ros::param::del("~KFDistWeight");
    lsd_slam::doKFReActivation = false;
    
    // Print essential info to observe in experiments
    fprintf(stderr, "KFUsageWeight = %f\n", lsd_slam::KFUsageWeight);
    fprintf(stderr, "KFDistWeight = %f\n", lsd_slam::KFDistWeight);
    fprintf(stderr, "allowNegativeIdepth = %d\n", lsd_slam::allowNegativeIdepths);
    printf("KFUsageWeight = %f\n", lsd_slam::KFUsageWeight);
    printf("KFDistWeight = %f\n", lsd_slam::KFDistWeight);
    printf("minUseGrad = %f\n", lsd_slam::minUseGrad);
    printf("allowNegativeIdepth = %d\n", lsd_slam::allowNegativeIdepths);
    printf("cameraPixelNoise2 = %f\n", lsd_slam::cameraPixelNoise2);
    
    // Settings for print or debug
    lsd_slam::consoleMode = true; // Default: false;

    if (lsd_slam::consoleMode)
    {
        lsd_slam::plotStereoImages = false; // Default: true
        lsd_slam::plotTracking = false;
        lsd_slam::plotTrackingIterationInfo = false;
        lsd_slam::printThreadingInfo = true;
        lsd_slam::displayDepthMap = false;
        systemGS->set_depthMapScreenshotFlag(true);
    }
    else
    {
        lsd_slam::plotStereoImages = false; // Default: true
        lsd_slam::plotTracking = true;
        lsd_slam::plotTrackingIterationInfo = true;
        lsd_slam::printThreadingInfo = true;
        lsd_slam::printRelocalizationInfo = true;
        lsd_slam::displayDepthMap = true;
    }

    XInitThreads();
    
    //--------------------------------------------------------------------------
    // Running Global-Shutter SLAM system
    //--------------------------------------------------------------------------
    runningIDX = 0;
    fakeTimeStamp = 0;
    for (int frameIdx=0;
         frameIdx<lsd_slam::initializationPhaseCount &&
         frameIdx<(int)imageFiles.size();
         frameIdx++)
    {
        
        cv::Mat imageDist = cv::imread(imageFiles[frameIdx], CV_LOAD_IMAGE_GRAYSCALE);
        
        if (imageDist.rows != height_input || imageDist.cols != width_input)
        {
            
            if(imageDist.rows * imageDist.cols == 0)
            {
                
                printf("failed to load image %s! skipping.\n",
                       imageFiles[frameIdx].c_str());
                
            }
            else
            {
                
                printf("image %s has wrong dimensions "
                       "- expecting %d x %d, found %d x %d. Skipping.\n",
                       imageFiles[frameIdx].c_str(),
                       width_output,height_output,imageDist.cols, imageDist.rows);
                
            }
            
            continue;
            
        }
        assert(imageDist.type() == CV_8U);
        
        undistorter->undistort(imageDist, imageInput);
        assert(imageInput.type() == CV_8U);

        printf("runningIDX %d, fakeTimeStamp %f\n", runningIDX, fakeTimeStamp);

        if(runningIDX == 0)
            systemGS->randomInit(imageInput.data, fakeTimeStamp, runningIDX);
        else
            systemGS->trackFrame(imageInput.data, runningIDX, hz == 0, fakeTimeStamp);
    
        // Increase frame index and time stamp
        runningIDX++;
        fakeTimeStamp += 1.0/hz; // seconds
        
        if(hz != 0)
            r.sleep();
        
        // For ROS messages and services
        ros::spinOnce();
        if(!ros::ok())
            break;
        
        // Finish if running index reaches the end of image frame
        if (runningIDX == (int)imageFiles.size())
            break;


    }
    
    //--------------------------------------------------------------------------
    // TODO: Transfer GS to RS depth
    
    // Check the depth is available
    bool transferDepthFrom_systemGS = false;
    lsd_slam::Frame *lastKeyframeGS;
    int startIDX = 0;
    if ((lsd_slam::initializationPhaseCount > 0) &&
        (systemGS->getKeyframesAll().size() > 0))
    {
    
        // Obtain the keyframe and starting idx
        // Now we will use this keyframe as the first starting frame
        // in RS-System
        lastKeyframeGS = systemGS->getKeyframesAll().back();
        transferDepthFrom_systemGS = true;
        startIDX = lastKeyframeGS->id();
        
    }

    
    //--------------------------------------------------------------------------
    // Running Rolling-Shutter SLAM sytem
    //--------------------------------------------------------------------------
    
    // Rolling Shutter SLAM System
    int width_bsplineSegImg = width_input;
    int height_bsplineSegImg = height_input;
    lsd_slam::SlamSystemForRollingShutter* systemRS =
    new lsd_slam::SlamSystemForRollingShutter(width_output,
                                              height_output,
                                              K,
                                              undistorter,
                                              width_bsplineSegImg,
                                              height_bsplineSegImg,
                                              numRowsPerSegment,
                                              methodRollingShutter,
                                              lsd_slam::doSlam,
                                              useGTmotion);
    
    systemRS->setVisualization(outputWrapper);
    if (lsd_slam::consoleMode == true)
    {
        systemRS->set_depthMapScreenshotFlag(true);
    }
    
    fprintf(stderr, "B-spline order k =  %d\n", systemRS->getTracker()->spline_k);
    printf("B-spline order k =  %d\n", systemRS->getTracker()->spline_k);
    
    // Read images and process them
    runningIDX = 0;
    fakeTimeStamp = 0;
    int startRowImage = 0; // index of the image row
    int startRowSeg = 0; // index of the segment row
    for(int frameIdx=startIDX; frameIdx<(int)imageFiles.size(); frameIdx++)
    {

        //----------------------------------------------------------------------
        // Read an image from the video file
        cv::Mat imageDist = cv::imread(imageFiles[frameIdx],
                                       CV_LOAD_IMAGE_GRAYSCALE);
        cv::Mat imageDistRGB = cv::imread(imageFiles[frameIdx],
                                          CV_LOAD_IMAGE_COLOR);
        
        // If the distortion emulation mode, then we distort the image
        if (emulDistort != 0)
        {
            
            cv::Mat emulImageDist = cv::Mat(width_input, height_input, CV_8UC1);
            undistorter->distortGray(imageDist, emulImageDist);
            
            // Save
            cv::Mat emulImageDistRGB = cv::Mat(width_input,
                                               height_input,
                                               CV_8UC3);
            undistorter->distortRGB(imageDistRGB, emulImageDistRGB);
            char buf[500];
            sprintf(buf, "%06d_emulDist.png", frameIdx);
            cv::imwrite(buf, emulImageDistRGB);
            
            // Replace with the emulated distorted image
            imageDist = emulImageDist;
            imageDistRGB = emulImageDistRGB;
            
        }

        if (imageDist.rows != height_input || imageDist.cols != width_input)
        {

            if(imageDist.rows * imageDist.cols == 0)
            {

                printf("failed to load image %s! skipping.\n",
                       imageFiles[frameIdx].c_str());

            }
            else
            {

                printf("image %s has wrong dimensions "
                       "- expecting %d x %d, found %d x %d. Skipping.\n",
                        imageFiles[frameIdx].c_str(),
                        width_output,height_output,imageDist.cols, imageDist.rows);

            }

            continue;

        }
        assert(imageDist.type() == CV_8U);
        //----------------------------------------------------------------------

        // Undistort the image
        if (enableUndistort == true)
        {
        
            undistorter->undistort(imageDist, imageInput);
            undistorter->undistortRGB(imageDistRGB, imageInputRGB);
            assert(imageInput.type() == CV_8U);
            
            // Undistort image mask for first frame (used for the rest)
            if (useImageMask == true)
            {
                
                undistorter->undistort(imageMaskDist, imageMask);
                assert(imageMask.type() == CV_8U);
                
            }
        
        }
        else // Already rectified input image
        {
            
            imageInput = imageDist;
            imageInputRGB = imageDistRGB;
        
        }
        char buf2[500];
        sprintf(buf2, "%06d_undist.jpg", frameIdx);
        cv::imwrite(buf2, imageInputRGB);

        //----------------------------------------------------------------------
        // Tracking
        printf("runningIDX = %d\n", runningIDX);
        if (runningIDX == 0)
        {

            if (useGTdepth == false)
            {

                if (transferDepthFrom_systemGS == true)
                {
                    
                    // Initialize id and timestamp of lastKeyframeGS
                    lastKeyframeGS->setID(runningIDX);
                    lastKeyframeGS->setTimestamp(fakeTimeStamp);
                    lastKeyframeGS->setImageMask(imageMask.data);

                    // Transfer initial depth obtained GS SLAM system
                    // to RS SLAM system
                    systemRS->depthTransferFrom(lastKeyframeGS);
                    
                }
                else
                {
                    
                
                    // Random depth initialization
                    systemRS->randomInitDepthRange(imageInput,
                                                   imageMask,
                                                   fakeTimeStamp,
                                                   runningIDX,
                                                   minInitDepth,
                                                   maxInitDepth);
                    
                }

            }
            else
            {

                
                // Use Ground truth depth init
                float *depthInput = (float*)malloc(sizeof(float)*width_output*height_output);
                
                
                // Select initial ground-truth depth loading method
                // based on validty of undistorter
                // Or distortion emulation mode uses an undistorted depth
                if (undistorter->isNoRectification() == true ||
                    emulDistort != 0)
                {
                    
                    // Load init GT
                    fprintf(stderr,
                            "Loading initial GT without radial distortion\n");
                    printf("Loading initial GT without radial distortion\n");
                    useGroundTruthDepth(depthInput,
                                        depthFiles,
                                        runningIDX,
                                        width_output, fx, fy, cx, cy,
                                        depthFormat,
                                        noiseOnInitGT);
                    
                }
                else
                {
                
                    // Load init GT considering radial distortion
                    fprintf(stderr,
                            "Loading initial GT with radial distortion\n");
                    printf("Loading initial GT with radial distortion\n");
                    useGroundTruthDepthRadial(undistorter,
                                              depthInput,
                                              depthFiles,
                                              runningIDX,
                                              width_output, height_output,
                                              fx,
                                              fy,
                                              cx,
                                              cy,
                                              depthFormat,
                                              noiseOnInitGT);
                    
                }

                if (useGTmotion == true)
                {

                    systemRS->gtDepthInitWithName(imageInput.data,
                                                depthInput,
                                                depthFiles[runningIDX],
                                                motionFiles[runningIDX],
                                                fakeTimeStamp,
                                                runningIDX);
                }
                else
                {

                    if (noiseOnInitGT == 0.0) // no noise on GT depth
                    {
                        systemRS->gtDepthInit(imageInput.data,
                                            imageMask.data,
                                            depthInput,
                                            fakeTimeStamp,
                                            runningIDX);

                    }
                    else
                    {
                        // Noise on GT depth
                        systemRS->gtDepthInitWithNoise(imageInput.data,
                                                     depthInput,
                                                     fakeTimeStamp,
                                                     runningIDX);
                        
                    }
                }

                delete depthInput;

            } // end if useGTdepth

            printf("runningIDX %d, startRowImage %d, startRowSeg %d\n",
                   runningIDX, startRowImage, startRowSeg);

            // Tracking
            systemRS->trackFrameRS_BN3(imageInput,
                                     imageMask,
                                     runningIDX,
                                     true,
                                     fakeTimeStamp,
                                     startRowImage,
                                     startRowSeg,
                                     rowPerSegment,
                                     framePerSegment,
                                     height_input,
                                     imageFiles,
                                     depthFiles,
                                     motionFiles);

        } // end if runningIDX == 0
        else
        {
            
            printf("runningIDX %d, startRowImage %d, startRowSeg %d\n",
                   runningIDX, startRowImage, startRowSeg);

            // Tracking frame for rolling shutter
            systemRS->trackFrameRS_BN3(imageInput,
                                     imageMask,
                                     runningIDX,
                                     true,
                                     fakeTimeStamp,
                                     startRowImage,
                                     startRowSeg,
                                     rowPerSegment,
                                     framePerSegment,
                                     height_input,
                                     imageFiles,
                                     depthFiles,
                                     motionFiles);
            
        }

        printf("reminder = %d\n", startRowImage % (int)rowPerSegment);

        // Update row index for image
        startRowImage = startRowImage + height_input + delayRows;

        // Update row index for segment
        if (startRowImage % (int)rowPerSegment == 0)
        {

            startRowSeg = startRowImage;

        }

        // Increase frame index and time stamp
        runningIDX++;
        fakeTimeStamp += 1.0/hz; // seconds

        // For ROS messages and services
        ros::spinOnce();
        if(!ros::ok())
            break;
        
        // Finish if running index reaches the end of image frame
        if (frameIdx == (int)imageFiles.size())
            break;

    } // end of for

    systemRS->finalize();
    delete systemRS;
    
    systemGS->finalize();
    delete systemGS;
    
    delete undistorter;
    delete outputWrapper;
    return 0;

}

