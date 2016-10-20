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
#include "Tracking/WarpResidual.h"
#include "util/globalFuncs.h"
#include "IOWrapper/ImageDisplay.h"

namespace lsd_slam
{

// Allocation of static warping pixel counters
int lsd_slam::WarpResidual::goodWarpPixelCountInPyramidImage_ = 0;
int lsd_slam::WarpResidual::badWarpPixelCountInPyramidImage_ = 0;
float lsd_slam::WarpResidual::usageCount_ = 0.0;
int lsd_slam::WarpResidual::lastGoodCount_ = 0;
int lsd_slam::WarpResidual::lastBadCount_ = 0;
int lsd_slam::WarpResidual::numWarpedPixelOutOfImage  = 0;

// Allocation of static debug images
//cv::Mat lsd_slam::WarpResidual::debugImageResiduals
//    = cv::Mat(480, 640, CV_8UC3);
//cv::Mat lsd_slam::WarpResidual::debugImageWeights
//    = cv::Mat(480, 640, CV_8UC3);
//cv::Mat lsd_slam::WarpResidual::debugImageOldImageSource
//    = cv::Mat(480, 640, CV_8UC3);
//cv::Mat lsd_slam::WarpResidual::debugImageOldImageWarped
//    = cv::Mat(480, 640, CV_8UC3);

#if 0
void WarpResidual::debugStart() const
{
    if (lsd_slam::plotTrackingIterationInfo ||
            lsd_slam::saveAllTrackingStagesInternal)
    {
        int other = lsd_slam::saveAllTrackingStagesInternal ? 255 : 0;
        fillCvMat(&debugImageResiduals,cv::Vec3b(other,other,255));
        fillCvMat(&debugImageWeights,cv::Vec3b(other,other,255));
        fillCvMat(&debugImageOldImageSource,cv::Vec3b(other,other,255));
        fillCvMat(&debugImageOldImageWarped,cv::Vec3b(other,other,255));
    }
}

void WarpResidual::debugFinish() const
{

    if (lsd_slam::plotTrackingIterationInfo == true)
    {
        lsd_slam::Util::displayImage( "!Weights",
            lsd_slam::SE3TrackerRollingShutter::debugImageWeights);
        lsd_slam::Util::displayImage( "!Intensities of second_frame "
                            "at transformed positions",
                            lsd_slam::WarpResidual::debugImageOldImageSource );
        lsd_slam::Util::displayImage( "!Intensities of second_frame "
                            "at pointcloud in first_frame",
                            lsd_slam::WarpResidual::debugImageOldImageWarped );
        lsd_slam::Util::displayImage( "!Residuals",
                            lsd_slam::WarpResidual::debugImageResiduals );


        // wait for key and handle it
        bool looping = true;
        while (looping)
        {
            int k = lsd_slam::Util::waitKey(1);
            if(k == -1)
            {
                if (lsd_slam::autoRunWithinFrame == true)
                    break;
                else
                    continue;
            }

            char key = k;
            if (key == ' ')
                looping = false;
            else
                handleKey(k);
        }
    }


}
#endif

} // end of namespace
