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

#pragma once
#include "util/settings.h"
#include "util/EigenCoreInclude.h"
#include "boost/thread/mutex.hpp"
#include <boost/thread/shared_mutex.hpp>
#include "DataStructures/Frame.h"
#include <algorithm>

namespace lsd_slam
{

class Frame;
class DepthMapPixelHypothesis;
class KeyFrameGraph;

/**
 * Point cloud used to track frame poses.
 * 
 * Basically this stores a point cloud generated from known frames. It is used to
 * track a new frame by finding a projection of the point cloud which makes it
 * look as much like the new frame as possible.
 * 
 * It is intended to use more than one old frame as source for the point cloud.
 * Also other data like Kinect depth data could be imported.
 * 
 * ATTENTION: as the level zero point cloud is not used for tracking, it is not
 * fully calculated. Only the weights are valid on this level!
 */
class TrackingReference
{
public:
	/** Creates an empty TrackingReference with optional preallocation per level. */
	TrackingReference();
	~TrackingReference();


    // Copy constructor (Jae-Hak Kim, 2015-04-19)
    TrackingReference (const TrackingReference &copyOfTrackingReference);

#if 0
    // Copy assignment (Jae-Hak Kim, 2015-04-19)
    TrackingReference& operator=(const TrackingReference & other)
    {

        //std::swap(keyframe, other.keyframe);
        this->importFrame(other.keyframe);
        std::swap(frameID, other.frameID);
        std::swap(posData, other.posData);
        std::swap(gradData, other.gradData);
        std::swap(colorAndVarData, other.colorAndVarData);
        std::swap(pointPosInXYGrid, other.pointPosInXYGrid);
        std::swap(numData, other.numData);

        return *this;

    }
#endif


    TrackingReference& operator=(const TrackingReference & other)
    {

        if (this != &other)
        {

            //------------------------------------------------------------------
            // Allocate new memory and copy
#if 1// USE copy constructor of Frame
            boost::unique_lock<boost::mutex> lock(accessMutex);
            Frame *new_keyframe = new Frame(*(other.keyframe));
            keyframeLock = other.keyframe->getActiveLock();
            frameID = other.frameID;
            wh_allocated = other.wh_allocated;
            lock.unlock();
#else // Use importFrame
            //this->importFrame(other.keyframe);
            boost::unique_lock<boost::mutex> lock(accessMutex);
            keyframeLock = other.keyframe->getActiveLock();
            keyframe = other.keyframe;
            frameID=keyframe->id();
            lock.unlock();
#endif

            Eigen::Vector3f* new_posData[PYRAMID_LEVELS];
            Eigen::Vector2f* new_gradData[PYRAMID_LEVELS];
            Eigen::Vector2f* new_colorAndVarData[PYRAMID_LEVELS];
            int* new_pointPosInXYGrid[PYRAMID_LEVELS];
            int new_numData[PYRAMID_LEVELS];

            for (int level=0;level<PYRAMID_LEVELS;level++)
            {

                int w = other.keyframe->width(level);
                int h = other.keyframe->height(level);

                new_posData[level] = new Eigen::Vector3f[w*h];
                new_gradData[level] = new Eigen::Vector2f[w*h];
                new_colorAndVarData[level] = new Eigen::Vector2f[w*h];
                new_pointPosInXYGrid[level] = new int[w*h];


                // Pointers of target
                Eigen::Vector3f* posDataPT = new_posData[level];
                Eigen::Vector2f* gradDataPT = new_gradData[level];
                Eigen::Vector2f* colorAndVarDataPT = new_colorAndVarData[level];
                int *idxPT = new_pointPosInXYGrid[level];

                // Pointers of copying source
                Eigen::Vector3f* copy_posDataPT =
                        other.posData[level];
                Eigen::Vector2f* copy_gradDataPT =
                        other.gradData[level];
                Eigen::Vector2f* copy_colorAndVarDataPT =
                        other.colorAndVarData[level];
                int *copy_idxPT = other.pointPosInXYGrid[level];

                // Copy data
                for (int i=0;i<other.numData[level];i++)
                {

                    // Copy source to target
                    *posDataPT = *copy_posDataPT;
                    *gradDataPT = *copy_gradDataPT;
                    *colorAndVarDataPT = *copy_colorAndVarDataPT;
                    *idxPT = *copy_idxPT;

                    // Increase pointers of target
                    posDataPT++;
                    gradDataPT++;
                    colorAndVarDataPT++;
                    idxPT++;

                    // Increase pointers of copying source
                    copy_posDataPT++;
                    copy_gradDataPT++;
                    copy_colorAndVarDataPT++;
                    copy_idxPT++;

                } // end of for numData

                new_numData[level] = other.numData[level];

            }


            //-----------------------------------------------------------------
            // Deallocate old memory
#if 1// USE copy constructor of Frame
            delete keyframe;
#endif
            for (int level=0; level<PYRAMID_LEVELS; level++)
            {
                delete[] posData[level];
                delete[] gradData[level];
                delete[] colorAndVarData[level];
                delete[] pointPosInXYGrid[level];
            }

            //------------------------------------------------------------------
            // Assign the new memory
#if 1 // Use copy constructor of Frame
            keyframe = new_keyframe;
#endif
            for (int level=0; level<PYRAMID_LEVELS; level++)
            {

                posData[level] = new_posData[level];
                gradData[level] = new_gradData[level];
                colorAndVarData[level] = new_colorAndVarData[level];
                pointPosInXYGrid[level] = new_pointPosInXYGrid[level];

                numData[level] = new_numData[level];

            }
#if 1 // Use copy constructor
#else// USE importFrame
            wh_allocated = other.wh_allocated;
#endif
        } // end if

        return *this;

    }


#if 0
TrackingReference& operator=(const TrackingReference & other)
{

    if (this != &other)
    {

            //------------------------------------------------------------------
            // Copy keyframe and frameID
            this->importFrame(other.keyframe);

            // Copy data for all leves except level 0
            for (int level=0;level<PYRAMID_LEVELS;level++)
            {

                // Allocate data
                int w = other.keyframe->width(level);
                int h = other.keyframe->height(level);
                if (posData[level] == nullptr)
                {

                    posData[level] = new Eigen::Vector3f[w*h];

                }
                if (pointPosInXYGrid[level] == nullptr)
                {

                    pointPosInXYGrid[level] = new int[w*h];

                }
                if (gradData[level] == nullptr)
                {
                    gradData[level] = new Eigen::Vector2f[w*h];

                }
                if (colorAndVarData[level] == nullptr)
                {

                    colorAndVarData[level] = new Eigen::Vector2f[w*h];

                }

                // Pointers of target
                Eigen::Vector3f* posDataPT = posData[level];
                int *idxPT = pointPosInXYGrid[level];
                Eigen::Vector2f* gradDataPT = gradData[level];
                Eigen::Vector2f* colorAndVarDataPT = colorAndVarData[level];

                // Pointers of copying source
                Eigen::Vector3f* copy_posDataPT =
                        other.posData[level];
                int *copy_idxPT = other.pointPosInXYGrid[level];
                Eigen::Vector2f* copy_gradDataPT =
                        other.gradData[level];
                Eigen::Vector2f* copy_colorAndVarDataPT =
                        other.colorAndVarData[level];

                // Copy data
                for (int i=0;i<other.numData[level];i++)
                {

                    // Copy source to target
                    *posDataPT = *copy_posDataPT;
                    *gradDataPT = *copy_gradDataPT;
                    *colorAndVarDataPT = *copy_colorAndVarDataPT;
                    *idxPT = *copy_idxPT;

                    // Increase pointers of target
                    posDataPT++;
                    gradDataPT++;
                    colorAndVarDataPT++;
                    idxPT++;

                    // Increase pointers of copying source
                    copy_posDataPT++;
                    copy_gradDataPT++;
                    copy_colorAndVarDataPT++;
                    copy_idxPT++;

                } // end of for numData

                numData[level] = other.numData[level];

            } // end of for pyramid level

        }

        return *this;

    }
#endif

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
	void importFrame(Frame* source);

	Frame* keyframe;
	boost::shared_lock<boost::shared_mutex> keyframeLock;
	int frameID;

	void makePointCloud(int level);
	void clearAll();
	void invalidate();
	Eigen::Vector3f* posData[PYRAMID_LEVELS];	// (x,y,z)
	Eigen::Vector2f* gradData[PYRAMID_LEVELS];	// (dx, dy)
	Eigen::Vector2f* colorAndVarData[PYRAMID_LEVELS];	// (I, Var)
	int* pointPosInXYGrid[PYRAMID_LEVELS];	// x + y*width
	int numData[PYRAMID_LEVELS];

private:
	int wh_allocated;
	boost::mutex accessMutex;
	void releaseAll();
};
}
