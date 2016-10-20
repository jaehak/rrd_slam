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
#include "util/SophusUtil.h"
#include "GlobalMapping/g2oTypeSim3Sophus.h"



namespace lsd_slam
{
class Frame;
class FramePoseStruct {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	FramePoseStruct(Frame* frame);
	virtual ~FramePoseStruct();

	// parent, the frame originally tracked on. never changes.
	FramePoseStruct* trackingParent;

	// set initially as tracking result (then it's a SE(3)),
	// and is changed only once, when the frame becomes a KF (->rescale).
	Sim3 thisToParent_raw;


	int frameID;
	Frame* frame;


	// whether this poseStruct is registered in the Graph. if true MEMORY WILL BE HANDLED BY GRAPH
	bool isRegisteredToGraph;

	// whether pose is optimized (true only for KF, after first applyPoseGraphOptResult())
	bool isOptimized;

	// true as soon as the vertex is added to the g2o graph.
	bool isInGraph;

	// graphVertex (if the frame has one, i.e. is a KF and has been added to the graph, otherwise 0).
	VertexSim3* graphVertex;

	void setPoseGraphOptResult(Sim3 camToWorld);
	void applyPoseGraphOptResult();
	Sim3 getCamToWorld(int recursionDepth = 0);
	void invalidateCache();

#if 1 // Copy assignment operator (Jae-Hak Kim, 2015-05-01)
    FramePoseStruct& operator=(const FramePoseStruct &other)
    {

        if (this != &other)
        {

            // Allocate a new memory
            trackingParent = other.trackingParent;
//            thisToParent_raw = Sim3(other.thisToParent_raw.rxso3(),
//                                    other.thisToParent_raw.translation());
            Sim3 tmp(other.thisToParent_raw.matrix());
            thisToParent_raw = tmp;
/*
            Sim3 *new_thisToParent_raw = new Sim3(
                        other.thisToParent_raw.rxso3(),
                        other.thisToParent_raw.translation());
*/
            frameID = other.frameID;
            frame = other.frame;
            isRegisteredToGraph = other.isRegisteredToGraph;
            isOptimized = other.isOptimized;
            graphVertex = other.graphVertex;
            cacheValidFor = other.cacheValidFor;
            cacheValidCounter = other.cacheValidCounter;
            camToWorld = other.camToWorld;
            camToWorld_new = other.camToWorld_new;
            hasUnmergedPose = other.hasUnmergedPose;
/*
            thisToParent_raw = *new_thisToParent_raw;
*/
        }

        return *this;

    }
#endif

#if 1 // Copy consturctor (Jae-Hak Kim, 2015-05-02)

    FramePoseStruct(const FramePoseStruct &copy)
    {
        
        trackingParent = copy.trackingParent;
        thisToParent_raw = copy.thisToParent_raw;
        frameID = copy.frameID;
        frame = copy.frame;
        isRegisteredToGraph = copy.isRegisteredToGraph;
        isOptimized = copy.isOptimized;
        graphVertex = copy.graphVertex;
        cacheValidFor = copy.cacheValidFor;
        cacheValidCounter = copy.cacheValidCounter;
        camToWorld = copy.camToWorld;
        camToWorld_new = copy.camToWorld_new;
        hasUnmergedPose = copy.hasUnmergedPose;
        
    }

#endif

private:
	int cacheValidFor;
	static int cacheValidCounter;

	// absolute position (camToWorld).
	// can change when optimization offset is merged.
	Sim3 camToWorld;

	// new, optimized absolute position. is added on mergeOptimization.
	Sim3 camToWorld_new;

	// whether camToWorld_new is newer than camToWorld
	bool hasUnmergedPose;
};

}
