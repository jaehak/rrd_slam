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
// Spline class for B-Spline curve of camera trajectory
//------------------------------------------------------------------------------
// Copyright@2014 University of Adelaide
//
// Author: Jae-Hak Kim (jaehak.kim@adelaide.edu.au)
//------------------------------------------------------------------------------
#include "Spline.h"

namespace lsd_slam
{

// Spline initializer with n control poses of T
#if 0
template <typename T>
void
Spline<T>::init(std::vector<Eigen::Matrix<T,4,4>> controlPoses_T)
{

    // Set the number of control poses
    this->n = controlPoses_T.size();

    // Check n
    if (this->n < 3)
    {

        printf("Error: Not enough number of control poses T\n");
        exit(1);

    }

    // Set the degree k - 1 of the B-Spline curve
    this->k = 4;

    // Set control poses
    this->controlposeTs = controlPoses_T;

    // Set the basis
    this->C << 6, 0, 0, 0,
               5, 3, -3, 1,
               1, 3, 3, -2,
               0, 0, 0, 1;
    this->C = (1.0/6.0)*this->C;

}

template <typename T>
Spline<T>::Spline(std::vector<Sophus::SE3Group<T>> Xi)
{

    std::vector<Eigen::Matrix<T,4,4>> controlPoses;
    for (int i=0; i<Xi.size(); i++)
    {

        Eigen::Matrix<T,4,4> Xi_mat = Xi[i].matrix();
        controlPoses.push_back(Xi_mat);

    }

    this->init(controlPoses);

}
#endif


} // end of namespace
