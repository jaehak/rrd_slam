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
#include "WarpResidual_BN3.h"

namespace lsd_slam
{

int lsd_slam::WarpResidual_BN3::goodWarpPixelCountInPyramidImage_ = 0;
int lsd_slam::WarpResidual_BN3::badWarpPixelCountInPyramidImage_ = 0;
float lsd_slam::WarpResidual_BN3::usageCount_ = 0.0;
int lsd_slam::WarpResidual_BN3::lastGoodCount_ = 0;
int lsd_slam::WarpResidual_BN3::lastBadCount_ = 0;

} // end of namespace
