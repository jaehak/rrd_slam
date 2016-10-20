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
// MyMats in C++
//
// Copyright@2015, University of Adelaide
// Author: Jae-Hak Kim <jaehak.kim@adelaide.edu.au>
//------------------------------------------------------------------------------
#ifndef __MYMATS__
#define __MYMATS__
#include <Eigen/Eigen>
#include <iostream>

namespace JHK
{

class MyMats
{
public:

    //--------------------------------------------------------------------------
    // Convert a vector/matrix into homogeneous coordinates
    //            x : [Eigen::MatrixX] a vector or matrix 
    // returns
    //            a vector or matrix in homogeneous form 
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
            
            if (x(m - 1, j) <= (1.0e-6) &&
                x(m - 1, j) >= (-1.0e-6))
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

};

} // end of namepsace
#endif
