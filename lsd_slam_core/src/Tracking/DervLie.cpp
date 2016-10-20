//
//  DervLie.cpp
//
//  Created by Jae-Hak Kim on 16/11/2015.
//
//
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

#include "DervLie.h"

namespace lsd_slam
{

//--------------------------------------------------------------------------
// Verified Jacobians and functions

//--------------------------------------------------------------------------
#if 0 // COMPACT VERSION
    
Eigen::Matrix<double,9,3>
DervLie::
dR_dw(Eigen::Matrix<double,3,1> w)
{
    
    Eigen::Matrix<double,9,3> M;
    
    Eigen::Matrix<double,3,1> e1, e2, e3;
    e1 << 1.0, 0.0, 0.0;
    e2 << 0.0, 1.0, 0.0;
    e3 << 0.0, 0.0, 1.0;
    
    double w1 = w[0];
    double w2 = w[1];
    double w3 = w[2];
    Eigen::Matrix<double,3,3> R = Sophus::SO3d::exp(w).matrix();
    double vnorm = w.norm();
    Eigen::Matrix<double,3,3> Id = Eigen::Matrix<double,3,3>::Identity();
    
    
    double theta = sqrt(w.transpose()*w);
    Eigen::Matrix<double,3,3> dR_dw1, dR_dw2, dR_dw3;
    if (theta < 1e-5)
    {
        dR_dw1 = skewsym(e1);
        dR_dw2 = skewsym(e2);
        dR_dw3 = skewsym(e3);
        
    }
    else
    {
        // Eq(9) in "A Compact Formula for the Derivative of a 3-D Rotation
        // in Exponential Coordinates" by G. Gallego and A. Yezzi, JMIV 2015
        dR_dw1 = (w1*skewsym(w) + skewsym(w.cross((Id - R)*e1)))*R/vnorm;
        dR_dw2 = (w2*skewsym(w) + skewsym(w.cross((Id - R)*e2)))*R/vnorm;
        dR_dw3 = (w3*skewsym(w) + skewsym(w.cross((Id - R)*e3)))*R/vnorm;
        
    }

    for(int j=0;j<3;j++)
    {
        
        if (j == 0)
        {
            
            M.block(0,j,3,1) = dR_dw1.block(0,0,3,1);
            M.block(3,j,3,1) = dR_dw1.block(0,1,3,1);
            M.block(6,j,3,1) = dR_dw1.block(0,2,3,1);
            
        }
        if (j == 1)
        {
            
            M.block(0,j,3,1) = dR_dw2.block(0,0,3,1);
            M.block(3,j,3,1) = dR_dw2.block(0,1,3,1);
            M.block(6,j,3,1) = dR_dw2.block(0,2,3,1);
            
        }
        if (j == 2)
        {
            
            M.block(0,j,3,1) = dR_dw3.block(0,0,3,1);
            M.block(3,j,3,1) = dR_dw3.block(0,1,3,1);
            M.block(6,j,3,1) = dR_dw3.block(0,2,3,1);
            
        }
        
    }
    
    return M;
    
}
    
#else
    
// dR_dw [3 x 3]: Verified by Ceres
Eigen::Matrix<double,9,3>
DervLie::
dR_dw(Eigen::Matrix<double,3,1> w)
{
    
    Eigen::Matrix<double,9,3> M;
    
    double theta = sqrt(w.transpose()*w);
    
    Eigen::Matrix<double,3,3> G4, G5, G6;
    G4 <<   0,  0,  0,
    0,  0, -1,
    0,  1,  0;
    G5 <<   0,  0,  1,
    0,  0,  0,
    -1,  0,  0,
    G6 <<   0, -1,  0,
    1,  0,  0,
    0,  0,  0;
    
    Eigen::Matrix<double,3,3> wx = skewsym(w);
    Eigen::Matrix<double,3,3> wx2 = wx*wx;
    
//    // Jacobian at identity
//    if (fabs(theta) < 1e-4)
//    {
//        
//        Eigen::Matrix<double,3,3> M1, M2, M3;
//        M << -skewsym(Eigen::Matrix<double,3,1>(1,0,0)),
//             -skewsym(Eigen::Matrix<double,3,1>(0,1,0)),
//             -skewsym(Eigen::Matrix<double,3,1>(0,0,1));
////        M.setZero();
////        M(5,0) =  1;
////        M(7,0) = -1;
////        M(2,1) = -1;
////        M(6,1) =  1;
////        M(1,2) =  1;
////        M(3,2) = -1;
//        
//    }
//    else // Jacobian at non-identity
    {
        
        double s1, s2, s3, s4;

        s1 = (theta*cos(theta) - sin(theta))/pow(theta,3.0);
        s2 = sin(theta)/theta;
        s3 = (theta*sin(theta) + 2*cos(theta) - 2.0)/pow(theta,4.0);
        s4 = (1.0 - cos(theta))/(theta*theta);
    
    // If theta is very close to zero, determine by Taylor series
    if (fabs(theta) < 1e-4)
    {
        
        // Taylor series of s1 = (theta*cos(theta) - sin(theta))/theta^3
        // is -1/3+theta^2/30-theta^4/840+theta^6/45360+O(theta^7)
        s1 = -1.0/3.0 + pow(theta,2.0)/30.0 - pow(theta,4.0)/840.0 +
        pow(theta,6.0)/45360.0;
        
        // Taylor series of s2 = sin(theta)/theta
        // is 1-theta^2/6+theta^4/120+O(theta^6)
        s2 = 1.0 - pow(theta,2.0)/6.0 + pow(theta,4.0)/120.0;
        
        // Taylor series of s3 = (theta*sin(theta) + 2*cos(theta) - 2)/theta^4
        // is -1/12+theta^2/180-theta^4/6720+O(theta^6)
        // Taylor series of s3 = (theta*sin(theta) + 2*cos(theta) - 2)/theta^3
        // is -theta/12+theta^3/180-theta^5/6720+O(theta^6)
        s3 = -theta/12.0 + pow(theta,3.0)/180.0 - pow(theta,5.0)/6720.0;
        
        // Taylor series of s4 = (1 - cos(theta))/theta^2
        // is 1/2-theta^2/24+theta^4/720-theta^6/40320+O(theta^7)
        s4 = 0.5 - pow(theta,2.0)/24.0 + pow(theta,4.0)/720.0
        - pow(theta,6.0)/40320.0;
        
    }
        
//        printf("[dR_dw] s1 = %f\n", s1);
//        printf("[dR_dw] s2 = %f\n", s2);
//        printf("[dR_dw] s3 = %f\n", s3);
//        printf("[dR_dw] s4 = %f\n", s4);

        double w1 = w(0);
        double w2 = w(1);
        double w3 = w(2);
        Eigen::Matrix<double,3,3> M1 =
        s1*w1*wx + s2*G4 + s3*w1*wx2 + s4*(G4*wx + wx*G4);
        Eigen::Matrix<double,3,3> M2 =
        s1*w2*wx + s2*G5 + s3*w2*wx2 + s4*(G5*wx + wx*G5);
        Eigen::Matrix<double,3,3> M3 =
        s1*w3*wx + s2*G6 + s3*w3*wx2 + s4*(G6*wx + wx*G6);

        M << M1(0,0), M2(0,0), M3(0,0),
        M1(1,0), M2(1,0), M3(1,0),
        M1(2,0), M2(2,0), M3(2,0),
        M1(0,1), M2(0,1), M3(0,1),
        M1(1,1), M2(1,1), M3(1,1),
        M1(2,1), M2(2,1), M3(2,1),
        M1(0,2), M2(0,2), M3(0,2),
        M1(1,2), M2(1,2), M3(1,2),
        M1(2,2), M2(2,2), M3(2,2);
        
    }
    
    return M;
    
}
    
#endif


//--------------------------------------------------------------------------
// dR_dw [9 x 3]: MATLAB version
Eigen::Matrix<double,9,3>
DervLie::
dR_dw_matlab(Eigen::Matrix<double,3,1> w)
{
    
    Eigen::Matrix<double,9,3> dR_dw;
    double w1 = w[0];
    double w2 = w[1];
    double w3 = w[2];
    
    double pfw2 = pow(fabs(w1),2.0)+pow(fabs(w2),2.0)+pow(fabs(w3),2.0);
    double pp54 = pow(pfw2,5.0/4.0);
    double spp12 = sin(pow(pfw2,1.0/4.0)*(1.0/2.0));
    double cpp12 = cos(pow(pfw2,1.0/4.0)*(1.0/2.0));
    double pp32 = pow(pfw2,3.0/2.0);
    double spfw2 = sqrt(pfw2);
    double cs12 = cpp12*spp12;
    double ps2 = pow(spp12,2.0);
    
    if ((fabs(w1) < 1e-4) &&
        (fabs(w2) < 1e-4) &&
        (fabs(w3) < 1e-4))
    {
        
        dR_dw.setZero();
        dR_dw(5,0) =  1.0;
        dR_dw(7,0) = -1.0;
        dR_dw(2,1) = -1.0;
        dR_dw(6,1) =  1.0;
        dR_dw(1,2) =  1.0;
        dR_dw(3,2) = -1.0;

        return dR_dw;
    }
    
    dR_dw(0,0) = (w2*w2)*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w1*1.0/pow(pfw2,2.0)*4.0+(w3*w3)*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w1*1.0/pow(pfw2,2.0)*4.0-(w2*w2)*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w1*1.0/pow(pfw2,3.0/2.0)*2.0-(w3*w3)*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w1*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(0,1) = (w2*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*-4.0)/(pfw2)+(w2*w2)*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w2*1.0/pow(pfw2,2.0)*4.0+(w3*w3)*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w2*1.0/pow(pfw2,2.0)*4.0-(w2*w2)*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w2*1.0/pow(pfw2,3.0/2.0)*2.0-(w3*w3)*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w2*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(0,2) = (w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*-4.0)/(pfw2)+(w2*w2)*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w3*1.0/pow(pfw2,2.0)*4.0+(w3*w3)*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w3*1.0/pow(pfw2,2.0)*4.0-(w2*w2)*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w3*1.0/pow(pfw2,3.0/2.0)*2.0-(w3*w3)*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w3*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(1,0) = (w2*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*2.0)/(pfw2)+(w3*pow(cos(sqrt(pfw2)*(1.0/2.0)),2.0)*w1)/(pfw2)-(w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w1)/(pfw2)-w1*w2*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w1*1.0/pow(pfw2,2.0)*4.0-w3*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w1*1.0/pow(pfw2,3.0/2.0)*2.0+w1*w2*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w1*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(1,1) = (w1*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*2.0)/(pfw2)+(w3*pow(cos(sqrt(pfw2)*(1.0/2.0)),2.0)*w2)/(pfw2)-(w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w2)/(pfw2)-w1*w2*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w2*1.0/pow(pfw2,2.0)*4.0-w3*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w2*1.0/pow(pfw2,3.0/2.0)*2.0+w1*w2*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w2*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(1,2) = cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*1.0/sqrt(pfw2)*2.0+(w3*pow(cos(sqrt(pfw2)*(1.0/2.0)),2.0)*w3)/(pfw2)-(w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w3)/(pfw2)-w1*w2*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w3*1.0/pow(pfw2,2.0)*4.0-w3*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w3*1.0/pow(pfw2,3.0/2.0)*2.0+w1*w2*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w3*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(2,0) = (w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*2.0)/(pfw2)-(w2*pow(cos(sqrt(pfw2)*(1.0/2.0)),2.0)*w1)/(pfw2)+(w2*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w1)/(pfw2)-w1*w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w1*1.0/pow(pfw2,2.0)*4.0+w2*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w1*1.0/pow(pfw2,3.0/2.0)*2.0+w1*w3*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w1*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(2,1) = cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*1.0/sqrt(pfw2)*-2.0-(w2*pow(cos(sqrt(pfw2)*(1.0/2.0)),2.0)*w2)/(pfw2)+(w2*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w2)/(pfw2)-w1*w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w2*1.0/pow(pfw2,2.0)*4.0+w2*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w2*1.0/pow(pfw2,3.0/2.0)*2.0+w1*w3*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w2*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(2,2) = (w1*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*2.0)/(pfw2)-(w2*pow(cos(sqrt(pfw2)*(1.0/2.0)),2.0)*w3)/(pfw2)+(w2*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w3)/(pfw2)-w1*w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w3*1.0/pow(pfw2,2.0)*4.0+w2*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w3*1.0/pow(pfw2,3.0/2.0)*2.0+w1*w3*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w3*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(3,0) = (w2*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*2.0)/(pfw2)-(w3*pow(cos(sqrt(pfw2)*(1.0/2.0)),2.0)*w1)/(pfw2)+(w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w1)/(pfw2)-w1*w2*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w1*1.0/pow(pfw2,2.0)*4.0+w3*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w1*1.0/pow(pfw2,3.0/2.0)*2.0+w1*w2*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w1*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(3,1) = (w1*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*2.0)/(pfw2)-(w3*pow(cos(sqrt(pfw2)*(1.0/2.0)),2.0)*w2)/(pfw2)+(w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w2)/(pfw2)-w1*w2*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w2*1.0/pow(pfw2,2.0)*4.0+w3*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w2*1.0/pow(pfw2,3.0/2.0)*2.0+w1*w2*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w2*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(3,2) = cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*1.0/sqrt(pfw2)*-2.0-(w3*pow(cos(sqrt(pfw2)*(1.0/2.0)),2.0)*w3)/(pfw2)+(w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w3)/(pfw2)-w1*w2*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w3*1.0/pow(pfw2,2.0)*4.0+w3*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w3*1.0/pow(pfw2,3.0/2.0)*2.0+w1*w2*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w3*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(4,0) = (w1*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*-4.0)/(pfw2)+(w1*w1)*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w1*1.0/pow(pfw2,2.0)*4.0+(w3*w3)*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w1*1.0/pow(pfw2,2.0)*4.0-(w1*w1)*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w1*1.0/pow(pfw2,3.0/2.0)*2.0-(w3*w3)*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w1*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(4,1) = (w1*w1)*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w2*1.0/pow(pfw2,2.0)*4.0+(w3*w3)*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w2*1.0/pow(pfw2,2.0)*4.0-(w1*w1)*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w2*1.0/pow(pfw2,3.0/2.0)*2.0-(w3*w3)*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w2*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(4,2) = (w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*-4.0)/(pfw2)+(w1*w1)*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w3*1.0/pow(pfw2,2.0)*4.0+(w3*w3)*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w3*1.0/pow(pfw2,2.0)*4.0-(w1*w1)*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w3*1.0/pow(pfw2,3.0/2.0)*2.0-(w3*w3)*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w3*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(5,0) = cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*1.0/sqrt(pfw2)*2.0+(w1*pow(cos(sqrt(pfw2)*(1.0/2.0)),2.0)*w1)/(pfw2)-(w1*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w1)/(pfw2)-w2*w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w1*1.0/pow(pfw2,2.0)*4.0-w1*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w1*1.0/pow(pfw2,3.0/2.0)*2.0+w2*w3*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w1*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(5,1) = (w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*2.0)/(pfw2)+(w1*pow(cos(sqrt(pfw2)*(1.0/2.0)),2.0)*w2)/(pfw2)-(w1*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w2)/(pfw2)-w2*w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w2*1.0/pow(pfw2,2.0)*4.0-w1*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w2*1.0/pow(pfw2,3.0/2.0)*2.0+w2*w3*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w2*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(5,2) = (w2*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*2.0)/(pfw2)+(w1*pow(cos(sqrt(pfw2)*(1.0/2.0)),2.0)*w3)/(pfw2)-(w1*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w3)/(pfw2)-w2*w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w3*1.0/pow(pfw2,2.0)*4.0-w1*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w3*1.0/pow(pfw2,3.0/2.0)*2.0+w2*w3*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w3*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(6,0) = (w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*2.0)/(pfw2)+(w2*pow(cos(sqrt(pfw2)*(1.0/2.0)),2.0)*w1)/(pfw2)-(w2*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w1)/(pfw2)-w1*w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w1*1.0/pow(pfw2,2.0)*4.0-w2*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w1*1.0/pow(pfw2,3.0/2.0)*2.0+w1*w3*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w1*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(6,1) = cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*1.0/sqrt(pfw2)*2.0+(w2*pow(cos(sqrt(pfw2)*(1.0/2.0)),2.0)*w2)/(pfw2)-(w2*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w2)/(pfw2)-w1*w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w2*1.0/pow(pfw2,2.0)*4.0-w2*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w2*1.0/pow(pfw2,3.0/2.0)*2.0+w1*w3*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w2*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(6,2) = (w1*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*2.0)/(pfw2)+(w2*pow(cos(sqrt(pfw2)*(1.0/2.0)),2.0)*w3)/(pfw2)-(w2*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w3)/(pfw2)-w1*w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w3*1.0/pow(pfw2,2.0)*4.0-w2*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w3*1.0/pow(pfw2,3.0/2.0)*2.0+w1*w3*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w3*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(7,0) = cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*1.0/sqrt(pfw2)*-2.0-(w1*pow(cos(sqrt(pfw2)*(1.0/2.0)),2.0)*w1)/(pfw2)+(w1*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w1)/(pfw2)-w2*w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w1*1.0/pow(pfw2,2.0)*4.0+w1*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w1*1.0/pow(pfw2,3.0/2.0)*2.0+w2*w3*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w1*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(7,1) = (w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*2.0)/(pfw2)-(w1*pow(cos(sqrt(pfw2)*(1.0/2.0)),2.0)*w2)/(pfw2)+(w1*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w2)/(pfw2)-w2*w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w2*1.0/pow(pfw2,2.0)*4.0+w1*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w2*1.0/pow(pfw2,3.0/2.0)*2.0+w2*w3*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w2*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(7,2) = (w2*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*2.0)/(pfw2)-(w1*pow(cos(sqrt(pfw2)*(1.0/2.0)),2.0)*w3)/(pfw2)+(w1*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w3)/(pfw2)-w2*w3*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w3*1.0/pow(pfw2,2.0)*4.0+w1*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w3*1.0/pow(pfw2,3.0/2.0)*2.0+w2*w3*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w3*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(8,0) = (w1*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*-4.0)/(pfw2)+(w1*w1)*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w1*1.0/pow(pfw2,2.0)*4.0+(w2*w2)*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w1*1.0/pow(pfw2,2.0)*4.0-(w1*w1)*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w1*1.0/pow(pfw2,3.0/2.0)*2.0-(w2*w2)*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w1*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(8,1) = (w2*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*-4.0)/(pfw2)+(w1*w1)*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w2*1.0/pow(pfw2,2.0)*4.0+(w2*w2)*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w2*1.0/pow(pfw2,2.0)*4.0-(w1*w1)*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w2*1.0/pow(pfw2,3.0/2.0)*2.0-(w2*w2)*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w2*1.0/pow(pfw2,3.0/2.0)*2.0;
    dR_dw(8,2) = (w1*w1)*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w3*1.0/pow(pfw2,2.0)*4.0+(w2*w2)*pow(sin(sqrt(pfw2)*(1.0/2.0)),2.0)*w3*1.0/pow(pfw2,2.0)*4.0-(w1*w1)*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w3*1.0/pow(pfw2,3.0/2.0)*2.0-(w2*w2)*cos(sqrt(pfw2)*(1.0/2.0))*sin(sqrt(pfw2)*(1.0/2.0))*w3*1.0/pow(pfw2,3.0/2.0)*2.0;
    
    return dR_dw;
    
}

//--------------------------------------------------------------------------
// dt_dw [3 x 3]: Verified by Ceres
Eigen::Matrix<double,3,3>
DervLie::
dt_dw(Eigen::Matrix<double,3,1> u,
      Eigen::Matrix<double,3,1> w)
{
    
    Eigen::Matrix<double,3,3> M;
    double theta = sqrt(w.transpose()*w);
    Eigen::Matrix<double,3,3> G4, G5, G6;
    G4 <<
    0,  0,  0,
    0,  0, -1,
    0,  1,  0;
    G5 <<
    0,  0,  1,
    0,  0,  0,
    -1,  0,  0,
    G6 <<
    0, -1,  0,
    1,  0,  0,
    0,  0,  0;
    Eigen::Matrix<double,3,3> wx = skewsym(w);
    Eigen::Matrix<double,3,3> wx2 = wx*wx;
    double s1, s2, s3, s4;
    
    // Use Jacobian at identity if theta = 0
    // Otherwise, use the approximated version for Jacobian at non-zero theta
    if (fabs(theta) < 1e-5)
    {
        //----------------------------------------------------------------------
        // Use Jacobian at identity, Eq (10.15) in J. Blanco's Tutorial
        M = skewsym(-u);
        return M;
        
    }
    else
    {
 

#if 0 // NO LONGER NEEDED SINCE WE HAVE A JACOBIAN AT IDENTITY
      // HOWEVER IT WOULD BE INTERESTING FINDING OUT A BUG IN THE BELOW CODE
      // I SUSPECT s4 times (G4*wx + wx*G4) has an issue in approximation
        //----------------------------------------------------------------------
        // Approximated version of Jacobian at non-zero theta
        double s1, s2, s3, s4;
        // If theta is very close to zero, determine by Taylor series
        if (fabs(theta) < 1e-4)
        {
            
            // Taylor series of s1 = (theta*sin(theta) - 2(1 - cos(theta)))/theta^4
            // is -1/12 + theta^2/180 - theta^4/6720 + theta^6/453600
            //    -theta^8/47900160 + theta^10/7264857600 - theta^12/1494484992000
            s1 = -1.0/12.0 + pow(theta,2.0)/180.0
                 -pow(theta,4.0)/6720.0;
            
            // Taylor series of s2 = (1.0 - cos(theta))/theta^2
            // is 1/2 - theta^2/24 + theta^4/720 - theta^6/40320 + theta^8/3628800
            //    -theta^10/479001600 + theta^12/87178291200
            s2 = 0.5 - theta*theta/24.0 + pow(theta,4.0)/720.0
            - pow(theta,6.0)/40320.0 + pow(theta,8.0)/3628800.0;
            
            // Taylor series of s3 = (3.0*sin(theta) - theta*cos(theta) - 2.0*theta)
            //                       /pow(theta,5.0)
            // is -1/60 + theta^2/1260 - theta^4/60480 + theta^6/4989600 -
            //    theta^8/622702080 + theta^10/108972864000 -
            //    theta^12/25406244864000
            s3 = -1.0/60.0 + theta*theta/1260.0 - pow(theta,4.0)/60480.0
            + pow(theta,6.0)/4989600 - pow(theta,8.0)/622702080.0;
            
            // Taylor series of s4 = (1 - sin(theta)/theta)/theta^2
            // is 1/6 - theta^2/120 + theta^4/5040 - theta^6/362880 +
            //    theta^8/39916800 - theta^10/6227020800 + theta^12/1307674368000;
            s4 = 1.0/6.0 - theta*theta/120.0 + pow(theta,4.0)/5040.0
            - pow(theta,6.0)/362880.0 + pow(theta,8.0)/39916800.0;
            
        }
        else
        {
            
            s1 = (theta*sin(theta) - 2.0*(1.0 - cos(theta)))/pow(theta,4.0);
            s2 = (1.0 - cos(theta))/(theta*theta);
            s3 = (3.0*sin(theta) - theta*cos(theta) - 2.0*theta)/pow(theta,5.0);
            s4 = (1.0 - sin(theta)/theta)/(theta*theta);
            
        }
#endif
        
        s1 = (theta*sin(theta) - 2.0*(1.0 - cos(theta)))/pow(theta,4.0);
        s2 = (1.0 - cos(theta))/(theta*theta);
        s3 = (3.0*sin(theta) - theta*cos(theta) - 2.0*theta)/pow(theta,5.0);
        s4 = (1.0 - sin(theta)/theta)/(theta*theta);
    
#if DEBUG
        printf("[dt_dw] s1 = %f\n", s1);
        printf("[dt_dw] s2 = %f\n", s2);
        printf("[dt_dw] s3 = %f\n", s3);
        printf("[dt_dw] s4 = %f\n", s4);
#endif
        
        double w1 = w(0);
        double w2 = w(1);
        double w3 = w(2);
        Eigen::Matrix<double,3,1> M1 =
        (s1*w1*wx + s2*G4 + s3*w1*wx2 + s4*(G4*wx + wx*G4))*u;
        Eigen::Matrix<double,3,1> M2 =
        (s1*w2*wx + s2*G5 + s3*w2*wx2 + s4*(G5*wx + wx*G5))*u;
        Eigen::Matrix<double,3,1> M3 =
        (s1*w3*wx + s2*G6 + s3*w3*wx2 + s4*(G6*wx + wx*G6))*u;
        
        M <<
        M1(0), M2(0), M3(0),
        M1(1), M2(1), M3(1),
        M1(2), M2(2), M3(2);
    
    }
   
    return M;
    
}
    
//--------------------------------------------------------------------------
// dt_dw (Right exp^{eps} multiplication) [3 x 3]: Verified by Ceres
Eigen::Matrix<double,3,3>
DervLie::
dt_dw_right(Eigen::Matrix<double,3,1> u,
      Eigen::Matrix<double,3,1> w)
{
    
    Eigen::Matrix<double,3,3> M;
    double theta = sqrt(w.transpose()*w);
    Eigen::Matrix<double,3,3> G4, G5, G6;
    G4 <<
    0,  0,  0,
    0,  0, -1,
    0,  1,  0;
    G5 <<
    0,  0,  1,
    0,  0,  0,
    -1,  0,  0,
    G6 <<
    0, -1,  0,
    1,  0,  0,
    0,  0,  0;
    Eigen::Matrix<double,3,3> wx = skewsym(w);
    Eigen::Matrix<double,3,3> wx2 = wx*wx;
    double s1, s2, s3, s4;
    
    // Use Jacobian at identity if theta = 0
    // Otherwise, use the approximated version for Jacobian at non-zero theta
    if (fabs(theta) < 1e-4)
    {
        //----------------------------------------------------------------------
        // Use Jacobian at identity, Eq (10.19) in J. Blanco's Tutorial
        M.setZero();
        return M;
        
    }
    else
    {
        
        
#if 0 // NO LONGER NEEDED SINCE WE HAVE A JACOBIAN AT IDENTITY
        // HOWEVER IT WOULD BE INTERESTING FINDING OUT A BUG IN THE BELOW CODE
        // I SUSPECT s4 times (G4*wx + wx*G4) has an issue in approximation
        //----------------------------------------------------------------------
        // Approximated version of Jacobian at non-zero theta
        double s1, s2, s3, s4;
        // If theta is very close to zero, determine by Taylor series
        if (fabs(theta) < 1e-5)
        {
            
            // Taylor series of s1 = (theta*sin(theta) - 2(1 - cos(theta)))/theta^4
            // is -1/12 + theta^2/180 - theta^4/6720 + theta^6/453600
            //    -theta^8/47900160 + theta^10/7264857600 - theta^12/1494484992000
            s1 = -1.0/12.0 + pow(theta,2.0)/180.0
            -pow(theta,4.0)/6720.0;
            
            // Taylor series of s2 = (1.0 - cos(theta))/theta^2
            // is 1/2 - theta^2/24 + theta^4/720 - theta^6/40320 + theta^8/3628800
            //    -theta^10/479001600 + theta^12/87178291200
            s2 = 0.5 - theta*theta/24.0 + pow(theta,4.0)/720.0
            - pow(theta,6.0)/40320.0 + pow(theta,8.0)/3628800.0;
            
            // Taylor series of s3 = (3.0*sin(theta) - theta*cos(theta) - 2.0*theta)
            //                       /pow(theta,5.0)
            // is -1/60 + theta^2/1260 - theta^4/60480 + theta^6/4989600 -
            //    theta^8/622702080 + theta^10/108972864000 -
            //    theta^12/25406244864000
            s3 = -1.0/60.0 + theta*theta/1260.0 - pow(theta,4.0)/60480.0
            + pow(theta,6.0)/4989600 - pow(theta,8.0)/622702080.0;
            
            // Taylor series of s4 = (1 - sin(theta)/theta)/theta^2
            // is 1/6 - theta^2/120 + theta^4/5040 - theta^6/362880 +
            //    theta^8/39916800 - theta^10/6227020800 + theta^12/1307674368000;
            s4 = 1.0/6.0 - theta*theta/120.0 + pow(theta,4.0)/5040.0
            - pow(theta,6.0)/362880.0 + pow(theta,8.0)/39916800.0;
            
        }
        else
        {
            
            s1 = (theta*sin(theta) - 2.0*(1.0 - cos(theta)))/pow(theta,4.0);
            s2 = (1.0 - cos(theta))/(theta*theta);
            s3 = (3.0*sin(theta) - theta*cos(theta) - 2.0*theta)/pow(theta,5.0);
            s4 = (1.0 - sin(theta)/theta)/(theta*theta);
            
        }
#endif
        
        s1 = (theta*sin(theta) - 2.0*(1.0 - cos(theta)))/pow(theta,4.0);
        s2 = (1.0 - cos(theta))/(theta*theta);
        s3 = (3.0*sin(theta) - theta*cos(theta) - 2.0*theta)/pow(theta,5.0);
        s4 = (1.0 - sin(theta)/theta)/(theta*theta);
        
        double w1 = w(0);
        double w2 = w(1);
        double w3 = w(2);
        Eigen::Matrix<double,3,1> M1 =
        (s1*w1*wx + s2*G4 + s3*w1*wx2 + s4*(G4*wx + wx*G4))*u;
        Eigen::Matrix<double,3,1> M2 =
        (s1*w2*wx + s2*G5 + s3*w2*wx2 + s4*(G5*wx + wx*G5))*u;
        Eigen::Matrix<double,3,1> M3 =
        (s1*w3*wx + s2*G6 + s3*w3*wx2 + s4*(G6*wx + wx*G6))*u;
        
        M <<
        M1(0), M2(0), M3(0),
        M1(1), M2(1), M3(1),
        M1(2), M2(2), M3(2);
        
    }
    
    return M;
    
}

//--------------------------------------------------------------------------
// dt_du [3 x 3] : Verified by Ceres
Eigen::Matrix<double,3,3>
DervLie::
dt_du(Eigen::Matrix<double,3,1> w)
{
    
    Eigen::Matrix<double,3,3> M;
    
    double theta = sqrt(w.transpose()*w);
    Eigen::Matrix<double,3,3> wx = skewsym(w);
    Eigen::Matrix<double,3,3> wx2 = wx*wx;
    
//    // Jacobian at identity
//    if (fabs(theta) < 1e-4)
//    {
//    
//        M = Eigen::MatrixXd::Identity(3,3);
//        
//    }
//    else
    {
        
        double A, B, C;
        // If theta is very close to zero, determine by Taylor series
        if (fabs(theta) < 1e-4)
        {
            
            // Taylor series of A = sin(theta)/theta
            // is 1 - theta^2/6 + theta^4/120 - theta^6/5040 + theta^8/362880
            A = 1.0 - (theta*theta)/6.0 + pow(theta,4.0)/120.0
            - pow(theta,6.0)/5040.0 + pow(theta,8.0)/362880.0;
            
            // Taylor series of B = (1.0 - cos(theta))/(theta^2)
            // is 1/2 - theta^2/24 + theta^4/720 - theta^6/40320
            //    + theta^8/3628800
            B = 0.5 - (theta*theta)/24.0 + pow(theta,4.0)/720.0
            - pow(theta,6.0)/40320.0 + pow(theta,8.0)/3628800.0;
            
            // Taylor series of C = (1.0 - sin(theta)/theta)/(theta*theta)
            // is 1/6 - theta^2/120 + theta^4/5040 - theta^6/362880
            //    + theta^8/39916800
            C = 1.0/6.0 - (theta*theta)/120.0 + pow(theta,4.0)/5040.0
            - pow(theta,6.0)/362880.0 + pow(theta,8.0)/39916800.0;
            
        }
        else
        {
            
            A = sin(theta)/theta;
            B = (1.0 - cos(theta))/(theta*theta);
            C = (1.0 - A)/(theta*theta);
            
        }
        
        Eigen::Matrix<double,3,3> V;
        V = Eigen::MatrixXd::Identity(3,3) + B*wx + C*wx2;
        M = V;
        
    }
    
    return M;
    
}
    
//--------------------------------------------------------------------------
// dt_du (Right exp^{eps} mult) [3 x 3] : Verified by Ceres
Eigen::Matrix<double,3,3>
DervLie::
dt_du_right(Eigen::Matrix<double,3,1> w)
{
    
    Eigen::Matrix<double,3,3> M;
    
    double theta = sqrt(w.transpose()*w);
    Eigen::Matrix<double,3,3> wx = skewsym(w);
    Eigen::Matrix<double,3,3> wx2 = wx*wx;
    
    // Jacobian at identity
    if (fabs(theta) < 1e-4)
    {
        
        M = Eigen::MatrixXd::Identity(3,3);
        
    }
    else
    {
        
//        double A, B, C;
//        //        // If theta is very close to zero, determine by Taylor series
//        //        if (fabs(theta) < 1e-6)
//        //        {
//        //
//        //            // Taylor series of A = sin(theta)/theta
//        //            // is 1 - theta^2/6 + theta^4/120 - theta^6/5040 + theta^8/362880
//        //            A = 1.0 - (theta*theta)/6.0 + pow(theta,4.0)/120.0
//        //            - pow(theta,6.0)/5040.0 + pow(theta,8.0)/362880.0;
//        //
//        //            // Taylor series of B = (1.0 - cos(theta))/(theta^2)
//        //            // is 1/2 - theta^2/24 + theta^4/720 - theta^6/40320
//        //            //    + theta^8/3628800
//        //            B = 0.5 - (theta*theta)/24.0 + pow(theta,4.0)/720.0
//        //            - pow(theta,6.0)/40320.0 + pow(theta,8.0)/3628800.0;
//        //
//        //            // Taylor series of C = (1.0 - sin(theta)/theta)/(theta*theta)
//        //            // is 1/6 - theta^2/120 + theta^4/5040 - theta^6/362880
//        //            //    + theta^8/39916800
//        //            C = 1.0/6.0 - (theta*theta)/120.0 + pow(theta,4.0)/5040.0
//        //            - pow(theta,6.0)/362880.0 + pow(theta,8.0)/39916800.0;
//        //
//        //        }
//        //        else
//        {
//            
//            A = sin(theta)/theta;
//            B = (1.0 - cos(theta))/(theta*theta);
//            C = (1.0 - A)/(theta*theta);
//            
//        }
//        
//        Eigen::Matrix<double,3,3> V;
//        V = Eigen::MatrixXd::Identity(3,3) + B*wx + C*wx2;
//        M = V;
        M = Sophus::SO3d::exp(w).matrix(); // Eq (10.19) in J. Blanco's tutorial
        
    }
    
    return M;
    
}

#if 0 // ORIGINAL
    
// dAB_dA [12 x 12] where A and B are in SE(3)
Eigen::Matrix<double,12,12>
DervLie::dAB_dA(Eigen::Matrix<double,4,4> A, Eigen::Matrix<double,4,4> B)
{

    Eigen::Matrix<double,12,12> M;
    M = Eigen::kroneckerProduct(B.transpose(), Eigen::Matrix3d::Identity());
#if DEBUG
    std::cout << "[dAB_dA] M = " << std::endl;
#endif
    std::cout << M << std::endl;
    
    return M;
    
}
    
#else // via Lie algebra
    
Eigen::Matrix<double,12,12>
DervLie::dAB_dA(Eigen::Matrix<double,4,4> A, Eigen::Matrix<double,4,4> B)
{
    
    Eigen::Matrix<double,12,12> M;
    Sophus::SE3d LieA = Sophus::SE3d(A);
    Sophus::SE3d LieB = Sophus::SE3d(B);
    
    Eigen::Matrix<double,12,6> dAB_dXiA = dBC_dXiB(LieA.log(), LieB.log());
    Eigen::Matrix<double,6,12> dXiA_dA = DervLie::dXiA_dA(A);

#if DEBUG
    std::cout << "[dAB_dA] dAB_dXiA = " << std::endl;
    std::cout << dAB_dXiA << std::endl;
    std::cout << "[dAB_dA] dXiA_dA = " << std::endl;
    std::cout << dXiA_dA << std::endl;
#endif
    
    M = dAB_dXiA*dXiA_dA;
//    std::cout << "[dAB_dA] M = " << std::endl;
//    std::cout << M << std::endl;
    
//// Or another equation
//    
//    Eigen::Matrix<double,12,6> dA_dXiA = DervLie::dA_dXiA(LieA.log());
//    M = DervLie::fast12x12_multJacobToMat(dA_dXiA*dXiA_dA, B);
    
    return M;
    
}
    
#endif
    
// dA_dXiA [12 x 6] : Verified by Ceres
Eigen::Matrix<double,12,6>
DervLie::
dA_dXiA(Eigen::Matrix<double,6,1> Xi)
{
    
    Eigen::Matrix<double,12,6> M;
    Eigen::Matrix<double,9,3> dRdu, dRdw;
    Eigen::Matrix<double,3,3> dtdu, dtdw;
    
    dRdu.setZero(9,3);
    Eigen::Matrix<double,3,1> u, w;
    u << Xi[0], Xi[1], Xi[2];
    w << Xi[3], Xi[4], Xi[5];
    dRdw = dR_dw_matlab(w);
    dtdu = dt_du(w);
    dtdw = dt_dw(u, w);
    
    M << dRdu, dRdw,
    dtdu, dtdw;
    
    return M;
    
}

// dA_dXiA at identity [12 x 6]
Eigen::Matrix<double,12,6>
DervLie::
dA_dXiA_at_Id(Eigen::Matrix<double,6,1> Xi)
{
    
    // Jacobian of exp^{eps} * D, Eq (10.15) in J. Blanco Tutorial
    
    Eigen::Matrix<double,12,6> M;
    Sophus::SE3d LieXi = Sophus::SE3d::exp(Xi);
    Eigen::Matrix<double,4,4> A = LieXi.matrix();
    Eigen::Matrix<double,3,3> R = A.block(0,0,3,3);
    Eigen::Matrix<double,3,1> t = A.block(0,3,3,1);

#if 0 // LEFT_EPS_MULT
    
    M.setZero();
    M.block(0,3,3,3) = -skewsym(R.block(0,0,3,1));
    M.block(3,3,3,3) = -skewsym(R.block(0,1,3,1));
    M.block(6,3,3,3) = -skewsym(R.block(0,2,3,1));
    

    M(9,0) = 1.0;
    M(10,1) = 1.0;
    M(11,2) = 1.0;
    M.block(9,3,3,3) = -skewsym(t);
    
#if 0 // PATCH ON
    
    Eigen::Matrix<double,3,1> w = Xi.tail(3);
    double theta = sqrt(w.transpose()*w);
    if (theta < 1e-5)
    {
        M.block(9,3,3,3) = -skewsym(t*1.5);
    }
    
#endif
    
#else // RIGHT_EPS_MULT
    
    Eigen::Matrix<double,3,1> dc1 = R.block(0,0,3,1);
    Eigen::Matrix<double,3,1> dc2 = R.block(0,1,3,1);
    Eigen::Matrix<double,3,1> dc3 = R.block(0,2,3,1);
    M.setZero();
    M.block(0,4,3,1) = -dc3;
    M.block(0,5,3,1) =  dc2;
    M.block(3,3,3,1) =  dc3;
    M.block(3,5,3,1) = -dc1;
    M.block(6,3,3,1) = -dc2;
    M.block(6,4,3,1) =  dc1;
    M.block(9,0,3,3) = R;
    
    
#if 0 // PATCH ON
    Eigen::Matrix<double,3,1> w = Xi.tail(3);
    double theta = sqrt(w.transpose()*w);
    if (theta < 1e-5)
    {
        
        M.block(9,3,3,3) = -skewsym(t)*0.5;
    
    }
#endif
    
#endif
    
    return M;
    
}

  
//--------------------------------------------------------------------------
Eigen::Matrix<double,12,6>
DervLie::
dA_dXiA_right(Eigen::Matrix<double,6,1> Xi) // Right exp^{eps} multiplication
{
    
    // Use Eq (10.19) in J. Blanco's tutorial
    
    Eigen::Matrix<double,12,6> M;
    Eigen::Matrix<double,9,3> dRdu, dRdw;
    Eigen::Matrix<double,3,3> dtdu, dtdw;
    
    dRdu.setZero(9,3);
    Eigen::Matrix<double,3,1> u, w;
    u << Xi[0], Xi[1], Xi[2];
    w << Xi[3], Xi[4], Xi[5];
    dRdw = dR_dw(w);
    dtdu = dt_du_right(w);
    dtdw = dt_dw_right(u, w);
    
    M << dRdu, dRdw,
    dtdu, dtdw;
    
    return M;
    
}


//--------------------------------------------------------------------------
// dw_dR [3 x 9] : Verified (not exact, homogeneous coords in rotation)
//Eigen::Matrix<double,3,9>
//DervLie::
//dw_dR(Eigen::Matrix<double,3,3> R)
//{
//
//    double r12 = R(0,1);
//    double r13 = R(0,2);
//    double r21 = R(1,0);
//    double r23 = R(1,2);
//    double r31 = R(2,0);
//    double r32 = R(2,1);
//
//    double alpha = (R.trace() - 1.0)/2.0;
//    double theta = acos(alpha);
//
//    printf("alpha = %f\n", alpha);
//    printf("theta = %f\n", theta);
//
//
//    double d02sin0_d0;
//    if ((theta < 1e-6) && (theta > -1e-6))
//    {
//
//        // If theta is close to zero, take the taylor expansion and
//        // approximate it.
//        d02sin0_d0 = theta/6.0 + 7.0*pow(theta,3)/180.0
//        + 31.0*pow(theta,5)/5040.0;
//
//    }
//    else
//    {
//
//        double s20 = sin(2.0*theta);
//        double c20 = cos(2.0*theta);
//        double s0 = sin(theta);
//        d02sin0_d0 = (-theta*s20*s0)/pow(c20 - 1.0,2.0) - s0/(c20 - 1.0);
//
//    }
//
//    double d0_dalpha;
//    if ((alpha < 1 + 1e-6) && (alpha > 1 - 1e-6)) // Near one
//    {
//
//        // Taylor series at alpha = 1
//        d0_dalpha = -1.0/(sqrt(2.0)*sqrt(alpha + 1.0))
//        -sqrt(alpha + 1.0)/(4.0*sqrt(2.0))
//        -3.0*pow((alpha + 1.0), 3.0/2.0)/(32.0*sqrt(2.0))
//        -5.0*pow((alpha + 1.0), 5.0/2.0)/(128.0*sqrt(2.0));
//
//    }
//    else
//    {
//
//        d0_dalpha = -1.0/sqrt(1.0 - alpha*alpha);
//
//    }
//
//    double d02sin0_rii = d02sin0_d0 * d0_dalpha * 0.5;
//    double s = d02sin0_rii;
//    double h;
//    if ((theta < 1e-6) && (theta > -1e-6))
//    {
//
//        // Taylor
//        h = 0.5 + (theta*theta)/12.0
//        + 7.0*pow(theta,4)/720.0
//        + (31.0*pow(theta,6))/30240.0;
//
//    }
//    else
//    {
//
//        h = theta/(2*sin(theta));
//
//    }
//
//    // dw_dR : Derivative of w (orientation) w.r.t. R (rotation)
//    Eigen::Matrix<double,3,9> dw_dR;
//    dw_dR << s*(r32-r23),  0,  0,  0, s*(r32-r23),  h,  0, -h, s*(r32-r23),
//    s*(r13-r31),  0, -h,  0, s*(r13-r31),  0,  h,  0, s*(r13-r31),
//    s*(r21-r12),  h,  0, -h, s*(r21-r12),  0,  0,  0, s*(r21-r12);
//
//    return dw_dR;
//
//}
    
#if 1 // DO_NOT_USE_dw_dR
    
Eigen::Matrix<double,3,9>
DervLie::
dw_dR_ver1(Eigen::Matrix<double,3,3> R)
{
    
    double r12 = R(0,1);
    double r13 = R(0,2);
    double r21 = R(1,0);
    double r23 = R(1,2);
    double r31 = R(2,0);
    double r32 = R(2,1);
    
    double alpha = (R.trace() - 1.0)/2.0;
    double theta = acos(alpha);
    
    //printf("alpha = %f\n", alpha);
    //printf("theta = %f\n", theta);
    
    double delta = 1e-6;
    double d02sin0_d0;
    if ((theta < delta) && (theta > -delta))
    {
        
        // If theta is close to zero, take the taylor expansion and
        // approximate it.
        d02sin0_d0 = theta/6.0 + 7.0*pow(theta,3.0)/180.0
        + 31.0*pow(theta,5.0)/5040.0;
        
    }
    else
    {
        
        //        d02sin0_d0 = (sin(theta)
        //        	- theta*cos(theta))/(2*sin(theta)*sin(theta));
        double s20 = sin(2.0*theta);
        double c20 = cos(2.0*theta);
        double s0 = sin(theta);
        d02sin0_d0 = (-theta*s20*s0)/pow(c20 - 1.0,2.0) - s0/(c20 - 1.0);
        
    }
    //printf("d02sin0_d0 = %f\n", d02sin0_d0);
    
    double d0_dalpha;
    //    if (((alpha < 1.0 + delta) && (alpha > 1.0 - delta)) // Near one
    if ((alpha < -1.0 + delta) && (alpha > -1.0 - delta)) // Near alpha = -1
    {
        
        // Taylor series at alpha = 1
        d0_dalpha = -1.0/(sqrt(2.0)*sqrt(alpha + 1.0))
        -sqrt(alpha + 1.0)/(4.0*sqrt(2.0))
        -3.0*pow((alpha + 1.0), 3.0/2.0)/(32.0*sqrt(2.0))
        -5.0*pow((alpha + 1.0), 5.0/2.0)/(128.0*sqrt(2.0));
        
    }
    else
    {
        
        d0_dalpha = -1.0/sqrt(1.0 - alpha*alpha);
        
    }
    //printf("d0_dalpha = %f\n", d0_dalpha);
    
    double d02sin0_rii = d02sin0_d0 * d0_dalpha * 0.5;
    //printf("d02sin0_rii = %f\n", d02sin0_rii);
    double s = d02sin0_rii;
    
    //    // NEW s derivation
    //    double x = R(0,0);
    //    double y = R(1,1);
    //    double z = R(2,2);
    //    d02sin0_d0 = -1.0/(4.0*(1 - 0.25*pow((x + y + z - 1.0),2.0)))
    //    - ((-x - y - z + 1.0)*acos(0.5*(x + y + z) - 1.0))
    //    / (8.0*pow(1.0 - 0.25*(x + y + z - 1)*(x + y + z -1), 3.0/2.0));
    //    printf("d02sin0_d0 = %f\n", d02sin0_d0);
    
    double h;
    if ((theta < delta) && (theta > -delta))
    {
        
        // Taylor
        h = 0.5 + (theta*theta)/12.0
        + 7.0*pow(theta,4.0)/720.0
        + (31.0*pow(theta,6.0))/30240.0;
        
    }
    else
    {
        
        h = theta/(2.0*sin(theta));
        
    }
    //printf("h = %f\n", h);
    
    double r11 = R(0,0);
    double r22 = R(1,1);
    double r33 = R(2,2);
    double alpha2 = (r11/2.0 + r22/2.0 + r33/2.0 - 0.5);
    double k = 1.0/(4.0*((alpha2*alpha2) - 1.0)) +
    (acos(alpha2)*(alpha2))/
    (4.0*pow((1.0 - (alpha2*alpha2)),(3.0/2.0)));
    //printf("s = %f\n", s);
    //printf("k = %f\n", k);
    
    double m = (theta*cos(theta) - sin(theta)) / (4.0*pow(sin(theta),3.0));
    
    
    if ((theta < delta) && (theta > -delta))
    {
        // m = (theta*cos(theta) - sin(theta)) / (4.0*pow(sin(theta),3.0))
        // Series expansion of m at theta = 0
        // is (-0.0833333 theta^3 - 0.0333333 theta^5 -
        //      0.00793651 theta^7 - 0.00148148 theta^9 -
        //      0.0002405 theta^11 + O(theta^13))/theta^3
        //
        // Taylor series of
        // m = -1/12 - theta^2/30 - theta^4/126 - theta^6/675 -
        //     theta^8/4158 - (691 theta^10)/19348875 - theta^12/200475
        //     -(7234 theta^14)/10854718875 - (43867 theta^16)/509533274250
        //     + O(theta^18)
        m = - 1.0/12.0
            - pow(theta,2.0)/30.0
            - pow(theta,4.0)/126.0
            - pow(theta,6.0)/675.0
            - pow(theta,8.0)/4158.0
            - (691.0*pow(theta,10.0))/19348875.0
            - pow(theta,12.0)/200475.0
            - (7234.0*pow(theta,14.0))/10854718875.0
            - (43867.0*pow(theta,16.0))/509533274250.0;
        
    }
    //printf("m m = %f\n", m);
    //printf("cos(theta) = %f\n", cos(theta));
    
    s = m;
    
    // dw_dR : Derivative of w (orientation) w.r.t. R (rotation)
    Eigen::Matrix<double,3,9> dw_dR;
    dw_dR << s*(r32-r23),  0,  0,  0, s*(r32-r23),  h,  0, -h, s*(r32-r23),
    s*(r13-r31),  0, -h,  0, s*(r13-r31),  0,  h,  0, s*(r13-r31),
    s*(r21-r12),  h,  0, -h, s*(r21-r12),  0,  0,  0, s*(r21-r12);
    
    return dw_dR;
    
}
    
#endif

////--------------------------------------------------------------------------
//// du_dw [3 x 3] : Verified
//Eigen::Matrix<double,3,3>
//DervLie::
//du_dw(Eigen::Matrix<double,3,1> u,
//      Eigen::Matrix<double,3,1> w)
//{
//
//    double w1 = w(0);
//    double w2 = w(1);
//    double w3 = w(2);
//    Eigen::Matrix<double,3,3> wx = skewsym(w);
//    Eigen::Matrix<double,3,3> wx2 = wx*wx;
//    Eigen::Matrix<double,9,3> dwx2_dw;
//
//    dwx2_dw <<   0, -2*w2, -2*w3,
//                w2,    w1,     0,
//                w3,     0,    w1,
//                w2,    w1,     0,
//             -2*w1,     0, -2*w3,
//                 0,    w3,    w2,
//                w3,     0,    w1,
//                 0,    w3,    w2,
//             -2*w1, -2*w2,     0;
//
//    double phi = sqrt(w1*w1 + w2*w2 + w3*w3);
//    Eigen::Matrix<double,9,3> dCwx2_dw;
//    double dD_dphi_phi;
//    double B, C, D;
//    if ((phi < 1e-6) && (phi > -1e-6))
//    {
//
//        B = 0.5 - phi*phi/24.0 + pow(phi,4)/720.0;
//        C = 1.0/6.0 - pow(phi,2)/120.0 + pow(phi,4)/5040.0;
//        dD_dphi_phi = phi/360.0 + pow(phi,3)/7560.0 +
//        pow(phi,5)/201600.0;
//        D = 1.0/12.0 + phi*phi/720.0 + pow(phi,4)/30240.0;
//
//    }
//    else
//    {
//
//        B = sin(phi)/phi;
//        C = (phi - sin(phi))/pow(phi,3);
//        dD_dphi_phi = (4.0 - 4*cos(phi) - phi*sin(phi) - phi*phi)/
//        (2*pow(phi,4)*(cos(phi) - 1.0));
//        double A1 = sin(phi)/phi;
//        double B1 = (1.0 - cos(phi))/(phi*phi);
//        D = (1.0/(phi*phi))*(1.0 - A1/(2.0*B1));
//
//    }
//
//    Eigen::Matrix<double,3,3> V = Eigen::MatrixXd::Identity(3,3)
//    + B*wx + C*wx2;
//    Eigen::Matrix<double,3,1> t = V*u;
//    double t1 = t(0);
//    double t2 = t(1);
//    double t3 = t(2);
//
//    // Equation (53)
//    Eigen::Matrix<double,3,3> M1 = skewsym(-t);
//    Eigen::Matrix<double,3,1> M2vec;
//    M2vec << -t1*w2*w2 + t2*w1*w2 - t1*w3*w3 + t3*w1*w3,
//    -t2*w1*w1 + t1*w2*w1 - t2*w3*w3 + t3*w2*w3,
//    -t3*w1*w1 + t1*w3*w1 - t3*w2*w2 + t2*w3*w2;
//    Eigen::Matrix<double,3,3> M2 = dD_dphi_phi*M2vec*w.transpose();
//    Eigen::Matrix<double,3,3> M3;
//    M3 << D*(t2*w2 + 1*t3*w3), D*(t2*w1 - 2*t1*w2), D*(t3*w1 - 2*t1*w3),
//    D*(t1*w2 - 2*t2*w1), D*(t1*w1 + 1*t3*w3), D*(t3*w2 - 2*t2*w3),
//    D*(t1*w3 - 2*t3*w1), D*(t2*w3 - 2*t3*w2), D*(t1*w1 + 1*t2*w2);
//
//    std::cout << "Eq 53" << std::endl;
//    std::cout << "M1 = " << std::endl << M1 << std::endl;
//    std::cout << "M2 = " << std::endl << M2 << std::endl;
//    std::cout << "M3 = " << std::endl << M3 << std::endl;
//
//    Eigen::Matrix<double,3,3> V_inv;
//    V_inv = Eigen::MatrixXd::Identity(3,3) - 0.5*wx + C*wx2;
//    Eigen::Matrix<double,3,3> du_dw;
////    du_dw = -0.5*M1 + M2 + C*M3 + V_inv*dt_dw(u, w); // What is the last term V_inv*dt_dw for?
//    du_dw = -0.5*M1 + M2 + C*M3;
//
//    return du_dw;
//
//}

//------------------------------------------------------------------------------
// Verified by Ceres
Eigen::Matrix<double,3,3>
DervLie::
du_dw(Eigen::Matrix<double,3,1> u, Eigen::Matrix<double,3,1> w)
{
    
    double w1 = w(0);
    double w2 = w(1);
    double w3 = w(2);
    Eigen::Matrix<double,3,3> wx = DervLie::skewsym(w);
    Eigen::Matrix<double,3,3> wx2 = wx*wx;
    Eigen::Matrix<double,9,3> dwx2_dw;
    
    dwx2_dw <<
         0, -2*w2, -2*w3,
        w2,    w1,     0,
        w3,     0,    w1,
        w2,    w1,     0,
     -2*w1,     0, -2*w3,
         0,    w3,    w2,
        w3,     0,    w1,
         0,    w3,    w2,
     -2*w1, -2*w2,     0;
    
    double phi = sqrt(w1*w1 + w2*w2 + w3*w3);
    Eigen::Matrix<double,9,3> dCwx2_dw;
    
    double A, B, C, G;
    // where A = sin(phi)/phi
    //       B = (1 - cos(phi))/phi^2
    //       C = (1 - A)/(phi^2)
    //       G = (1/phi^2)*(1 - A/(2*B))
    double dG_dphi_phi;
    // where dG_dphi = (4 - 4*cos(phi) - phi*sin(phi) - phi^2)
    //                 /(2*phi^3*(cos(phi) - 1))
    //       dG_dphi_phi = (4 - 4*cos(phi) - phi*sin(phi) - phi^2)
    //                     /(2*phi^4*(cos(phi) - 1))
    
    if ((phi < 1e-4) && (phi > -1e-4))
    {
        
        // Series expansion of A = sin(phi)/phi at phi = 0
        // is 1 - phi^2/6 + phi^4/120 + O(theta^6)
        A = 1.0 - phi*phi/6.0 + pow(phi,4.0)/120.0;
        
        // Series expansion of B = (1 - cos(phi))/phi^2 at phi = 0
        // is 1/2 - phi^2/24 + phi^4/720 + O(phi^6)
        B = 0.5 - phi*phi/24.0 + pow(phi,4.0)/720.0;
        
        // Series expansion of C = (1/phi^2)*(1 - A) at phi = 0
        // is 1/6-phi^2/120+phi^4/5040+O(phi^6)
        C = 1.0/6.0 - phi*phi/120.0 + pow(phi,4.0)/5040.0;
        
        // Series expansion of G = (1/phi^2)*(1 - A/(2*B))
        // is 1/12 + phi^2/720 + phi^4/30240 + O(phi^6)
        G = 1.0/12.0 + phi*phi/720.0 + pow(phi,4.0)/30240.0;
        
        // Series expansion of dC/dphi*(1/phi)
        // is 1/360 + phi^2/7560 + phi^4/201600 + O(phi^6) at phi = 0
        dG_dphi_phi = 1.0/360.0 + pow(phi,2.0)/7560.0 +
        pow(phi,4.0)/201600.0;
        
        
    }
    else
    {
        
        A = sin(phi)/phi;
        B = (1.0 - cos(phi))/(phi*phi);
        C = (1.0/(phi*phi))*(1.0 - A);
        G = (1.0/(phi*phi))*(1.0 - A/(2.0*B));
        
        dG_dphi_phi = (4.0 - 4.0*cos(phi) - phi*sin(phi) - (phi*phi))
        /(2.0*pow(phi,4.0)*(cos(phi) - 1.0));
        
    }
    
#if 0
    printf("phi = %f\n", phi);
    printf("sin(phi)/phi = %f\n", A);
    printf("(1.0 - cos(phi))/(phi*phi) = %f\n", B);
    printf("(1.0/(phi*phi))*(1.0 - A) = %f\n", C);
    printf("dC_dphi_phi = %f\n", dG_dphi_phi);
#endif
    
    Eigen::Matrix<double,3,3> V = Eigen::MatrixXd::Identity(3,3) + B*wx + C*wx2;
    Eigen::Matrix<double,3,1> t = V*u;
    double t1 = t(0);
    double t2 = t(1);
    double t3 = t(2);
    
    // Equation (53)
    Eigen::Matrix<double,3,3> M1 = DervLie::skewsym(-t);
    Eigen::Matrix<double,3,1> M2vec;
    M2vec << -t1*w2*w2 + t2*w1*w2 - t1*w3*w3 + t3*w1*w3,
    -t2*w1*w1 + t1*w2*w1 - t2*w3*w3 + t3*w2*w3,
    -t3*w1*w1 + t1*w3*w1 - t3*w2*w2 + t2*w3*w2;
    Eigen::Matrix<double,3,3> M2 = dG_dphi_phi*M2vec*w.transpose();
    Eigen::Matrix<double,3,3> M3;
    M3 << (t2*w2 + 1*t3*w3), (t2*w1 - 2*t1*w2), (t3*w1 - 2*t1*w3),
    (t1*w2 - 2*t2*w1), (t1*w1 + 1*t3*w3), (t3*w2 - 2*t2*w3),
    (t1*w3 - 2*t3*w1), (t2*w3 - 2*t3*w2), (t1*w1 + 1*t2*w2);
    
    
    Eigen::Matrix<double,3,3> V_inv;
    V_inv = Eigen::MatrixXd::Identity(3,3) - 0.5*wx + G*wx2;
    Eigen::Matrix<double,3,3> du_dw;
    du_dw = -0.5*M1 + M2 + G*M3;
    
    return du_dw;
    
}

#if 1 // DO_NOT_USE_dXiA_dA
//--------------------------------------------------------------------------
// dXiA_dA [6 x 12]: Verified
Eigen::Matrix<double,6,12>
DervLie::
dXiA_dA_ver2(Eigen::Matrix<double,4,4> A)
{
    
    Eigen::Matrix<double,3,3> R = A.block(0,0,3,3);
    Eigen::Matrix<double,3,1> t = A.block(0,3,3,1);
    
    double r12 = R(0,1);
    double r13 = R(0,2);
    double r21 = R(1,0);
    double r23 = R(1,2);
    double r31 = R(2,0);
    double r32 = R(2,1);
    
    double alpha = (R.trace() - 1.0)/2.0;
    double theta = acos(alpha);
    
#if 0//TO_BE_DELETED
    double d02sin0_d0;
    if ((theta < 1e-6) && (theta > -1e-6))
    {
        
        // If theta is close to zero, take the taylor expansion and
        // approximate it.
        d02sin0_d0 = theta/6.0 + 7.0*pow(theta,3.0)/180.0
        + 31.0*pow(theta,5.0)/5040.0;
        
    }
    else
    {
        
        d02sin0_d0 = (sin(theta)
                      - theta*cos(theta))/(2*sin(theta)*sin(theta));
        
    }
    
    double d0_dalpha;
    if ((alpha < 1e-6) && (alpha > -1e-6))
    {
        
        // Taylor series at alpha = 1
        d0_dalpha = -1.0/(sqrt(2.0)*sqrt(alpha + 1.0))
        -sqrt(alpha + 1.0)/(4.0*sqrt(2.0))
        -3.0*pow((alpha + 1.0), 3.0/2.0)/(32.0*sqrt(2.0))
        -5.0*pow((alpha + 1.0), 5.0/2.0)/(128.0*sqrt(2.0));
        
    }
    else
    {
        
        d0_dalpha = -1.0/sqrt(1.0 - alpha*alpha);
        
    }
#endif
    
    
    double h;
    if ((theta < 1e-6) && (theta > -1e-6))
    {
        
        // Taylor
        h = 0.5 + theta*theta/12.0
        + 7.0*pow(theta,4.0)/720.0
        + 31.0*pow(theta,6.0)/30240.0;
        
    }
    else
    {
        
        h = theta/(2*sin(theta));
        
    }
    
    // dw_dR : Derivative of w (orientation) w.r.t. R (rotation)
    Eigen::Matrix<double,3,9> dw_dR_ = dw_dR_ver1(R);
    
    // dw_dt : Derivative of w (orientation) w.r.t. t (translation)
    Eigen::Matrix<double,3,3> dw_dt;
    dw_dt.setZero(3,3);
    
    Eigen::Matrix<double,3,1> w;
    w << h*(r32 - r23), h*(r13 - r31), h*(r21 - r12);
    
    Eigen::Matrix<double,3,3> wx = skewsym(w);
    Eigen::Matrix<double,3,3> wx2 = wx*wx;
    double A1, B, D;
    if (fabs(theta) < 1e-6)
    {
        A1 = 1.0 - (theta*theta)/6.0 + pow(theta,4.0)/120.0;
        B = 0.5 - (theta*theta)/24.0 + pow(theta,4.0)/720.0;
        D = 1.0/12.0 + (theta*theta)/720.0 + pow(theta,4.0)/30240.0;
    }
    else
    {
        
        A1 = sin(theta)/theta;
        B = (1.0 - cos(theta))/(theta*theta);
        D = (1.0/(theta*theta))*(1.0 - A1/(2*B));
    }
    
    Eigen::Matrix<double,3,3> V_inv;
    V_inv = Eigen::MatrixXd::Identity(3,3) - 0.5*wx + D*wx2;
    
    Eigen::Matrix<double,3,1> u = V_inv*t;
    Eigen::Matrix<double,3,3> du_dw_ = du_dw(u, w);
    
#if 0//DEBUG
    std::cout << "du_dw = " << std::endl;
    std::cout << du_dw_ << std::endl;
#endif
    // du_dR : Derivative of u (translation) w.r.t. R (rotation)
    Eigen::Matrix<double,3,9> du_dR = du_dw_*dw_dR_;
    
    // du_dt : Derivative of u (translation) w.r.t. t (translation)
    Eigen::Matrix<double,3,3> du_dt;
    du_dt = V_inv;
    
    //------------------------------------------------------------------
    // Set dXiA_dA from du_dR, du_dt, dw_dR and dw_dt,
    // where A = [R t].
    //------------------------------------------------------------------
    Eigen::Matrix<double,6,12> M;
    M << du_dR, du_dt,
    dw_dR_, dw_dt;
    
    return M;
    
}
#endif

//--------------------------------------------------------------------------
// Verified working (deprecated, use Jacob_x_Mat instead)
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>
DervLie::
Jacob_times_Vec(const Eigen::Ref<const
                Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> >& J,
                const Eigen::Matrix<double,Eigen::Dynamic,1> vec)
{
    
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> J_out;
    J_out.resize(J.rows()/vec.rows(), J.cols());
#if 0//DEBUG
    std::cout << "J_out size = " << J_out.rows()
    << "x" <<  J_out.cols() << std::endl;
#endif
    for (int j=0; j<J.cols(); j++)
    {
        
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> tmp_J;
        tmp_J.resize(J.rows()/vec.rows(),vec.rows());
        
        // Copy a column of J to a matrix
        int i=0;
        for (int q=0; q<tmp_J.cols(); q++)
        {
            
            for (int p=0; p<tmp_J.rows(); p++)
            {
                
                tmp_J(p,q) = J(i,j);
                i++;
                
            }
            
        }
        
#if 0 //DEBUG
        std::cout << tmp_J.rows() << "x" << tmp_J.cols() << " ";
        std::cout << "tmp_J = " << std::endl << tmp_J << std::endl;
        std::cout << vec.rows() << "x" << vec.cols() << " ";
        std::cout << "vec = " << std::endl << vec << std::endl;
#endif
        // Compute the right multilication
        tmp_J = tmp_J*vec;
#if 0 //DEBUG
        std::cout << tmp_J.rows() << "x" << tmp_J.cols() << " ";
        std::cout << "tmp_J = " << std::endl << tmp_J << std::endl;
        std::cout << std::endl;
#endif
        // Copy back the result to the jacobian column
        i=0;
        for (int q=0; q<tmp_J.cols(); q++)
        {
            
            for (int p=0; p<tmp_J.rows(); p++)
            {
                
                J_out(i,j) = tmp_J(p,q);
                i++;
                
            }
            
        }
        
    }
    
    return J_out;
    
}

//--------------------------------------------------------------------------
// d(BC)_dXiB [12 x 6]: Verified by Ceres
#if 0 // ORIGINAL
Eigen::Matrix<double,12,6>
DervLie::
dBC_dXiB(Eigen::Matrix<double,6,1> XiB,
         Eigen::Matrix<double,6,1> XiC)
{
    
    Eigen::Matrix<double,12,6> M;
    
    Sophus::SE3d LieC = Sophus::SE3d::exp(XiC);
#if 0// NOT_USED
    Sophus::SE3d LieB = Sophus::SE3d::exp(XiB);
    Eigen::Matrix<double,4,4> B = LieB.matrix();
#endif
    Eigen::Matrix<double,4,4> C = LieC.matrix();
    
#if 0// NOTUSED
    Eigen::Matrix<double,3,3> RB = B.block(0,0,3,3);
    Eigen::Matrix<double,3,1> tB = B.block(0,3,3,1);
#endif
    
    Eigen::Matrix<double,3,3> RC = C.block(0,0,3,3);
    Eigen::Matrix<double,3,1> tC = C.block(0,3,3,1);
    
    Eigen::Matrix<double,12,6> dB_dXiB = dA_dXiA(XiB);
    
    // Set dRA_duB
    Eigen::Matrix<double,9,3> dRB_duB = dB_dXiB.block(0,0,9,3);
    Eigen::Matrix<double,9,3> dRA_duB =
    fast9x3_Jacob9x3_x_Mat3x3(dRB_duB, RC);
#if 0 //DEBUG
    std::cout << "dRA_duB = " << std::endl;
    std::cout << dRA_duB << std::endl;
#endif
    // Set dRA_dwB
    Eigen::Matrix<double,9,3> dRB_dwB = dB_dXiB.block(0,3,9,3);
    Eigen::Matrix<double,9,3> dRA_dwB =
    fast9x3_Jacob9x3_x_Mat3x3(dRB_dwB, RC);
#if 0//DEBUG
    std::cout << "dRA_dwB = " << std::endl;
    std::cout << dRA_dwB << std::endl;
#endif
    // set dtA_duB
    Eigen::Matrix<double,3,3> dtB_duB = dB_dXiB.block(9,0,3,3);
    Eigen::Matrix<double,3,3> dtA_duB =
    fast3x3_Jacob9x3_x_Vec3x1(dRB_duB, tC) + dtB_duB;
#if 0//DEBUG
    std::cout << "dtA_duB = " << std::endl;
    std::cout << dtA_duB << std::endl;
#endif
    // Set dtA_dwB
    Eigen::Matrix<double,3,3> dtB_dwB = dB_dXiB.block(9,3,3,3);
    Eigen::Matrix<double,3,3> dtA_dwB =
    fast3x3_Jacob9x3_x_Vec3x1(dRB_dwB, tC) + dtB_dwB;
    
    M << dRA_duB, dRA_dwB,
    dtA_duB, dtA_dwB;
    
    return M;
    
}
#else // SHORT form
    
Eigen::Matrix<double,12,6>
DervLie::
dBC_dXiB(Eigen::Matrix<double,6,1> XiB,
         Eigen::Matrix<double,6,1> XiC)

{

    Eigen::Matrix<double,4,4> C = Sophus::SE3d::exp(XiC).matrix();
    Eigen::Matrix<double,12,6> dB_dXiB = DervLie::dA_dXiA(XiB);
    
    return fast12x6_multJacobToMat(dB_dXiB, C);
    
}

#endif

//--------------------------------------------------------------------------
// d(BC)_dXiC [12 x 6]: Verified
#if 0 // ORIGINAL
Eigen::Matrix<double,12,6>
DervLie::
dBC_dXiC(Eigen::Matrix<double,6,1> XiB,
         Eigen::Matrix<double,6,1> XiC)
{
    
    Eigen::Matrix<double,12,6> M;
    
    Sophus::SE3d LieB = Sophus::SE3d::exp(XiB);
    Eigen::Matrix<double,4,4> B = LieB.matrix();
#if 0//NOT USED
    Sophus::SE3d LieC = Sophus::SE3d::exp(XiC);
    Eigen::Matrix<double,4,4> C = LieC.matrix();
#endif
    
    Eigen::Matrix<double,3,3> RB = B.block(0,0,3,3);
#if 0// NOT_USED
    Eigen::Matrix<double,3,1> tB = B.block(0,3,3,1);
    Eigen::Matrix<double,3,3> RC = C.block(0,0,3,3);
    Eigen::Matrix<double,3,1> tC = C.block(0,3,3,1);
#endif
    Eigen::Matrix<double,12,6> dC_dXiC = dA_dXiA(XiC);
    
    // Set dRA_duC
    Eigen::Matrix<double,9,3> dRC_duC = dC_dXiC.block(0,0,9,3);
    Eigen::Matrix<double,9,3> dRA_duC =
    fast9x3_Mat3x3_x_Jacob9x3(RB, dRC_duC);
#if 0//DEBUG
    std::cout << "dRA_duC = " << std::endl;
    std::cout << dRA_duC << std::endl;
#endif
    // Set dRA_dwC
    Eigen::Matrix<double,9,3> dRC_dwC = dC_dXiC.block(0,3,9,3);
    Eigen::Matrix<double,9,3> dRA_dwC =
    fast9x3_Mat3x3_x_Jacob9x3(RB, dRC_dwC);
#if 0//DEBUG
    std::cout << "dRA_dwC = " << std::endl;
    std::cout << dRA_dwC << std::endl;
#endif
    // set dtA_duC
    Eigen::Matrix<double,3,3> dtC_duC = dC_dXiC.block(9,0,3,3);
    Eigen::Matrix<double,3,3> dtA_duC =
    fast3x3_Mat3x3_x_Jacob3x3(RB, dtC_duC);
#if 0//DEBUG
    std::cout << "dtA_duC = " << std::endl;
    std::cout << dtA_duC << std::endl;
#endif
    // Set dtA_dwC
    Eigen::Matrix<double,3,3> dtC_dwC = dC_dXiC.block(9,3,3,3);
    Eigen::Matrix<double,3,3> dtA_dwC =
    fast3x3_Mat3x3_x_Jacob3x3(RB, dtC_dwC);
    
    M << dRA_duC, dRA_dwC,
    dtA_duC, dtA_dwC;
    
    return M;
    
}
#else // SHORT form

Eigen::Matrix<double,12,6>
DervLie::
dBC_dXiC(Eigen::Matrix<double,6,1> XiB,
         Eigen::Matrix<double,6,1> XiC)
{
    
    Eigen::Matrix<double,4,4> B = Sophus::SE3d::exp(XiB).matrix();
    Eigen::Matrix<double,12,6> dC_dXiC = DervLie::dA_dXiA(XiC);
    
    return fast12x6_multMatToJacob(B, dC_dXiC);
    
}

#endif

//--------------------------------------------------------------------------
// This is a special multiplication of a matrix to a Jacobian matrix, which
// is reshaped in order to have each column of the Jacobian matrix
// may be correspond to the size of the matrix left multiplying
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>
DervLie::
Mat_x_Jacob(const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> M,
            const Eigen::Ref<const
            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> >& J)
{
    
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> J_out;
    J_out.resize(M.rows()*(J.rows()/M.cols()) , J.cols());
#if 0//DEBUG
    std::cout << "J_out size = " << J_out.rows()
    << "x" <<  J_out.cols() << std::endl;
#endif
    // Determine the size of a matrix which is reshaped from a column
    // of the Jacobian matrix. This matrix will be multiplied with the
    // input matrix M
    int new_row_tmpJ = M.cols();
    int new_col_tmpJ = J.rows()/M.cols();
    
    for (int j=0; j<J.cols(); j++)
    {
        
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> tmp_J;
        tmp_J.resize(new_row_tmpJ, new_col_tmpJ);
        
        // Copy a column of J to a matrix
        int i=0;
        for (int q=0; q<tmp_J.cols(); q++)
        {
            
            for (int p=0; p<tmp_J.rows(); p++)
            {
                
                tmp_J(p,q) = J(i,j);
                i++;
                
            }
            
        }
        
#if 0 //DEBUG
        std::cout << tmp_J.rows() << "x" << tmp_J.cols() << " ";
        std::cout << "tmp_J = " << std::endl << tmp_J << std::endl;
        std::cout << vec.rows() << "x" << vec.cols() << " ";
        std::cout << "vec = " << std::endl << vec << std::endl;
#endif
        // Compute the multilication
        tmp_J = M*tmp_J;
#if 0 //DEBUG
        std::cout << tmp_J.rows() << "x" << tmp_J.cols() << " ";
        std::cout << "tmp_J = " << std::endl << tmp_J << std::endl;
        std::cout << std::endl;
#endif
        // Copy back the result to the jacobian column
        i=0;
        for (int q=0; q<tmp_J.cols(); q++)
        {
            
            for (int p=0; p<tmp_J.rows(); p++)
            {
                
                J_out(i,j) = tmp_J(p,q);
                i++;
                
            }
            
        }
        
    }
    
    return J_out;
    
}

//--------------------------------------------------------------------------
// This is a special multiplication of a matrix to a Jacobian matrix, which
// is reshaped in order to have each column of the Jacobian matrix
// may be correspond to the size of the matrix left multiplying
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>
DervLie::
Jacob_x_Mat(const Eigen::Ref<const
            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> >& J,
            const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> M)
{
    
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> J_out;
    J_out.resize((J.rows()/M.rows())*M.cols() , J.cols());
#if 0//DEBUG
    std::cout << "J_out size = " << J_out.rows()
    << "x" <<  J_out.cols() << std::endl;
#endif
    // Determine the size of a matrix which is reshaped from a column
    // of the Jacobian matrix. This matrix will be multiplied with the
    // input matrix M
    int new_row_tmpJ = J.rows()/M.rows();
    int new_col_tmpJ = M.rows();
    
    for (int j=0; j<J.cols(); j++)
    {
        
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> tmp_J;
        tmp_J.resize(new_row_tmpJ, new_col_tmpJ);
        
        // Copy a column of J to a matrix
        int i=0;
        for (int q=0; q<tmp_J.cols(); q++)
        {
            
            for (int p=0; p<tmp_J.rows(); p++)
            {
                
                tmp_J(p,q) = J(i,j);
                i++;
                
            }
            
        }
        
#if 0 //DEBUG
        std::cout << tmp_J.rows() << "x" << tmp_J.cols() << " ";
        std::cout << "tmp_J = " << std::endl << tmp_J << std::endl;
        std::cout << vec.rows() << "x" << vec.cols() << " ";
        std::cout << "vec = " << std::endl << vec << std::endl;
#endif
        // Compute the multilication
        tmp_J = tmp_J*M;
#if 0 //DEBUG
        std::cout << tmp_J.rows() << "x" << tmp_J.cols() << " ";
        std::cout << "tmp_J = " << std::endl << tmp_J << std::endl;
        std::cout << std::endl;
#endif
        // Copy back the result to the jacobian column
        i=0;
        for (int q=0; q<tmp_J.cols(); q++)
        {
            
            for (int p=0; p<tmp_J.rows(); p++)
            {
                
                J_out(i,j) = tmp_J(p,q);
                i++;
                
            }
            
        }
        
    }
    
    return J_out;
    
}


//--------------------------------------------------------------------------
// Valid
Eigen::Matrix<double,3,3>
DervLie::
skewsym(Eigen::Matrix<double,3,1> x)
{
    
    Eigen::Matrix<double,3,3> M;
    M <<     0, -x(2),  x(1),
    x(2),    0,  -x(0),
    -x(1),  x(0),    0;
    
    return M;
}

//--------------------------------------------------------------------------
// Verified by Ceres
Eigen::Matrix<double,9,3>
DervLie::
dRinv_dw(Eigen::Matrix<double,3,1> w)
{
    
    Eigen::Matrix<double,9,3> M;
    
    double theta = sqrt(w.transpose()*w);
    
    Eigen::Matrix<double,3,3> G4, G5, G6;
    G4 <<   0,  0,  0,
    0,  0, -1,
    0,  1,  0;
    G5 <<   0,  0,  1,
    0,  0,  0,
    -1,  0,  0,
    G6 <<   0, -1,  0,
    1,  0,  0,
    0,  0,  0;
    
    Eigen::Matrix<double,3,3> wx = skewsym(w);
    Eigen::Matrix<double,3,3> wx2 = wx*wx;
    
    // Jacobian at identity
    if (fabs(theta) < 1e-6)
    {
        
        M << -skewsym(Eigen::Matrix<double,3,1>(1,0,0)),
             -skewsym(Eigen::Matrix<double,3,1>(0,1,0)),
             -skewsym(Eigen::Matrix<double,3,1>(0,0,1));
        
    }
    else // Jacobian at non-identity
    {
        
        double s1, s2, s3, s4;
//        // If theta is very close to zero, determine by Taylor series
//        if (fabs(theta) < 1e-6)
//        {
//            
//            s1 = -1.0/3.0 + pow(theta,2.0)/30.0 + pow(theta,4.0)/840.0;
//            s2 = 1 - pow(theta,2.0)/6.0 + pow(theta,4.0)/120.0;
//            s3 = -1.0/12.0 + pow(theta,2.0)/180.0 + pow(theta,4.0)/6720.0;
//            s4 = 0.5 - pow(theta,2.0)/24.0 + pow(theta,4.0)/720.0;
//            
//        }
//        else
//        {
        
            s1 = (theta*cos(theta) - sin(theta))/pow(theta,3.0);
            s2 = sin(theta)/theta;
            s3 = (theta*sin(theta) + 2*cos(theta) - 2.0)/pow(theta,4.0);
            s4 = (1.0 - cos(theta))/(theta*theta);
            
//        }
        
        double w1 = w(0);
        double w2 = w(1);
        double w3 = w(2);
        Eigen::Matrix<double,3,3> M1 =
        -s1*w1*wx - s2*G4 + s3*w1*wx2 + s4*(G4*wx + wx*G4);
        Eigen::Matrix<double,3,3> M2 =
        -s1*w2*wx - s2*G5 + s3*w2*wx2 + s4*(G5*wx + wx*G5);
        Eigen::Matrix<double,3,3> M3 =
        -s1*w3*wx - s2*G6 + s3*w3*wx2 + s4*(G6*wx + wx*G6);
        
        M << M1(0,0), M2(0,0), M3(0,0),
        M1(1,0), M2(1,0), M3(1,0),
        M1(2,0), M2(2,0), M3(2,0),
        M1(0,1), M2(0,1), M3(0,1),
        M1(1,1), M2(1,1), M3(1,1),
        M1(2,1), M2(2,1), M3(2,1),
        M1(0,2), M2(0,2), M3(0,2),
        M1(1,2), M2(1,2), M3(1,2),
        M1(2,2), M2(2,2), M3(2,2);

    }
    
    return M;
    
}

//--------------------------------------------------------------------------
// Verified by Ceres
#if 0 // ORIGINAL
Eigen::Matrix<double,12,6>
DervLie::
dAinv_dXiA(Eigen::Matrix<double,6,1> Xi)
{
    
    Eigen::Matrix<double,12,6> M;
    Eigen::Matrix<double,9,3> dRinvdu, dRinvdw;
    Eigen::Matrix<double,3,3> dtdu, dtdw;
    
    dRinvdu.setZero(9,3);
    Eigen::Matrix<double,3,1> u, w;
    u << Xi[0], Xi[1], Xi[2];
    w << Xi[3], Xi[4], Xi[5];
    dRinvdw = dRinv_dw(w);
    dtdu = dt_du(w);
    dtdw = dt_dw(u, w);
    
    Sophus::SE3d LieXi = Sophus::SE3d::exp(Xi);
    Eigen::Matrix<double,4,4> A = LieXi.matrix();
    Eigen::Matrix<double,3,3> Rinv = A.block(0,0,3,3).transpose();
    Eigen::Matrix<double,3,1> t = A.block(0,3,3,1);
    
    M <<  dRinvdu, dRinvdw,
    -fast3x3_Jacob9x3_x_Vec3x1(dRinvdu, t) -
    fast3x3_Mat3x3_x_Jacob3x3(Rinv, dtdu),
    -fast3x3_Jacob9x3_x_Vec3x1(dRinvdw, t) -
    fast3x3_Mat3x3_x_Jacob3x3(Rinv, dtdw);
    
    return M;
    
}
#elseif 1 // SECOND ORINGIAL
Eigen::Matrix<double,12,6>
DervLie::
dAinv_dXiA(Eigen::Matrix<double,6,1> Xi)
{
    
    Eigen::Matrix<double,12,6> M;
    Sophus::SE3d LieXi = Sophus::SE3d::exp(Xi);
    Eigen::Matrix<double,4,4> A = LieXi.matrix();
    Eigen::Matrix<double,3,3> R = A.block(0,0,3,3);
    Eigen::Matrix<double,3,1> t = A.block(0,3,3,1);
    
    // d(A^{-1})/dA is
    //
    //   [ d(R^{-1})_/dR       d(R^{-1})/dt ]
    //   [ d(-R^{-1}*t)/dR  d(-R^{-1}*t)/dt ]
    //
    // becoming
    //
    //   [ E_{9x9}    0_{9x3} ]
    //   [ -E_{9x9}*t -R^{-1} ]
    //
    // where E_{9x9} is a sparse jacobian
    //
    // E_{9x9} =
    //
    //   [1 0 0 0 0 0 0 0 0 ]
    //   [0 0 0 1 0 0 0 0 0 ]
    //   [0 0 0 0 0 0 1 0 0 ]
    //   [0 1 0 0 0 0 0 0 0 ]
    //   [0 0 0 0 1 0 0 0 0 ]
    //   [0 0 0 0 0 0 0 1 0 ]
    //   [0 0 1 0 0 0 0 0 0 ]
    //   [0 0 0 0 0 1 0 0 0 ]
    //   [0 0 0 0 0 0 0 0 1 ]
    //
    Eigen::Matrix<double,9,9> E;
    E.setZero();
    E(0,0) = 1.0;
    E(1,3) = 1.0;
    E(3,1) = 1.0;
    E(2,6) = 1.0;
    E(6,2) = 1.0;
    E(4,4) = 1.0;
    E(5,7) = 1.0;
    E(7,5) = 1.0;
    E(8,8) = 1.0;
    Eigen::Matrix<double,12,12> dAinv_dA;
    Eigen::Matrix<double,9,9> dRinv_dR = E;
    Eigen::Matrix<double,9,3> dRinv_dt = Eigen::MatrixXd::Zero(9,3);
    Eigen::Matrix<double,9,9> dnegRinv_dR = -E;
    Eigen::Matrix<double,3,9> dnegRinv_t_dR
        = fast3x9_Jacob9x9_x_Mat3x1(dnegRinv_dR, t);
    Eigen::Matrix<double,3,3> dnegRinv_t_dt = -R.transpose();
    
    dAinv_dA <<
    dRinv_dR, dRinv_dt,
    dnegRinv_t_dR, dnegRinv_t_dt;
    
    // Out matrix M is
    // d(A^{-1})/dA * dA/dXiA
    M = dAinv_dA*dA_dXiA(Xi);
    
    return M;
    
}
#else // MORE COMPACT EXPRESSION
    
Eigen::Matrix<double,12,6>
DervLie::
dAinv_dXiA_right(Eigen::Matrix<double,6,1> Xi)
{
    
    Eigen::Matrix<double,12,6> dAinv_dXiAinv = DervLie::dA_dXiA_right(-Xi);
    Eigen::Matrix<double,6,6> dXiAinv_dXiA =
        -Eigen::Matrix<double,6,6>::Identity();
    return dAinv_dXiAinv*dXiAinv_dXiA;
    
}
    
Eigen::Matrix<double,12,6>
DervLie::
dAinv_dXiA(Eigen::Matrix<double,6,1> Xi)
{
    
    Eigen::Matrix<double,12,6> dAinv_dXiAinv = DervLie::dA_dXiA(-Xi);
    Eigen::Matrix<double,6,6> dXiAinv_dXiA =
    -Eigen::Matrix<double,6,6>::Identity();
    return dAinv_dXiAinv*dXiAinv_dXiA;

}
    
Eigen::Matrix<double,12,6>
DervLie::
dAinv_dXiA_at_Id(Eigen::Matrix<double,6,1> Xi)
{
    
    Eigen::Matrix<double,12,6> dAinv_dXiAinv = DervLie::dA_dXiA_at_Id(-Xi);
    Eigen::Matrix<double,6,6> dXiAinv_dXiA =
    -Eigen::Matrix<double,6,6>::Identity();
    return dAinv_dXiAinv*dXiAinv_dXiA;
    
}
    
#endif

//--------------------------------------------------------------------------
// Verified by Ceres
#if 0 // ORIGINAL
Eigen::Matrix<double,12,6>
DervLie::
dBinvC_dXiB(Eigen::Matrix<double,6,1> XiB,
            Eigen::Matrix<double,6,1> XiC)
{
    
    Eigen::Matrix<double,12,6> M;
    
    Sophus::SE3d LieC = Sophus::SE3d::exp(XiC);
#if 0//NOT_USED
    Sophus::SE3d LieB = Sophus::SE3d::exp(XiB);
    Eigen::Matrix<double,4,4> B = LieB.matrix();
#endif
    Eigen::Matrix<double,4,4> C = LieC.matrix();
    
#if 0//NOT USED
    Eigen::Matrix<double,3,3> RB = B.block(0,0,3,3);
    Eigen::Matrix<double,3,1> tB = B.block(0,3,3,1);
#endif
    Eigen::Matrix<double,3,3> RC = C.block(0,0,3,3);
    Eigen::Matrix<double,3,1> tC = C.block(0,3,3,1);
    
    Eigen::Matrix<double,12,6> dBinv_dXiB = dAinv_dXiA(XiB);
    
    // Set dRA_duB = d(RBinv_RC)/duB
    Eigen::Matrix<double,9,3> dRB_duB = dBinv_dXiB.block(0,0,9,3);
    Eigen::Matrix<double,9,3> dRA_duB =
    fast9x3_Jacob9x3_x_Mat3x3(dRB_duB, RC);
#if 0//DEBUG
    std::cout << "dRA_duB = " << std::endl;
    std::cout << dRA_duB << std::endl;
#endif
    // Set dRA_dwB
    Eigen::Matrix<double,9,3> dRB_dwB = dBinv_dXiB.block(0,3,9,3);
    Eigen::Matrix<double,9,3> dRA_dwB =
    fast9x3_Jacob9x3_x_Mat3x3(dRB_dwB, RC);
#if 0//DEBUG
    std::cout << "dRA_dwB = " << std::endl;
    std::cout << dRA_dwB << std::endl;
#endif
    // set dtA_duB
    Eigen::Matrix<double,3,3> dtB_duB = dBinv_dXiB.block(9,0,3,3);
    Eigen::Matrix<double,3,3> dtA_duB =
    fast3x3_Jacob9x3_x_Vec3x1(dRB_duB, tC) + dtB_duB;
#if 0//DEBUG
    std::cout << "dtA_duB = " << std::endl;
    std::cout << dtA_duB << std::endl;
#endif
    // Set dtA_dwB
    Eigen::Matrix<double,3,3> dtB_dwB = dBinv_dXiB.block(9,3,3,3);
    Eigen::Matrix<double,3,3> dtA_dwB =
    fast3x3_Jacob9x3_x_Vec3x1(dRB_dwB, tC) + dtB_dwB;
    
    M << dRA_duB, dRA_dwB,
    dtA_duB, dtA_dwB;
    
    return M;
    
}
#else // SHORT form
    
Eigen::Matrix<double,12,6>
DervLie::
dBinvC_dXiB(Eigen::Matrix<double,6,1> XiB,
            Eigen::Matrix<double,6,1> XiC)
{
    
    Eigen::Matrix<double,4,4> C = Sophus::SE3d::exp(XiC).matrix();
    Eigen::Matrix<double,12,6> dBinv_dXiB = DervLie::dAinv_dXiA(XiB);
    
    return fast12x6_multJacobToMat(dBinv_dXiB, C);
    
}

Eigen::Matrix<double,12,6>
DervLie::
dBinvC_dXiB_at_Id(Eigen::Matrix<double,6,1> XiB,
            Eigen::Matrix<double,6,1> XiC)
{
    
    Eigen::Matrix<double,4,4> C = Sophus::SE3d::exp(XiC).matrix();
    Eigen::Matrix<double,4,4> B = Sophus::SE3d::exp(XiB).matrix();

    Eigen::Matrix<double,12,12> dAinvC_dAinv
        = Eigen::kroneckerProduct(C.transpose(), Eigen::Matrix3d::Identity());
    
    Eigen::Matrix<double,12,12> dAinv_dA = DervLie::dAinv_dA(B);
    
#if 0 // LEFT_EPS_MULT
    
    Eigen::Matrix<double,12,6> depsB_deps = DervLie::dEps_plus_D_dXiEps(B);
    return dAinvC_dAinv*dAinv_dA*depsB_deps;
    
#else // RIGHT_EPS_MULT

    Eigen::Matrix<double,12,6> dBeps_deps = DervLie::dD_plus_Eps_dXiEps(B);
    return dAinvC_dAinv*dAinv_dA*dBeps_deps;

#endif
    
    
}

Eigen::Matrix<double,12,6>
DervLie::
dAinv_negEps_D_deps_at_Id(Eigen::Matrix<double,6,1> XiA,
                  Eigen::Matrix<double,6,1> XiD)
{
    
    // Variant of Equation (10.28) in Blanco's tutorial
    // Note: We use -XiA for A inverse. Then, exp(-XiA) becomes -exp(XiA).
    //       d(A^{-1}*exp(-eps)*D)/d(eps) = -d(A^{-1}*exp(eps)*D)/d(eps)
    //       This is essentially same as calling dBinvC_dXiB_at_Id(XiA, XiD).
    
    Eigen::Matrix<double,12,6> M;
    Eigen::Matrix3d RA = Sophus::SE3d::exp(-XiA).matrix().block(0,0,3,3);
    Eigen::Vector3d ta = Sophus::SE3d::exp(-XiA).matrix().block(0,3,3,1);
    Sophus::SE3d LieD = Sophus::SE3d::exp(XiD);
    Eigen::Matrix4d D = LieD.matrix();
    Eigen::Vector3d dc1 = D.block(0,0,3,1);
    Eigen::Vector3d dc2 = D.block(0,1,3,1);
    Eigen::Vector3d dc3 = D.block(0,2,3,1);
    Eigen::Vector3d dc4 = D.block(0,3,3,1);

#if 0 // LEFT_EPS_MULT
    
#if 0 // PATCH ON
    
    Eigen::Vector3d w = Sophus::SE3d::exp(-XiA).log().tail(3);
    double theta = sqrt(w.transpose()*w);
    if (theta < 1e-6)
    {
        
        dc4 = -0.5*ta + dc4;
        
    }
    
#endif
    
    M.setZero();
    M.block(0,3,3,3) = -RA*skewsym(dc1);
    M.block(3,3,3,3) = -RA*skewsym(dc2);
    M.block(6,3,3,3) = -RA*skewsym(dc3);
    M.block(9,0,3,3) = RA;
    M.block(9,3,3,3) = -RA*skewsym(dc4);
    
    return -M;
    
#else // RIGHT_EPS_MULT
    
//    Eigen::Matrix<double,3,3> RD = D.block(0,0,3,3);
//    M.setZero();
//    M.block(0,4,3,1) = -RA*dc3;
//    M.block(0,5,3,1) =  RA*dc2;
//    M.block(3,3,3,1) =  RA*dc3;
//    M.block(3,5,3,1) = -RA*dc1;
//    M.block(6,3,3,1) = -RA*dc2;
//    M.block(6,4,4,1) =  RA*dc1;
//    M.block(9,0,3,3) =  RA*RD;
//    
//    return -M;
    
    Eigen::Matrix<double,4,4> Ainv = Sophus::SE3d::exp(-XiA).matrix();
    Eigen::Matrix<double,4,4> B = Sophus::SE3d::exp(XiD).matrix();
    M = -dEps_plus_D_dXiEps(Ainv*B);
    
    return M;
    
#endif
    
}
    

#endif
    
#if 0 // ORIGINAL
    
Eigen::Matrix<double,12,12>
DervLie::
dAinv_dA(Eigen::Matrix<double,4,4> A)
{
    
    std::cout << "[dAinv_A] A = " << std::endl;
    std::cout << A << std::endl;
    
    Eigen::Matrix<double,12,12> M;
    
    // T_33, Tranpose permutation matrix of size 9 x 9
    Eigen::Matrix<double,9,9> T_33;
    T_33.setZero();
#if 0 // ORIGINAL as tutorial
    
    T_33(0,0) = 1.0;
    T_33(1,3) = 1.0;
    T_33(3,1) = 1.0;
    T_33(2,6) = 1.0;
    T_33(6,2) = 1.0;
    T_33(4,4) = 1.0;
    T_33(5,7) = 1.0;
    T_33(7,5) = 1.0;
    T_33(8,8) = 1.0;
    
#else // from ceres
//    0         1         2         3         4         5         6         7         8
//    ---------------------------------------------------------------------------------
//    0         0         0         0         0         0         0         0         0
//    0   -0.5000         0    0.5000         0         0         0         0         0
//    0         0   -0.5000         0         0         0    0.5000         0         0
//    0    0.5000         0   -0.5000         0         0         0         0         0
//    0         0         0         0         0         0         0         0         0
//    0         0         0         0         0   -0.5000         0    0.5000         0
//    0         0    0.5000         0         0         0   -0.5000         0         0
//    0         0         0         0         0    0.5000         0   -0.5000         0
//    0         0         0         0         0         0         0         0         0
    
    T_33(1,1) = -0.5;
    T_33(1,3) =  0.5;
    T_33(2,2) = -0.5;
    T_33(2,6) =  0.5;
    T_33(3,1) =  0.5;
    T_33(3,3) = -0.5;
    T_33(5,5) = -0.5;
    T_33(5,7) =  0.5;
    T_33(6,2) =  0.5;
    T_33(6,6) = -0.5;
    T_33(7,5) =  0.5;
    T_33(7,7) = -0.5;
    
#endif
    
    
    Eigen::Matrix<double,3,3> R = A.block(0,0,3,3);
    Eigen::Matrix<double,3,1> t = A.block(0,3,3,1);
    
    M.setZero();
    M.block(0,0,9,9) = T_33;
    M.block(9,0,3,9)
        = Eigen::kroneckerProduct(Eigen::Matrix3d::Identity(), -t.transpose());
    M.block(9,9,3,3) = -R.transpose();
    
    std::cout << "[dAinv_dA] M = " << std::endl;
    std::cout << M << std::endl;
    
    return M;
    
}

#else // NEW VERSION via Quaternion
    
//------------------------------------------------------------------------------
// dAinv_dA via quaternion
// d(A^{-1})/dA = [d(R^{-1}
    
Eigen::Matrix<double,12,12>
DervLie::
dAinv_dA(Eigen::Matrix<double,4,4> A)
{
    Eigen::Matrix<double,12,12> M;
    
    Eigen::Matrix<double,6,1> XiA = Sophus::SE3d(A).log();
    Eigen::Matrix<double,12,6> dAinv_dXiA = DervLie::dAinv_dXiA(XiA);
    Eigen::Matrix<double,6,12> dXiA_dA = DervLie::dXiA_dA(A);
    
    return dAinv_dXiA*dXiA_dA;
    
}
    
#endif
    

//--------------------------------------------------------------------------
// Verified by Ceres
#if 0 // ORIGINAL
Eigen::Matrix<double,12,6>
DervLie::
dBinvC_dXiC(Eigen::Matrix<double,6,1> XiB,
            Eigen::Matrix<double,6,1> XiC)
{
    
    Eigen::Matrix<double,12,6> M;
    
    Sophus::SE3d LieB = Sophus::SE3d::exp(XiB);
    Eigen::Matrix<double,4,4> B = LieB.matrix();
#if 0// NOT_USED
    Sophus::SE3d LieC = Sophus::SE3d::exp(XiC);
    Eigen::Matrix<double,4,4> C = LieC.matrix();
#endif
    
    Eigen::Matrix<double,3,3> RB = B.block(0,0,3,3);
#if 0// NOT_USED
    Eigen::Matrix<double,3,1> tB = B.block(0,3,3,1);
    Eigen::Matrix<double,3,3> RC = C.block(0,0,3,3);
    Eigen::Matrix<double,3,1> tC = C.block(0,3,3,1);
#endif
    Eigen::Matrix<double,3,3> RBinv = RB.inverse();
    
    Eigen::Matrix<double,12,6> dC_dXiC = dA_dXiA(XiC);
    
    
    // Set dRA_duC
    Eigen::Matrix<double,9,3> dRC_duC = dC_dXiC.block(0,0,9,3);
    Eigen::Matrix<double,9,3> dRA_duC =
    fast9x3_Mat3x3_x_Jacob9x3(RBinv, dRC_duC);
#if 0//DEBUG
    std::cout << "dRA_duC = " << std::endl;
    std::cout << dRA_duC << std::endl;
#endif
    // Set dRA_dwC
    Eigen::Matrix<double,9,3> dRC_dwC = dC_dXiC.block(0,3,9,3);
    Eigen::Matrix<double,9,3> dRA_dwC =
    fast9x3_Mat3x3_x_Jacob9x3(RBinv, dRC_dwC);
#if 0//DEBUG
    std::cout << "dRA_dwC = " << std::endl;
    std::cout << dRA_dwC << std::endl;
#endif
    // set dtA_duC
    Eigen::Matrix<double,3,3> dtC_duC = dC_dXiC.block(9,0,3,3);
    Eigen::Matrix<double,3,3> dtA_duC =
    fast3x3_Mat3x3_x_Jacob3x3(RBinv, dtC_duC);
#if 0//DEBUG
    std::cout << "dtA_duC = " << std::endl;
    std::cout << dtA_duC << std::endl;
#endif
    // Set dtA_dwC
    Eigen::Matrix<double,3,3> dtC_dwC = dC_dXiC.block(9,3,3,3);
    Eigen::Matrix<double,3,3> dtA_dwC =
    fast3x3_Mat3x3_x_Jacob3x3(RBinv, dtC_dwC);
    
    M << dRA_duC, dRA_dwC,
    dtA_duC, dtA_dwC;
    
    return M;
    
}
#else
Eigen::Matrix<double,12,6>
DervLie::
dBinvC_dXiC(Eigen::Matrix<double,6,1> XiB,
            Eigen::Matrix<double,6,1> XiC)
{

    return DervLie::dBC_dXiC(-XiB, XiC);
    
}
    
Eigen::Matrix<double,12,6>
DervLie::
dBinvC_dXiC_at_Id(Eigen::Matrix<double,6,1> XiB,
                  Eigen::Matrix<double,6,1> XiC)
{
    
    // Variant of Equation (10.28) in Blanco's tutorial
    // Note: We use -XiA for A inverse. Then, exp(-XiA) becomes -exp(XiA).
    //       d(A^{-1}*exp(eps)*D)/d(eps)
    
    Eigen::Matrix<double,12,6> M;
    Eigen::Matrix3d RA = Sophus::SE3d::exp(-XiB).matrix().block(0,0,3,3);
    Eigen::Vector3d ta = Sophus::SE3d::exp(-XiB).matrix().block(0,3,3,1);
    Eigen::Matrix4d D = Sophus::SE3d::exp(XiC).matrix();
    Eigen::Vector3d dc1 = D.block(0,0,3,1);
    Eigen::Vector3d dc2 = D.block(0,1,3,1);
    Eigen::Vector3d dc3 = D.block(0,2,3,1);
    Eigen::Vector3d dc4 = D.block(0,3,3,1);
 
#if 0 // LEFT_EPS_MULT
    
#if 0 // Patch ON
    
    Eigen::Vector3d w = XiC.tail(3);
    double theta = sqrt(w.transpose()*w);
    if (theta < 1e-5)
    {
        
        dc4 = 1.5*dc4;
        
    }
    
#endif
    
    M.setZero();
    M.block(0,3,3,3) = -RA*skewsym(dc1);
    M.block(3,3,3,3) = -RA*skewsym(dc2);
    M.block(6,3,3,3) = -RA*skewsym(dc3);
    M.block(9,0,3,3) = RA;
    M.block(9,3,3,3) = -RA*skewsym(dc4);
    
    return M;
    
#else // RIGHT_EPS_MULT
    
    Eigen::Matrix<double,3,3> RD = D.block(0,0,3,3);
    M.setZero();
    M.block(0,4,3,1) = -RA*dc3;
    M.block(0,5,3,1) =  RA*dc2;
    M.block(3,3,3,1) =  RA*dc3;
    M.block(3,5,3,1) = -RA*dc1;
    M.block(6,3,3,1) = -RA*dc2;
    M.block(6,4,4,1) =  RA*dc1;
    M.block(9,0,3,3) =  RA*RD;
    
    return M;
    
#endif
    
}

#endif

//--------------------------------------------------------------------------
// Verified (but not exact)
// d(exp( k * (XiA^{-1} concat XiB)) / d XiA
Eigen::Matrix<double,12,6>
DervLie::
dExp_k_XiAinv_Concat_XiB_dXiA(double k,
                              Eigen::Matrix<double,6,1> XiA,
                              Eigen::Matrix<double,6,1> XiB)
{
    
    Eigen::Matrix<double,12,6> M;
    
    Sophus::SE3d LieA = Sophus::SE3d::exp(XiA);
    Sophus::SE3d LieB = Sophus::SE3d::exp(XiB);
    Sophus::SE3d LieAinvB = LieA.inverse()*LieB;
    Eigen::Matrix<double,6,1> XiAinvXiB = LieAinvB.log();
    Eigen::Matrix<double,4,4> AinvB = LieAinvB.matrix();
    
    Eigen::Matrix<double,12,6> dExpKXiAinvConcatXiB_dKXiAinvConcatXiB;
    dExpKXiAinvConcatXiB_dKXiAinvConcatXiB = dA_dXiA(k*XiAinvXiB);
    
    Eigen::Matrix<double,6,6> dKXiAinvConcatXiB_dXiAinvXiB;
    dKXiAinvConcatXiB_dXiAinvXiB = k*Eigen::MatrixXd::Identity(6,6);
    
    Eigen::Matrix<double,6,12> dLogAinvB_dAinvB;
    dLogAinvB_dAinvB = lsd_slam::DervLie::dXiA_dA(AinvB); // The reason of inexact
    
    Eigen::Matrix<double,12,6> dAinvB_dXiA = dBinvC_dXiB(XiA, XiB);
    
    // Chain rule
    M = dExpKXiAinvConcatXiB_dKXiAinvConcatXiB *
    dKXiAinvConcatXiB_dXiAinvXiB *
    dLogAinvB_dAinvB *
    dAinvB_dXiA;

    return M;
    
}

Eigen::Matrix<double,12,12>
DervLie::
dExp_k_XiAinv_Concat_XiB_dA(double k,
                            Eigen::Matrix<double,6,1> XiA,
                            Eigen::Matrix<double,6,1> XiB)
{
    
    Eigen::Matrix<double,12,12> M;
    
    Sophus::SE3d LieA = Sophus::SE3d::exp(XiA);
    Sophus::SE3d LieB = Sophus::SE3d::exp(XiB);
    Sophus::SE3d LieAinvB = LieA.inverse()*LieB;
    Eigen::Matrix<double,6,1> XiAinvXiB = LieAinvB.log();
    Eigen::Matrix<double,4,4> AinvB = LieAinvB.matrix();
    
    Eigen::Matrix<double,12,6> dExpKXiAinvConcatXiB_dKXiAinvConcatXiB;
    dExpKXiAinvConcatXiB_dKXiAinvConcatXiB = dA_dXiA(k*XiAinvXiB);
    
    std::cout << "k*XiAinvXiB = " << std::endl;
    Eigen::Matrix<double,6,1> k_XiAinvXiB = k*XiAinvXiB;
    std::cout << k_XiAinvXiB << std::endl;
    
    Eigen::Matrix<double,6,6> dKXiAinvConcatXiB_dXiAinvXiB;
    dKXiAinvConcatXiB_dXiAinvXiB = k*Eigen::MatrixXd::Identity(6,6);
    
    Eigen::Matrix<double,6,12> dLogAinvB_dAinvB;
    dLogAinvB_dAinvB = lsd_slam::DervLie::dXiA_dA(AinvB);
    
    std::cout << "AinvB = " << std::endl;
    std::cout << AinvB << std::endl;
    
    Eigen::Matrix<double,12,12> dAinvB_dA =
    dAB_dA(LieA.inverse().matrix(), LieB.matrix())*dAinv_dA(LieA.matrix());
    
    Eigen::Matrix<double,12,12> M_tmp;
    M_tmp = dExpKXiAinvConcatXiB_dKXiAinvConcatXiB * k*dLogAinvB_dAinvB;
    std::cout << "M_tmp = " << std::endl;
    std::cout << M_tmp << std::endl;
    
    // Chain rule
    M = dExpKXiAinvConcatXiB_dKXiAinvConcatXiB *
    //dKXiAinvConcatXiB_dXiAinvXiB *
    k*dLogAinvB_dAinvB *
    dAinvB_dA;
    
    return M;
    
}

//------------------------------------------------------------------------------
Eigen::Matrix<double,12,6>
DervLie::
dExp_k_XiAinv_Concat_XiB_dXiA_at_Id(double k,
                              Eigen::Matrix<double,6,1> XiA,
                              Eigen::Matrix<double,6,1> XiB)
{
    
    Eigen::Matrix<double,12,6> M;
    
    Sophus::SE3d LieA = Sophus::SE3d::exp(XiA);
    Sophus::SE3d LieB = Sophus::SE3d::exp(XiB);
    Sophus::SE3d LieAinvB = LieA.inverse()*LieB;
    Eigen::Matrix<double,6,1> XiAinvXiB = LieAinvB.log();
    Eigen::Matrix<double,4,4> AinvB = LieAinvB.matrix();
    
    Eigen::Matrix<double,12,6> dExpKXiAinvConcatXiB_dKXiAinvConcatXiB;
    dExpKXiAinvConcatXiB_dKXiAinvConcatXiB = dA_dXiA(k*XiAinvXiB);
    
    Eigen::Matrix<double,6,6> dKXiAinvConcatXiB_dXiAinvXiB;
    dKXiAinvConcatXiB_dXiAinvXiB = k*Eigen::MatrixXd::Identity(6,6);
    
    Eigen::Matrix<double,6,12> dLogAinvB_dAinvB;
    dLogAinvB_dAinvB = lsd_slam::DervLie::dXiA_dA(AinvB);

#if 0 // USE BinvC_dXiB_at_Id which does not handling zero angle
    
    Eigen::Matrix<double,12,6> dAinvB_dXiA = dBinvC_dXiB_at_Id(XiA, XiB);
    
#else // USE dAinv_negEps_D_deps_at_Id handling zero angle
    
    Eigen::Matrix<double,12,6> dAinvB_dXiA
        = dAinv_negEps_D_deps_at_Id(XiA, XiB);

#endif
    
    // Chain rule
    M = dExpKXiAinvConcatXiB_dKXiAinvConcatXiB *
    dKXiAinvConcatXiB_dXiAinvXiB *
    dLogAinvB_dAinvB *
    dAinvB_dXiA;
    
    return M;
    
}

//--------------------------------------------------------------------------
// Verified (but not exact)
// d(exp( k * (XiA^{-1} concat XiB)) / d XiB
Eigen::Matrix<double,12,6>
DervLie::
dExp_k_XiAinv_Concat_XiB_dXiB(double k,
                              Eigen::Matrix<double,6,1> XiA,
                              Eigen::Matrix<double,6,1> XiB)
{
    
    Eigen::Matrix<double,12,6> M;
    
    Sophus::SE3d LieA = Sophus::SE3d::exp(XiA);
    Sophus::SE3d LieB = Sophus::SE3d::exp(XiB);
    Sophus::SE3d LieAinvB = LieA.inverse()*LieB;
    Eigen::Matrix<double,6,1> XiAinvXiB = LieAinvB.log();
    Eigen::Matrix<double,4,4> AinvB = LieAinvB.matrix();
    
    Eigen::Matrix<double,12,6> dExpKXiAinvConcatXiB_dKXiAinvConcatXiB;
    dExpKXiAinvConcatXiB_dKXiAinvConcatXiB = dA_dXiA(k*XiAinvXiB);
    
    Eigen::Matrix<double,6,6> dKXiAinvConcatXiB_dXiAinvXiB;
    dKXiAinvConcatXiB_dXiAinvXiB = k*Eigen::MatrixXd::Identity(6,6);
    
    Eigen::Matrix<double,6,12> dLogAinvB_dAinvB;
    dLogAinvB_dAinvB = dXiA_dA(AinvB); // The reason of inexact
    
    Eigen::Matrix<double,12,6> dAinvB_dXiB = dBinvC_dXiC(XiA, XiB);
    
    // Chain rule
    M = dExpKXiAinvConcatXiB_dKXiAinvConcatXiB *
    dKXiAinvConcatXiB_dXiAinvXiB *
    dLogAinvB_dAinvB *
    dAinvB_dXiB;
    
    return M;
    
}

//--------------------------------------------------------------------------
// d(exp( k * (XiA^{-1} concat XiB)) / d XiB
Eigen::Matrix<double,12,6>
DervLie::
dExp_k_XiAinv_Concat_XiB_dXiB_at_Id(double k,
                              Eigen::Matrix<double,6,1> XiA,
                              Eigen::Matrix<double,6,1> XiB)
{
    
    Eigen::Matrix<double,12,6> M;
    
    Sophus::SE3d LieA = Sophus::SE3d::exp(XiA);
    Sophus::SE3d LieB = Sophus::SE3d::exp(XiB);
    Sophus::SE3d LieAinvB = LieA.inverse()*LieB;
    Eigen::Matrix<double,6,1> XiAinvXiB = LieAinvB.log();
    Eigen::Matrix<double,4,4> AinvB = LieAinvB.matrix();
    
    Eigen::Matrix<double,12,6> dExpKXiAinvConcatXiB_dKXiAinvConcatXiB;
    dExpKXiAinvConcatXiB_dKXiAinvConcatXiB = dA_dXiA(k*XiAinvXiB);
    
    Eigen::Matrix<double,6,6> dKXiAinvConcatXiB_dXiAinvXiB;
    dKXiAinvConcatXiB_dXiAinvXiB = k*Eigen::MatrixXd::Identity(6,6);
    
    Eigen::Matrix<double,6,12> dLogAinvB_dAinvB;
    dLogAinvB_dAinvB = lsd_slam::DervLie::dXiA_dA(AinvB);
    
    //Eigen::Matrix<double,12,6> dAinvB_dXiB = dBinvC_dXiC_at_Id(XiA, XiB);
    Eigen::Matrix<double,12,6> dAinvB_dXiB = dD_plus_Eps_dXiEps(AinvB);
    
    // Chain rule
    M = dExpKXiAinvConcatXiB_dKXiAinvConcatXiB *
    dKXiAinvConcatXiB_dXiAinvXiB *
    dLogAinvB_dAinvB *
    dAinvB_dXiB;
    
    return M;
    
}

//--------------------------------------------------------------------------
// Verified by Ceres (exactly)
// d(exp( k * (XiA^{-1} concat XiB)) / d k
Eigen::Matrix<double,12,1>
DervLie::
dExp_k_XiAinv_Concat_XiB_dK(double k,
                            Eigen::Matrix<double,6,1> XiA,
                            Eigen::Matrix<double,6,1> XiB)
{
    
    Eigen::Matrix<double,12,1> M;
    
    Sophus::SE3d LieA = Sophus::SE3d::exp(XiA);
    Sophus::SE3d LieB = Sophus::SE3d::exp(XiB);
    Sophus::SE3d LieAinvB = LieA.inverse()*LieB;
    Eigen::Matrix<double,6,1> XiAinvXiB = LieAinvB.log();
#if 0//NOT_USED
    Eigen::Matrix<double,4,4> AinvB = LieAinvB.matrix();
#endif
    Eigen::Matrix<double,12,6> dExpKXiAinvConcatXiB_dKXiAinvConcatXiB;
    dExpKXiAinvConcatXiB_dKXiAinvConcatXiB = dA_dXiA(k*XiAinvXiB);
    
    Eigen::Matrix<double,6,1> dKXiAinvConcatXiB_dK;
    dKXiAinvConcatXiB_dK = XiAinvXiB;
    
    // Chain rule
    M = dExpKXiAinvConcatXiB_dKXiAinvConcatXiB *
    dKXiAinvConcatXiB_dK;
    
    return M;
    
}

//--------------------------------------------------------------------------
// 3x3 Fast Jacob_x_Vec
// Verified by Ceres
Eigen::Matrix<double,3,3>
DervLie::
fast3x3_Jacob9x3_x_Vec3x1(const Eigen::Ref<const
                          Eigen::Matrix<double,9,3> >& J,
                          const Eigen::Matrix<double,3,1> V)
{
    
    Eigen::Matrix<double,3,3> Z;
    
    for (int i=0; i<3; i++)
    {
        
        Z(0,i) = J(0,i)*V(0) + J(3,i)*V(1) + J(6,i)*V(2);
        Z(1,i) = J(1,i)*V(0) + J(4,i)*V(1) + J(7,i)*V(2);
        Z(2,i) = J(2,i)*V(0) + J(5,i)*V(1) + J(8,i)*V(2);
        
    }
    
    return Z;
    
}

//--------------------------------------------------------------------------
// Verified by Ceres
Eigen::Matrix<double,3,6>
DervLie::
dRXt_dXi(Eigen::Matrix<double,3,1> X,
         Eigen::Matrix<double,6,1> Xi)
{
   
    // dy/dXi where y = R*X + t and exp(Xi) = [R t; 0 0 0 1]
    // X is the input point with respect to source (keyframe)
    // y is the warped point with respect to target (curframe)
    // Xi is a transformation from source to target (key to cur)
    
    Eigen::Matrix<double,3,6> M;
    
#if 0 // ORIGINAL

    Eigen::Matrix<double,3,1> u, w;
    u << Xi(0), Xi(1), Xi(2); // upsilon translation
    w << Xi(3), Xi(4), Xi(5); // omega rotation
    
    Eigen::Matrix<double,9,3> dRdu;
    dRdu.setZero(9,3);
    
    Eigen::Matrix<double,3,3> dydu;
    Eigen::Matrix<double,3,3> J1;
    J1 = fast3x3_Jacob9x3_x_Vec3x1(dRdu, X);
    Eigen::Matrix<double,3,3> dtdu = dt_du(w);
    dydu = J1 + dtdu;
    
    Eigen::Matrix<double,3,3> dydw;
    Eigen::Matrix<double,9,3> dRdw = dR_dw(w);
    Eigen::Matrix<double,3,3> dtdw = dt_dw(u, w);
    dydw = fast3x3_Jacob9x3_x_Vec3x1(dRdw, X) + dtdw;
    
    M << dydu, dydw;
    
#else
    
    Eigen::Matrix<double,1,4> p;
    p << X.transpose(), 1;
    Eigen::Matrix<double,3,12> dAp_dA
        = Eigen::kroneckerProduct(p, Eigen::Matrix3d::Identity());
    
    M = dAp_dA*dA_dXiA(Xi);
    
#endif
    
    return M;
    
}
    
Eigen::Matrix<double,3,12>
DervLie::dAx_dA(Eigen::Matrix<double,4,4> A, Eigen::Matrix<double,4,1> x)
{
    // d(A*x)/d(A) where A in SE(3) and x is a 4 x 1 homogeneous vector
    return Eigen::kroneckerProduct(x.transpose(),
                                   Eigen::Matrix3d::Identity());
        
}

//--------------------------------------------------------------------------
Eigen::Matrix<double,12,12>
DervLie::
fast12x12_multMatToJacob(const Eigen::Matrix<double,4,4> M,
                         const Eigen::Ref<const
                         Eigen::Matrix<double,12,12> >& J)
{

    Eigen::Matrix<double,12,12> Z;
    
    for (int j=0;j<J.cols();j++)
    {
        
        // Equation Z from MATLAB
        // Z =
        //      j1*m1_1 + j2*m1_2 + j3*m1_3
        //      j1*m2_1 + j2*m2_2 + j3*m2_3
        //      j1*m3_1 + j2*m3_2 + j3*m3_3
        //      j4*m1_1 + j5*m1_2 + j6*m1_3
        //      j4*m2_1 + j5*m2_2 + j6*m2_3
        //      j4*m3_1 + j5*m3_2 + j6*m3_3
        //      j7*m1_1 + j8*m1_2 + j9*m1_3
        //      j7*m2_1 + j8*m2_2 + j9*m2_3
        //      j7*m3_1 + j8*m3_2 + j9*m3_3
        //      j10*m1_1 + j11*m1_2 + j12*m1_3
        //      j10*m2_1 + j11*m2_2 + j12*m2_3
        //      j10*m3_1 + j11*m3_2 + j12*m3_3
        Z(0,j) =  J(0,j)*M(0,0) + J(1,j)*M(0,1) + J(2,j)*M(0,2);
        Z(1,j) =  J(0,j)*M(1,0) + J(1,j)*M(1,1) + J(2,j)*M(1,2);
        Z(2,j) =  J(0,j)*M(2,0) + J(1,j)*M(2,1) + J(2,j)*M(2,2);
        Z(3,j) =  J(3,j)*M(0,0) + J(4,j)*M(0,1) + J(5,j)*M(0,2);
        Z(4,j) =  J(3,j)*M(1,0) + J(4,j)*M(1,1) + J(5,j)*M(1,2);
        Z(5,j) =  J(3,j)*M(2,0) + J(4,j)*M(2,1) + J(5,j)*M(2,2);
        Z(6,j) =  J(6,j)*M(0,0) + J(7,j)*M(0,1) + J(8,j)*M(0,2);
        Z(7,j) =  J(6,j)*M(1,0) + J(7,j)*M(1,1) + J(8,j)*M(1,2);
        Z(8,j) =  J(6,j)*M(2,0) + J(7,j)*M(2,1) + J(8,j)*M(2,2);
        Z(9,j) =  J(9,j)*M(0,0) + J(10,j)*M(0,1) + J(11,j)*M(0,2);
        Z(10,j) = J(9,j)*M(1,0) + J(10,j)*M(1,1) + J(11,j)*M(1,2);
        Z(11,j) = J(9,j)*M(2,0) + J(10,j)*M(2,1) + J(11,j)*M(2,2);
        
    }
    
    return Z;

}

//--------------------------------------------------------------------------
Eigen::Matrix<double,12,12>
DervLie::
fast12x12_multJacobToMat(const Eigen::Ref<const
                         Eigen::Matrix<double,12,12> >& J,
                         const Eigen::Matrix<double,4,4> M)
{
    
    Eigen::Matrix<double,12,12> Z;
    
    for (int j=0;j<J.cols();j++)
    {
        
        // Equation Z from MATLAB
        // Z =
        //       j1*m1_1 + j4*m2_1 + j7*m3_1
        //       j2*m1_1 + j5*m2_1 + j8*m3_1
        //       j3*m1_1 + j6*m2_1 + j9*m3_1
        //       j1*m1_2 + j4*m2_2 + j7*m3_2
        //       j2*m1_2 + j5*m2_2 + j8*m3_2
        //       j3*m1_2 + j6*m2_2 + j9*m3_2
        //       j1*m1_3 + j4*m2_3 + j7*m3_3
        //       j2*m1_3 + j5*m2_3 + j8*m3_3
        //       j3*m1_3 + j6*m2_3 + j9*m3_3
        // j10 + j1*m1_4 + j4*m2_4 + j7*m3_4
        // j11 + j2*m1_4 + j5*m2_4 + j8*m3_4
        // j12 + j3*m1_4 + j6*m2_4 + j9*m3_4
        
        Z(0,j) =            J(0,j)*M(0,0) + J(3,j)*M(1,0) + J(6,j)*M(2,0);
        Z(1,j) =            J(1,j)*M(0,0) + J(4,j)*M(1,0) + J(7,j)*M(2,0);
        Z(2,j) =            J(2,j)*M(0,0) + J(5,j)*M(1,0) + J(8,j)*M(2,0);
        Z(3,j) =            J(0,j)*M(0,1) + J(3,j)*M(1,1) + J(6,j)*M(2,1);
        Z(4,j) =            J(1,j)*M(0,1) + J(4,j)*M(1,1) + J(7,j)*M(2,1);
        Z(5,j) =            J(2,j)*M(0,1) + J(5,j)*M(1,1) + J(8,j)*M(2,1);
        Z(6,j) =            J(0,j)*M(0,2) + J(3,j)*M(1,2) + J(6,j)*M(2,2);
        Z(7,j) =            J(1,j)*M(0,2) + J(4,j)*M(1,2) + J(7,j)*M(2,2);
        Z(8,j) =            J(2,j)*M(0,2) + J(5,j)*M(1,2) + J(8,j)*M(2,2);
        Z(9,j) =   J(9,j) + J(0,j)*M(0,3) + J(3,j)*M(1,3) + J(6,j)*M(2,3);
        Z(10,j) = J(10,j) + J(1,j)*M(0,3) + J(4,j)*M(1,3) + J(7,j)*M(2,3);
        Z(11,j) = J(11,j) + J(2,j)*M(0,3) + J(5,j)*M(1,3) + J(8,j)*M(2,3);
        
    }
    
    return Z;
    
}

    
//--------------------------------------------------------------------------
// Verified
Eigen::Matrix<double,12,6>
DervLie::
fast12x6_multJacobToMat(const Eigen::Ref<const
                        Eigen::Matrix<double,12,6> >& J,
                        const Eigen::Matrix<double,4,4> M)
{
    Eigen::Matrix<double,12,6> Z;
    for (int j=0;j<J.cols();j++)
    {
        
        // Equation Z from MATLAB
        // Z =
        //       j1*m1_1 + j4*m2_1 + j7*m3_1
        //       j2*m1_1 + j5*m2_1 + j8*m3_1
        //       j3*m1_1 + j6*m2_1 + j9*m3_1
        //       j1*m1_2 + j4*m2_2 + j7*m3_2
        //       j2*m1_2 + j5*m2_2 + j8*m3_2
        //       j3*m1_2 + j6*m2_2 + j9*m3_2
        //       j1*m1_3 + j4*m2_3 + j7*m3_3
        //       j2*m1_3 + j5*m2_3 + j8*m3_3
        //       j3*m1_3 + j6*m2_3 + j9*m3_3
        // j10 + j1*m1_4 + j4*m2_4 + j7*m3_4
        // j11 + j2*m1_4 + j5*m2_4 + j8*m3_4
        // j12 + j3*m1_4 + j6*m2_4 + j9*m3_4
        
        Z(0,j) =            J(0,j)*M(0,0) + J(3,j)*M(1,0) + J(6,j)*M(2,0);
        Z(1,j) =            J(1,j)*M(0,0) + J(4,j)*M(1,0) + J(7,j)*M(2,0);
        Z(2,j) =            J(2,j)*M(0,0) + J(5,j)*M(1,0) + J(8,j)*M(2,0);
        Z(3,j) =            J(0,j)*M(0,1) + J(3,j)*M(1,1) + J(6,j)*M(2,1);
        Z(4,j) =            J(1,j)*M(0,1) + J(4,j)*M(1,1) + J(7,j)*M(2,1);
        Z(5,j) =            J(2,j)*M(0,1) + J(5,j)*M(1,1) + J(8,j)*M(2,1);
        Z(6,j) =            J(0,j)*M(0,2) + J(3,j)*M(1,2) + J(6,j)*M(2,2);
        Z(7,j) =            J(1,j)*M(0,2) + J(4,j)*M(1,2) + J(7,j)*M(2,2);
        Z(8,j) =            J(2,j)*M(0,2) + J(5,j)*M(1,2) + J(8,j)*M(2,2);
        Z(9,j) =   J(9,j) + J(0,j)*M(0,3) + J(3,j)*M(1,3) + J(6,j)*M(2,3);
        Z(10,j) = J(10,j) + J(1,j)*M(0,3) + J(4,j)*M(1,3) + J(7,j)*M(2,3);
        Z(11,j) = J(11,j) + J(2,j)*M(0,3) + J(5,j)*M(1,3) + J(8,j)*M(2,3);
        
    }
    
    return Z;
    
}

//--------------------------------------------------------------------------
// Verified
Eigen::Matrix<double,12,1>
DervLie::
fast12x1_multJacobToMat(const Eigen::Ref<const
                        Eigen::Matrix<double,12,1> >& J,
                        const Eigen::Matrix<double,4,4> M)
{
    Eigen::Matrix<double,12,1> Z;
    //for (int j=0;j<J.cols();j++)
    int j = 0;
    {
        
        Z(0,j) =            J(0,j)*M(0,0) + J(3,j)*M(1,0) + J(6,j)*M(2,0);
        Z(1,j) =            J(1,j)*M(0,0) + J(4,j)*M(1,0) + J(7,j)*M(2,0);
        Z(2,j) =            J(2,j)*M(0,0) + J(5,j)*M(1,0) + J(8,j)*M(2,0);
        Z(3,j) =            J(0,j)*M(0,1) + J(3,j)*M(1,1) + J(6,j)*M(2,1);
        Z(4,j) =            J(1,j)*M(0,1) + J(4,j)*M(1,1) + J(7,j)*M(2,1);
        Z(5,j) =            J(2,j)*M(0,1) + J(5,j)*M(1,1) + J(8,j)*M(2,1);
        Z(6,j) =            J(0,j)*M(0,2) + J(3,j)*M(1,2) + J(6,j)*M(2,2);
        Z(7,j) =            J(1,j)*M(0,2) + J(4,j)*M(1,2) + J(7,j)*M(2,2);
        Z(8,j) =            J(2,j)*M(0,2) + J(5,j)*M(1,2) + J(8,j)*M(2,2);
        Z(9,j) =   J(9,j) + J(0,j)*M(0,3) + J(3,j)*M(1,3) + J(6,j)*M(2,3);
        Z(10,j) = J(10,j) + J(1,j)*M(0,3) + J(4,j)*M(1,3) + J(7,j)*M(2,3);
        Z(11,j) = J(11,j) + J(2,j)*M(0,3) + J(5,j)*M(1,3) + J(8,j)*M(2,3);
        
    }
    
    return Z;
    
}

//--------------------------------------------------------------------------
// Verified
Eigen::Matrix<double,12,6>
DervLie::
fast12x6_multMatToJacob(const Eigen::Matrix<double,4,4> M,
                        const Eigen::Ref<const
                        Eigen::Matrix<double,12,6> >& J)
{
    
    Eigen::Matrix<double,12,6> Z;
    
    for (int j=0;j<J.cols();j++)
    {
        
        // Equation Z from MATLAB
        // Z =
        //      j1*m1_1 + j2*m1_2 + j3*m1_3
        //      j1*m2_1 + j2*m2_2 + j3*m2_3
        //      j1*m3_1 + j2*m3_2 + j3*m3_3
        //      j4*m1_1 + j5*m1_2 + j6*m1_3
        //      j4*m2_1 + j5*m2_2 + j6*m2_3
        //      j4*m3_1 + j5*m3_2 + j6*m3_3
        //      j7*m1_1 + j8*m1_2 + j9*m1_3
        //      j7*m2_1 + j8*m2_2 + j9*m2_3
        //      j7*m3_1 + j8*m3_2 + j9*m3_3
        //      j10*m1_1 + j11*m1_2 + j12*m1_3
        //      j10*m2_1 + j11*m2_2 + j12*m2_3
        //      j10*m3_1 + j11*m3_2 + j12*m3_3
        Z(0,j) =  J(0,j)*M(0,0) + J(1,j)*M(0,1) + J(2,j)*M(0,2);
        Z(1,j) =  J(0,j)*M(1,0) + J(1,j)*M(1,1) + J(2,j)*M(1,2);
        Z(2,j) =  J(0,j)*M(2,0) + J(1,j)*M(2,1) + J(2,j)*M(2,2);
        Z(3,j) =  J(3,j)*M(0,0) + J(4,j)*M(0,1) + J(5,j)*M(0,2);
        Z(4,j) =  J(3,j)*M(1,0) + J(4,j)*M(1,1) + J(5,j)*M(1,2);
        Z(5,j) =  J(3,j)*M(2,0) + J(4,j)*M(2,1) + J(5,j)*M(2,2);
        Z(6,j) =  J(6,j)*M(0,0) + J(7,j)*M(0,1) + J(8,j)*M(0,2);
        Z(7,j) =  J(6,j)*M(1,0) + J(7,j)*M(1,1) + J(8,j)*M(1,2);
        Z(8,j) =  J(6,j)*M(2,0) + J(7,j)*M(2,1) + J(8,j)*M(2,2);
        Z(9,j) =  J(9,j)*M(0,0) + J(10,j)*M(0,1) + J(11,j)*M(0,2);
        Z(10,j) = J(9,j)*M(1,0) + J(10,j)*M(1,1) + J(11,j)*M(1,2);
        Z(11,j) = J(9,j)*M(2,0) + J(10,j)*M(2,1) + J(11,j)*M(2,2);
        
    }
    
    return Z;
    
}

// Verified
Eigen::Matrix<double,12,1>
DervLie::
fast12x1_multMatToJacob(const Eigen::Matrix<double,4,4> M,
                        const Eigen::Ref<const
                        Eigen::Matrix<double,12,1> >& J)
{
    
    Eigen::Matrix<double,12,1> Z;
    
    //for (int j=0;j<J.cols();j++)
    int j = 0;
    {
        
        Z(0,j) =  J(0,j)*M(0,0) + J(1,j)*M(0,1) + J(2,j)*M(0,2);
        Z(1,j) =  J(0,j)*M(1,0) + J(1,j)*M(1,1) + J(2,j)*M(1,2);
        Z(2,j) =  J(0,j)*M(2,0) + J(1,j)*M(2,1) + J(2,j)*M(2,2);
        Z(3,j) =  J(3,j)*M(0,0) + J(4,j)*M(0,1) + J(5,j)*M(0,2);
        Z(4,j) =  J(3,j)*M(1,0) + J(4,j)*M(1,1) + J(5,j)*M(1,2);
        Z(5,j) =  J(3,j)*M(2,0) + J(4,j)*M(2,1) + J(5,j)*M(2,2);
        Z(6,j) =  J(6,j)*M(0,0) + J(7,j)*M(0,1) + J(8,j)*M(0,2);
        Z(7,j) =  J(6,j)*M(1,0) + J(7,j)*M(1,1) + J(8,j)*M(1,2);
        Z(8,j) =  J(6,j)*M(2,0) + J(7,j)*M(2,1) + J(8,j)*M(2,2);
        Z(9,j) =  J(9,j)*M(0,0) + J(10,j)*M(0,1) + J(11,j)*M(0,2);
        Z(10,j) = J(9,j)*M(1,0) + J(10,j)*M(1,1) + J(11,j)*M(1,2);
        Z(11,j) = J(9,j)*M(2,0) + J(10,j)*M(2,1) + J(11,j)*M(2,2);
        
    }
    
    return Z;
    
}

//--------------------------------------------------------------------------
// 9x3 Fast Mat_x_Jacob
// Verified
Eigen::Matrix<double,9,3>
DervLie::
fast9x3_Mat3x3_x_Jacob9x3(const Eigen::Matrix<double,3,3> M,
                          const Eigen::Ref<const
                          Eigen::Matrix<double,9,3> >& J)
{
    
    Eigen::Matrix<double,9,3> Z;
    
    for (int i=0; i<3; i++)
    {
        
        Z(0,i) = J(0,i)*M(0,0) + J(1,i)*M(0,1) + J(2,i)*M(0,2);
        Z(1,i) = J(0,i)*M(1,0) + J(1,i)*M(1,1) + J(2,i)*M(1,2);
        Z(2,i) = J(0,i)*M(2,0) + J(1,i)*M(2,1) + J(2,i)*M(2,2);
        Z(3,i) = J(3,i)*M(0,0) + J(4,i)*M(0,1) + J(5,i)*M(0,2);
        Z(4,i) = J(3,i)*M(1,0) + J(4,i)*M(1,1) + J(5,i)*M(1,2);
        Z(5,i) = J(3,i)*M(2,0) + J(4,i)*M(2,1) + J(5,i)*M(2,2);
        Z(6,i) = J(6,i)*M(0,0) + J(7,i)*M(0,1) + J(8,i)*M(0,2);
        Z(7,i) = J(6,i)*M(1,0) + J(7,i)*M(1,1) + J(8,i)*M(1,2);
        Z(8,i) = J(6,i)*M(2,0) + J(7,i)*M(2,1) + J(8,i)*M(2,2);
        
    }
    
    return Z;
    
}

//--------------------------------------------------------------------------
// 9x3 Fast Jacob_x_Mat
// Verified
Eigen::Matrix<double,9,3>
DervLie::
fast9x3_Jacob9x3_x_Mat3x3(const Eigen::Ref<const
                          Eigen::Matrix<double,9,3> >& J,
                          const Eigen::Matrix<double,3,3> M)
{
    
    Eigen::Matrix<double,9,3> Z;
    
    for (int i=0; i<3; i++)
    {
        
        Z(0,i) = J(0,i)*M(0,0) + J(3,i)*M(1,0) + J(6,i)*M(2,0);
        Z(1,i) = J(1,i)*M(0,0) + J(4,i)*M(1,0) + J(7,i)*M(2,0);
        Z(2,i) = J(2,i)*M(0,0) + J(5,i)*M(1,0) + J(8,i)*M(2,0);
        Z(3,i) = J(0,i)*M(0,1) + J(3,i)*M(1,1) + J(6,i)*M(2,1);
        Z(4,i) = J(1,i)*M(0,1) + J(4,i)*M(1,1) + J(7,i)*M(2,1);
        Z(5,i) = J(2,i)*M(0,1) + J(5,i)*M(1,1) + J(8,i)*M(2,1);
        Z(6,i) = J(0,i)*M(0,2) + J(3,i)*M(1,2) + J(6,i)*M(2,2);
        Z(7,i) = J(1,i)*M(0,2) + J(4,i)*M(1,2) + J(7,i)*M(2,2);
        Z(8,i) = J(2,i)*M(0,2) + J(5,i)*M(1,2) + J(8,i)*M(2,2);
        
    }
    
    return Z;
    
}

//--------------------------------------------------------------------------
// 3x3 Fast Mat_x_Jacob
// Verified
Eigen::Matrix<double,3,3>
DervLie::
fast3x3_Mat3x3_x_Jacob3x3(const Eigen::Matrix<double,3,3> M,
                          const Eigen::Ref<const
                          Eigen::Matrix<double,3,3> >& J)
{
    
    Eigen::Matrix<double,3,3> Z;
    
    // Eq from MATLAB
    // Z =
    // [ j1_1*m1_1 + j2_1*m1_2 + j3_1*m1_3, j1_2*m1_1 + j2_2*m1_2 + j3_2*m1_3, j1_3*m1_1 + j2_3*m1_2 + j3_3*m1_3]
    // [ j1_1*m2_1 + j2_1*m2_2 + j3_1*m2_3, j1_2*m2_1 + j2_2*m2_2 + j3_2*m2_3, j1_3*m2_1 + j2_3*m2_2 + j3_3*m2_3]
    // [ j1_1*m3_1 + j2_1*m3_2 + j3_1*m3_3, j1_2*m3_1 + j2_2*m3_2 + j3_2*m3_3, j1_3*m3_1 + j2_3*m3_2 + j3_3*m3_3]
    for (int i=0; i<3; i++)
    {
        
        Z(0,i) = J(0,i)*M(0,0) + J(1,i)*M(0,1) + J(2,i)*M(0,2);
        Z(1,i) = J(0,i)*M(1,0) + J(1,i)*M(1,1) + J(2,i)*M(1,2);
        Z(2,i) = J(0,i)*M(2,0) + J(1,i)*M(2,1) + J(2,i)*M(2,2);
        
    }
    
    return Z;
    
}

//--------------------------------------------------------------------------
// UN-VERIFIED


Eigen::Matrix<double,12,6>
DervLie::
addJacobToMat(Eigen::Matrix<double,12,6> J,
              Eigen::Matrix<double,4,4> M)
{
    
    for (int j=0; j<J.cols(); j++)
    {
        
        Eigen::Matrix<double,4,4> tmp_J;
        
        // Copy a column of J to a matrix
        int i=0;
        for (int q=0; q<tmp_J.cols(); q++)
        {
            
            for (int p=0; p<tmp_J.rows(); p++)
            {
                
                tmp_J(p,q) = J(i,j);
                i++;
                
            }
            
        }
        tmp_J(3,0) = 0;
        tmp_J(3,1) = 0;
        tmp_J(3,2) = 0;
        tmp_J(3,3) = 1;
        
        // Compute the right addition
        tmp_J = tmp_J + M;
        
        // Copy back the result to the jacobian column
        i=0;
        for (int q=0; q<tmp_J.cols(); q++)
        {
            
            for (int p=0; p<tmp_J.rows(); p++)
            {
                
                J(i,j) = tmp_J(p,q);
                i++;
                
            }
            
        }
        
    }
    
    return J;
    
}

Eigen::Matrix<double,12,Eigen::Dynamic>
DervLie::
multJacobToMat(const Eigen::Ref<const
               Eigen::Matrix<double,12,Eigen::Dynamic> >& J,
               const Eigen::Matrix<double,4,4> M)
{
    
    Eigen::Matrix<double,12,Eigen::Dynamic> J_out = J;
    
    for (int j=0; j<J.cols(); j++)
    {
        
        Eigen::Matrix<double,4,4> tmp_J;
        
        // Copy a column of J to a matrix
        int i=0;
        for (int q=0; q<4; q++)
        {
            
            for (int p=0; p<3; p++)
            {
                
                tmp_J(p,q) = J(i,j);
                i++;
                
            }
            
        }
        tmp_J(3,0) = 0;
        tmp_J(3,1) = 0;
        tmp_J(3,2) = 0;
        tmp_J(3,3) = 0;
        
        // Compute the right multilication
        tmp_J = tmp_J*M;
        
        // Copy back the result to the jacobian column
        i=0;
        for (int q=0; q<4; q++)
        {
            
            for (int p=0; p<3; p++)
            {
                
                J_out(i,j) = tmp_J(p,q);
                i++;
                
            }
            
        }
        
    }
    
    return J_out;
    
}


Eigen::Matrix<double,12,Eigen::Dynamic>
DervLie::
multMatToJacob(const Eigen::Matrix<double,4,4> M,
               const Eigen::Ref<const
               Eigen::Matrix<double,12,Eigen::Dynamic> >& J)
{
    
    Eigen::Matrix<double,12,Eigen::Dynamic> J_out = J;
    
    for (int j=0; j<J.cols(); j++)
    {
        
        Eigen::Matrix<double,4,4> tmp_J;
        
        // Copy a column of J to a matrix
        int i=0;
        for (int q=0; q<4; q++)
        {
            
            for (int p=0; p<3; p++)
            {
                
                tmp_J(p,q) = J(i,j);
                i++;
                
            }
            
        }
        tmp_J(3,0) = 0;
        tmp_J(3,1) = 0;
        tmp_J(3,2) = 0;
        tmp_J(3,3) = 0;
        
        // Compute the left multilication
        tmp_J = M*tmp_J;
        
        // Copy back the result to the jacobian column
        i=0;
        for (int q=0; q<4; q++)
        {
            
            for (int p=0; p<3; p++)
            {
                
                J_out(i,j) = tmp_J(p,q);
                i++;
                
            }
            
        }
        
    }
    
    return J_out;
    
}



//--------------------------------------------------------------------------

Eigen::Matrix<double,12,6>
DervLie::
dA_dXiB(Eigen::Matrix<double,4,4> A)
{
    
    Eigen::Matrix<double,12,6> M;
    M << 0, 0, 0,      0,   A(2,0), -A(1,0),
    0, 0, 0, -A(2,0),       0,  A(0,0),
    0, 0, 0,  A(1,0), -A(0,0),       0,
    //
    0, 0, 0,      0,   A(2,1), -A(1,1),
    0, 0, 0, -A(2,1),       0,  A(0,1),
    0, 0, 0,  A(1,1), -A(0,1),       0,
    //
    0, 0, 0,      0,   A(2,2), -A(1,2),
    0, 0, 0, -A(2,2),       0,  A(0,2),
    0, 0, 0,  A(1,2), -A(0,2),       0,
    //
    1, 0, 0,      0,   A(2,3), -A(1,3),
    0, 1, 0, -A(2,3),       0,  A(0,3),
    0, 0, 1,  A(1,3), -A(0,3),       0;
    
    return M;
    
}

Eigen::Matrix<double,12,6>
DervLie::
dA_dXiC(Eigen::Matrix<double,4,4> A,
        Eigen::Matrix<double,4,4> B)
{
    
    Eigen::Matrix<double,3,1> R_A1 = A.block(0,0,3,1);
    Eigen::Matrix<double,3,1> R_A2 = A.block(0,1,3,1);
    Eigen::Matrix<double,3,1> R_A3 = A.block(0,2,3,1);
    Eigen::Matrix<double,3,1> R_B1 = B.block(0,0,3,1);
    Eigen::Matrix<double,3,1> R_B2 = B.block(0,1,3,1);
    Eigen::Matrix<double,3,1> R_B3 = B.block(0,2,3,1);
    Eigen::Matrix<double,3,3> R_B = B.block(0,0,3,3);
    Eigen::Matrix<double,3,1> t_A = A.block(0,3,3,1);
    Eigen::Matrix<double,3,1> t_B = B.block(0,3,3,1);
    
    Eigen::Matrix<double,3,3> tB_cross_R_B = skewsym(t_B)*R_B;
    
    Eigen::Matrix<double,12,6> M;
    M << Eigen::MatrixXd::Zero(3,3), R_B1.cross(R_A1),
    R_B2.cross(R_A1), R_B3.cross(R_A1),
    Eigen::MatrixXd::Zero(3,3), R_B1.cross(R_A2),
    R_B2.cross(R_A2), R_B3.cross(R_A2),
    Eigen::MatrixXd::Zero(3,3), R_B1.cross(R_A3),
    R_B2.cross(R_A3), R_B3.cross(R_A3),
    R_B,
    skewsym(R_B1)*t_A + tB_cross_R_B.col(0),
    skewsym(R_B2)*t_A + tB_cross_R_B.col(1),
    skewsym(R_B3)*t_A + tB_cross_R_B.col(2);
    
    return M;
    
}

//Eigen::Matrix<double,6,12>
//DervLie::
//dXiA_dA(Eigen::Matrix<double,4,4> A)
//{
//
//    Eigen::Matrix<double,3,3> R = A.block(0,0,3,3);
//    Eigen::Matrix<double,3,1> t = A.block(0,3,3,1);
//
//    double r12 = R(0,1);
//    double r13 = R(0,2);
//    double r21 = R(1,0);
//    double r23 = R(1,2);
//    double r31 = R(2,0);
//    double r32 = R(2,1);
//
//    double alpha = (R.trace() - 1.0)/2.0;
//    double theta = acos(alpha);
//    double d02sin0_d0;
//
//    if ((theta < 1e-6) && (theta > -1e-6))
//    {
//        // If theta is close to zero, take the taylor expansion and
//        // approximate it.
//        d02sin0_d0 = theta/6.0;
//
//    }
//    else
//    {
//
//        d02sin0_d0 = (sin(theta)
//                      - theta*cos(theta))/(2*sin(theta)*sin(theta));
//
//    }
//
//    double d0_dalpha = -1.0/sqrt(1.0-alpha*alpha);
//    if ((theta < 1e-6) && (theta > -1e-6))
//    {
//
//        d0_dalpha = 0.0;
//
//    }
//    else
//    {
//
//        d0_dalpha = -1.0/sqrt(1.0 - alpha*alpha);
//
//    }
//    double d02sin0_rii = d02sin0_d0 * d0_dalpha * 0.5;
//
//    double s = d02sin0_rii;
//    double h;
//    if ((theta < 1e-6) && (theta > -1e-6))
//    {
//
//        h = 0.5 + theta*theta/12.0;
//
//    }
//    else
//    {
//
//        h = theta/(2*sin(theta));
//
//    }
//
//    // Derivative of w (orientation) w.r.t. R (rotation)
//    Eigen::Matrix<double,3,9> dw_dR;
//    dw_dR << s*(r32-r23),  0, 0, 0, s*(r32-r23), -h,  0, h, s*(r32-r23),
//    s*(r13-r31),  0, h, 0, s*(r13-r31),  0, -h, 0, s*(r13-r31),
//    s*(r21-r12), -h, 0, h, s*(r21-r12),  0,  0, 0, s*(r21-r12);
//
//    // Derivative of w (orientation) w.r.t. t (translation)
//    Eigen::Matrix<double,3,3> dw_dt;
//    dw_dt.setZero(3,3);
//
//    Eigen::Matrix<double,3,1> w;
//    w << h*(r32 - r23), h*(r13 - r31), h*(r21 - r12);
//    double w1 = w(0);
//    double w2 = w(1);
//    double w3 = w(2);
//    Eigen::Matrix<double,3,3> wx = skewsym(w);
//    Eigen::Matrix<double,3,3> wx_copy = wx;
//    Eigen::Matrix<double,3,3> wx2 = wx*wx_copy;
//    Eigen::Matrix<double,9,3> dwx2_dw;
//    dwx2_dw <<   0, -2*w2, -2*w3,
//    w2,    w1,     0,
//    w3,     0,    w1,
//    w2,    w1,     0,
//    -2*w1,     0, -2*w3,
//    0,    w3,    w2,
//    w3,     0,    w1,
//    0,    w3,    w2,
//    -2*w1, -2*w2,     0;
//
//    double phi = sqrt(w1*w1 + w2*w2 + w3*w3);
//    Eigen::Matrix<double,9,3> dCwx2_dw;
//    double C;
//    if ((phi < 1e-6) && (phi > -1e-6))
//    {
//
//        dCwx2_dw.setZero();
//        C = 0.0;
//
//    }
//    else
//    {
//
//        double dC_dphi = -(4*cos(phi) + phi*sin(phi) + phi*phi - 4)/
//        (2*phi*phi*phi*(cos(phi) - 1.0));
//        Eigen::Matrix<double,1,3> dphi_dw;
//        dphi_dw << w1/sqrt(w1*w1 + w2*w2 + w3*w3),
//        w2/sqrt(w1*w1 + w2*w2 + w3*w3),
//        w3/sqrt(w1*w1 + w2*w2 + w3*w3);
//
//        Eigen::Matrix<double,1,3> dC_dw = dC_dphi*dphi_dw;
//        Eigen::Matrix<double,9,3> dC_dw_x_wx2;
//        dC_dw_x_wx2 << wx2(0,0)*dC_dw(0), wx2(0,0)*dC_dw(1), wx2(0,0)*dC_dw(2),
//        wx2(1,0)*dC_dw(0), wx2(1,0)*dC_dw(1), wx2(1,0)*dC_dw(2),
//        wx2(2,0)*dC_dw(0), wx2(2,0)*dC_dw(1), wx2(2,0)*dC_dw(2),
//        wx2(0,1)*dC_dw(0), wx2(0,1)*dC_dw(1), wx2(0,1)*dC_dw(2),
//        wx2(1,1)*dC_dw(0), wx2(1,1)*dC_dw(1), wx2(1,1)*dC_dw(2),
//        wx2(2,1)*dC_dw(0), wx2(2,1)*dC_dw(1), wx2(2,1)*dC_dw(2),
//        wx2(0,2)*dC_dw(0), wx2(0,2)*dC_dw(1), wx2(0,2)*dC_dw(2),
//        wx2(1,2)*dC_dw(0), wx2(1,2)*dC_dw(1), wx2(1,2)*dC_dw(2),
//        wx2(2,2)*dC_dw(0), wx2(2,2)*dC_dw(1), wx2(2,2)*dC_dw(2);
//
//        double A1 = sin(phi)/phi;
//        double B1 = (1.0 - cos(phi))/(phi*phi);
//        C = (1/(phi*phi))*(1 - A1/(2*B1));
//        dCwx2_dw = dC_dw_x_wx2 + C*dwx2_dw;
//
//    }
//
//    Eigen::Matrix<double,9,3> dwx_dw;
//    dwx_dw << 0,  0,  0,
//    0,  0,  1,
//    0, -1,  0,
//    0,  0, -1,
//    0,  0,  0,
//    1,  0,  0,
//    0,  1,  0,
//    -1,  0,  0,
//    0,  0,  0;
//
//    Eigen::Matrix<double,9,1> du_dw_vec;
//    du_dw_vec = (-0.5*dwx_dw + dCwx2_dw)*t;
//
//    Eigen::Matrix<double,3,3> du_dw;
//    du_dw << du_dw_vec(0), du_dw_vec(3), du_dw_vec(6),
//    du_dw_vec(1), du_dw_vec(4), du_dw_vec(7),
//    du_dw_vec(2), du_dw_vec(5), du_dw_vec(8);
//
//    // Derivative of u (translation) w.r.t. R (rotation)
//    Eigen::Matrix<double,3,9> du_dR = du_dw*dw_dR;
//
//    // Derivative of u (translation) w.r.t. t (translation)
//    Eigen::Matrix<double,3,3> du_dt;
//    Eigen::Matrix<double,3,3> V_inv;
//    V_inv = Eigen::MatrixXd::Identity(3,3) - 0.5*wx + C*wx2;
//    du_dt = V_inv;
//
//    //----------------------------------------------------------------------
//    // Set dXiA_dA from du_dR, du_dt, dw_dR and dw_dt,
//    // where A = [R t].
//    //----------------------------------------------------------------------
//    Eigen::Matrix<double,6,12> M;
//    M << du_dR, du_dt,
//    dw_dR, dw_dt;
//
//    return M;
//
//}

Eigen::Matrix<double,6,12>
DervLie::
dXiA_dA(Eigen::Matrix<double,4,4> A)
{
    
    Eigen::Matrix<double,3,3,Eigen::ColMajor> R = A.block(0,0,3,3);
    Eigen::Matrix<double,3,1> t = A.block(0,3,3,1);
    
    double alpha = (R.trace() - 1.0)/2.0;
    //printf("alpha = %f, acos(alpha) = %f\n", alpha, acos(alpha));
    
    // dw_dR : Derivative of w (orientation) w.r.t. R (rotation)
    Eigen::Matrix<double,3,9> dw_dR_ = dw_dR_matlab(R);
    
    // dw_dt : Derivative of w (orientation) w.r.t. t (translation)
    Eigen::Matrix<double,3,3> dw_dt;
    dw_dt.setZero(3,3);
    
    Sophus::SO3d lieOmega = Sophus::SO3d(R);
    Eigen::Matrix<double,3,1> w = lieOmega.log();
    
    Eigen::Matrix<double,3,3> wx = skewsym(w);
    Eigen::Matrix<double,3,3> wx2 = wx*wx;
    double theta = sqrt(w.transpose()*w);
    
//    printf("theta = %f\n", theta);
    
    Eigen::Matrix<double,3,3> V_inv;
    if (fabs(theta) < 1e-04)
    {
        
        V_inv = Eigen::MatrixXd::Identity(3,3);
        
    }
    else
    {
        
        double A1, B, G;
//        if (fabs(theta) < 1e-05)
//        {
//            A1 = 1.0 - (theta*theta)/6.0 + pow(theta,4.0)/120.0;
//            B = 0.5 - (theta*theta)/24.0 + pow(theta,4.0)/720.0;
//            G = 1.0/12.0 + (theta*theta)/720.0 + pow(theta,4.0)/30240.0;
//        }
//        else
//        {
//            
            A1 = sin(theta)/theta;
            B = (1.0 - cos(theta))/(theta*theta);
            //G = (1.0/(theta*theta))*(1.0 - A1/(2.0*B));
            G = (1.0/(theta*theta))*(1.0 - theta/(2.0*tan(theta/2.0))); // As Sophus
        
//        }
    
        V_inv = Eigen::MatrixXd::Identity(3,3) - 0.5*wx + G*wx2;
        
    }
    
    
    Eigen::Matrix<double,3,1> u = V_inv*t;
    Eigen::Matrix<double,3,3> du_dw_ = du_dw(u, w);
    
    // du_dR : Derivative of u (translation) w.r.t. R (rotation)
    Eigen::Matrix<double,3,9> du_dR = du_dw_*dw_dR_;
    
    // du_dt : Derivative of u (translation) w.r.t. t (translation)
    Eigen::Matrix<double,3,3> du_dt;
    du_dt = V_inv;
    
    //----------------------------------------------------------------------
    // Set dXiA_dA from du_dR, du_dt, dw_dR and dw_dt,
    // where A = [R t].
    //----------------------------------------------------------------------
    Eigen::Matrix<double,6,12> M;
    M << du_dR, du_dt,
    dw_dR_, dw_dt;
    
    return M;
    
}

//Eigen::Matrix<double,6,12>
//DervLie::
//dXiA_dA_at_Id(Eigen::Matrix<double,4,4> A)
//{
//    
//    return DervLie::dXiA_dA(A)*dEps_plus_D_dXiEps(A);
//    
//}
    
//------------------------------------------------------------------------------
Eigen::Matrix<double,6,7>
DervLie::
dXiA_dXiA7(Eigen::Matrix<double,6,1> XiA)
{
    
    Eigen::Matrix<double,6,7> M;
    Eigen::Matrix<double,3,4> dw_dq;
    Sophus::SE3d LieXiA = Sophus::SE3d::exp(XiA);
    double qw = LieXiA.so3().unit_quaternion().w();
    double qx = LieXiA.so3().unit_quaternion().x();
    double qy = LieXiA.so3().unit_quaternion().y();
    double qz = LieXiA.so3().unit_quaternion().z();
    
    dw_dq(0,0) = (1.0/(qw*qw)*qx*-2.0)/(1.0/(qw*qw)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))+1.0);
    
    dw_dq(0,1) = atan(sqrt(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))/qw)*1.0/sqrt(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))*2.0-qx*fabs(qx)*((qx/fabs(qx)))*atan(sqrt(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))/qw)*1.0/pow(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0),3.0/2.0)*2.0+(qx*fabs(qx)*((qx/fabs(qx)))*2.0)/(qw*(1.0/(qw*qw)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))+1.0)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0)));
    
    dw_dq(0,2) = qx*fabs(qy)*((qy/fabs(qy)))*atan(sqrt(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))/qw)*1.0/pow(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0),3.0/2.0)*-2.0+(qx*fabs(qy)*((qy/fabs(qy)))*2.0)/(qw*(1.0/(qw*qw)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))+1.0)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0)));
    
    dw_dq(0,3) = qx*fabs(qz)*((qz/fabs(qz)))*atan(sqrt(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))/qw)*1.0/pow(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0),3.0/2.0)*-2.0+(qx*fabs(qz)*((qz/fabs(qz)))*2.0)/(qw*(1.0/(qw*qw)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))+1.0)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0)));
    
    dw_dq(1,0) = (1.0/(qw*qw)*qy*-2.0)/(1.0/(qw*qw)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))+1.0);
    
    dw_dq(1,1) = qy*fabs(qx)*((qx/fabs(qx)))*atan(sqrt(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))/qw)*1.0/pow(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0),3.0/2.0)*-2.0+(qy*fabs(qx)*((qx/fabs(qx)))*2.0)/(qw*(1.0/(qw*qw)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))+1.0)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0)));
    
    dw_dq(1,2) = atan(sqrt(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))/qw)*1.0/sqrt(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))*2.0-qy*fabs(qy)*((qy/fabs(qy)))*atan(sqrt(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))/qw)*1.0/pow(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0),3.0/2.0)*2.0+(qy*fabs(qy)*((qy/fabs(qy)))*2.0)/(qw*(1.0/(qw*qw)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))+1.0)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0)));
    
    dw_dq(1,3) = qy*fabs(qz)*((qz/fabs(qz)))*atan(sqrt(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))/qw)*1.0/pow(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0),3.0/2.0)*-2.0+(qy*fabs(qz)*((qz/fabs(qz)))*2.0)/(qw*(1.0/(qw*qw)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))+1.0)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0)));
    
    dw_dq(2,0) = (1.0/(qw*qw)*qz*-2.0)/(1.0/(qw*qw)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))+1.0);
    
    dw_dq(2,1) = qz*fabs(qx)*((qx/fabs(qx)))*atan(sqrt(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))/qw)*1.0/pow(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0),3.0/2.0)*-2.0+(qz*fabs(qx)*((qx/fabs(qx)))*2.0)/(qw*(1.0/(qw*qw)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))+1.0)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0)));
    
    dw_dq(2,2) = qz*fabs(qy)*((qy/fabs(qy)))*atan(sqrt(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))/qw)*1.0/pow(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0),3.0/2.0)*-2.0+(qz*fabs(qy)*((qy/fabs(qy)))*2.0)/(qw*(1.0/(qw*qw)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))+1.0)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0)));
    
    dw_dq(2,3) = atan(sqrt(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))/qw)*1.0/sqrt(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))*2.0-qz*fabs(qz)*((qz/fabs(qz)))*atan(sqrt(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))/qw)*1.0/pow(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0),3.0/2.0)*2.0+(qz*fabs(qz)*((qz/fabs(qz)))*2.0)/(qw*(1.0/(qw*qw)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0))+1.0)*(pow(fabs(qx),2.0)+pow(fabs(qy),2.0)+pow(fabs(qz),2.0)));
    
    M.block(0,0,3,3) = Eigen::Matrix<double,3,3>::Identity();
    M.block(3,0,3,4) = dw_dq;
    
    return M;
    
}

//------------------------------------------------------------------------------
Eigen::Matrix<double,7,6>
DervLie::
dXiA7_dXiA(Eigen::Matrix<double,6,1> XiA)
{

    Eigen::Matrix<double,7,6> M;
    double w1 = XiA(3);
    double w2 = XiA(4);
    double w3 = XiA(5);
    Eigen::Matrix<double,4,3> dq_dw;
    
    double pw2 = pow(fabs(w1),2.0)+pow(fabs(w2),2.0)+pow(fabs(w3),2.0);
    double sw12 = sqrt(pw2)*(1.0/2.0);
    
    dq_dw(0,0) = sin(sw12)*w1*1.0/sqrt(pw2)*(-1.0/2.0);
    dq_dw(0,1) = sin(sw12)*w2*1.0/sqrt(pw2)*(-1.0/2.0);
    dq_dw(0,2) = sin(sw12)*w3*1.0/sqrt(pw2)*(-1.0/2.0);
    dq_dw(1,0) = sin(sw12)*1.0/sqrt(pw2)+(w1*cos(sw12)*w1*(1.0/2.0))/(pw2)-w1*sin(sw12)*w1*1.0/pow(pw2,3.0/2.0);
    dq_dw(1,1) = (w1*cos(sw12)*w2*(1.0/2.0))/(pw2)-w1*sin(sw12)*w2*1.0/pow(pw2,3.0/2.0);
    dq_dw(1,2) = (w1*cos(sw12)*w3*(1.0/2.0))/(pw2)-w1*sin(sw12)*w3*1.0/pow(pw2,3.0/2.0);
    dq_dw(2,0) = (w2*cos(sw12)*w1*(1.0/2.0))/(pw2)-w2*sin(sw12)*w1*1.0/pow(pw2,3.0/2.0);
    dq_dw(2,1) = sin(sw12)*1.0/sqrt(pw2)+(w2*cos(sw12)*w2*(1.0/2.0))/(pw2)-w2*sin(sw12)*w2*1.0/pow(pw2,3.0/2.0);
    dq_dw(2,2) = (w2*cos(sw12)*w3*(1.0/2.0))/(pw2)-w2*sin(sw12)*w3*1.0/pow(pw2,3.0/2.0);
    dq_dw(3,0) = (w3*cos(sw12)*w1*(1.0/2.0))/(pw2)-w3*sin(sw12)*w1*1.0/pow(pw2,3.0/2.0);
    dq_dw(3,1) = (w3*cos(sw12)*w2*(1.0/2.0))/(pw2)-w3*sin(sw12)*w2*1.0/pow(pw2,3.0/2.0);
    dq_dw(3,2) = sin(sw12)*1.0/sqrt(pw2)+(w3*cos(sw12)*w3*(1.0/2.0))/(pw2)-w3*sin(sw12)*w3*1.0/pow(pw2,3.0/2.0);
    
    M.block(0,0,3,3) = Eigen::Matrix<double,3,3>::Identity();
    M.block(3,0,4,3) = dq_dw;
    
    return M;
    
}
    
//------------------------------------------------------------------------------
Eigen::Matrix<double,12,6>
DervLie::
dD_plus_Eps_dXiEps(Eigen::Matrix<double,4,4> D)
{
    // Jacobian of a pose D in SE3 w.r.t. the lie algebra taken at XiEps = 0
    
    Eigen::Matrix<double,3,3> R = D.block(0, 0, 3, 3);
    Eigen::Matrix<double,3,1> dc1 = D.block(0, 0, 3, 1);
    Eigen::Matrix<double,3,1> dc2 = D.block(0, 1, 3, 1);
    Eigen::Matrix<double,3,1> dc3 = D.block(0, 2, 3, 1);
    
    // From Eq (10.19) in A tutorial of SE3 transformation.. by J. Blanco
    Eigen::Matrix<double,12,6> M;
    M <<      0,      0,      0,       0, -dc3(0),  dc2(0),
              0,      0,      0,       0, -dc3(1),  dc2(1),
              0,      0,      0,       0, -dc3(2),  dc2(2),
              0,      0,      0,  dc3(0),       0, -dc1(0),
              0,      0,      0,  dc3(1),       0, -dc1(1),
              0,      0,      0,  dc3(2),       0, -dc1(2),
              0,      0,      0, -dc2(0),  dc1(0),       0,
              0,      0,      0, -dc2(1),  dc1(1),       0,
              0,      0,      0, -dc2(2),  dc1(2),       0,
         R(0,0), R(0,1), R(0,2),       0,       0,       0,
         R(1,0), R(1,1), R(1,2),       0,       0,       0,
         R(2,0), R(2,1), R(2,2),       0,       0,       0;
    
    return M;
    
}
    
Eigen::Matrix<double,12,6>
DervLie::
dEps_plus_D_dXiEps(Eigen::Matrix<double,4,4> D)
{
    // Jacobian of a pose D in SE3 w.r.t. the lie algebra taken at XiEps = 0
    
    Eigen::Matrix<double,3,1> dc1 = D.block(0, 0, 3, 1);
    Eigen::Matrix<double,3,1> dc2 = D.block(0, 1, 3, 1);
    Eigen::Matrix<double,3,1> dc3 = D.block(0, 2, 3, 1);
    Eigen::Matrix<double,3,1> dt = D.block(0, 3, 3, 1);
    Eigen::Matrix<double,3,3> skew_dc1 = skewsym(dc1);
    Eigen::Matrix<double,3,3> skew_dc2 = skewsym(dc2);
    Eigen::Matrix<double,3,3> skew_dc3 = skewsym(dc3);
    Eigen::Matrix<double,3,3> skew_dt = skewsym(dt);
    
    // From Eq (10.15) in A tutorial of SE3 transformation.. by J. Blanco
    Eigen::Matrix<double,12,6> M;
    M <<
    0,      0,      0,  -skew_dc1(0,0),  -skew_dc1(0,1),  -skew_dc1(0,2),
    0,      0,      0,  -skew_dc1(1,0),  -skew_dc1(1,1),  -skew_dc1(1,2),
    0,      0,      0,  -skew_dc1(2,0),  -skew_dc1(2,1),  -skew_dc1(2,2),
    
    0,      0,      0,  -skew_dc2(0,0),  -skew_dc2(0,1),  -skew_dc2(0,2),
    0,      0,      0,  -skew_dc2(1,0),  -skew_dc2(1,1),  -skew_dc2(1,2),
    0,      0,      0,  -skew_dc2(2,0),  -skew_dc2(2,1),  -skew_dc2(2,2),
    
    0,      0,      0,  -skew_dc3(0,0),  -skew_dc3(0,1),  -skew_dc3(0,2),
    0,      0,      0,  -skew_dc3(1,0),  -skew_dc3(1,1),  -skew_dc3(1,2),
    0,      0,      0,  -skew_dc3(2,0),  -skew_dc3(2,1),  -skew_dc3(2,2),
    
    1,      0,      0,   -skew_dt(0,0),   -skew_dt(0,1),   -skew_dt(0,2),
    0,      1,      0,   -skew_dt(1,0),   -skew_dt(1,1),   -skew_dt(1,2),
    0,      0,      1,   -skew_dt(2,0),   -skew_dt(2,1),   -skew_dt(2,2);
    
    return M;
    
}
    
//------------------------------------------------------------------------------
// dw_dR (Derivative of lie algebra so(3) wrt rotation matrix
//        through quaternion as Sophus library implementation)
//------------------------------------------------------------------------------
Eigen::Matrix<double,3,9>
DervLie::
dw_dR_matlab(Eigen::Matrix<double,3,3> R)
{
    
    double r11, r12, r13;
    double r21, r22, r23;
    double r31, r32, r33;
    r11=R(0,0); r12=R(0,1); r13=R(0,2);
    r21=R(1,0); r22=R(1,1); r23=R(1,2);
    r31=R(2,0); r32=R(2,1); r33=R(2,2);
    double t = sqrt(r11+r22+r33+1.0);
    double f_r12_r21 = fabs(r12-r21);
    double f_r13_r31 = fabs(r13-r31);
    double f_r23_r32 = fabs(r23-r32);
    
    // Handling divided by zero
    double r12_r21_over_f_r12_r21;
    if (f_r12_r21 < 1e-6)
    {
        r12_r21_over_f_r12_r21 = 1.0;
    }
    else
    {
        r12_r21_over_f_r12_r21 = (r12-r21)/f_r12_r21;
    }
    double r13_r31_over_f_r13_r31;
    if (f_r13_r31 < 1e-6)
    {
        r13_r31_over_f_r13_r31 = 1.0;
    }
    else
    {
        r13_r31_over_f_r13_r31 = (r13-r31)/f_r13_r31;
    }
    double r23_r32_over_f_r23_r32;
    if (f_r23_r32 < 1e-6)
    {
        r23_r32_over_f_r23_r32 = 1.0;
    }
    else
    {
        r23_r32_over_f_r23_r32 = (r23-r32)/f_r23_r32;
    }
    
    // dw1_dR
    Eigen::Matrix<double,1,9> dw1_dR;
    
    // Handling divided by zero when r12==r21, r13==r31, r23==r32
    double alpha = (R.trace() - 1.0)/2.0;
    double theta = acos(alpha);
    
#if DEBUG
    printf("[dw_dR_matlab] theta = %f\n", theta);
    printf("[dw_dR_matlab] f_r12_r21 = %f\n", f_r12_r21);
    printf("[dw_dR_matlab] f_r13_r31 = %f\n", f_r13_r31);
    printf("[dw_dR_matlab] f_r23_r32 = %f\n", f_r23_r32);
#endif
    
    if ((f_r12_r21 < 1e-6) && (f_r13_r31 < 1e-6) && (f_r23_r32 < 1e-6)
        && fabs(theta) < 1e-6)
    {
        
        Eigen::Matrix<double,3,9> M;
        M.setZero();
        M(2,1) = 0.5;
        M(1,2) = -0.5;
        M(2,3) = -0.5;
        M(0,5) = 0.5;
        M(1,6) = 0.5;
        M(0,7) = -0.5;
        
        return M;
        
    }
    
    dw1_dR(0,0) = 1.0/pow(t,3.0)*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*(r23-r32)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*(1.0/2.0)+((sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/pow(r11+r22+r33+1.0,3.0/2.0)+1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*((pow(f_r12_r21,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r13_r31,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r23_r32,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t)*1.0/t)*(r23-r32)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0))-(atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*(r23-r32)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*((pow(f_r12_r21,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r13_r31,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r23_r32,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t)*(1.0/2.0))/t;
    dw1_dR(0,1) = (f_r12_r21*((r12_r21_over_f_r12_r21))*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/pow(t,2.0)*(r23-r32)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*(-1.0/4.0))/t+(f_r12_r21*((r12_r21_over_f_r12_r21))*1.0/pow(t,2.0)*(r23-r32)*1.0/t*(1.0/2.0))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0)*(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)));
    dw1_dR(0,2) = (f_r13_r31*((r13_r31_over_f_r13_r31))*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/pow(t,2.0)*(r23-r32)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*(-1.0/4.0))/t+(f_r13_r31*((r13_r31_over_f_r13_r31))*1.0/pow(t,2.0)*(r23-r32)*1.0/t*(1.0/2.0))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0)*(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)));
    dw1_dR(0,3) = (f_r12_r21*((r12_r21_over_f_r12_r21))*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/pow(t,2.0)*(r23-r32)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*(1.0/4.0))/t-(f_r12_r21*((r12_r21_over_f_r12_r21))*1.0/pow(t,2.0)*(r23-r32)*1.0/t*(1.0/2.0))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0)*(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)));
    dw1_dR(0,4) = 1.0/pow(t,3.0)*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*(r23-r32)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*(1.0/2.0)+((sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/pow(r11+r22+r33+1.0,3.0/2.0)+1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*((pow(f_r12_r21,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r13_r31,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r23_r32,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t)*1.0/t)*(r23-r32)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0))-(atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*(r23-r32)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*((pow(f_r12_r21,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r13_r31,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r23_r32,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t)*(1.0/2.0))/t;
    dw1_dR(0,5) = (atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)))/t-(f_r23_r32*((r23_r32_over_f_r23_r32))*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/pow(t,2.0)*(r23-r32)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*(1.0/4.0))/t+(f_r23_r32*((r23_r32_over_f_r23_r32))*1.0/pow(t,2.0)*(r23-r32)*1.0/t*(1.0/2.0))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0)*(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)));
    dw1_dR(0,6) = (f_r13_r31*((r13_r31_over_f_r13_r31))*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/pow(t,2.0)*(r23-r32)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*(1.0/4.0))/t-(f_r13_r31*((r13_r31_over_f_r13_r31))*1.0/pow(t,2.0)*(r23-r32)*1.0/t*(1.0/2.0))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0)*(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)));
    dw1_dR(0,7) = -(atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)))/t+(f_r23_r32*((r23_r32_over_f_r23_r32))*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/pow(t,2.0)*(r23-r32)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*(1.0/4.0))/t-(f_r23_r32*((r23_r32_over_f_r23_r32))*1.0/pow(t,2.0)*(r23-r32)*1.0/t*(1.0/2.0))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0)*(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)));
    dw1_dR(0,8) = 1.0/pow(t,3.0)*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*(r23-r32)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*(1.0/2.0)+((sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/pow(r11+r22+r33+1.0,3.0/2.0)+1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*((pow(f_r12_r21,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r13_r31,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r23_r32,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t)*1.0/t)*(r23-r32)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0))-(atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*(r23-r32)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*((pow(f_r12_r21,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r13_r31,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r23_r32,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t)*(1.0/2.0))/t;
    
    // dw2_dR
    Eigen::Matrix<double,1,9> dw2_dR;
    dw2_dR(0,0) = 1.0/pow(t,3.0)*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*(r13-r31)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*(-1.0/2.0)-((sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/pow(r11+r22+r33+1.0,3.0/2.0)+1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*((pow(f_r12_r21,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r13_r31,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r23_r32,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t)*1.0/t)*(r13-r31)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0))+(atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*(r13-r31)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*((pow(f_r12_r21,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r13_r31,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r23_r32,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t)*(1.0/2.0))/t;
    dw2_dR(0,1) = (f_r12_r21*((r12_r21_over_f_r12_r21))*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/pow(t,2.0)*(r13-r31)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*(1.0/4.0))/t-(f_r12_r21*((r12_r21_over_f_r12_r21))*1.0/pow(t,2.0)*(r13-r31)*1.0/t*(1.0/2.0))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0)*(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)));
    dw2_dR(0,2) = -(atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)))/t+(f_r13_r31*((r13_r31_over_f_r13_r31))*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/pow(t,2.0)*(r13-r31)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*(1.0/4.0))/t-(f_r13_r31*((r13_r31_over_f_r13_r31))*1.0/pow(t,2.0)*(r13-r31)*1.0/t*(1.0/2.0))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0)*(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)));
    dw2_dR(0,3) = (f_r12_r21*((r12_r21_over_f_r12_r21))*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/pow(t,2.0)*(r13-r31)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*(-1.0/4.0))/t+(f_r12_r21*((r12_r21_over_f_r12_r21))*1.0/pow(t,2.0)*(r13-r31)*1.0/t*(1.0/2.0))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0)*(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)));
    dw2_dR(0,4) = 1.0/pow(t,3.0)*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*(r13-r31)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*(-1.0/2.0)-((sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/pow(r11+r22+r33+1.0,3.0/2.0)+1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*((pow(f_r12_r21,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r13_r31,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r23_r32,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t)*1.0/t)*(r13-r31)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0))+(atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*(r13-r31)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*((pow(f_r12_r21,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r13_r31,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r23_r32,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t)*(1.0/2.0))/t;
    dw2_dR(0,5) = (f_r23_r32*((r23_r32_over_f_r23_r32))*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/pow(t,2.0)*(r13-r31)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*(1.0/4.0))/t-(f_r23_r32*((r23_r32_over_f_r23_r32))*1.0/pow(t,2.0)*(r13-r31)*1.0/t*(1.0/2.0))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0)*(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)));
    dw2_dR(0,6) = (atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)))/t-(f_r13_r31*((r13_r31_over_f_r13_r31))*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/pow(t,2.0)*(r13-r31)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*(1.0/4.0))/t+(f_r13_r31*((r13_r31_over_f_r13_r31))*1.0/pow(t,2.0)*(r13-r31)*1.0/t*(1.0/2.0))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0)*(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)));
    dw2_dR(0,7) = (f_r23_r32*((r23_r32_over_f_r23_r32))*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/pow(t,2.0)*(r13-r31)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*(-1.0/4.0))/t+(f_r23_r32*((r23_r32_over_f_r23_r32))*1.0/pow(t,2.0)*(r13-r31)*1.0/t*(1.0/2.0))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0)*(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)));
    dw2_dR(0,8) = 1.0/pow(t,3.0)*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*(r13-r31)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*(-1.0/2.0)-((sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/pow(r11+r22+r33+1.0,3.0/2.0)+1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*((pow(f_r12_r21,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r13_r31,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r23_r32,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t)*1.0/t)*(r13-r31)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0))+(atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*(r13-r31)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*((pow(f_r12_r21,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r13_r31,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r23_r32,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t)*(1.0/2.0))/t;
    
    // dw3_dR
    Eigen::Matrix<double,1,9> dw3_dR;
    dw3_dR(0,0) = 1.0/pow(t,3.0)*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*(r12-r21)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*(1.0/2.0)+((sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/pow(r11+r22+r33+1.0,3.0/2.0)+1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*((pow(f_r12_r21,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r13_r31,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r23_r32,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t)*1.0/t)*(r12-r21)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0))-(atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*(r12-r21)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*((pow(f_r12_r21,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r13_r31,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r23_r32,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t)*(1.0/2.0))/t;
    dw3_dR(0,1) = (atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)))/t-(f_r12_r21*((r12_r21_over_f_r12_r21))*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/pow(t,2.0)*(r12-r21)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*(1.0/4.0))/t+(f_r12_r21*((r12_r21_over_f_r12_r21))*1.0/pow(t,2.0)*(r12-r21)*1.0/t*(1.0/2.0))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0)*(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)));
    dw3_dR(0,2) = (f_r13_r31*((r13_r31_over_f_r13_r31))*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/pow(t,2.0)*(r12-r21)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*(-1.0/4.0))/t+(f_r13_r31*((r13_r31_over_f_r13_r31))*1.0/pow(t,2.0)*(r12-r21)*1.0/t*(1.0/2.0))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0)*(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)));
    dw3_dR(0,3) = -(atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)))/t+(f_r12_r21*((r12_r21_over_f_r12_r21))*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/pow(t,2.0)*(r12-r21)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*(1.0/4.0))/t-(f_r12_r21*((r12_r21_over_f_r12_r21))*1.0/pow(t,2.0)*(r12-r21)*1.0/t*(1.0/2.0))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0)*(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)));
    dw3_dR(0,4) = 1.0/pow(t,3.0)*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*(r12-r21)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*(1.0/2.0)+((sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/pow(r11+r22+r33+1.0,3.0/2.0)+1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*((pow(f_r12_r21,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r13_r31,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r23_r32,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t)*1.0/t)*(r12-r21)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0))-(atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*(r12-r21)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*((pow(f_r12_r21,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r13_r31,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r23_r32,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t)*(1.0/2.0))/t;
    dw3_dR(0,5) = (f_r23_r32*((r23_r32_over_f_r23_r32))*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/pow(t,2.0)*(r12-r21)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*(-1.0/4.0))/t+(f_r23_r32*((r23_r32_over_f_r23_r32))*1.0/pow(t,2.0)*(r12-r21)*1.0/t*(1.0/2.0))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0)*(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)));
    dw3_dR(0,6) = (f_r13_r31*((r13_r31_over_f_r13_r31))*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/pow(t,2.0)*(r12-r21)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*(1.0/4.0))/t-(f_r13_r31*((r13_r31_over_f_r13_r31))*1.0/pow(t,2.0)*(r12-r21)*1.0/t*(1.0/2.0))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0)*(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)));
    dw3_dR(0,7) = (f_r23_r32*((r23_r32_over_f_r23_r32))*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*1.0/pow(t,2.0)*(r12-r21)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*(1.0/4.0))/t-(f_r23_r32*((r23_r32_over_f_r23_r32))*1.0/pow(t,2.0)*(r12-r21)*1.0/t*(1.0/2.0))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0)*(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)));
    dw3_dR(0,8) = 1.0/pow(t,3.0)*atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*(r12-r21)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*(1.0/2.0)+((sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/pow(r11+r22+r33+1.0,3.0/2.0)+1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*((pow(f_r12_r21,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r13_r31,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r23_r32,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t)*1.0/t)*(r12-r21)*1.0/sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0)))/(t*((pow(f_r12_r21,2.0)*1.0/pow(t,2.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0))/(r11+r22+r33+1.0)+1.0))-(atan(sqrt(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0))*1.0/t*2.0)*(r12-r21)*1.0/pow(pow(f_r12_r21,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r13_r31,2.0)*1.0/pow(t,2.0)*(1.0/4.0)+pow(f_r23_r32,2.0)*1.0/pow(t,2.0)*(1.0/4.0),3.0/2.0)*((pow(f_r12_r21,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r13_r31,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t+(pow(f_r23_r32,2.0)*1.0/pow(t,3.0)*(1.0/4.0))/t)*(1.0/2.0))/t;
    
    Eigen::Matrix<double,3,9> M;
    M.block(0,0,1,9) = dw1_dR;
    M.block(1,0,1,9) = dw2_dR;
    M.block(2,0,1,9) = dw3_dR;
    
    return M;
    
}
    
Eigen::Matrix<double,12,6>
DervLie::
dsBC_dXiB(double s,
          Eigen::Matrix<double,6,1> XiB,
          Eigen::Matrix<double,6,1> XiC)
{
    
    Eigen::Matrix<double,12,6> M;
    
    M = s*dBC_dXiB(XiB, XiC);
    
    return M;
    
}
    
//--------------------------------------------------------------------------
// 3x9 Fast Jacob_x_Mat
// Verified
Eigen::Matrix<double,3,9>
DervLie::
fast3x9_Jacob9x9_x_Mat3x1(const Eigen::Ref<const
                          Eigen::Matrix<double,9,9> >& J,
                          const Eigen::Matrix<double,3,1> M)
{
    
    Eigen::Matrix<double,3,9> Z;
    
    // From matlab
    // reshape(J(:,1),3,3)*M
    // ans =
    // j1_1*m1 + j4_1*m2 + j7_1*m3
    // j2_1*m1 + j5_1*m2 + j8_1*m3
    // j3_1*m1 + j6_1*m2 + j9_1*m3
    for (int i=0; i<9; i++)
    {
        
        Z(0,i) = J(0,i)*M(0,0) + J(3,i)*M(1,0) + J(6,i)*M(2,0);
        Z(1,i) = J(1,i)*M(0,0) + J(4,i)*M(1,0) + J(7,i)*M(2,0);
        Z(2,i) = J(2,i)*M(0,0) + J(5,i)*M(1,0) + J(8,i)*M(2,0);
        
    }
    
    return Z;
    
}

} // end of namespace lsd_slam