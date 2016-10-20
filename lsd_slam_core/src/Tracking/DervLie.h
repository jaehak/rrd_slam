//
//  DervLie.h
//
//  Created by Jae-Hak Kim on 16/11/2015.
//
//  Derivative of Lie Algebra and Lie Group
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

#ifndef __lsd_slam_core__DervLie__
#define __lsd_slam_core__DervLie__

#include <iostream>
#include <stdio.h>
#include "sophus/se3.hpp"
#include <unsupported/Eigen/KroneckerProduct>

namespace lsd_slam
{

class DervLie
{
    
public:
        
    //--------------------------------------------------------------------------
    // Verified Jacobians and functions

    //--------------------------------------------------------------------------
    // dR_dw [9 x 3]: Verified by Ceres
    static
    Eigen::Matrix<double,9,3> dR_dw(Eigen::Matrix<double,3,1> w) ;
    
    //--------------------------------------------------------------------------
    // dR_dw [9 x 3] matlab version
    static
    Eigen::Matrix<double,9,3>
    dR_dw_matlab(Eigen::Matrix<double,3,1> w);

    //--------------------------------------------------------------------------
    // dt_dw [3 x 3]: Verified by Ceres
    static
    Eigen::Matrix<double,3,3> dt_dw(Eigen::Matrix<double,3,1> u,
                                    Eigen::Matrix<double,3,1> w) ;
    
    //--------------------------------------------------------------------------
    // dt_dw (Right exp^{eps} multiplication) [3 x 3]: Verified by Ceres
    static
    Eigen::Matrix<double,3,3> dt_dw_right(Eigen::Matrix<double,3,1> u,
                                         Eigen::Matrix<double,3,1> w) ;
    //--------------------------------------------------------------------------
    // dt_du [3 x 3] : Verified by Ceres
    static
    Eigen::Matrix<double,3,3> dt_du(Eigen::Matrix<double,3,1> w) ;

    //--------------------------------------------------------------------------
    // dt_du (Right exp^{eps} mult) [3 x 3]
    static
    Eigen::Matrix<double,3,3> dt_du_right(Eigen::Matrix<double,3,1> w);
    
    //--------------------------------------------------------------------------
    // dAB_dA [12 x 12]
    static
    Eigen::Matrix<double,12,12>
    dAB_dA(Eigen::Matrix<double,4,4> A, Eigen::Matrix<double,4,4> B);

    //--------------------------------------------------------------------------
    // dA_dXiA [12 x 6] : Verified by Ceres
    static
    Eigen::Matrix<double,12,6> dA_dXiA(Eigen::Matrix<double,6,1> Xi) ;
    
    //--------------------------------------------------------------------------
    // dA_dXiA at identity [12 x 6]
    static Eigen::Matrix<double,12,6>
    dA_dXiA_at_Id(Eigen::Matrix<double,6,1> Xi);
    
    //--------------------------------------------------------------------------
    // dA_dXiA (Right exp^{eps} multiplication) [12 x 6]
    static
    Eigen::Matrix<double,12,6>
    dA_dXiA_right(Eigen::Matrix<double,6,1> Xi); // Left matrix multiplication
    
    //--------------------------------------------------------------------------
    // d(A^{-1}*(-eps)*D)/d(eps) Equation (10.28) in J. Blanco's tutorial
    static
    Eigen::Matrix<double,12,6>
    dAinv_negEps_D_deps_at_Id(Eigen::Matrix<double,6,1> XiA,
                              Eigen::Matrix<double,6,1> XiD);
    
    //--------------------------------------------------------------------------
    // d(A*x)/d(A)
    static
    Eigen::Matrix<double,3,12>
    dAx_dA(Eigen::Matrix<double,4,4> A, Eigen::Matrix<double,4,1> x);
    
    
    #if 1 // DO_NOT_USE_dw_dR
    //--------------------------------------------------------------------------
    // dw_dR [3 x 9] : Verified (not exact, homogeneous coords in rotation)
    static
    Eigen::Matrix<double,3,9> dw_dR_ver1(Eigen::Matrix<double,3,3> R) ;
    #endif

    //--------------------------------------------------------------------------
    // du_dw [3 x 3] : Verified
    static
    Eigen::Matrix<double,3,3> du_dw(Eigen::Matrix<double,3,1> u,
                                    Eigen::Matrix<double,3,1> w) ;

    #if 1 // DO_NOT_USE_dXiA_dA

    //--------------------------------------------------------------------------
    // dXiA_dA [6 x 12]: Verified
    static
    Eigen::Matrix<double,6,12> dXiA_dA_ver2(Eigen::Matrix<double,4,4> A);
        
    #endif

    //--------------------------------------------------------------------------
    // Verified working (deprecated, use Jacob_x_Mat instead)
    static
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>
    Jacob_times_Vec(const Eigen::Ref<const
                    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> >& J,
                    const Eigen::Matrix<double,Eigen::Dynamic,1> vec);

    //--------------------------------------------------------------------------
    // d(BC)_dXiB [12 x 6]: Verified by Ceres
    static
    Eigen::Matrix<double,12,6> dBC_dXiB(Eigen::Matrix<double,6,1> XiB,
                                        Eigen::Matrix<double,6,1> XiC);

    //--------------------------------------------------------------------------
    // d(BC)_dXiC [12 x 6]: Verified
    static
    Eigen::Matrix<double,12,6> dBC_dXiC(Eigen::Matrix<double,6,1> XiB,
                                        Eigen::Matrix<double,6,1> XiC);

    //--------------------------------------------------------------------------
    // This is a special multiplication of a matrix to a Jacobian matrix, which
    // is reshaped in order to have each column of the Jacobian matrix
    // may be correspond to the size of the matrix left multiplying
    static
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>
    Mat_x_Jacob(const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> M,
                const Eigen::Ref<const
                Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> >& J) ;

    //--------------------------------------------------------------------------
    // This is a special multiplication of a matrix to a Jacobian matrix, which
    // is reshaped in order to have each column of the Jacobian matrix
    // may be correspond to the size of the matrix left multiplying
    static
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>
    Jacob_x_Mat(
        const Eigen::Ref<const
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> >& J,
        const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> M);

    //--------------------------------------------------------------------------
    // Valid
    static
    Eigen::Matrix<double,3,3> skewsym(Eigen::Matrix<double,3,1> x);

    //--------------------------------------------------------------------------
    // Verified by Ceres
    static
    Eigen::Matrix<double,9,3> dRinv_dw(Eigen::Matrix<double,3,1> w);

    //--------------------------------------------------------------------------
    // Verified by Ceres
    static
    Eigen::Matrix<double,12,6> dAinv_dXiA(Eigen::Matrix<double,6,1> Xi);

    //--------------------------------------------------------------------------
    static
    Eigen::Matrix<double,12,6> dAinv_dXiA_right(Eigen::Matrix<double,6,1> Xi);

    //--------------------------------------------------------------------------
    static
    Eigen::Matrix<double,12,6> dAinv_dXiA_at_Id(Eigen::Matrix<double,6,1> Xi);

    //--------------------------------------------------------------------------
    // dAinv_dA (12 x 12)
    static
    Eigen::Matrix<double,12,12> dAinv_dA(Eigen::Matrix<double,4,4> A);
    
    //--------------------------------------------------------------------------
    // Verified by Ceres
    static
    Eigen::Matrix<double,12,6> dBinvC_dXiB(Eigen::Matrix<double,6,1> XiB,
                                           Eigen::Matrix<double,6,1> XiC);

    //--------------------------------------------------------------------------
    // Verified by Ceres
    static
    Eigen::Matrix<double,12,6> dBinvC_dXiB_at_Id(Eigen::Matrix<double,6,1> XiB,
                                           Eigen::Matrix<double,6,1> XiC);
    
    //--------------------------------------------------------------------------
    // Verified by Ceres
    static
    Eigen::Matrix<double,12,6> dBinvC_dXiC(Eigen::Matrix<double,6,1> XiB,
                                           Eigen::Matrix<double,6,1> XiC);

    //--------------------------------------------------------------------------
    static
    Eigen::Matrix<double,12,6> dBinvC_dXiC_at_Id(Eigen::Matrix<double,6,1> XiB,
                                           Eigen::Matrix<double,6,1> XiC);
    
    //--------------------------------------------------------------------------
    static
    Eigen::Matrix<double,12,12>
    dExp_k_XiAinv_Concat_XiB_dA(double k,
                                Eigen::Matrix<double,6,1> XiA,
                                Eigen::Matrix<double,6,1> XiB);
    
    //--------------------------------------------------------------------------
    // Verified (but not exact)
    // d(exp( k * (XiA^{-1} concat XiB)) / d XiA
    static
    Eigen::Matrix<double,12,6>
    dExp_k_XiAinv_Concat_XiB_dXiA(double k,
                                  Eigen::Matrix<double,6,1> XiA,
                                  Eigen::Matrix<double,6,1> XiB);

    //--------------------------------------------------------------------------
    static
    Eigen::Matrix<double,12,6>
    dExp_k_XiAinv_Concat_XiB_dXiA_at_Id(double k,
                                  Eigen::Matrix<double,6,1> XiA,
                                  Eigen::Matrix<double,6,1> XiB);

    //--------------------------------------------------------------------------
    // Verified (but not exact)
    // d(exp( k * (XiA^{-1} concat XiB)) / d XiB
    static
    Eigen::Matrix<double,12,6>
    dExp_k_XiAinv_Concat_XiB_dXiB(double k,
                                  Eigen::Matrix<double,6,1> XiA,
                                  Eigen::Matrix<double,6,1> XiB);

    //--------------------------------------------------------------------------
    static
    Eigen::Matrix<double,12,6>
    dExp_k_XiAinv_Concat_XiB_dXiB_at_Id(double k,
                                  Eigen::Matrix<double,6,1> XiA,
                                  Eigen::Matrix<double,6,1> XiB);
    
    //--------------------------------------------------------------------------
    // Verified by Ceres (exactly)
    // d(exp( k * (XiA^{-1} concat XiB)) / d k
    static
    Eigen::Matrix<double,12,1>
    dExp_k_XiAinv_Concat_XiB_dK(double k,
                                Eigen::Matrix<double,6,1> XiA,
                                Eigen::Matrix<double,6,1> XiB);

    //--------------------------------------------------------------------------
    // 3x3 Fast Jacob_x_Vec
    // Verified by Ceres
    static
    Eigen::Matrix<double,3,3>
    fast3x3_Jacob9x3_x_Vec3x1(const Eigen::Ref<const
                              Eigen::Matrix<double,9,3> >& J,
                              const Eigen::Matrix<double,3,1> V);

    //--------------------------------------------------------------------------
    // Verified by Ceres
    static
    Eigen::Matrix<double,3,6> dRXt_dXi(Eigen::Matrix<double,3,1> X,
                                       Eigen::Matrix<double,6,1> Xi);
    
    //--------------------------------------------------------------------------
    static
    Eigen::Matrix<double,12,6> dsBC_dXiB(double s,
                                         Eigen::Matrix<double,6,1> XiB,
                                         Eigen::Matrix<double,6,1> XiC);

    //--------------------------------------------------------------------------
    static
    Eigen::Matrix<double,12,12>
    fast12x12_multMatToJacob(const Eigen::Matrix<double,4,4> M,
                             const Eigen::Ref<const
                             Eigen::Matrix<double,12,12> >& J);

    //--------------------------------------------------------------------------
    static
    Eigen::Matrix<double,12,12>
    fast12x12_multJacobToMat(const Eigen::Ref<const
                             Eigen::Matrix<double,12,12> >& J,
                             const Eigen::Matrix<double,4,4> M);
    
    //--------------------------------------------------------------------------
    // Verified
    static
    Eigen::Matrix<double,12,6>
    fast12x6_multJacobToMat(const Eigen::Ref<const
                            Eigen::Matrix<double,12,6> >& J,
                            const Eigen::Matrix<double,4,4> M);

    //--------------------------------------------------------------------------
    // Verified
    static
    Eigen::Matrix<double,12,1>
    fast12x1_multJacobToMat(const Eigen::Ref<const
                            Eigen::Matrix<double,12,1> >& J,
                            const Eigen::Matrix<double,4,4> M);

    //--------------------------------------------------------------------------
    // Verified
    static
    Eigen::Matrix<double,12,6>
    fast12x6_multMatToJacob(const Eigen::Matrix<double,4,4> M,
                            const Eigen::Ref<const
                            Eigen::Matrix<double,12,6> >& J);

    // Verified
    static
    Eigen::Matrix<double,12,1>
    fast12x1_multMatToJacob(const Eigen::Matrix<double,4,4> M,
                            const Eigen::Ref<const
                            Eigen::Matrix<double,12,1> >& J);

    //--------------------------------------------------------------------------
    // 9x3 Fast Mat_x_Jacob
    // Verified
    static
    Eigen::Matrix<double,9,3>
    fast9x3_Mat3x3_x_Jacob9x3(const Eigen::Matrix<double,3,3> M,
                              const Eigen::Ref<const
                              Eigen::Matrix<double,9,3> >& J);

    //--------------------------------------------------------------------------
    // 9x3 Fast Jacob_x_Mat
    // Verified
    static
    Eigen::Matrix<double,9,3>
    fast9x3_Jacob9x3_x_Mat3x3(const Eigen::Ref<const
                              Eigen::Matrix<double,9,3> >& J,
                              const Eigen::Matrix<double,3,3> M);

    //--------------------------------------------------------------------------
    // 3x3 Fast Mat_x_Jacob
    // Verified
    static
    Eigen::Matrix<double,3,3>
    fast3x3_Mat3x3_x_Jacob3x3(const Eigen::Matrix<double,3,3> M,
                              const Eigen::Ref<const
                              Eigen::Matrix<double,3,3> >& J);

    //--------------------------------------------------------------------------
    // 3x9 Fast Jacob_x_Mat
    // Verified
    static
    Eigen::Matrix<double,3,9>
    fast3x9_Jacob9x9_x_Mat3x1(const Eigen::Ref<const
                              Eigen::Matrix<double,9,9> >& J,
                              const Eigen::Matrix<double,3,1> M);
    
    //--------------------------------------------------------------------------
    // Verified
    static
    Eigen::Matrix<double,6,12> dXiA_dA(Eigen::Matrix<double,4,4> A);

    //--------------------------------------------------------------------------
    static
    Eigen::Matrix<double,6,7> dXiA_dXiA7(Eigen::Matrix<double,6,1> XiA);
    
    //--------------------------------------------------------------------------
    static
    Eigen::Matrix<double,7,6> dXiA7_dXiA(Eigen::Matrix<double,6,1> XiA);
    
//    //--------------------------------------------------------------------------
//    static
//    Eigen::Matrix<double,6,12> dXiA_dA_at_Id(Eigen::Matrix<double,4,4> A);
    
#if 1// UNVERIFIED

    static
    Eigen::Matrix<double,12,6>
    addJacobToMat(Eigen::Matrix<double,12,6> J,
                  Eigen::Matrix<double,4,4> M);

    static
    Eigen::Matrix<double,12,Eigen::Dynamic>
    multJacobToMat(const Eigen::Ref<const
                   Eigen::Matrix<double,12,Eigen::Dynamic> >& J,
                   const Eigen::Matrix<double,4,4> M);

    static
    Eigen::Matrix<double,12,Eigen::Dynamic>
    multMatToJacob(const Eigen::Matrix<double,4,4> M,
                   const Eigen::Ref<const
                   Eigen::Matrix<double,12,Eigen::Dynamic> >& J);

    //--------------------------------------------------------------------------

    static
    Eigen::Matrix<double,12,6> dA_dXiB(Eigen::Matrix<double,4,4> A);

    static
    Eigen::Matrix<double,12,6> dA_dXiC(Eigen::Matrix<double,4,4> A,
                                       Eigen::Matrix<double,4,4> B);
    
    static
    Eigen::Matrix<double,12,6>
    dD_plus_Eps_dXiEps(Eigen::Matrix<double,4,4> D);

    static
    Eigen::Matrix<double,12,6>
    dEps_plus_D_dXiEps(Eigen::Matrix<double,4,4> D);

    static
    Eigen::Matrix<double,3,9>
    dw_dR_matlab(Eigen::Matrix<double,3,3> R);

#endif
    
    //--------------------------------------------------------------------------
    // UN-VERIFIED but VALID
    
    static
    Sophus::SE3d
    concat(Sophus::SE3d a, Sophus::SE3d b)
    {
        
        return Sophus::SE3d(Sophus::SE3d::exp(a.log())*
                            Sophus::SE3d::exp(b.log()));
        
    }
    
};
    
} // end of namespace lsd_slam


#endif /* defined(__lsd_slam_core__DervLie__) */
