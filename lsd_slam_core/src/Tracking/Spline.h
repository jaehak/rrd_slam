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
// HISTORY
//------------------------------------------------------------------------------
// 2014-11-10 Writing started by Jae-Hak Kim
//            porting from MATLAB code, +JHK/@Spline/Spline.m
//------------------------------------------------------------------------------
#ifndef SPLINE_H
#define SPLINE_H
#include <vector>
#include "sophus/se3.hpp"
//#include "util/SophusUtil.h"
#include "ceres/ceres.h"
#include "Tracking/JetOps.h"
#include <stdio.h>
#define EIGEN_USE_NEW_STDVECTOR
#include <Eigen/StdVector>

namespace lsd_slam
{

//------------------------------------------------------------------------------
// B-Spline class for camera trajectory
template <typename T>
class Spline
{

    typedef typename Sophus::SE3Group<T> SE3T;

protected:

    int n = 0; // Number of control poses
    int k = 0; // k - 1 degree of the B-spline curve

    // n control poses for 4x4 transformation
    //EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix<T,4,4>);
    std::vector<Eigen::Matrix<T,4,4>,
    Eigen::aligned_allocator<Eigen::Matrix<T,4,4>> > controlPoseTs;

    // A cummulative basis in matrix representation
    Eigen::Matrix<T,4,4> C4; // Basis for k = 4
    Eigen::Matrix<T,3,3> C3; // Basis for k = 3
    Eigen::Matrix<T,2,2> C2; // Basis for k = 2

public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    //--------------------------------------------------------------------------
    // Constructor with n control pose Ts
    //     Xis : [nx1 SE3] n 4x4 matrices in SE(3)
    Spline(std::vector<Sophus::SE3Group<T> > Xi, int k)
    {

        //EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix<T,4,4>);
        
        std::vector<Eigen::Matrix<T,4,4>,
        Eigen::aligned_allocator<Eigen::Matrix<T,4,4>> > input_controlPoses;
        for (int i=0; i<(int)Xi.size(); i++)
        {

            Eigen::Matrix<T,4,4> Xi_mat = Xi[i].matrix();
            input_controlPoses.push_back(Xi_mat);

        }
        this->initWithTandK(input_controlPoses, k);

    }
    
    // Copy constructor
    Spline(const Spline &other) = delete; // Disabled
//    {
//        this->n = other->n;
//        this->k = other->k;
//        for (int i=0;i<other->controlPoseTs.size();i++)
//        {
//            
//            this->controlPoseTs->push_back(other->controlPoseTs[i]);
//            
//        }
//        this->C4 = other->C4;
//        this->C3 = other->C3;
//        this->C2 = other->C2;
//
//    }
    
    // Copy assignment operator
    Spline& operator= (const Spline& other) = delete; // Disabled
//    {
//        
//        Spline tmp(other); // Re-use copy constructor
//        *this = std::move(tmp); // re-use move-assignment
//        return *this;
//        
//    }

    // Move assignment operator
    Spline& operator= (Spline&& other) noexcept = delete; // Disabled
//    {
//        // Delete data
//        this->controlPoseTs.clear();
//        for(int i=0;i<other.controlPoseTs.size();i++)
//        {
//            
//            this->controlPoseTs->push_back(other.controlPoseTs[i]);
//            
//        }
//        other->controlPosesTs.clear();
//        return *this;
//    }

    // Get number of control poses
    int numControlPoses() { return n; }
    
    // Get order
    int order() { return k; }

    // Transform spline control poses
    void transform(Sophus::SE3Group<T> Xi)
    {

        for (int i=0; i<this->controlposeTs.size(); i++)
        {

            this->controlposeTs[i] = Xi * this->controlposeTs[i];

        }

    }

    // Convenient constructor
    Spline(std::vector<Sophus::SE3Group<T> > Xi) : Spline(Xi, 3)
    {

    }
    
    // Destructor
    ~Spline()
    {
        
        controlPoseTs.clear();
        
    }
    
    // control poses
    std::vector<Eigen::Matrix<T,4,4>, Eigen::aligned_allocator<Eigen::Matrix<T,4,4>> > *getControlPoses()
    {
        return &controlPoseTs;
    }


private:

    // Similar to designated constructor
    //void init(std::vector<Eigen::Matrix4d> controlposeTs);

    void initWithTandK(
            std::vector<Eigen::Matrix<T,4,4>,
                       Eigen::aligned_allocator<Eigen::Matrix<T,4,4>>> input_controlPoseTs,
            int k)
    {

        // Set the number of control poses
        this->n = input_controlPoseTs.size();

        // Check n
        if (this->n < k)
        {

            printf("Error: Not enough number of control poses T\n");
            exit(1);

        }

        // Set the degree k - 1 of the B-Spline curve
        this->k = k;

        // Set control poses
        for (int i=0;i<(int)input_controlPoseTs.size();i++)
        {
            
            this->controlPoseTs.push_back(input_controlPoseTs[i]);
            
        }

        if (k == 4)
        {
            // Set the basis
            Eigen::Matrix<T,4,4> C;
            C << T(6.0), T(0.0),  T(0.0),  T(0.0),
                 T(5.0), T(3.0), T(-3.0),  T(1.0),
                 T(1.0), T(3.0),  T(3.0), T(-2.0),
                 T(0.0), T(0.0),  T(0.0),  T(1.0);
            this->C4 = T(1.0/6.0)*C;
        }

        if (k == 3)
        {
            // Set the basis
            Eigen::Matrix<T,3,3> C;
            C << T(2.0), T(1.0),  T(0.0),
                 T(0.0), T(2.0), T(0.0),
                 T(0.0), T(-1.0),  T(1.0);
            this->C3 = T(1.0/2.0)*C;
        }

        if (k == 2)
        {
            // Set the basis
            this->C2 << T(1.0), T(0.0),
                        T(0.0), T(1.0);
        }

    }

public:
    
    //-----------------------------------------------------------------
    // Update control poses with a new control poses
    void updateControlPoses(std::vector<Sophus::SE3Group<T> > Xi)
    {
        
        // Clear existing control poses
        if (this->n > 0)
        {
            
            this->controlPoseTs.clear();
            
        }
        
        // Assign new control poses to a vector of matrices
        //std::vector<Eigen::Matrix<T,4,4> > controlPoses;
        for (int i=0; i<(int)Xi.size(); i++)
        {
            
            Eigen::Matrix<T,4,4> Xi_mat = Xi[i].matrix();
            this->controlPoseTs.push_back(Xi_mat);
            
        }
        
        // Check number of control poses
        if ((int)this->controlPoseTs.size() < k)
        {
            
            printf("Error: Not enough number of control poses T\n");
            exit(1);
            
        }
        
        // Set the number of control poses
        this->n = this->controlPoseTs.size();
        
        // Upddate control poses
        //this->controlPoseTs = controlPoses;
        
    }
    
    //-----------------------------------------------------------------
    Eigen::Matrix<T,4,4> getPoseOnSpline(T u, int i)
    {
        Eigen::Matrix<T,4,4> M;
        M = getPoseOnSpline_CtrlToRef(u, i); // Default getPoseOneSpline
        //M = getPoseOnSpline_RefToCtrl(u, i);
        return M;
        
    }
    
    //-----------------------------------------------------------------
    // Compute pose on the spline trajectory
    //     u : [double] A parameter of the curve such that 0.0 <= u
    //                  <= 1.0 within the i-th interval. The interval
    //                  is a segment between two control poses
    //                  at i and i + 1
    //     i : [int] index of the interval where u(t) is defined
    //               (0-based index starting from 0, 1, ...)
    // returns
    //
    //     [4x4 double] the pose on the spline
    //
    // CtrlToRef version
    Eigen::Matrix<T,4,4> getPoseOnSpline_CtrlToRef(T u, int i)
    {

        typedef Sophus::SE3Group<T> SE3T;
        // Refer to Eq (4) in paper, Spline Fusion: A continuous-time
        // representation for visual-inertial fusion with application
        // to rolling shutter cameras

        // Check u
        if ((u < T(0.0)) || (u > T(1.0)))
        {

            printf("Error: u should be real numbers in [0 .. 1].\n");
            exit(1);

        }

        // Check index if it is out of range or not integer
        if ((i < 0) || (i > this->n - (this->k - 1) + 1) || ((i % 1) != 0))
        {

            printf("i=%d, this->n=%d, "
                   "this->k = %d, (i mod 1)=%d\n", i, this->n, this->k, (i%1));
            printf("Error: the given interval i (%d) is out of range. "
                   "It should be an integer value between 0 and %d\n",
                   i,
                   this->n - (this->k - 1) + 1);
            exit(1);

        }
        
        // Check if spline is initialized
        if (n == 0 || k == 0)
        {
            
            printf("Error: Spline object is not initialized.\n");
            exit(1);
            
        }

        // (i - 1)-th control pose
        Eigen::Matrix<T,4,4> T_w_im1 = this->controlPoseTs[i];
#if 0//DEBUG
        printf("T_w_im1 = Ts[%d]\n", i);
        printf("T_w_im1=\n");
        for (int p=0; p<4; p++)
        {
            for (int q=0; q<4; q++)
            {
                printf("%.2f ", ceres::JetOps<T>::GetScalar(T_w_im1(p,q)));
            }
            printf("\n");
        }
#endif
        // tilde(B(u))_j
        Eigen::Matrix<T, 4, 1> Bu;
        if (this->k == 4)
        {
            Bu = this->C4*Eigen::Matrix<T,4,1>(T(1.0), T(u), T(u*u), T(u*u*u));
        }
        if (this->k == 3)
        {
            Eigen::Matrix<T, 3, 1> tmp =
                    this->C3*Eigen::Matrix<T,3,1>(T(1.0), T(u), T(u*u));
            Bu(0,0) = tmp[0];
            Bu(1,0) = tmp[1];
            Bu(2,0) = tmp[2];
            Bu(3,0) = T(0);
        }
        if (this->k == 2)
        {
            Eigen::Matrix<T, 2, 1> tmp =
                    this->C2*Eigen::Matrix<T,2,1>(T(1.0), T(u));
            Bu(0,0) = tmp[0];
            Bu(1,0) = tmp[1];
            Bu(2,0) = T(0);
            Bu(3,0) = T(0);
        }
#if 0//DEBUG
        printf("\nBu = \n");
        for (int p=0; p<4; p++)
            printf("%f ",ceres::JetOps<T>::GetScalar(Bu(p)));
        printf("\n");
#endif
        // Multiplying
        Eigen::Matrix<T,4,4> mulTerm;
        for (int k=0; k<4; k++)
        {

            for (int m=0; m<4; m++)
            {

                if (m == k)
                {

                    mulTerm(k,m) = T(1.0);

                }
                else
                {

                    mulTerm(k,m) = T(0.0);

                }

            }

        }

//        for (int k=0; k<4; k++)
//        {

//            mulTerm(k,k) = T(1.0);

//        }
        //for (int j=1; j<=3; j++)
        for (int j=1; j<=this->k - 1; j++)
        {

            // InvT(i + j - 2 + 1)
            Eigen::Matrix<T,4,4> invT_ipjm2p1 =
                    this->controlPoseTs[i + j - 2 + 1].inverse();

            // T(i + j -1 + 1)
            Eigen::Matrix<T,4,4> T_ipjm1p1 =
                    this->controlPoseTs[i + j - 1 + 1];

            // Omega_{i + j}
            SE3T LieOmega_ipj = SE3T(invT_ipjm2p1*T_ipjm1p1);
            Eigen::Matrix<T,6,1> Omega_ipj = LieOmega_ipj.log();
#if 0//DEBUG
            printf("Omega_ipj = Ts[%d].inverse() * Ts[%d]\n",
                   i + j - 2 + 1, i + j - 1 + 1);
            printf("\nOmega_ipj = \n");
            for (int p=0; p<6; p++)
                printf("%.2f ", ceres::JetOps<T>::GetScalar(Omega_ipj(p)));
            printf("\n");
#endif

            // Multiplication over j = 1 to 3
            Eigen::Matrix<T,4,4> term = SE3T::exp(Bu(j)*Omega_ipj).matrix();

#if 0//DEBUG
            printf("\nterm = \n");
            for (int p=0; p<4; p++)
            {
                for (int q=0; q<4; q++)
                {
                    printf("%.2f ", ceres::JetOps<T>::GetScalar(term(p,q)));
                }
                printf("\n");
            }
#endif

            Eigen::Matrix<T,4,4> tmp = mulTerm * term;
            mulTerm = tmp;

        }

#if 0//DEBUG
        printf("\nmulTerm = \n");
        for (int p=0; p<4; p++)
        {
            for (int q=0; q<4; q++)
            {
                printf("%.2f ", ceres::JetOps<T>::GetScalar(mulTerm(p,q)));
            }
            printf("\n");
        }
#endif

        // The pose on the spline trajectory at the parameter u defined on
        // the i-ths interval
        Eigen::Matrix<T,4,4> pose = T_w_im1 * mulTerm;

#if 0//DEBUG

        Eigen::Matrix<T,4,4> OneMat;
        OneMat.setOnes();
        Eigen::Matrix<T,4,4> tmp = OneMat - pose;
        double sqsum = 0.0;
        for (int p=0; p<4; p++)
        {
            for (int q=0; q<4; q++)
            {
                double val = ceres::JetOps<T>::GetScalar(tmp(p,q));
                sqsum = sqsum + val;
            }
        }
        if (sqsum > 0.0)
        {
            printf("\nComputed pose from spline is BIG = \n");
            for (int p=0; p<4; p++)
            {
                for (int q=0; q<4; q++)
                {
                    printf("%.2f ", ceres::JetOps<T>::GetScalar(pose(p,q)));
                }
                printf("\n");
            }
        }

#endif

        return pose;

    }
    
    //-----------------------------------------------------------------
    // Compute pose on the spline trajectory
    // This Inv version assumes each T is represented by a pose of
    // the reference frame with respect to the control point
    //
    //     u : [double] A parameter of the curve such that 0.0 <= u
    //                  <= 1.0 within the i-th interval. The interval
    //                  is a segment between two control poses
    //                  at i and i + 1
    //     i : [int] index of the interval where u(t) is defined
    //               (0-based index starting from 0, 1, ...)
    // returns
    //
    //     [4x4 double] the pose on the spline
    //
    // RefToCtrl version
    Eigen::Matrix<T,4,4> getPoseOnSpline_RefToCtrl(T u, int i)
    {
        
        typedef Sophus::SE3Group<T> SE3T;
        // Refer to Eq (4) in paper, Spline Fusion: A continuous-time
        // representation for visual-inertial fusion with application
        // to rolling shutter cameras
        
        // Check u
        if ((u < T(0.0)) || (u > T(1.0)))
        {
            
            printf("Error: u should be real numbers in [0 .. 1].\n");
            exit(1);
            
        }
        
        // Check index if it is out of range or not integer
        if ((i < 0) || (i > this->n - (this->k - 1) + 1) || ((i % 1) != 0))
        {
            
            printf("Error: the given interval i (%d) is out of range. "
                   "It should be an integer value between 0 and %d\n",
                   i,
                   this->n - (this->k - 1) + 1);
            exit(1);
            
        }
        
        // (i - 1)-th control pose
        Eigen::Matrix<T,4,4> T_w_im1 = this->controlposeTs[i];
#if 0//DEBUG
        printf("T_w_im1 = Ts[%d]\n", i);
        printf("T_w_im1=\n");
        for (int p=0; p<4; p++)
        {
            for (int q=0; q<4; q++)
            {
                printf("%.2f ", ceres::JetOps<T>::GetScalar(T_w_im1(p,q)));
            }
            printf("\n");
        }
#endif
        // tilde(B(u))_j
        Eigen::Matrix<T, 4, 1> Bu;
        if (this->k == 4)
        {
            Bu = this->C4*Eigen::Matrix<T,4,1>(T(1.0), T(u), T(u*u), T(u*u*u));
        }
        if (this->k == 3)
        {
            Eigen::Matrix<T, 3, 1> tmp =
            this->C3*Eigen::Matrix<T,3,1>(T(1.0), T(u), T(u*u));
            Bu(0,0) = tmp[0];
            Bu(1,0) = tmp[1];
            Bu(2,0) = tmp[2];
            Bu(3,0) = T(0);
        }
        if (this->k == 2)
        {
            Eigen::Matrix<T, 2, 1> tmp =
            this->C2*Eigen::Matrix<T,2,1>(T(1.0), T(u));
            Bu(0,0) = tmp[0];
            Bu(1,0) = tmp[1];
            Bu(2,0) = T(0);
            Bu(3,0) = T(0);
        }
#if 0//DEBUG
        printf("\nBu = \n");
        for (int p=0; p<4; p++)
            printf("%f ",ceres::JetOps<T>::GetScalar(Bu(p)));
        printf("\n");
#endif
        // Multiplying
        Eigen::Matrix<T,4,4> mulTerm;
        for (int k=0; k<4; k++)
        {
            
            for (int m=0; m<4; m++)
            {
                
                if (m == k)
                {
                    
                    mulTerm(k,m) = T(1.0);
                    
                }
                else
                {
                    
                    mulTerm(k,m) = T(0.0);
                    
                }
                
            }
            
        }
        
        //        for (int k=0; k<4; k++)
        //        {
        
        //            mulTerm(k,k) = T(1.0);
        
        //        }
        //for (int j=1; j<=3; j++)
        for (int j=1; j<=this->k - 1; j++)
        {
            
            // T(i + j - 1 + 1)^{-1}
            Eigen::Matrix<T,4,4> invT_ipjm1p1 =
            this->controlposeTs[i + j - 1 + 1].inverse();
            
            // T(i + j - 2 + 1)
            Eigen::Matrix<T,4,4> T_ipjm2p1 =
            this->controlposeTs[i + j - 2 + 1];
            
            // Omega_{i + j}
            SE3T LieOmega_ipj = SE3T(T_ipjm2p1*invT_ipjm1p1);
            Eigen::Matrix<T,6,1> Omega_ipj = LieOmega_ipj.log();
#if 0//DEBUG
            printf("Omega_ipj = Ts[%d].inverse() * Ts[%d]\n",
                   i + j - 2 + 1, i + j - 1 + 1);
            printf("\nOmega_ipj = \n");
            for (int p=0; p<6; p++)
                printf("%.2f ", ceres::JetOps<T>::GetScalar(Omega_ipj(p)));
            printf("\n");
#endif
            
            // Multiplication over j = 1 to 3
            Eigen::Matrix<T,4,4> term = SE3T::exp(Bu(j)*Omega_ipj).matrix();
            
#if 0//DEBUG
            printf("\nterm = \n");
            for (int p=0; p<4; p++)
            {
                for (int q=0; q<4; q++)
                {
                    printf("%.2f ", ceres::JetOps<T>::GetScalar(term(p,q)));
                }
                printf("\n");
            }
#endif
            
            Eigen::Matrix<T,4,4> tmp = term.inverse() * mulTerm;
            mulTerm = tmp;
            
        }
        
#if 0//DEBUG
        printf("\nmulTerm = \n");
        for (int p=0; p<4; p++)
        {
            for (int q=0; q<4; q++)
            {
                printf("%.2f ", ceres::JetOps<T>::GetScalar(mulTerm(p,q)));
            }
            printf("\n");
        }
#endif
        
        // The pose on the spline trajectory at the parameter u defined on
        // the i-ths interval
        Eigen::Matrix<T,4,4> pose = mulTerm * T_w_im1;
        
#if 0//DEBUG
        
        Eigen::Matrix<T,4,4> OneMat;
        OneMat.setOnes();
        Eigen::Matrix<T,4,4> tmp = OneMat - pose;
        double sqsum = 0.0;
        for (int p=0; p<4; p++)
        {
            for (int q=0; q<4; q++)
            {
                double val = ceres::JetOps<T>::GetScalar(tmp(p,q));
                sqsum = sqsum + val;
            }
        }
        if (sqsum > 0.0)
        {
            printf("\nComputed pose from spline is BIG = \n");
            for (int p=0; p<4; p++)
            {
                for (int q=0; q<4; q++)
                {
                    printf("%.2f ", ceres::JetOps<T>::GetScalar(pose(p,q)));
                }
                printf("\n");
            }
        }
        
#endif
        
        return pose;
        
    }
    
}; // end of class

} // end of namespace

#endif // SPLINE_H
