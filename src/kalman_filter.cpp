#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

#define PI 3.14159265358979323846

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

// Not used
/*
void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    Q_ = Q_in;
    H_ = H_in;
    R_ = R_in;
}*/

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
    // x_{k+1) = F * x_k + w, w~N(0,Q)
    // this updates x and P to one step ahead
    x_ = F_*x_;
    Eigen::MatrixXd FT = F_.transpose();
    P_ = F_*P_*F_*FT + Q_;
    
    //std::cout << x_ << std::endl;
    //std::cout << P_ << std::endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
    // need the full kalman filter algo here
    VectorXd innovation = z - H_ * x_;
    Eigen::MatrixXd HT = H_.transpose();
    MatrixXd innovationCov = H_*P_*HT + R_;
    Eigen::MatrixXd innoCovInv = innovationCov.inverse();
    MatrixXd kalmanGain = P_*HT*innoCovInv;
    
    x_ = x_ + kalmanGain*innovation;
    P_ = (Id_-kalmanGain*H_)*P_;
    
    if (x_.sum()>1e10)
        throw("we have an out of range problem");
    
    //std::cout << x_ << std::endl;
    //std::cout << P_ << std::endl;
    //std::cout << std::endl;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
    // similar to the above but we need H as a function and we require the jacobian to find the innovation
    VectorXd hx(3);
    hx(0) = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
    hx(1) = atan2(x_(1), x_(0));
    hx(2) = (x_(0)*x_(2) + x_(1)*x_(3))/hx(0);
    
    VectorXd innovation = z - hx;
    
    //normalizing phi
    while (innovation[1] > (2*PI)) {
        innovation[1] -= (2*PI);
    }
    while (innovation[1] < -(2*PI)) {
        innovation[1] += 2*PI;
    }
    
    MatrixXd J_h = Tools().CalculateJacobian(x_);
    
    MatrixXd innovationCov = R_ + J_h*P_*J_h.transpose();
    MatrixXd kalmanGain = P_*J_h.transpose()*innovationCov.inverse();
    
    x_ = x_ + kalmanGain*innovation;
    P_ = (Id_-kalmanGain*J_h)*P_;
    
    if (x_.sum()>1e10)
        throw("we have an out of range problem");
    
    //std::cout << x_ << std::endl;
    //std::cout << P_ << std::endl;
    //std::cout << std::endl;
    
}
