#include "kalman_filter.h"
#include<iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_* Ht * Si;

  // new estimate
  x_ = x_ + K * y;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K*H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];

  VectorXd z_pred(3);
  z_pred(0) = sqrt(px*px+py*py);
  /*
  Before and while calculating the Jacobian matrix Hj,
  make sure your code avoids dividing by zero.
  For example, both the x and y values might be zero
  or px*px + py*py might be close to zero.
  */
  if(z_pred(0)<0.000001){
    if(px>0.0){
        px += 0.001;
    }
    else{
        px -= 0.001;
    }
    if(py>0.0){
        py += 0.001;
    }
    else{
        py -= 0.001;
    }
  }

  z_pred(1) = atan2(py,px);
  z_pred(2) = (px*vx+py*vy)/z_pred(0);

  VectorXd y = z - z_pred;
  /*
  In C++, atan2() returns values between -pi and pi.
  When calculating phi in y = z - h(x) for radar measurements,
  the resulting angle phi in the y vector should be adjusted so that it is between -pi and pi.
  The Kalman filter is expecting small angle values between the range -pi and pi.
  HINT: when working in radians, you can add 2π or subtract 2π until the angle is within the desired range.
  */
  y[1] -= (2 * M_PI) * floor((y[1] + M_PI) / (2 * M_PI));

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

  // new estimate
  x_ = x_ + K * y;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K*H_) * P_;
}
