#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // Check the validity of the following inputs:
  // * the estimation vector size should not be zero
  if(estimations.size() == 0){
    cout << "Error! The estimation vector size should not be zero.\n";
    return rmse;
  }
  // * the estimation vector size should equal ground truth vector size
  if(estimations.size() != ground_truth.size()){
    cout << "Error! The estimation vector size should equal ground truth vector size.\n";
    return rmse;
  }

  // accumulate squared residuals
  for(int i=0;i<estimations.size();++i){
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse /= estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
  Hj << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;

  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // define the same calculation repeatedly
  float den = px*px + py*py;

  // check division by zero
  if(den == 0){
    cout<<"CalculateJacobian() - Error - Division by Zero\n";
    return Hj;
  }


  // define the same calculation repeatedly - continue
  float den_sqrt = sqrt(den);
  float den_3_2 = den*den_sqrt;

  // compute the Jacobian matrix
  Hj << px/den_sqrt,                py/den_sqrt,                0,             0,
        -py/den,                    px/den,                     0,             0,
        py*(vx*py - vy*px)/den_3_2, px*(vy*px - vx*py)/den_3_2, px / den_sqrt, py / den_sqrt;

  return Hj;
}
