#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}


VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
{

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size() != ground_truth.size() || estimations.size() == 0)
  {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  //accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); i++)
  {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  //calculate the mean
  rmse = rmse / estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) 
{
  MatrixXd Hj(3,4);
  
  //recover state parameters
  float pos_x = x_state(0);
  float pos_y = x_state(1);
  float vel_x = x_state(2);
  float vel_y = x_state(3);

  //pre-compute a set of terms to avoid repeated calculation
  float c1 = (pos_x * pos_x) + (pos_y * pos_y);
  float c2 = sqrt(c1);
  float c3 = (c1 * c2);

  //check division by zero
  if(fabs(c1) < 0.0001)
  {
    cout << "CalculateJacobian () - Error - Division by Zero" << endl;
    return Hj;	
  }
  //compute the Jacobian matrix 3X3
  //    1.______________________________________________2.______________________________________________3.____________4.___________
  Hj << (pos_x / c2),                                   (pos_y / c2),                                   (0),          (0),
       -(pos_y / c1),                                   (pos_x / c1),                                   (0),          (0),
        (pos_y * (vel_x * pos_y - vel_y * pos_x) / c3), (pos_x * (pos_x * vel_y - pos_y * vel_x) / c3), (pos_x / c2), (pos_y / c2);

  return Hj;
}

















