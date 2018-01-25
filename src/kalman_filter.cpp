#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) 
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() 
{

  //predict the state

  // Update the state vector 
  x_ = F_ * x_; //(+ u = 0)

  //transpose of transition matrix
  MatrixXd F_trans = F_.transpose();

  //calculate state covariance matrix
  P_ = F_ * P_ * F_trans + Q_ ;

}

void KalmanFilter::Update(const VectorXd &z, bool EKF) 
{
  VectorXd z_pred;
  VectorXd y;
  if(EKF)
  {

    // update the state
    VectorXd z_radar;
    z_radar = VectorXd(3);

    //recover state parameters
    float pos_x = x_(0);
    float pos_y = x_(1);
    float vel_x = x_(2);
    float vel_y = x_(3);

    //check division by zero
    z_radar(0) = sqrt( pow(pos_x, 2) + pow(pos_y, 2) );
    z_radar(1) = atan2(pos_y,pos_x);

    if( sqrt( pow(pos_x, 2) + pow(pos_y, 2) ) < .0001 )
      z_radar(2) = ( (pos_x * vel_x) + (pos_y * vel_y) ) / 0.0001;
    else
      z_radar(2) = ( (pos_x * vel_x) + (pos_y * vel_y) ) / sqrt( ( pow(pos_x, 2) + pow(pos_y, 2) ) );

    z_pred = z_radar;
    y = z - z_pred;
  
    bool in_pi = false;
    while (in_pi == false) 
    {
      if (y(1) > 3.14159) 
        y(1) = y(1) - 6.2831;
      else if (y(1) < -3.14159) 
        y(1) = y(1) + 6.2831;
      else 
        in_pi = true;
      
    }

  }
  else
  {
    z_pred = H_ * x_;
    y = z - z_pred;
  }
  //update the state by using Kalman Filter equations
  MatrixXd H_trans       = H_.transpose();
  MatrixXd S             = H_ * P_ * H_trans + R_;
  MatrixXd S_inverse     = S.inverse();
  MatrixXd PHt           = P_ * H_trans;
  MatrixXd K             = PHt * S_inverse;

  //new estimate P`
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * H_) * P_;
 
}
