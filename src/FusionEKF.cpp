#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
Tools tools;
MatrixXd CalculateJacobian(const VectorXd& x_state);

//Constructor!
FusionEKF::FusionEKF() 
{
  is_initialized_ = false;

  previous_timestamp_ = 0;

  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_      = MatrixXd(3, 4);


  ekf_.Q_  = MatrixXd(4, 4);
  ekf_.H_  = MatrixXd(4, 4);


  R_laser_ << 0.0225, 0,
              0,      0.0225;

  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

/*  Hj_      << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1, 0;*/

  ekf_.H_  << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1, 0,
              0, 0, 0, 1;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) 
{
  //Initialization Step:
  if (!is_initialized_) 
  {
    //first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
    {      

      //convert radar from polar to cartesian coordinates and initialize state.
      //parameters come in from measurement_pack[rho, theta, d_rho]      
      float rho   = measurement_pack.raw_measurements_[0];
      float theta = measurement_pack.raw_measurements_[	1];
      float rho_dot = measurement_pack.raw_measurements_[2];

      ekf_.x_(0) = rho * cos(theta);
      ekf_.x_(1) = rho * sin(theta);
      ekf_.x_(2) = rho_dot * cos(theta);
      ekf_.x_(3) = rho_dot * sin(theta);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) 
    {

      //laser readings come in X and Y only     
      ekf_.x_(0) = measurement_pack.raw_measurements_(0);
      ekf_.x_(1) = measurement_pack.raw_measurements_(1);
    }

  ekf_.P_  = MatrixXd(4, 4);
  ekf_.P_  << 1, 0, 0,    0,
              0, 1, 0,    0,
              0, 0, 1000, 0,
              0, 0, 0,    1000;

    //init timestamp
    previous_timestamp_ = measurement_pack.timestamp_;

    //end of initialization and first loop.
    is_initialized_ = true;

    return;
  }

  //Prediction Step:
  
  //calculate delta t and update the state transition matrix
  float delta_t = (measurement_pack.timestamp_ - previous_timestamp_ ) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  ekf_.F_  = MatrixXd(4, 4);
  ekf_.F_  << 1, 0, delta_t, 0,
              0, 1, 0,  delta_t,
              0, 0, 1,  0,
              0, 0, 0,  1;
            
  //set the noise
  noise_ax = 9.0;
  noise_ay = 9.0;

  //process covariance matrix Q
            //1.____________________________________2.____________________________________3.______________________________________4.______________________________
  ekf_.Q_ <<  ( pow(delta_t, 4) / (4 * noise_ax) ), (0),                                  ( pow(delta_t, 3) / (2 * noise_ax) ),   (0),
              (0),                                  ( pow(delta_t, 4) / (4 * noise_ay) ), (0),                                    ( pow(delta_t, 3) / (2 * noise_ay) ),
              ( pow(delta_t, 3) / (2 * noise_ax) ), (0),                                  ( pow(delta_t,2) * (noise_ax) ),        (0),
              (0),                                  ( pow(delta_t, 3) / (2 * noise_ay) ), (0),                                    ( pow(delta_t, 2) * (noise_ay) );

  ekf_.Predict();

  //Update Step:
  
  //check if incoming data is from radar or lidar and process accordingly.
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
  {
    //radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.Update(measurement_pack.raw_measurements_,true);
  } 
  else 
  {
    //laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_,false);
  }

  //print output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
