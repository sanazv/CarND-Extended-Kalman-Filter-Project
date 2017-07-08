#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
    

  //measurement covariance matrix - laser
  H_laser_ << 1.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0;
  R_laser_ << 0.0225, 0.0,
        0.0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0., 0.,
        0., 0.0009, 0.,
        0., 0., 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
    Q_ = Eigen::MatrixXd(4,4);

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    //ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
        double rho = measurement_pack.raw_measurements_(0);
        double phi = measurement_pack.raw_measurements_(1);
        double dphidt = measurement_pack.raw_measurements_(2);
        ekf_.x_(0) = rho * sin(phi);
        ekf_.x_(1) = rho * cos(phi);
        ekf_.x_(2) = dphidt * sin(phi);
        ekf_.x_(3) = dphidt * cos(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
        ekf_.x_(0) = measurement_pack.raw_measurements_(0);
        ekf_.x_(1) = measurement_pack.raw_measurements_(1);
        // with the laser measurment we don't have any info in velocity, we start with the
        // assumption that it is not moving
        ekf_.x_(2) = 0.0;
        ekf_.x_(3) = 0.0;
    }
      // set the previous timestamp to allow for time updating
      previous_timestamp_ = measurement_pack.timestamp_;
      
      ekf_.P_ = Eigen::MatrixXd::Zero(4, 4);
      ekf_.Id_ = Eigen::MatrixXd::Identity(4, 4);
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
    
    double noise_ax = 9.0; // \sig_{ax}^2
    double noise_ay = 9.0; // \sig_{ay}^2
    long long currentTime = measurement_pack.timestamp_;
    double del_t = (currentTime-previous_timestamp_)/1000000.0;
    double del4 = pow(del_t, 4.0);
    double del3 = pow(del_t, 3.0);
    double del2 = pow(del_t, 2.0);
    Q_ = Eigen::MatrixXd(4,4);
    Q_ << 0.25*del4*noise_ax , 0.0, 0.5*del3*noise_ax, 0.0,
            0.0, 0.25*del4*noise_ay , 0.0, 0.5*del3*noise_ay,
            0.5*del3*noise_ax , 0.0, del2*noise_ax, 0.0,
            0.0, 0.5*del3*noise_ay , 0.0, del2*noise_ay;
    ekf_.Q_ = Q_; // set the process covariance
    
    Eigen::MatrixXd F = Eigen::MatrixXd(4,4);
    F << 1.0, 0.0, del_t, 0.0,
        0.0, 1.0, 0.0, del_t,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0;
    ekf_.F_ = F;
    ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
      ekf_.R_ = R_laser_;
      ekf_.H_ = H_laser_;
      ekf_.Update(measurement_pack.raw_measurements_);
  }

    // update the timestamp
    previous_timestamp_ = currentTime;
    
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
    
    double x = 0.;
}
