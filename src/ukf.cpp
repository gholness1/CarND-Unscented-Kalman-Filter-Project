/*************************************************
 * Gary Holness
 * Unscented Kalman Filter Project
 *
 * Note:  To printout NIS, go to debug.h and set
 *        #define PRINT_NIS true  
 * 
 *        To make silent, 
 ************************************************/

#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "debug.h"

#define SMALLNUM 0.0001

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  NIS_laser_ = 0.0;
  NIS_radar_ = 0.0;
  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_pred_ = MatrixXd(5, 5);
  P_aug_ = MatrixXd(7,7);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  //std_a_ = 30;
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 30;
  std_yawdd_ = M_PI/6.0; //0.2
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15; //0.01

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15; //0.15

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  n_x_ = 5;
  n_aug_= 7;

  lambda_= 3 - n_aug_;

  NUM_SIGPTS_= 2*n_aug_ + 1;

  weights_ = VectorXd(NUM_SIGPTS_);
  weights_.setZero();

  for (int i= 0; i < NUM_SIGPTS_; i++) {
    if (i==0) {
      weights_(i) = ((double)lambda_)/(((double)lambda_) + ((double)n_aug_));
    } else {
      weights_(i) = 0.5/(lambda_ + n_aug_);
    }
  }

  Xsig_pred_= MatrixXd(n_x_, NUM_SIGPTS_);
  Xsig_gen_ = MatrixXd(n_x_, 2*n_x_ + 1);
  Xsig_aug_ = MatrixXd(n_aug_, NUM_SIGPTS_);
  Xsig_pred_aug_ = MatrixXd(n_aug_, NUM_SIGPTS_);

  x_.setZero();

  P_ << std_laspx_,  0,  0,  0,  0,
        0,  std_laspy_,  0,  0,  0,
        0,  0,  1,  0,  0,
	0,  0,  0,  1,  0,
	0,  0,  0,  0,  1;

  Xsig_pred_.setZero();
  Xsig_gen_.setZero();
  Xsig_aug_.setZero();
  Xsig_pred_aug_.setZero();

  ///* Radar measurement noise covariance matrix
  R_radar_ = MatrixXd(3,3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
	      0, 0, std_radrd_*std_radrd_;

  ///* LIDAR measurement noise covariance matrix
  R_laser_ = MatrixXd(2,2);
  R_laser_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;


  is_initialized_= false;

}

UKF::~UKF() {
  std::cout << "UKF::~UKF: called" << std::endl;
}


/****
 * AugmentedSigmaPoints
 *
 */
void UKF::AugmentedSigmaPoints() {
     if (DEBUG_AUGMENTED_SIGMA_POINTS) {
       std::cout << "UKF::AugmentedSigmaPoints" << std::endl;
     }

     MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
     P_aug.fill(0.0);

     MatrixXd P_aug_sqrt= MatrixXd(n_aug_,n_aug_);
     P_aug_sqrt.fill(0.0);

     MatrixXd Xsig_aug= MatrixXd(n_aug_, NUM_SIGPTS_);
     Xsig_aug.fill(0.0);

     VectorXd x_aug = VectorXd(n_aug_); 
     x_aug.setZero();

     x_aug << x_(0), x_(1), x_(2), x_(3), x_(4), 0, 0;

     P_aug.topLeftCorner(n_x_,n_x_) = P_;

     //populate lower right corner with Q matrix
     P_aug(5,5) = std_a_*std_a_;
     P_aug(6,6) = std_yawdd_*std_yawdd_;

     if (DEBUG_AUGMENTED_SIGMA_POINTS) {
        std::cout << "UKF:AugmentedSigmaPoints:  P_aug = " << P_aug << std::endl;
     }      

     P_aug_sqrt = P_aug.llt().matrixL();

     Xsig_aug_.col(0) = x_aug;

     MatrixXd scaledDev = MatrixXd(n_aug_,n_aug_);


     scaledDev= sqrt(lambda_ + n_aug_) * P_aug_sqrt;


     for (int i= 0; i < n_aug_; i++) {
       Xsig_aug_.col(i+1) = x_aug + scaledDev.col(i);
       Xsig_aug_.col(i+1 + n_aug_) = x_aug - scaledDev.col(i);
     }

     if (DEBUG_AUGMENTED_SIGMA_POINTS) {
       std::cout << "UKF:AugmentedSigmaPoints:  Xsig_aug= " << Xsig_aug_ << std::endl;
     }

}


/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  double rho= 0.0;
  double rhodot = 0.0;
  double phi = 0.0;

  double px = 0.0;
  double py = 0.0;
  //double vx = 0.0;
  //double vy = 0.0;
 
  double v= 0.0;
  double psi= 0.0;
  double psidot= 0.0;

  if (!is_initialized_) {
    if (DEBUG_PROCESS_MEASUREMENT) {
      std::cout << "UKF::ProcessMeasurement:  is_initialized==false, initializing..." << std::endl;
    }
  
    /***
     * If initial sensor measurement is RADAR, then convert from polar
     * coordinates to Cartesian coordinates because the state vector
     * is maintained in Cartesian coordinates.   Set velocity to zero.
     *
     */

    if ( (meas_package.sensor_type_ == MeasurementPackage::RADAR) && use_radar_) {
      rho= meas_package.raw_measurements_(0);
      phi= meas_package.raw_measurements_(1);
      rhodot= meas_package.raw_measurements_(2);

      px= rho*cos(phi);
      py= rho*sin(phi);
      v= rhodot;
      psi= phi;
      psidot= 0.0;

    } else if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && use_laser_) {
     /***
      * Initialize state from LASER
      */
      px= meas_package.raw_measurements_(0);
      py= meas_package.raw_measurements_(1);
      v= 0.0;
      psi= 0.0;
      psidot= 0.0;
   }


   /***
    * Numerical precision adustment
    * Snap px to reasonable small number (1/10,000th currenty)
    * Snap py to reasonable small number
    */
    if (fabs(px) < SMALLNUM ) {
      if (px < 0)
        px = -SMALLNUM;
      else
        px = SMALLNUM;
    }

    if (fabs(py) < SMALLNUM ) {
      if (py < 0)
        py = -SMALLNUM;
      else
        py = SMALLNUM;
    }

#if NORMALIZE_ANGLES
    while (psi > M_PI) {
      psi= psi - 2.*M_PI;
    }
    while (psi < -M_PI) {
      psi= psi + 2.*M_PI;
    }
#endif

    x_ << px, py, v, psi, psidot;

    previous_timestamp_ = meas_package.timestamp_;


    is_initialized_ = true;


    if (DEBUG_PROCESS_MEASUREMENT) {
      std::cout << "UKF::ProcessMeasurement:  x_= " << x_ << std::endl;
      std::cout << "UKF::ProcessMeasurement:  initializing...DONE" << std::endl;
    }

    return;
  }


  double dt= (meas_package.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = meas_package.timestamp_;

  if (DEBUG_PROCESS_MEASUREMENT) {
    std::cout << "UKF::ProcessMeasurement:  dt= " << dt << std::endl;
  }

  
  Prediction(dt);

  if ( (meas_package.sensor_type_ == MeasurementPackage::RADAR) && use_radar_) {
    UpdateRadar(meas_package);
  } else if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && use_laser_) {
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  if (DEBUG_PREDICTION) {
    std::cout << "UKF::Prediction" << std::endl;
  }


  /***
   * generate the sigma points
   */
  Xsig_aug_.fill(0.0);

  AugmentedSigmaPoints();


  /***
   * Predict Sigma points
   */
  Xsig_pred_.setZero();

  for (int i= 0; i < NUM_SIGPTS_; i++) {
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  if (DEBUG_PREDICTION) {
    std::cout << "UKF::Prediction:  sigma points predicted" << std::endl;
    std::cout << "UKF::Prediciton:  Xsig_pred_= " << Xsig_pred_ << std::endl;
  }




  /***
   * Now that predicted sigma points are stored, calculate the
   * predicted state mean and the predicted state transition
   * covariance matrix
   */

  if (DEBUG_PREDICTION) {
    std::cout << "UKF::Prediction:  calculating state mean..." << std::endl;
  }

  x_.setZero(); 

 /***
  * predict the state mean vector
  */
  for (int i=0; i < NUM_SIGPTS_; i++) {
    x_= x_ + weights_(i) * Xsig_pred_.col(i);
  } 

  if (DEBUG_PREDICTION) {
    std::cout << "UKF::Prediction:  calculating state mean...DONE" << std::endl;
    std::cout << "UKF::Prediciton:  x_ " << x_ << std::endl;   
  }


 /***
  * predict the state covariance matrix
  */

  if (DEBUG_PREDICTION) {
    std::cout << "UKF::Prediction:  calculating state covariance matrix..." << std::endl;
  }

  VectorXd diffVec = VectorXd(n_x_);

  P_.setZero();

  for (int i=0; i < NUM_SIGPTS_; i++) {
    diffVec.setZero();
    diffVec= Xsig_pred_.col(i) - x_;  

   /***
    * normalize the angle psi to range
    * [ -pi, pi]
    */

    if (false) { //debug statement
      std::cout << "UKF::Prediction:  diffVec(3)= " << diffVec(3) << std::endl;
    }

#if NORMALIZE_ANGLES
    while (diffVec(3) > M_PI) diffVec(3)-=2.*M_PI;
    while (diffVec(3)< -M_PI) diffVec(3)+=2.*M_PI;
#endif

    P_ = P_ + weights_(i) * (diffVec * diffVec.transpose());
  }

  if (DEBUG_PREDICTION) {
    std::cout << "UKF::Prediction:  calculating state covariance matrix...DONE" << std::endl;
    std::cout << "UKF::Prediction:  P_ = " << P_ << std::endl;
  }


}



/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  int n_z= 2;
 
  //transform sigma points into measurement space

  double px= 0;
  double py= 0;

  MatrixXd Zsig = MatrixXd(n_z, NUM_SIGPTS_);
  Zsig.setZero();

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.setZero();

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.setZero();

  MatrixXd R= MatrixXd(n_z,n_z);
  R.setZero();

  for (int i= 0; i < NUM_SIGPTS_; i++) {
    px= Xsig_pred_(0,i);
    py= Xsig_pred_(1,i);

    if (abs(px) < 0.0001) {
      if (px < 0)
        px = -0.0001;
      else
        px = 0.0001;
    }

    if (abs(py) < 0.0001) {
      if (py < 0)
        py = -0.0001;
      else
        py = 0.0001;
    }

    Zsig(0,i) = px;
    Zsig(1,i) = py;
  }

  for (int i= 0; i < NUM_SIGPTS_; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }


  VectorXd diffVec = VectorXd(n_z);

  for (int i= 0; i < NUM_SIGPTS_; i++) {
    diffVec = Zsig.col(i) - z_pred;

    S =  S + weights_(i) * diffVec * diffVec.transpose();
  }

  S= S + R_laser_;


  if (DEBUG_UPDATE_LIDAR) {
    std::cout << "UKF::UpdateLidar:  z_pred = " << z_pred << std::endl;
    std::cout << "UKF::UpdateLidar:  S = " << S << std::endl;
  }



  /***
   * Update the measurement state
   */

  VectorXd z = VectorXd(n_z);

  /***
   * Depending on sensor type, get measurment
   * LASER:  px, py
   */
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);


  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.setZero();

  VectorXd stateDiff = VectorXd(n_x_);
  VectorXd measDiff = VectorXd(n_z);


  for (int i= 0; i < NUM_SIGPTS_; i++) {
     stateDiff = Xsig_pred_.col(i) - x_;
     measDiff = Zsig.col(i) - z_pred;

     Tc = Tc + weights_(i) * (stateDiff* measDiff.transpose());
  }

  //calculate Kalman gain K;
  MatrixXd K= MatrixXd(n_x_,n_z);

  K.setZero();

  K= Tc * S.inverse();


  /***
   * update state mean and covariance matrix
   *
   */

  VectorXd z_diff = VectorXd(n_z);
  z_diff= z - z_pred;


  NIS_laser_ = z_diff.transpose()*S.inverse() * z_diff;

  if (PRINT_NIS) {
    std::cout << "UKF::UpdateLidar: NIS_laser_= " << NIS_laser_ << std::endl;
  }

  x_ = x_ + K * z_diff;

  P_ = P_ - K * S * K.transpose();
}


/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */


  int n_z= 3;

  double px= 0;
  double py= 0;
  double v1= 0;
  double v2= 0;
  double v= 0;
  double psi= 0;
  //double psidot= 0;
  double rho_tx= 0.0;
  double phi_tx= 0.0;
  double rhodot_tx= 0.0;


  MatrixXd Zsig = MatrixXd(n_z, NUM_SIGPTS_);
  Zsig.setZero();

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.setZero();
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.setZero();
 

  /***
   * Transform predicted sigma points into measurment
   * space of radar
   *
   */
  for (int i= 0; i < NUM_SIGPTS_; i++) {
    px= Xsig_pred_(0,i);
    py= Xsig_pred_(1,i);
    v=  Xsig_pred_(2,i);
    psi=Xsig_pred_(3,i);

    //psidot= Xsig_pred_(4,i);
    v1 = v* cos(psi);
    v2 = v* sin(psi);
    
    rho_tx= sqrt(px*px + py*py);
    phi_tx= atan2(py,px);

    rhodot_tx=  (px*v1 + py*v2)/rho_tx;

#if 0
    if ( ((px < 0.001) && (py < 0.001))  ||
         ((px < 0.001) && (py > 0.001)) 
       ) {
      phi_tx= 0.0;
      rhodot_tx= 0.0;
    }

#if NORMALIZE_ANGLES
    while (phi_tx < -M_PI) {
      phi_tx += 2.*M_PI;
    }
    while (phi_tx > M_PI) {
      phi_tx -= 2.*M_PI;    
    }
#endif
#endif
 
    Zsig(0,i) = rho_tx;
    Zsig(1,i) = phi_tx;
    Zsig(2,i) = rhodot_tx;
    
  }
  
  //calculate mean predicted measurement
 
  z_pred.setZero();
 
  for (int i= 0; i < NUM_SIGPTS_; i++) {
    z_pred = z_pred + weights_(i)*Zsig.col(i);
  }

  
  VectorXd diffVec = VectorXd(n_z);
  
  for (int i= 0; i < NUM_SIGPTS_; i++) {
     diffVec = Zsig.col(i) - z_pred;

#if NORMALIZE_ANGLES
     while (diffVec(1) > M_PI) diffVec(1) -= 2.*M_PI;
     while (diffVec(1) < -M_PI) diffVec(1) += 2.*M_PI;
#endif
   
     S = S +  weights_(i)* diffVec * diffVec.transpose();
  }
  
  S= S + R_radar_;

  if (DEBUG_UPDATE_RADAR) {
    std::cout << "UKF::UpdateRadar:  z_pred_ = " << z_pred_ << std::endl;
    std::cout << "UKF::UpdateRadar:  S = " << S << std::endl;
  }




  /***
   * Update the Measurement State
   */
 
 
  VectorXd z = VectorXd(n_z);

  /***
   * Depending on sensor type, get measurment
   *
   * RADAR:  rho, phi, rhodot
   */
   z(0) = meas_package.raw_measurements_(0);
   z(1) = meas_package.raw_measurements_(1);
   z(2) = meas_package.raw_measurements_(2);

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.setZero();

  VectorXd stateDiff = VectorXd(n_x_);
  VectorXd measDiff = VectorXd(n_z);


  for (int i= 0; i < NUM_SIGPTS_; i++) {
     stateDiff = Xsig_pred_.col(i) - x_;
     while (stateDiff(3) > M_PI)  stateDiff(3) -= 2.*M_PI;
     while (stateDiff(3) < -M_PI) stateDiff(3) += 2.*M_PI;

     measDiff = Zsig.col(i) - z_pred;
     while (measDiff(1) > M_PI)  measDiff(1) -= 2.*M_PI;
     while (measDiff(1) < -M_PI) measDiff(1) += 2.*M_PI;

     Tc = Tc + weights_(i) * (stateDiff* measDiff.transpose());
  }


  //calculate Kalman gain K;
  MatrixXd K= MatrixXd(n_x_,n_z);

  K.setZero();

  K= Tc * S.inverse();

  /***
   * update state mean and covariance matrix
   *
   */

  VectorXd z_diff = VectorXd(n_z);
  z_diff= z - z_pred;

#if NORMALIZE_ANGLES
  while (z_diff(1) > M_PI)  z_diff(1) -= 2.*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;
#endif

  NIS_radar_ = z_diff.transpose()*S.inverse() * z_diff;

  if (PRINT_NIS) {
    std::cout << "UKF::UpdateRadar: NIS_radar_= " << NIS_radar_ << std::endl;
  }

  x_ = x_ + K * z_diff;

  P_ = P_ - K * S * K.transpose();
}
