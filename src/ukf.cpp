#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

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

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

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


  is_initialized_ = false;
  n_x_ = 5; 
  n_z_ = 3;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  // Initialize weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.0);


  // initializing laser Matrixs R, H
  R_laser_ = MatrixXd(2, 2);
  H_laser_ = MatrixXd(2, 5);

  //define spreading parameter
  lambda_ = 3 - n_aug_;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);
  Zsig_pred_ = MatrixXd(n_z_, 2 * n_aug_ + 1);
  Zsig_pred_.fill(0.0);

  
  //measurement covariance matrix - laser
  R_laser_ << std_laspx_ * std_laspx_, 0,
	  0, std_laspx_ * std_laspx_;

  H_laser_ << 1, 0, 0, 0, 0,
	  0, 1, 0, 0, 0;

  x_ << 1, 1, 1, 1, 1;
  P_ << 1, 0, 0, 0, 0,
	  0, 1, 0, 0, 0,
	  0, 0, 1, 0, 0,
	  0, 0, 0, 1, 0,
	  0, 0, 0, 0, 1;
}

UKF::~UKF() {}

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
	if (!is_initialized_) {
		x_.fill(0.0);
		if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			x_[0] = meas_package.raw_measurements_[0];
			x_[1] = meas_package.raw_measurements_[1];
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) 
		{
			float rho = meas_package.raw_measurements_[0];
			float phi = meas_package.raw_measurements_[1];
			float rho_dot = meas_package.raw_measurements_[2];
			x_[0] = rho * cos(phi);
			x_[1] = rho * sin(phi);
		}
		is_initialized_ = true;
		time_us_ = meas_package.timestamp_;
		return;
	}
	double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
	time_us_ = meas_package.timestamp_;
	Prediction(dt);

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		UpdateRadar(meas_package);
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
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
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	Xsig_aug.fill(0.0);

	/* Sigma Points Logic */
	//create augmented mean vector
	VectorXd x_aug = VectorXd(7);
	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(7, 7);
	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	P_aug.fill(0.0);
	P_aug.topLeftCorner(5, 5) = P_;
	P_aug(5, 5) = std_a_ * std_a_;
	P_aug(6, 6) = std_yawdd_ * std_yawdd_;

	MatrixXd A = P_aug.llt().matrixL();

	/*  create augmented sigma points  */
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i< n_aug_; i++)
	{
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
	}


	/* Sigma Points prediction   */
	for (int i = 0; i < (2 * n_aug_ + 1); i++) {
		VectorXd input_x = Xsig_aug.col(i);
		float px = input_x[0];
		float py = input_x[1];
		float v = input_x[2];
		float psi = input_x[3];
		float psi_dot = input_x[4];
		float mu_a = input_x[5];
		float mu_psi_dot_dot = input_x[6];

		VectorXd term2 = VectorXd(5);
		VectorXd term3 = VectorXd(5);
		VectorXd result = VectorXd(5);
		if (psi_dot < 0.001) {
			term2 << v * cos(psi) * delta_t, v * sin(psi) * delta_t, 0, psi_dot * delta_t, 0;
			term3 << 0.5 * delta_t*delta_t * cos(psi) * mu_a,
				0.5 * delta_t*delta_t * sin(psi) * mu_a,
				delta_t * mu_a,
				0.5 * delta_t*delta_t * mu_psi_dot_dot,
				delta_t * mu_psi_dot_dot;
			result = Xsig_aug.col(i).head(5) + term2 + term3;
		}
		else {
			term2 << (v / psi_dot) * (sin(psi + psi_dot * delta_t) - sin(psi)),
				(v / psi_dot) * (-cos(psi + psi_dot * delta_t) + cos(psi)),
				0,
				psi_dot * delta_t,
				0;

			term3 << 0.5 * delta_t*delta_t * cos(psi) * mu_a,
				0.5 * delta_t*delta_t * sin(psi) * mu_a,
				delta_t * mu_a,
				0.5 * delta_t*delta_t * mu_psi_dot_dot,
				delta_t * mu_psi_dot_dot;
			result = Xsig_aug.col(i).head(5) + term2 + term3;
		}

		Xsig_pred_.col(i) = result;
	}

	// Predict Mean and Convariance
	x_.fill(0.0);
	for (int i = 0; i<2 * n_aug_ + 1; i++) {
		x_ = x_ + weights_[i] * Xsig_pred_.col(i);
	}
	P_.fill(0.0);
	for (int i = 0; i<2 * n_aug_ + 1; i++) {
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		while (x_diff[3]> M_PI)
			x_diff[3] -= 2.*M_PI;
		while (x_diff[3] <-M_PI)
			x_diff[3] += 2.*M_PI;
		P_ = P_ + weights_[i] * x_diff * x_diff.transpose();
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
	float px = meas_package.raw_measurements_[0];
	float py = meas_package.raw_measurements_[1];
	VectorXd z = VectorXd(2);
	z << px, py;

	VectorXd z_pred = H_laser_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_laser_.transpose();
	MatrixXd PHt = P_ * Ht;
	MatrixXd S = H_laser_ * PHt + R_laser_;
	MatrixXd Si = S.inverse();
	MatrixXd K = PHt * Si;

	//Update estimate
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	x_ = x_ + (K * y);
	P_ = (I - K * H_laser_) * P_;

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


	float rho = meas_package.raw_measurements_[0];
	float phi = meas_package.raw_measurements_[1];
	float rho_dot = meas_package.raw_measurements_[2];

	if (fabs(phi) > M_PI) {
		phi -= round(phi / (2. * M_PI)) * (2. * M_PI);
	}

	VectorXd z = VectorXd(n_z_);
	z << rho, phi, rho_dot;

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z_);
	MatrixXd S = MatrixXd(n_z_, n_z_);
	PredictRadarMeasurement(&z_pred, &S);

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z_);

	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)  //2n+1 simga points
	{

		VectorXd z_diff = Zsig_pred_.col(i) - z_pred;
		//angle normalization
		if (fabs(z_diff(1)) > M_PI) {
			z_diff(1) -= round(z_diff(1) / (2. * M_PI)) * (2. * M_PI);
		}

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		//angle normalization
		if (fabs(x_diff(3)) > M_PI) {
			x_diff(3) -= round(x_diff(3) / (2. * M_PI)) * (2. * M_PI);
		}

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();

	}

	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//residual
	VectorXd z_diff = z - z_pred;

	//angle normalization
	if (fabs(z_diff(1)) > M_PI) {
		z_diff(1) -= round(z_diff(1) / (2. * M_PI)) * (2. * M_PI);
	}

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K * S*K.transpose();

	//print result
	std::cout << "Updated state x: " << std::endl << x_ << std::endl;
	std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;


}


void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out) {
	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

												// extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// Check Divide by Zero and reset p_x & p_y
		if (fabs(p_x) < 0.0001 && fabs(p_y) < 0.0001) {
			p_x = 0.0001;
			p_y = 0.0001;
		}

		// measurement model
		Zsig_pred_(0, i) = sqrt(p_x*p_x + p_y * p_y);
		Zsig_pred_(1, i) = atan2(p_y, p_x);
		Zsig_pred_(2, i) = (p_x*v1 + p_y * v2) / sqrt(p_x*p_x + p_y * p_y);
	}

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z_);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights_(i) * Zsig_pred_.col(i);
	}

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z_, n_z_);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
												//residual
		VectorXd z_diff = Zsig_pred_.col(i) - z_pred;

		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z_, n_z_);
	R << std_radr_ * std_radr_, 0, 0,
		0, std_radphi_*std_radphi_, 0,
		0, 0, std_radrd_*std_radrd_;
	S = S + R;

	//print result
	std::cout << "z_pred: " << std::endl << z_pred << std::endl;
	std::cout << "S: " << std::endl << S << std::endl;

	//write result
	*z_out = z_pred;
	*S_out = S;
}
