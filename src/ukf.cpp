#include "ukf.h"
#include "Eigen/Dense"
#include "Eigen"
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

	float ra = 7;

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 0.05*5;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.01*8;

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
	use_laser_ = true;
	use_radar_ = true;

	n_x_ = 5;
	n_aug_ = 7;
	lambda_ = 3 - n_aug_;
	Xsig_pred_ = MatrixXd(n_x_, (2*n_aug_+1));
	Xsig_pred_.fill(0.0);

	weights_ = VectorXd(2*n_aug_+1);
	weights_(0) = lambda_/(lambda_ + n_aug_);
	for(int i=1; i<(2*n_aug_+1);i++)
	{
		weights_(i) = 0.5/(lambda_ + n_aug_);
	}
	P_ = MatrixXd::Identity(5, 5)*0.004;
	x_<<0.0,0.0,0.0,0.0,0.0;
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
	if (!is_initialized_)
	{
	    cout<<"UKF:\n";
		time_us_ = meas_package.timestamp_;
				if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			/**
      Convert radar from polar to cartesian coordinates and initialize state.
			 */
			float rho = meas_package.raw_measurements_(0);
			float theta = meas_package.raw_measurements_(1);
			float r_dot = meas_package.raw_measurements_(2);
			x_ << rho*cos(theta), rho*sin(theta),
					r_dot*cos(theta), r_dot*sin(theta),0;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			/**
      Initialize state.
			 */
			x_ <<meas_package.raw_measurements_(0),
					meas_package.raw_measurements_(1), 0, 0,0;
		}
		is_initialized_ = true;
		return;
	}

	if((meas_package.sensor_type_== MeasurementPackage::LASER)&use_laser_)
	{
		//cout<<"laser"<<endl;
		double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
		time_us_ = meas_package.timestamp_;
		Prediction(dt);
		//cout<<"laser update"<<endl;
		UpdateLidar(meas_package);
	}
	else if((meas_package.sensor_type_== MeasurementPackage::RADAR)&use_radar_)
	{
		//cout<<"radar"<<endl;
		double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
		time_us_ = meas_package.timestamp_;
		Prediction(dt);

		UpdateRadar(meas_package);
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

	//---------------parepare sigma points-----------------
	// augmented state
	VectorXd x_aug_(n_aug_);
	x_aug_.fill(0.0);
	x_aug_.head(n_x_) = x_;

	// augmented covariance matrix
	MatrixXd P_aug(n_aug_, n_aug_);
	P_aug.fill(0.0);
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug(n_x_, n_x_) = std_a_*std_a_;
	P_aug(n_x_+1, n_x_+1) = std_yawdd_*std_yawdd_;

	//calculate square root of P
	MatrixXd A = P_aug.llt().matrixL();

	double a=std::sqrt(lambda_+n_aug_);

	// generate augmented sigma points
	MatrixXd Xsig_aug_(n_aug_, (2*n_aug_+1));
	Xsig_aug_.col(0) = x_aug_;
	for(int i=0;i<n_aug_;i++)
	{
		Xsig_aug_.col(1+i) = x_aug_+a*A.col(i);
		Xsig_aug_.col(1+n_aug_+i) = x_aug_-a*A.col(i);
	}

	//--------------sigma points prediction-----------------
	VectorXd x0(n_aug_);
	VectorXd x1(n_x_);
	VectorXd x2(n_x_);
	x1.fill(0.0);
	x2.fill(0.0);
	float d2 = delta_t*delta_t;
	using std::sin;
	using std::cos;

	for(int i=0;i<(2*n_aug_+1);i++)
	{
		x0 = Xsig_aug_.col(i);
		if(x0(4)!=0)
		{
			x1(0) = x0(2)/x0(4)*( sin(x0(3)+x0(4)*delta_t) - sin(x0(3)) );
			x1(1) = x0(2)/x0(4)*(-cos(x0(3)+x0(4)*delta_t) + cos(x0(3)) );
			x1(3) = x0(4)*delta_t;
		}
		else
		{
			x1(0) = x0(2)*cos(x0(3))*delta_t;
			x1(1) = x0(2)*sin(x0(3))*delta_t;
			x1(3) = x0(4)*delta_t;
		}
		x2(0) = 0.5*d2*cos(x0(3))*x0(5);
		x2(1) = 0.5*d2*sin(x0(3))*x0(5);
		x2(2) = delta_t*x0(5);
		x2(3) = 0.5*d2*x0(6);
		x2(4) = delta_t*x0(6);
		Xsig_pred_.col(i) = x0.head(n_x_)+x1+x2;
	}

	// predict mean and covariance of state
	x_.fill(0.0);
	P_.fill(0.0);
	for(int i=0;i<(2*n_aug_+1);i++)
	{
		x_ += weights_(i)*Xsig_pred_.col(i);
	}

	for(int i=0;i<(2*n_aug_+1);i++)
	{
		VectorXd tmp1 = Xsig_pred_.col(i) - x_;
		P_ += weights_(i)*tmp1*tmp1.transpose();
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

	VectorXd z(2);
	z = meas_package.raw_measurements_;

	MatrixXd Zsig(2, (2*n_aug_+1));
	Zsig.fill(0.0);
	Zsig.block(0,0,2, (2*n_aug_+1)) = Xsig_pred_.block(0,0,2, (2*n_aug_+1));

	VectorXd zpred(2);
	zpred.fill(0.0);
	for(int i=0;i<(2*n_aug_+1);i++)
	{
		zpred += weights_(i)*Zsig.col(i);
	}

	MatrixXd S(2,2);
	S.fill(0.0);
	for(int i=0;i<(2*n_aug_+1);i++)
	{
		VectorXd dif = Zsig.col(i) - zpred;
		S += weights_(i)*dif*dif.transpose();
	}
	S(0, 0) += std_laspx_*std_laspx_;
	S(1, 1) += std_laspy_*std_laspy_;

	MatrixXd Tc(n_x_, 2);
	Tc.fill(0.0);
	for(int i=0;i<(2*n_aug_+1);i++)
	{
		VectorXd dif = Zsig.col(i) - zpred;
		Tc += weights_(i)*(Xsig_pred_.col(i) - x_)*dif.transpose();
	}

	MatrixXd K = Tc*S.inverse();
	x_ += K*(z - zpred);
	P_ += -K*S*K.transpose();

	VectorXd zdif = z - zpred;
	double nis = zdif.transpose()*S.inverse()*zdif;
	cout<<"laser update nis = "<<nis<<endl;
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

	// measurement data
	VectorXd z(3);
	z = meas_package.raw_measurements_;

	// create matrix for sigma points in measurement space
	MatrixXd Zsig(3, (2*n_aug_+1));
	Zsig.fill(0.0);

	// measurement covariance matrix
	MatrixXd S(3,3);
	S.fill(0.0);

	// mean predicted measurement
	VectorXd zpred(3);
	zpred.fill(0.0);

	using std::sqrt;
	using std::atan2;
	using std::sin;
	using std::cos;

	// transfer sigma points into measurement space
	for(int i=0;i<(2*n_aug_+1);i++)
	{
		VectorXd x0 = Xsig_pred_.col(i);
		double rho = sqrt(x0(0)*x0(0) + x0(1)*x0(1));
		double phi = atan2(x0(1), x0(0));
		if (rho>1e-9)
		{
		double rd = (x0(0)*cos(x0(3))+x0(1)*sin(x0(3)))*x0(2)/rho;
		Zsig.col(i) << rho, phi, rd;
		}
		else
		{
			cout<<"r = 0===================================="<<endl;
			return;
		}
	}

	// calculate predicted mean measurement
	for(int i=0;i<(2*n_aug_+1);i++)
	{
		zpred += weights_(i)*Zsig.col(i);
	}

	// calculate innovation covariance matrix S
	for(int i=0;i<(2*n_aug_+1);i++)
	{
		VectorXd dif = Zsig.col(i) - zpred;
		S += weights_(i)*dif*dif.transpose();
	}

	S(0, 0) += std_radr_*std_radr_;
	S(1, 1) += std_radphi_*std_radphi_;
	S(2, 2) = std_radrd_*std_radrd_;

	// matrix for cross correlation
	MatrixXd Tc(n_x_, 3);
	Tc.fill(0.0);
	// calculate cross correlation
	for(int i=0;i<(2*n_aug_+1);i++)
	{
		VectorXd dif = Zsig.col(i) - zpred;
		Tc += weights_(i)*(Xsig_pred_.col(i) - x_)*dif.transpose();
	}

	// Kalman gain
	MatrixXd K = Tc*S.inverse();

	VectorXd r = (z - zpred);
	while (r(1)> PI)
	{
		r(1) -= 2*PI;
	}
	while (r(1)< -PI)
	{
		r(1)+= 2*PI;
	}
	// update state mean and cov matrix
	x_ += K*r;
	P_ += -K*S*K.transpose();

	VectorXd zdif = z - zpred;
	double nis = zdif.transpose()*S.inverse()*zdif;
	cout<<"radar update nis = "<<nis<<endl;
}
