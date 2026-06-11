#pragma once
#include "linalg.h"

class KalmanFilter {
private:
	Vector9 m_state;						// state vector (position, velocity, etc.)
	Matrix9x9 m_covariance;					// post-update covariance P_post
	Matrix9x9 m_covariance_pred;			// predicted covariance P⁻
	Matrix9x9 m_processNoise;				// process noise covariance
	Matrix9x9 m_measurementNoise;			// measurement noise covariance
	Matrix9x9 m_measurementMatrix;			// measurement matrix (maps state to measurement space)
public:
	KalmanFilter(const Vector9& initialState, const Matrix9x9& covariance, const Matrix9x9& covariance_pred, const Matrix9x9& processNoise, const Matrix9x9& measurementNoise, const Matrix9x9& measurementMatrix)
		: m_state{ initialState }, m_covariance{ covariance }, m_covariance_pred{ covariance_pred }, m_processNoise{ processNoise }, m_measurementNoise{ measurementNoise }, m_measurementMatrix{ measurementMatrix } {
	}

	Matrix9x9 systemsDynamics() const
	{
		// Define the system dynamics matrix 
		Vector9 F1{ 0,0,0, 1,0,0, 0,0,0 }; // px += vx
		Vector9 F2{ 0,0,0, 0,1,0, 0,0,0 }; // py += vy
		Vector9 F3{ 0,0,0, 0,0,1, 0,0,0 }; // pz += vz
		Vector9 F4{ 0,0,0, 0,0,0, 1,0,0 }; // vx += ax
		Vector9 F5{ 0,0,0, 0,0,0, 0,1,0 }; // vy += ay
		Vector9 F6{ 0,0,0, 0,0,0, 0,0,1 }; // vz += az
		Vector9 F7{ 0,0,0, 0,0,0, 0,0,0 }; // ax = random walk
		Vector9 F8{ 0,0,0, 0,0,0, 0,0,0 };
		Vector9 F9{ 0,0,0, 0,0,0, 0,0,0 };
		return (Matrix9x9(F1, F2, F3, F4, F5, F6, F7, F8, F9));
	}

	Matrix9x9 gainMatrix()
	{
		// compute Kalman gain using the formula K = P * H^T * (H * P * H^T + R)^-1
		return{ m_covariance * m_measurementMatrix.trans() * (m_measurementMatrix * m_covariance * m_measurementMatrix.trans() + m_measurementNoise).inv() };
	}

	Matrix9x9 CovarianceBeforeError(const Matrix9x9& gain)
	{
		// compute new error covariance using the formula P' = (I - K * H) * P
		Matrix9x9 identity{};
		return{ (identity - gain * m_measurementMatrix) * m_covariance_pred };
	}

	Matrix9x9 DiscreteTransition(const double& sampTime)
	{
		Matrix9x9 identity{};
		Matrix9x9 systemDynamics{ systemsDynamics() };
		return{ identity + systemDynamics * sampTime };
	}

	Matrix9x9 CovarianceAfterError(const double& sampTime)
	{
		// compute covariance matrix after update using the formula M = (I + F * Ts) * P' * (I + F * Ts)^T + Q, where F is the system dynamics matrix, Ts is the sample time, P' is the covariance before error update, and Q is the process noise covariance
		return{ (DiscreteTransition(sampTime) * m_covariance * DiscreteTransition(sampTime).trans() + m_processNoise) };
	}

	void predict(double dt)
	{
		Matrix9x9 phi{ DiscreteTransition(dt) };
		m_covariance_pred = phi * m_covariance * phi.trans() + m_processNoise; // P⁻
		m_state = phi * m_state; // x⁻
	}

	void update(const Vector9& measurement)
	{
		Matrix9x9 gain{ m_covariance_pred * m_measurementMatrix.trans() *
						(m_measurementMatrix * m_covariance_pred * m_measurementMatrix.trans() + m_measurementNoise).inv() };
		m_state = m_state + gain * (measurement - m_measurementMatrix * m_state);
		Matrix9x9 identity{};
		m_covariance = (identity - gain * m_measurementMatrix) * m_covariance_pred; // P_post
	}

	Vector9 getState() const { return m_state; }
};