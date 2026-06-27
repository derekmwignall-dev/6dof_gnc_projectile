#pragma once
#include "linalg.h"
#include "linalg_mxn.h"

class KalmanFilter {
private:
	Vector9   m_state;						// state vector (position, velocity, etc.)
	Matrix9x9 m_covariance;					// post-update covariance P_post
	Matrix9x9 m_covariance_pred;			// predicted covariance P⁻
	Matrix9x9 m_processNoise;				// process noise covariance
	Matrix3x3 m_measurementNoise;			// measurement noise covariance
	Matrix3x9 m_measurementMatrix;			// measurement matrix (maps state to measurement space)

	// Innovation diagnostics (populated each update step)
	Vector3   m_innovation;					// pre-fit residual: y = z - H*x⁻
	Matrix3x3 m_innovationCov;				// innovation covariance: S = H*P⁻*Hᵀ + R
public:
	KalmanFilter(const Vector9& initialState, const Matrix9x9& covariance, const Matrix9x9& covariance_pred, const Matrix9x9& processNoise, const Matrix3x3& measurementNoise, const Matrix3x9& measurementMatrix)
		: m_state{ initialState }, m_covariance{ covariance }, m_covariance_pred{ covariance_pred }, m_processNoise{ processNoise }, m_measurementNoise{ measurementNoise }, m_measurementMatrix{ measurementMatrix } {
	}
	void setR(Matrix3x3& R) { m_measurementNoise = R; }

	Matrix9x9 systemsDynamics(const double& dt) const
	{
		// Define the system dynamics matrix 
		Vector9 F1{ 1,0,0, dt,0,0, 0.5 * dt * dt,0,0 }; // px += p0x + (vx * dt) + (0.5 * ax * dt^2)
		Vector9 F2{ 0,1,0, 0,dt,0, 0,0.5 * dt * dt,0 }; // py += p0y + (vy * dt) + (0.5 * ay * dt^2)
		Vector9 F3{ 0,0,1, 0,0,dt, 0,0,0.5 * dt * dt }; // pz += p0z + (vz * dt) + (0.5 * az * dt^2)
		Vector9 F4{ 0,0,0, 1,0,0,  dt,0,0 };            // vx += v0x + (ax * dt)
		Vector9 F5{ 0,0,0, 0,1,0,  0,dt,0 };            // vy += v0y + (ay * dt)
		Vector9 F6{ 0,0,0, 0,0,1,  0,0,dt };			// vz += v0z + (az * dt)
		Vector9 F7{ 0,0,0, 0,0,0,  1,0,0 };				// ax += a0x
		Vector9 F8{ 0,0,0, 0,0,0,  0,1,0 };				// ay += a0y
		Vector9 F9{ 0,0,0, 0,0,0,  0,0,1 };				// az += a0z
		return (Matrix9x9(F1, F2, F3, F4, F5, F6, F7, F8, F9));
	}

	void predict(double dt)
	{
		Matrix9x9 phi{ systemsDynamics(dt) };
		m_covariance_pred = phi * m_covariance * phi.trans() + m_processNoise; // P⁻
		m_state = phi * m_state; // x⁻
	}

	void update(const Vector3& measurement)
	{
		// Pre-fit residual (innovation) — computed against the prior x⁻, before update
		m_innovation = measurement - m_measurementMatrix * m_state;  // y = z - H*x⁻

		// Innovation covariance S = H*P⁻*Hᵀ + R
		m_innovationCov = m_measurementMatrix * m_covariance_pred * m_measurementMatrix.trans() + m_measurementNoise;

		// Gate: reject if innovation is too large
		// Chi-squared threshold for 3 DOF at 99.9% = 16.27
		Matrix3x3 S_inv{ m_innovationCov.inv() };
		double mahal{ (m_innovation * S_inv).dotP(m_innovation) }; // eᵀS⁻¹e
		if (mahal > 16.27) return; // skip this measurement

		// Kalman gain K = P⁻*Hᵀ * S⁻¹
		Matrix9x3 gain{ m_covariance_pred * m_measurementMatrix.trans() * m_innovationCov.inv() };

		// State update
		m_state = m_state + gain * m_innovation;

		// Joseph-form covariance update (pre-existing, numerically stable)
		Matrix9x9 identity{};
		Matrix9x9 P1{ identity - gain * m_measurementMatrix };
		Matrix9x9 S{ gain * m_measurementNoise * gain.trans() };
		m_covariance = P1 * m_covariance_pred * P1.trans() + S; // P_post
	}

	Vector9   getState()          const { return m_state; }
	Vector3   getInnovation()     const { return m_innovation; }     // y = z - H*x⁻  [ft]
	Matrix3x3 getInnovationCov()  const { return m_innovationCov; }  // S = H*P⁻*Hᵀ + R

	double computeNEES(const Vector9& true_state) // Normalized estimated error square (NEES) 
	{
		Vector9 error{ true_state - m_state };
		
		return { (error * m_covariance.inv()).dotP(error) };
	}
};