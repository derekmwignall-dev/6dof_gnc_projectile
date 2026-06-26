#pragma once

#include "linalg.h"
#include "motion.h"
#include "integrator.h"
#include "filters.h"
#include <random>

class Target {
private:
	Vector3 m_position{};
	Vector3 m_velocity{};
	Vector3 m_acceleration{};
	Quaternion m_orientation{};

public:
	Target(const Vector3& position, const Vector3& velocity, const Vector3& acceleration, const Quaternion& orientation)
		: m_position{ position }, m_velocity{ velocity }, m_acceleration{ acceleration }, m_orientation{ orientation } {}

	Vector3 getPosition() const { return m_position; }
	Vector3 getVelocity() const { return m_velocity; }
	Vector3 getAcceleration() const { return m_acceleration; }
	Quaternion getOrientation() const { return m_orientation; }

	Vector3 update(double dt)
	{
		m_position = m_position + (m_velocity * dt) + (m_acceleration * (0.5 * dt * dt));
		m_velocity = m_velocity + (m_acceleration * dt);

		// align orientation with velocity
		if (m_velocity.magnitude() > 1e-10)
		{
			Vector3 vel_dir{ m_velocity.normalize() };
			Vector3 body{ 1.0, 0.0, 0.0 };
			Vector3 axis{ body.crossP(vel_dir) };
			double  angle{ std::acos(body.dotP(vel_dir)) };
			if (axis.magnitude() > 1e-10)
				m_orientation = Quaternion::fromAxisAngle(
					axis.normalize(),
					angle * (180.0 / Constants::pi));
		}
		return m_position;
	}
	
	double computeLatAngle(const ProjectileMotion& p) const // computes the line-of-sight angle between the projectile and the target, which is used in the guidance system to determine the direction of the acceleration command
	{
		Vector3 r_TM{ m_position - p.getPosition() }; // relative position vector from projectile to target
		return std::atan2(r_TM.getY(), r_TM.getX()); // return the angle in radians 
	}
};

class Guidance {
private: 
	const double m_N{ 3.0 }; // navigation constant for proportional navigation
	Vector3 m_aT_smoothed{};
public:
	Guidance(const double N) : m_N{ N } {} // constructor to set navigation constant, with default value of 3.0
	Guidance() = default; // default constructor to allow for default navigation constant

	double time2go(const ProjectileMotion& p, const Target& t) // estimates time to go until intercept based on current relative position and velocity of projectile and target
	{
		Vector3 r_TM{ t.getPosition() - p.getPosition() }; // relative position vector from projectile to target
		Vector3 v_TM{ t.getVelocity() - p.getVelocity() }; // relative velocity vector from projectile to target
		double Vc{ -r_TM.normalize().dotP(v_TM) }; // closing velocity

		if (Vc <= 0.0)
			return std::numeric_limits<double>::infinity();
		return r_TM.magnitude() / Vc;
	}

	// computes the acceleration command for the projectile to intercept the target using ZEM-based augmented proportional navigation guidance law
	Vector3 computeAcceleration(const Vector9& filteredState,
		const Vector3& missilePos,
		const Vector3& missileVel)
	{
		Vector3 tgt_pos{ filteredState.get(0), filteredState.get(1), filteredState.get(2) };
		Vector3 tgt_vel{ filteredState.get(3), filteredState.get(4), filteredState.get(5) };
		Vector3 a_T{ filteredState.get(6), filteredState.get(7), filteredState.get(8) };

		Vector3 r{ tgt_pos - missilePos };
		Vector3 v{ tgt_vel - missileVel };

		double Vc{ -r.normalize().dotP(v) };
		if (Vc <= 0.0) return Vector3{ 0.0, 0.0, 0.0 };

		double t_go{ r.magnitude() / Vc };
		//if (t_go < 0.3) return Vector3{ 0.0, 0.0, 0.0 };

		// Gravity vector in inertial frame
		Vector3 gravity{ 0.0, 0.0, -Constants::gravity };
		double alpha = 0.05;  // smoothing factor (small = slow/stable)

		m_aT_smoothed = m_aT_smoothed * (1.0 - alpha) + a_T * alpha;
		Vector3 ZEM{ r + v * t_go + m_aT_smoothed * (0.5 * t_go * t_go) 
					- gravity * (0.5 * t_go * t_go) };

		return ZEM * (m_N / (t_go * t_go));
	}

	double missDistance(const ProjectileMotion& p, const Target& t)
	{
		Vector3 r_TM{ t.getPosition() - p.getPosition() }; // relative position vector from projectile to target
		return r_TM.magnitude(); // return the magnitude of the relative position vector as the miss distance
	}
};

// Class controls the thrust of the projectile, which determines how to achieve a desired accleration and direction through thrust and orientation
class Autopilot {
private:
	double m_Kp_att{};			// attitude control gain
	double m_Kp_spd{};				// speed control gain
	double m_Kd_att{};					// proportional gain
	double m_V_cmd{};				// commanded speed
	double m_thrust_max{};

public:
	Autopilot(double Kp_att, double Kp_spd, double Kd_att, double V_cmd, double thrust_max)
		: m_Kp_att{ Kp_att }, m_Kp_spd{ Kp_spd }, m_Kd_att{ Kd_att }, m_V_cmd{ V_cmd }, m_thrust_max{ thrust_max } {}

	Vector3 attitudeCommand(const ProjectileMotion& p, const Vector3& nC)
	{
		// nC already contains gravity compensation from guidance
		// do NOT add gravity here — it causes double-correction
		Vector3 desired_dir{ nC.normalize() };
		if (nC.magnitude() < 1e-10)
			desired_dir = p.getVelocity().normalize();

		Vector3 body_axis{ 1.0, 0.0, 0.0 };
		Vector3 current_dir{ p.getOrientation().rotate(body_axis) };

		Vector3 error{ current_dir.crossP(desired_dir) };
		Vector3 omega{ p.getAngularVelocity() };

		return error * m_Kp_att - omega * m_Kd_att;
	}

	double speedCommand(const Forces& force, const ProjectileMotion& p) // How the projectile should accomplish what we want (i.e. required speed)
	{
		double V_current{ p.getVelocity().magnitude() };
		double V_error{ m_V_cmd - V_current };

		double sin_gamma{ -p.getVelocity().getZ() / V_current };

		// base thrust: overcome drag + gravity component + speed error
		double T{ m_Kp_spd * V_error
				  + force.dragForce(p.getPosition(), p.getVelocity()).magnitude()
				  + force.gravForce().magnitude() * sin_gamma };

		if (T >= m_thrust_max) return m_thrust_max;
		if (T <= 0.0) return 0.0;
		return T;
	}

	Vector3 thrustForce(const ProjectileMotion& p, const Forces& forces) // Translate to thrust force and directional change
	{
		double T{ speedCommand(forces, p) };
		Vector3 body_axis{ 1.0, 0.0, 0.0 };
		Vector3 thrust_dir{ p.getOrientation().rotate(body_axis) };
		return thrust_dir * T;
	}
};