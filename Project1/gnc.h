#pragma once
#include "linalg.h"
#include "motion.h"
#include "integrator.h"

class Target {
private:
	Vector3 m_position{};
	Vector3 m_velocity{};
	Vector3 m_acceleration{};

public:
	Target(const Vector3& position, const Vector3& veloctity, const Vector3& acceleration) : m_position{ position }, m_velocity{ veloctity }, m_acceleration{ acceleration } {}

	Vector3 getPosition() const { return m_position; }
	Vector3 getVelocity() const { return m_velocity; }
	Vector3 getAcceleration() const { return m_acceleration; }

	Vector3 update(double dt) // updates target's position and velocity based on its current velocity and acceleration, using simple kinematic equations
	{
		m_position = m_position + (m_velocity * dt) + (m_acceleration * (0.5 * dt * dt)); // update position using s = ut + 0.5at^2
		m_velocity = m_velocity + (m_acceleration * dt); // update velocity using v = u + at
		return m_position; // return updated position for use in guidance system
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
public:
	Guidance(const double N) : m_N{ N } {} // constructor to set navigation constant, with default value of 3.0
	Guidance() = default; // default constructor to allow for default navigation constant

	Vector3 time2go(const ProjectileMotion& p, const Target& t) // estimates time to go until intercept based on current relative position and velocity of projectile and target
	{
		Vector3 r_TM{ t.getPosition() - p.getPosition() }; // relative position vector from projectile to target
		Vector3 v_TM{ t.getVelocity() - p.getVelocity() }; // relative velocity vector from projectile to target
		double closing_velocity{ -r_TM.normalize().dotP(v_TM) }; // closing velocity, which is the rate at which the projectile is closing the distance to the target along the line of sight
		if (closing_velocity <= 0.0)
			return Vector3{ std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() }; // if closing velocity is zero or negative, return infinity to indicate that intercept is not possible
		return r_TM * (1.0 / closing_velocity); // return time to go as distance divided by closing velocity
	}

	// computes the acceleration command for the projectile to intercept the target using proportional navigation guidance law
	Vector3 computeAcceleration(const ProjectileMotion& p, const Target& t)
	{
		Vector3 r_TM{ t.getPosition() - p.getPosition() };
		Vector3 v_TM{ t.getVelocity() - p.getVelocity() };

		double Vc{ -r_TM.normalize().dotP(v_TM) };

		Vector3 los_unit{ r_TM.normalize() };
		Vector3 los_rate_vec{ r_TM.crossP(v_TM) * (1.0 / (r_TM.magnitude() * r_TM.magnitude())) };

		// Cross back with LOS unit vector to get lateral acceleration in the guidance plane
		Vector3 accel{ los_rate_vec.crossP(los_unit) * (m_N * Vc) };
		return accel;
	}

	double missDistance(const ProjectileMotion& p, const Target& t)
	{
		Vector3 r_TM{ t.getPosition() - p.getPosition() }; // relative position vector from projectile to target
		return r_TM.magnitude(); // return the magnitude of the relative position vector as the miss distance
	}
};

class Navigation {
private:
	Guidance m_guidance; // instance of Guidance class to compute acceleration commands for the projectile based on the target's position and velocity
public:
	Navigation(double N) : m_guidance{ N } {}
	// computes the acceleration command for the projectile to intercept the target using the Guidance class's computeAcceleration function, which implements the proportional navigation guidance law
	Vector3 computeCommand(const ProjectileMotion& p, const Target& t) 
	{
		return m_guidance.computeAcceleration(p, t);
	}

};