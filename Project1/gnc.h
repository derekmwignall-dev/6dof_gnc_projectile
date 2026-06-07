#pragma once
#include "linalg.h"
#include "motion.h"
#include "integrator.h"

class Target {
private:
	Vector3 m_position{};
public:
	Target(const Vector3& position) : m_position{ position } {}
	Vector3 getPosition() const { return m_position; }

};

class Guidance {
private: 
	const double m_N{ 3.0 }; // navigation constant for proportional navigation
public:
	Guidance(const double N) : m_N{ N } {} // constructor to set navigation constant, with default value of 3.0

	// computes the acceleration command for the projectile to intercept the target using proportional navigation guidance law
	Vector3 computeAcceleration(const ProjectileMotion& p, const Target& t)
	{
		Vector3 r{ t.getPosition() - p.getPosition() }; // relative position vector from projectile to target
		Vector3 v{ p.getVelocity() }; // velocity of projectile
		Vector3 a{ (m_N * v.magnitude() * r.crossP(v).normalize()) * (1/ r.magnitude()) }; // acceleration command based on proportional navigation formula
		return a;
	}
};

class Navigation {
public:
	Vector3 F_cmd(const ProjectileMotion& p, const Target& t)
	{
		Guidance guidance;
		Vector3 a_cmd{ guidance.computeAcceleration(p, t) }; // compute acceleration command from guidance system
		return a_cmd; // return acceleration command for use in integrator to update projectile's velocity and position
	}
};