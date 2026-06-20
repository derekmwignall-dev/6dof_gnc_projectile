#pragma once
#include "linalg.h"
#include "constants.h"
#include <cmath>
#include <iostream>

// Class that performs calculations related to the motion of a projectile, including its position, velocity, orientation, and angular velocity.
// It also includes functions to update these values based on the forces applied to the projectile and its trajectory.
class ProjectileMotion {
private:
	Vector3 m_position {};
	Vector3 m_velocity {};
	Vector3 m_angular_velocity{};
	Quaternion m_orientation{};

public:
	ProjectileMotion(Vector3 position, Vector3 velocity, Quaternion orientation, Vector3 angular_velocity)
		: m_position{ position }, m_velocity{ velocity }, m_orientation{ orientation }, m_angular_velocity{ angular_velocity } {
	}

	// extracts values from class for use in main()
	Vector3 getPosition() const { return m_position; }
	Vector3 getVelocity() const { return m_velocity; }
	Vector3 getAngularVelocity() const { return m_angular_velocity; }
	Quaternion getOrientation() const { return m_orientation; }

	// sets values in class for use in main() or for updating position and velocity based on current forces applied or trajectory
	void setPosition(const Vector3& position) { m_position = position; }
	void setVelocity(const Vector3& velocity) { m_velocity = velocity; }
	void setOrientation(const Quaternion& orientation) { m_orientation = orientation; }
	void setAngularVelocity(const Vector3& angular_velocity) { m_angular_velocity = angular_velocity; }

	// updates position and velocity based on current forces applied or trajectory, which is done by adding the change in position and velocity (delta) to the current position and velocity. 
	// The delta is calculated as v*dt for position and a*dt for velocity, where dt is the time step of the simulation.
	void updatePosition(const Vector3& delta) { m_position = m_position + delta; } // passes v*dt as delta
	void updateVelocity(const Vector3& delta) { m_velocity = m_velocity + delta; } // passes a*dt as delta
	void updateOrientation(const Quaternion& delta) { m_orientation = (delta * m_orientation).normalize(); } // passes a quaternion representing the change in orientation as delta, which is calculated from the angular velocity and time step
	void updateAngularVelocity(const Vector3& delta) { m_angular_velocity = m_angular_velocity + delta; } // passes alpha*dt as delta
};

// Calculates forces on object 
class Forces {
private:
	double m_cross_sec {};
	double m_weight {}; 
	double m_length {};
	double m_diameter {};
	double m_drag_coefficient {};	// can be calculated with Cd = D / (A * .5 * r * V^2), but currently just set as an input
	double m_charge { 0 };
	Matrix3x3 m_inertia_tensor{};	// can be calculated based on the shape and mass distribution of the object, but currently just set as an input
public:
	Forces(double weight, double cross_sec, double length, double diameter, double drag_coefficient, double charge, Matrix3x3 inertia_tensor)
		: m_weight{ weight }, m_cross_sec{ cross_sec }, m_length{ length }, m_diameter{ diameter }, m_drag_coefficient{ drag_coefficient }, m_charge{ charge }, m_inertia_tensor {
		inertia_tensor
	} {
	}

	// extracts values from class for use in main()
	double getMass() const { return m_weight/Constants::gravity; }
	double getWeight() const { return m_weight; }
	double getCrossSec() const { return m_cross_sec; }
	double getDragCoefficient() const { return m_drag_coefficient; }
	double getCharge() const { return m_charge; }

	Matrix3x3 getInertiaTensor() const { return m_inertia_tensor; }

	// Find the air density based on an approximation at different altitudes
	double getAirDensity(const Vector3& position) const
	{
		double altitude{ position.getZ() };
		altitude = std::max(altitude, 0.0); // prevent negative altitude
		if (altitude < 36089.0)
			return 0.002377 * std::pow((1.0 - altitude / 145442.0), 4.255876);
		else
			return 0.0007087 * std::exp(-(altitude - 36089.0) / 20806.0);
	}

	// Function for gravitational forces
	Vector3 gravForce() const
	{
		return Vector3{ 0.0, 0.0, -1 } * m_weight; // only acts in -z-direction (downwards)
	}

	// Function for drag forces
	Vector3 dragForce(const Vector3& position, const Vector3& velocity) const
	{
		double coefficient { 0.5 * getAirDensity(position) * m_cross_sec * m_drag_coefficient }; // calculate scalar value of drag force
		double speed { velocity.magnitude() };
		return (-coefficient * speed * speed) * velocity.normalize();
	}
	
	// Function for lift force
	Vector3 liftForce(const Vector3& position, const Vector3& velocity, const Quaternion& orientation) const
	{
		Vector3 body_axis{ 1.0, 0.0, 0.0 };
		Vector3 current_dir{ orientation.rotate(body_axis) };

		double S_plan{ m_length*m_diameter }; // Planform area (approximated for a cylinder)
		double S_ref{ Constants::pi * m_diameter * m_diameter / 4.0 }; // Reference area (cylinder)
		
		double angle_of_attack{ velocity.angleBetween(current_dir)}; // angle between current orientation direction and current velocity direction

		double aoa_clamped{ std::min(angle_of_attack, 0.3) }; // max ~17 degrees
		double lift_coefficient{ (2.0 * aoa_clamped) +
								 (1.5 * S_plan * aoa_clamped * aoa_clamped / S_ref) };

		double dynamicP{ 0.5 * getAirDensity(position) * velocity.magnitude() * velocity.magnitude() };

		Vector3 velDir = velocity.normalize();
		Vector3 rejection = current_dir - velDir * velDir.dotP(current_dir);

		if (rejection.magnitude() < 1e-10 || velocity.magnitude() < 1e-10)
			return Vector3{ 0.0, 0.0, 0.0 };

		Vector3 liftDir = rejection.normalize(); // Induce force perpendicularly to the current relative airflow (i.e. upwards)
		Vector3 result{ (lift_coefficient * dynamicP * S_ref) * liftDir }; 
		
		return result;
	}

	// Function for Lorentz Forces
	Vector3 lorentzForce(const Vector3& velocity, const Vector3& E, const Vector3& B) const
	{
		return m_charge*(E + (velocity.crossP(B)));
	}

	// Determine torque in missile
	Vector3 torque(const Vector3& r, const Vector3& F) const
	{
		return r.crossP(F);
	}

	// calculates angular acceleration based on the inverse of the inertia tensor and the torque, which can be used to update the angular velocity
	Vector3 angularAccel(const Vector3& torque) const
	{
		return m_inertia_tensor.inv() * torque; 
	}
};


namespace Inertia {
	inline Matrix3x3 sphere(double mass, double radius) {
		double I{ 0.4 * mass * radius * radius };
		return Matrix3x3{ Vector3{I,0,0}, Vector3{0,I,0}, Vector3{0,0,I} };
	}
	inline Matrix3x3 cylinder(double mass, double radius, double height)
	{
		double Iaxial{ 0.5 * mass * radius * radius };
		double Iperp{ 0.25 * mass * radius * radius + (1.0 / 12.0) * mass * height * height };
		return Matrix3x3{ Vector3{Iperp,0,0}, Vector3{0,Iperp,0}, Vector3{0,0,Iaxial} };
	}
}