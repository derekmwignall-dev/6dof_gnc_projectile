#pragma once
#include "motion.h"

void rk4Step(ProjectileMotion& projectile, const Forces& forces,
	const Vector3& E, const Vector3& B, double dt,
	const Vector3& r, const Vector3& f_total,
	const Vector3& f_thrust, const Vector3& alpha_cmd)
{
	// Compute k1 values
	Vector3 k1_pos{ projectile.getPosition() };
	Vector3 k1_vel{ projectile.getVelocity() };
	Vector3 f1{ forces.gravForce() + forces.dragForce(k1_pos, k1_vel) + forces.lorentzForce(k1_vel, E, B) + forces.liftForce(k1_pos, k1_vel, projectile.getOrientation()) + f_total};
	Vector3 k1_acc{ f1 * (1.0/forces.getMass())};
	Vector3 k1_ang_vel{ projectile.getAngularVelocity() };
	Vector3 k1_ang_acc{ forces.angularAccel(forces.torque(r, f_thrust)) + alpha_cmd };

	// Compute k2 values
	Vector3 k2_pos{ projectile.getPosition() + (k1_vel * (dt / 2.0)) };
	Vector3 k2_vel{ projectile.getVelocity() + (k1_acc * (dt / 2.0)) };
	Vector3 f2{ forces.gravForce() + forces.dragForce(k2_pos, k2_vel) + forces.lorentzForce(k2_vel, E, B) + forces.liftForce(k2_pos, k2_vel, projectile.getOrientation()) + f_total };
	Vector3 k2_acc{ f2 * (1.0 / forces.getMass()) };
	Vector3 k2_ang_vel{ projectile.getAngularVelocity() + (k1_ang_acc * (dt / 2.0)) };
	Vector3 k2_ang_acc{ forces.angularAccel(forces.torque(r, f_thrust)) + alpha_cmd };

	// Compute k3 values
	Vector3 k3_pos{ projectile.getPosition() + (k2_vel * (dt / 2.0)) };
	Vector3 k3_vel{ projectile.getVelocity() + (k2_acc * (dt / 2.0)) };
	Vector3 f3{ forces.gravForce() + forces.dragForce(k3_pos, k3_vel) + forces.lorentzForce(k3_vel, E, B) + forces.liftForce(k3_pos, k3_vel, projectile.getOrientation()) + f_total };
	Vector3 k3_acc{ f3 * (1.0 / forces.getMass()) };
	Vector3 k3_ang_vel{ projectile.getAngularVelocity() + (k2_ang_acc * (dt / 2.0)) };
	Vector3 k3_ang_acc{ forces.angularAccel(forces.torque(r, f_thrust)) + alpha_cmd };

	// Compute k4 values
	Vector3 k4_pos{ projectile.getPosition() + (k3_vel * dt) };
	Vector3 k4_vel{ projectile.getVelocity() + (k3_acc * dt) };
	Vector3 f4{ forces.gravForce() + forces.dragForce(k4_pos, k4_vel) + forces.lorentzForce(k4_vel, E, B) + forces.liftForce(k4_pos, k4_vel, projectile.getOrientation()) + f_total };
	Vector3 k4_acc{ f4 * (1.0 / forces.getMass()) };
	Vector3 k4_ang_vel{ projectile.getAngularVelocity() + (k3_ang_acc * dt) };
	Vector3 k4_ang_acc{ forces.angularAccel(forces.torque(r, f_thrust)) + alpha_cmd };

	// Update position and velocity using weighted average of slopes
	Vector3 delta_pos{ (k1_vel + (k2_vel * 2.0) + (k3_vel * 2.0) + k4_vel) * (dt / 6.0) };
	Vector3 delta_vel{ (k1_acc + (k2_acc * 2.0) + (k3_acc * 2.0) + k4_acc) * (dt / 6.0) };
	Vector3 delta_ang_vel{ (k1_ang_acc + (k2_ang_acc * 2.0) + (k3_ang_acc * 2.0) + k4_ang_acc) * (dt / 6.0) };
	Quaternion delta_orient{};
	if (delta_ang_vel.magnitude() > 1e-10)
		delta_orient = Quaternion::fromAxisAngle(delta_ang_vel.normalize(), delta_ang_vel.magnitude());
	// else leave as identity quaternion

	projectile.updatePosition(delta_pos);
	projectile.updateVelocity(delta_vel);
	projectile.updateAngularVelocity(delta_ang_vel);
	projectile.updateOrientation(delta_orient);
}
	