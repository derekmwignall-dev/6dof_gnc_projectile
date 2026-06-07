#include <iostream>
#include <fstream>
#include "motion.h"
#include "linalg.h"
#include "integrator.h"

int main()
{
    // Simulation parameters
    double dt{ 0.01 };
    int numSteps{ 1000 };
	Vector3 r{ 0.0, 0.0, 0.5 }; // position vector from center of mass to point of force application (for torque calculations)

    // Define E and B fields (zero for first test)
    Vector3 E{ 1.0, 0.0, 10.0 };
    Vector3 B{ 0.0, 0.0, 1.0 };

    // Inertia tensor for a sphere
    Matrix3x3 I{ Inertia::sphere(1.0, 0.1) };

	// Initial orientation (identity quaternion for no rotation)
	Quaternion initialOrientation{};

	// Instantiate forces object with mass, cross-sectional area, drag coefficient, charge, and inertia tensor
    Forces forces{ 1.0, 0.01, 0.3, 1.0, I };

	// Instantiate projectile motion -- (position, velocity, orientation, angular_velocity)
    ProjectileMotion projectile{ Vector3{ 0.0, 0.0, 0.0}, Vector3{ 50.0, 0.0, 50.0}, initialOrientation, Vector3{ 10.0, 0.0, 10.0 } };

    // Open output file
    std::ofstream outFile{ "trajectory.csv" };
    outFile << "t,x,y,z,qw,qx,qy,qz\n"; // CSV header

    // Simulation loop
    for (int step{ 0 }; step < numSteps; ++step)
    {
		double t{ step * dt }; // calculate current time based on step and dt

		rk4Step(projectile, forces, E, B, dt, r); // update projectile's position and velocity using RK4 integrator
		
        // Write current time and position to file
		Vector3 pos{ projectile.getPosition() };
		Quaternion orient{ projectile.getOrientation() };

		// Write current time, position, and orientation to file in CSV format
		outFile << t << "," << pos.getX() << "," << pos.getY() << "," << pos.getZ() << "," << orient.getW() << "," << orient.getX() << "," << orient.getY() << "," << orient.getZ() << '\n';
    }

    outFile.close();
    return 0;
}