#include <iostream>
#include <fstream>
#include "motion.h"
#include "linalg.h"
#include "integrator.h"
#include "gnc.h"

int main()
{
    // Simulation parameters
    double dt{ 0.01 };
    int numSteps{ 10000 }; // maximum number of simulation steps before forced exit

    // Lever arm: position vector from center of mass to point of force application (for torque calculations)
    Vector3 r{ 0.0, 0.0, 0.5 };

    // Define E and B fields (zero for first test)
    Vector3 E{ 0.0, 0.0, 0.0 };
    Vector3 B{ 0.0, 0.0, 0.0 };

    // Inertia tensor for a sphere
    Matrix3x3 I{ Inertia::sphere(1.0, 0.1) };

    // Initial orientation (identity quaternion for no rotation)
    Quaternion initialOrientation{};

    // Instantiate forces object with mass, cross-sectional area, drag coefficient, charge, and inertia tensor
    Forces forces{ 1.0, 0.01, 0.3, 1.0, I };

    // Instantiate projectile motion -- (position, velocity, orientation, angular_velocity)
    ProjectileMotion projectile{ Vector3{ 0.0, 0.0, 1000 }, Vector3{ 50.0, 0.0, 1.0 }, initialOrientation, Vector3{ 10.0, 0.0, 10.0 } };

    // Instantiate navigation (wraps guidance with N = 3.0)
    Navigation navigation{ 3.0 };

    // Instantiate target with initial position, velocity, and acceleration
    Target target{ Vector3{ 100.0, 0.0, 1000 }, Vector3{ 0.0, 10.0, 0.0 }, Vector3{ 0.0, 0.0, 0.0 } };

    // Open output file
    std::ofstream outFile{ "trajectory.csv" };
    outFile << "t,x,y,z,qw,qx,qy,qz,tx,ty,tz\n"; // CSV header

    for (int step{ 0 }; step < numSteps; ++step)
    {
        double t{ step * dt }; // current simulation time in seconds

        // 1. Compute guidance command via Navigation (not Guidance directly)
        Vector3 nC{ navigation.computeCommand(projectile, target) };
        Vector3 f_guidance{ nC * forces.getMass() };

        // 2. Integrate projectile state forward one timestep
        rk4Step(projectile, forces, E, B, dt, r, f_guidance);

        // 3. Update target state
        target.update(dt);

        // 4. Log current state
        Vector3 pos{ projectile.getPosition() };
        Vector3 tpos{ target.getPosition() };
        Quaternion orient{ projectile.getOrientation() };

        outFile << t << "," << pos.getX() << "," << pos.getY() << "," << pos.getZ() << ","
            << orient.getW() << "," << orient.getX() << "," << orient.getY() << "," << orient.getZ() << ","
            << tpos.getX() << "," << tpos.getY() << "," << tpos.getZ() << '\n';

        // 5. Check termination conditions
        if ((target.getPosition() - projectile.getPosition()).magnitude() < 1.0)
        {
            std::cout << "Intercept at t = " << t << "s\n";
            break;
        }
        else if (projectile.getPosition().getZ() < 0.0)
        {
            std::cout << "Ground impact at t = " << t << "s\n";
            break;
        }
    }

    outFile.close();
    return 0;
}