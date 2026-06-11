#include <iostream>
#include <fstream>
#include "motion.h"
#include "linalg.h"
#include "integrator.h"
#include "gnc.h"
#include "filters.h"
#include "settings.h"

int main()
{
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

    // Instantiate measurement matrix to identity
    Matrix9x9 H{ };

    // Extract initials from settings file
    SimConfig cfg{ SimConfig::load("settings.cfg") };

    ProjectileMotion projectile{ cfg.missilePos, cfg.missileVel, initialOrientation, Vector3{} };
    Target           target{ cfg.targetPos,  cfg.targetVel,  cfg.targetAcc };
    Guidance         guidance{ cfg.navConstant };
    Sensor           sensor{ cfg.sensorSigmaPos, cfg.sensorSigmaVel, cfg.sensorSigmaAcc };
    KalmanFilter     kalman{ cfg.kalmanX0, cfg.kalmanP0, cfg.kalmanP0, cfg.kalmanQ, cfg.kalmanR, H };
    double           dt{ cfg.dt };
    int              numSteps{ cfg.numSteps };

    // Open output file
    std::ofstream outFile{ "trajectory.csv" };
    outFile << "t,x,y,z,qw,qx,qy,qz,tx,ty,tz\n"; // CSV header

    for (int step{ 0 }; step < numSteps; ++step)
    {
        double t{ step * dt };

        // 1. Sense
        Vector9 measurement{ sensor.measure(projectile, target) };

        // 2. Filter
        kalman.predict(dt);
        kalman.update(measurement);
        Vector9 x_hat{ kalman.getState() };

        // 3. Guidance
        Vector3 nC{ guidance.computeAcceleration(x_hat) };
        Vector3 f_guidance{ nC * forces.getMass() };

        // 4. Dynamics
        rk4Step(projectile, forces, E, B, dt, r, f_guidance);

        // 5. Target update
        target.update(dt);

        // 6. Termination
        Vector3 relPos{ target.getPosition() - projectile.getPosition() };
        if (relPos.magnitude() < 1.0) {
            std::cout << "Intercept at t = " << t << "s, miss distance = "
                << relPos.magnitude() << "m\n";
            break;
        }

        // 7. Log
        Vector3 pos{ projectile.getPosition() };
        Quaternion orient{ projectile.getOrientation() };
        Vector3 tar{ target.getPosition() };
        outFile << t << "," << pos.getX() << "," << pos.getY() << "," << pos.getZ() << ","
                            << orient.getW() << "," << orient.getX() << "," << orient.getY() << "," << orient.getZ() << "," 
                            << tar.getX() << "," << tar.getY() << "," << tar.getZ() << "\n";
    }

    outFile.close();
    return 0;
}