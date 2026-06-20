#pragma once
#include <iostream>
#include <fstream>
#include "motion.h"
#include "linalg.h"
#include "integrator.h"
#include "gnc.h"
#include "filters.h"
#include "settings.h"
#include "radar.h"

struct SimResult {
    double missDistance{};
    double interceptTime{};
    bool hit{ false };
	Vector3 initialPos{};  // perturbed initial position
	Vector3 initialVel{};  // perturbed initial velocity
};

SimResult runSimulation(const SimConfig& cfg, bool logToFile = false)
{
    SimResult result{};

    Vector3 E{ 0.0, 0.0, 0.0 };
    Vector3 B{ 0.0, 0.0, 0.0 };
    Matrix9x9 H{};
    double maxNC{ 0.0 };
    double prevRange{ std::numeric_limits<double>::max() };

    // Initial orientation
    Vector3 vel_dir{ cfg.missileVel.normalize() };
    Vector3 body{ 1.0, 0.0, 0.0 };
    Vector3 axis{ body.crossP(vel_dir) };
    double angle{ std::acos(body.dotP(vel_dir)) };
    Quaternion initialOrientation{};
    if (axis.magnitude() > 1e-10)
        initialOrientation = Quaternion::fromAxisAngle(
            axis.normalize(),
            angle * (180.0 / Constants::pi));

    Forces           forces{ cfg.weight, cfg.crossSec, cfg.length, cfg.diameter,
                             cfg.dragCoeff, cfg.charge, cfg.inertiaTensor };
    ProjectileMotion projectile{ cfg.missilePos, cfg.missileVel, initialOrientation, Vector3{} };
    Target           target{ cfg.targetPos, cfg.targetVel, cfg.targetAcc, cfg.targetOrientation };
    Guidance         guidance{ cfg.navConstant };
    RCSModel         rcs{ cfg.radarWavelength, cfg.diameter / 2.0, cfg.length };
    Seeker           seeker{ cfg.radarPower, cfg.radarGain, cfg.radarWavelength,
                             cfg.radarBandwidth, cfg.radarNoiseFig, cfg.noiseTemp, rcs };
    KalmanFilter     kalman{ cfg.kalmanX0, cfg.kalmanP0, cfg.kalmanP0, cfg.kalmanQ, cfg.kalmanR, H };
    Autopilot        autopilot{ cfg.Kp_att, cfg.Kp_spd, cfg.Kd_att, cfg.V_cmd, cfg.thrust_max };

    double dt0{ cfg.dt };
    double dt{ dt0 };
    double t{};

    // Only open file if logging enabled (main run, not Monte Carlo)
    std::ofstream outFile{};
    if (logToFile)
    {
        outFile.open("trajectory.csv");
        outFile << "t,x,y,z,qw,qx,qy,qz,tx,ty,tz,Vc,range,nCx,nCy,nCz,ox,oy,oz\n";
    }
    std::cout << "Starting simulation, numSteps=" << cfg.numSteps << " dt=" << dt << "\n";
    for (int step{ 0 }; step < cfg.numSteps; ++step)
    {
        t += dt;

        Vector9 measurement{ seeker.measure(projectile, target) };
        kalman.predict(dt);
        kalman.update(measurement);
        Vector9 x_hat{ kalman.getState() };

        Vector3 nC{ guidance.computeAcceleration(x_hat) };
        Vector3 alpha_cmd{ autopilot.attitudeCommand(projectile, nC) };
        Vector3 f_thrust{ autopilot.thrustForce(projectile, forces) };
        Vector3 f_total{ f_thrust };
        Vector3 r_world{ projectile.getOrientation().rotate(cfg.leverArm) };

        rk4Step(projectile, forces, E, B, dt, r_world, f_total, f_thrust, alpha_cmd);
        target.update(dt);

        if (nC.magnitude() > maxNC) maxNC = nC.magnitude();

        Vector3 r_TM{ target.getPosition() - projectile.getPosition() };
        Vector3 v_TM{ target.getVelocity() - projectile.getVelocity() };
        double  Vc{ -r_TM.normalize().dotP(v_TM) };
        double  range{ r_TM.magnitude() };

        // Closest approach detection
        if (range > prevRange)
        {
            result.missDistance = prevRange;
            result.interceptTime = t;
            result.hit = prevRange < 10.0;
            break;
        }
        prevRange = range;

        if (range < 100) dt = dt0 / 10;
        else             dt = dt0;

        // Log only if enabled
        if (logToFile)
        {
            Vector3    pos{ projectile.getPosition() };
            Quaternion orient{ projectile.getOrientation() };
            Vector3    current_dir{ orient.rotate(Vector3{1,0,0}) };
            Vector3    tar{ target.getPosition() };
            outFile << t << "," << pos.getX() << "," << pos.getY() << "," << pos.getZ() << ","
                << orient.getW() << "," << orient.getX() << "," << orient.getY() << "," << orient.getZ() << ","
                << tar.getX() << "," << tar.getY() << "," << tar.getZ() << ","
                << Vc << "," << range << ","
                << nC.getX() << "," << nC.getY() << "," << nC.getZ() << ","
                << current_dir.getX() << "," << current_dir.getY() << "," << current_dir.getZ() << "\n";
        }
    }

    if (logToFile) outFile.close();
    return result;
}
