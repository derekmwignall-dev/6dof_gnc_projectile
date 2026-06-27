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
#include "linalg_mxn.h"

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
    Matrix3x9 H{ Vector9(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                 Vector9(0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                 Vector9(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) };
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
        outFile << "t,x,y,z,qw,qx,qy,qz,tx,ty,tz,ktx,kty,ktz,Vc,range,nCx,nCy,nCz,tvx,tvy,tvz,ktvx,ktvy,ktvz,ox,oy,oz,innov_x,innov_y,innov_z,nees\n";
    }
    std::cout << "Starting simulation, numSteps=" << cfg.numSteps << " dt=" << dt << "\n";
    for (int step{ 0 }; step < cfg.numSteps; ++step)
    {
        t += dt;

        SeekerMeasurement meas{ seeker.measure(projectile, target) };

        double sp{ meas.sigma };
        Matrix3x3 R_k{ Vector3(sp * sp, 0, 0),
                       Vector3(0, sp * sp, 0),
                       Vector3(0, 0, sp * sp) };
        kalman.setR(R_k);

        kalman.predict(dt);
        kalman.update(meas.position);

        Vector9 x_hat{ kalman.getState() };
        Vector3 innov{ kalman.getInnovation() };  // y = z - H*x⁻  (pre-fit residual)

        Vector3 nC{ guidance.computeAcceleration(x_hat,
                                          projectile.getPosition(),
                                          projectile.getVelocity()) };

        double maxG = 30.0 * Constants::gravity;  // structural limit
        if (nC.magnitude() > maxG)
            nC = nC.normalize() * maxG;

        // After computing nC, decompose into:
        // - Axial component (along missile velocity): handled by thrust as-is
        // - Lateral component: apply directly as a lateral force
        Vector3 vel_hat = projectile.getVelocity().normalize();
        Vector3 nC_axial = vel_hat * vel_hat.dotP(nC);
        Vector3 nC_lateral = nC - nC_axial;

        Vector3 alpha_cmd{ autopilot.attitudeCommand(projectile, nC) };
        Vector3 f_thrust{ autopilot.thrustForce(projectile, forces) };
        Vector3 f_lateral = nC_lateral * forces.getMass();
        Vector3 f_total = f_thrust + f_lateral;  // was just f_thrust
        Vector3 r_world{ projectile.getOrientation().rotate(cfg.leverArm) };

        rk4Step(projectile, forces, E, B, dt, r_world, f_total, f_thrust, alpha_cmd);
        target.update(dt);

        if (nC.magnitude() > maxNC) maxNC = nC.magnitude();

        Vector3 r_TM{ target.getPosition() - projectile.getPosition() };
        Vector3 v_TM{ target.getVelocity() - projectile.getVelocity() };
        double  Vc{ -r_TM.normalize().dotP(v_TM) };
        double  range{ r_TM.magnitude() };

        // Closest approach detection
        if (range > prevRange) {
            result.missDistance = prevRange;
            result.interceptTime = t;
            result.hit = prevRange < 10.0;
            break;
        }
        prevRange = range;


        // Log only if enabled
        if (logToFile)
        {
            Vector3    pos{ projectile.getPosition() };
            Quaternion orient{ projectile.getOrientation() };
            Vector3    current_dir{ orient.rotate(Vector3{1,0,0}) };
            Vector3    tar{ target.getPosition() };
            Vector3    tacc{ target.getAcceleration() };
            Vector3    kTar{ x_hat.getX1(), x_hat.getX2(), x_hat.getX3() };Vector3    tvel{ target.getVelocity() };
            Vector3    kTvel{ x_hat.getX4(), x_hat.getX5(), x_hat.getX6() };
            
            Vector9    tar_state{ tar.getX(), tar.getY(), tar.getZ(), tvel.getX(), tvel.getY(), tvel.getZ(), tacc.getX(), tacc.getY(), tacc.getZ() };
            double     nees{ kalman.computeNEES(tar_state) };

            outFile << t << "," << pos.getX() << "," << pos.getY() << "," << pos.getZ() << ","
                << orient.getW() << "," << orient.getX() << "," << orient.getY() << "," << orient.getZ() << ","
                << tar.getX() << "," << tar.getY() << "," << tar.getZ() << ","
                << kTar.getX() << "," << kTar.getY() << "," << kTar.getZ() << ","
                << Vc << "," << range << ","
                << nC.getX() << "," << nC.getY() << "," << nC.getZ() << ","
                << tvel.getX() << "," << tvel.getY() << "," << tvel.getZ() << ","
                << kTvel.getX() << "," << kTvel.getY() << "," << kTvel.getZ() << ","
                << current_dir.getX() << "," << current_dir.getY() << "," << current_dir.getZ() << ","
                << innov.getX() << "," << innov.getY() << "," << innov.getZ() << ","
                << nees << "\n";
        }
    }

    if (logToFile) outFile.close();
    return result;
}