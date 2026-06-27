#pragma once
#include "linalg.h"
#include "motion.h"

#include <fstream>
#include <string>

struct SimConfig {
    // Raw inputs
    Vector3 missilePos;
    Vector3 missileVel;
    Vector3 targetPos;
    Vector3 targetVel;
    Vector3 targetAcc;
    double  navConstant;
    double  dt;
    int     numSteps;
    double  weight{};
    double  dragCoeff{};
    double  charge{};
    double  thrust_max{};
    double  Kd_att{};
    double  V_cmd{};
    double  length{};
	double  diameter{};
	double  Kp_spd{};
	double  Kp_att{};
	double  mc_sigma_vel{};
	double  mc_sigma_pos{};
	int     mc_runs{};
	double  radarPower{};
	double  radarGain{};
	double  radarWavelength{};
	double  radarBandwidth{};
	double  radarNoiseFig{};
	double  noiseTemp{};
    Quaternion targetOrientation{};  // default identity

    // Derived — computed automatically from the above
    Vector3 relPos;       // targetPos - missilePos
    Vector3 relVel;       // targetVel - missileVel
    Vector9 kalmanX0;     // initial Kalman state [relPos, relVel, targetAcc]
    Matrix9x9 kalmanP0;   // initial covariance
    Matrix9x9 kalmanQ;    // process noise
    Matrix3x3 kalmanR;    // measurement noise
    double sensorSigmaPos;
    double sensorSigmaVel;
    double sensorSigmaAcc;
    Matrix3x3 inertiaTensor;
    Vector3   leverArm;
	double  crossSec;

    // Extracts from settings file for use in main.cpp
    static SimConfig load(const std::string& path)
    {
        std::ifstream file{ path };
        if (!file) throw std::runtime_error("Could not open config: " + path);

        SimConfig cfg{};
        std::string key;

        while (file >> key) {
            if (key == "missile_pos")	   file >> cfg.missilePos;
            else if (key == "missile_vel") file >> cfg.missileVel;
            else if (key == "target_pos")  file >> cfg.targetPos;
            else if (key == "target_vel")  file >> cfg.targetVel;
            else if (key == "target_acc")  file >> cfg.targetAcc;
            else if (key == "nav_constant")file >> cfg.navConstant;
            else if (key == "dt")          file >> cfg.dt;
            else if (key == "num_steps")   file >> cfg.numSteps;
            else if (key == "sigma_pos")   file >> cfg.sensorSigmaPos;
            else if (key == "sigma_vel")   file >> cfg.sensorSigmaVel;
            else if (key == "sigma_acc")   file >> cfg.sensorSigmaAcc;
            else if (key == "weight")      file >> cfg.weight;
            else if (key == "drag_coeff")  file >> cfg.dragCoeff;
            else if (key == "charge")      file >> cfg.charge;
            else if (key == "T_max")       file >> cfg.thrust_max;
            else if (key == "Kd_att")      file >> cfg.Kd_att;
            else if (key == "V_cmd")       file >> cfg.V_cmd;
			else if (key == "length")	   file >> cfg.length;
			else if (key == "diameter")	   file >> cfg.diameter;
			else if (key == "Kp_att")      file >> cfg.Kp_att;
			else if (key == "Kp_spd")	   file >> cfg.Kp_spd;
			else if (key == "mc_runs")	   file >> cfg.mc_runs;
			else if (key == "mc_sigma_pos")file >> cfg.mc_sigma_pos;
			else if (key == "mc_sigma_vel")file >> cfg.mc_sigma_vel;
			else if (key == "radar_power") file >> cfg.radarPower;
			else if (key == "radar_gain")  file >> cfg.radarGain;
			else if (key == "radar_wavelength")
										   file >> cfg.radarWavelength;
			else if (key == "radar_bandwidth")
										   file >> cfg.radarBandwidth;
			else if (key == "radar_noise_fig")
										   file >> cfg.radarNoiseFig;
			else if (key == "noise_temp")  file >> cfg.noiseTemp;

        }

        cfg.derive();
        return cfg;
    }
    void derive()
    {
		crossSec = diameter * diameter / 4 * Constants::pi;
        
        kalmanX0 = Vector9{
            targetPos.getX(), targetPos.getY(), targetPos.getZ(),
            targetVel.getX(), targetVel.getY(), targetVel.getZ(),
            targetAcc.getX(), targetAcc.getY(), targetAcc.getZ()
        };

        // Covariance seeded from sensor sigmas
        double sp2{ sensorSigmaPos * sensorSigmaPos };
        double sv2{ sensorSigmaVel * sensorSigmaVel };
        double sa2{ sensorSigmaAcc * sensorSigmaAcc };
        //double compQ{ sa2 * dt };
        // Set initial Kalman filter matrices
        kalmanP0 = Matrix9x9::diagonal(sp2, sp2, sp2, sv2, sv2, sv2, sa2, sa2, sa2);
        kalmanR = Matrix3x3{ Vector3(sp2, 0, 0),
                             Vector3(0, sp2, 0),
                             Vector3(0, 0, sp2) };

        double Sa{ 100.0 };  // ft²/s⁴ 
        double dt5{ dt * dt * dt * dt * dt / 20.0 };
        double dt3v{ dt * dt * dt / 3.0 };
        double dt4{ dt * dt * dt * dt / 8.0 };
        double dt3{ dt * dt * dt / 6.0 };
        double dt2{ dt * dt / 2.0 };

        kalmanQ = { Matrix9x9{Vector9(dt5, 0.0, 0.0,  dt4, 0.0, 0.0,  dt3, 0.0, 0.0),
                              Vector9(0.0, dt5, 0.0,  0.0, dt4, 0.0,  0.0, dt3, 0.0),
                              Vector9(0.0, 0.0, dt5,  0.0, 0.0, dt4,  0.0, 0.0, dt3),
                              Vector9(dt4, 0.0, 0.0,  dt3v,0.0, 0.0,  dt2, 0.0, 0.0),
                              Vector9(0.0, dt4, 0.0,  0.0, dt3v,0.0,  0.0, dt2, 0.0),
                              Vector9(0.0, 0.0, dt4,  0.0, 0.0, dt3v, 0.0, 0.0, dt2),
                              Vector9(dt3, 0.0, 0.0,  dt2, 0.0, 0.0,  dt,  0.0, 0.0),
                              Vector9(0.0, dt3, 0.0,  0.0, dt2, 0.0,  0.0, dt,  0.0),
                              Vector9(0.0, 0.0, dt3,  0.0, 0.0, dt2,  0.0, 0.0, dt) } * Sa };

        inertiaTensor = Inertia::cylinder(weight / Constants::gravity, diameter / 2, length);
        leverArm = Vector3{ length/2.0, 0.0, 0.0 };

        Vector3 tvel_dir{ targetVel.normalize() };
        Vector3 tbody{ 1.0, 0.0, 0.0 };
        Vector3 taxis{ tbody.crossP(tvel_dir) };
        double  tangle{ std::acos(tbody.dotP(tvel_dir)) };
        if (taxis.magnitude() > 1e-10)
            targetOrientation = Quaternion::fromAxisAngle(
                taxis.normalize(),
                tangle * (180.0 / Constants::pi));
    }
};