#pragma once
#include "linalg.h"

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

    // Derived — computed automatically from the above
    Vector3 relPos;       // targetPos - missilePos
    Vector3 relVel;       // targetVel - missileVel
    Vector9 kalmanX0;     // initial Kalman state [relPos, relVel, targetAcc]
    Matrix9x9 kalmanP0;   // initial covariance
    Matrix9x9 kalmanQ;    // process noise
    Matrix9x9 kalmanR;    // measurement noise
    double sensorSigmaPos;
    double sensorSigmaVel;
    double sensorSigmaAcc;

    // Extracts from settings file for use in main.cpp
    static SimConfig load(const std::string& path)
    {
        std::ifstream file{ path };
        if (!file) throw std::runtime_error("Could not open config: " + path);

        SimConfig cfg{};
        std::string key;

        while (file >> key) {
            if (key == "missile_pos") file >> cfg.missilePos;
            else if (key == "missile_vel") file >> cfg.missileVel;
            else if (key == "target_pos")  file >> cfg.targetPos;
            else if (key == "target_vel")  file >> cfg.targetVel;
            else if (key == "target_acc")  file >> cfg.targetAcc;
            else if (key == "nav_constant") file >> cfg.navConstant;
            else if (key == "dt")          file >> cfg.dt;
            else if (key == "num_steps")   file >> cfg.numSteps;
            else if (key == "sigma_pos")   file >> cfg.sensorSigmaPos;
            else if (key == "sigma_vel")   file >> cfg.sensorSigmaVel;
            else if (key == "sigma_acc")   file >> cfg.sensorSigmaAcc;
        }

        cfg.derive();
        return cfg;
    }

private:
    void derive()
    {
        relPos = targetPos - missilePos;
        relVel = targetVel - missileVel;
        kalmanX0 = Vector9{
            relPos.getX(), relPos.getY(), relPos.getZ(),
            relVel.getX(), relVel.getY(), relVel.getZ(),
            targetAcc.getX(), targetAcc.getY(), targetAcc.getZ()
        };

        // Covariance seeded from sensor sigmas
        double sp2{ sensorSigmaPos * sensorSigmaPos };
        double sv2{ sensorSigmaVel * sensorSigmaVel };
        double sa2{ sensorSigmaAcc * sensorSigmaAcc };

        // Set initial Kalman filter matrices
        kalmanP0 = Matrix9x9::diagonal(sp2, sp2, sp2, sv2, sv2, sv2, sa2, sa2, sa2);
        kalmanQ = Matrix9x9::diagonal(0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 1.0, 1.0, 1.0);
        kalmanR = Matrix9x9::diagonal(sp2, sp2, sp2, sv2, sv2, sv2, sa2, sa2, sa2);
    }
};