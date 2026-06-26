#pragma once
#include <iostream>
#include <cmath>
#include <random>
#include "constants.h"
#include "motion.h"
#include "gnc.h"

class RCSModel {
private:
	double m_wavelength{};
	double m_radius{};
	double m_length{};
	mutable std::default_random_engine m_generator;
public:
	RCSModel(double wavelength, double radius, double length)
		: m_wavelength{ wavelength }, m_radius{ radius }, m_length{ length } {
	}


	double computeRCScylinder(double aspect_angle, bool average) const // Compute baseline cross section for a cylinder
	{
		double rcs_broad{ 2*Constants::pi * m_radius * m_length * m_length / m_wavelength };
		double rcs_cap{ 4 * Constants::pi * Constants::pi * Constants::pi * m_radius * m_radius * m_radius * m_radius / m_wavelength  / m_wavelength };

		if (average) return{ rcs_broad/2 + rcs_cap/2};
		else return{ rcs_broad * std::sin(aspect_angle) * std::sin(aspect_angle) + rcs_cap * std::cos(aspect_angle) * std::cos(aspect_angle) };
	}

	double computeRCSsphere() const // Compute baseline cross section for a sphere
	{
		return{ Constants::pi * m_radius * m_radius};
	}

    double applySwerling2(double mean_rcs) const
    {
        if (mean_rcs < 1e-10) return 1e-10; // minimum RCS floor
        std::uniform_real_distribution<double> dist(1e-10, 1.0);
        return -mean_rcs * std::log(dist(m_generator));
    }

	double compute(const Vector3& missilePos, const Vector3& targetPos,
		const Quaternion& targetOrientation) const
	{
		// LOS vector from target to missile
		Vector3 los{ (missilePos - targetPos).normalize() };

		// target body axis (assuming target points along x)
		Vector3 target_axis{ targetOrientation.rotate(Vector3{1,0,0}) };

		// aspect angle between LOS and target body axis
		double aspect_angle{ los.angleBetween(target_axis) };

		// geometric RCS at this aspect angle
		double geometric_rcs{ computeRCScylinder(aspect_angle, false) };

		// apply Swerling 2 fluctuation
		return applySwerling2(geometric_rcs);
	}

    double computeAverageRCS() const
    {
        return computeRCScylinder(0.0, true); // average flag = true
    }

};


// Class generates a noise value for position, velocity, and acceleration based on radar cross-section (RCS) of the target, then adds it to the current pos/vel/acc
class Seeker {
private:
    double m_power{};        // transmit power (watts)
    double m_gain{};         // antenna gain (dimensionless)
    double m_wavelength{};   // radar wavelength (ft)
    double m_bandwidth{};    // receiver bandwidth (Hz)
    double m_noise_fig{};    // noise figure (linear, not dB)
    double m_noise_temp{};   // system noise temperature (K)
    RCSModel m_rcs;
    mutable std::default_random_engine m_generator;

public:
    Seeker(double power, double gain, double wavelength, double bandwidth,
        double noise_fig, double noise_temp, const RCSModel& rcs)
        : m_power{ power }, m_gain{ gain }, m_wavelength{ wavelength },
        m_bandwidth{ bandwidth }, m_noise_fig{ noise_fig },
        m_noise_temp{ noise_temp }, m_rcs{ rcs } {
    }

    double computeSNR(double rcs, double range) const
    {
        if (range < 1e-10) return 1e10; // avoid division by zero

        double P_n{ Constants::k * m_noise_temp * m_bandwidth * m_noise_fig };
        double P_r{ m_power * m_gain * m_gain * m_wavelength * m_wavelength * rcs /
                       (std::pow(4.0 * Constants::pi, 3.0) * std::pow(range, 4.0)) };
        return P_r / P_n;
    }

    double computeNoiseSigma(double snr, double range) const
    {
        double snr_clamped{ std::max(snr, 0.01) }; // minimum SNR floor
        return range / std::sqrt(2.0 * snr_clamped);
    }

    Vector3 measure(const ProjectileMotion& p, const Target& t) const
    {
        // compute RCS at current aspect angle
        double rcs{ m_rcs.compute(p.getPosition(), t.getPosition(),
                                     t.getOrientation()) };

        // compute range
        Vector3 r_TM{ t.getPosition() - p.getPosition() };
        double  range{ r_TM.magnitude() };

        // compute SNR and noise sigma
        double snr{ computeSNR(rcs, range) };
        double sigma_pos{ computeNoiseSigma(snr, range) };

        // guard against invalid sigma values
        sigma_pos = std::max(sigma_pos, 0.01);  // minimum 0.01 ft

        // generate noisy measurement
        std::normal_distribution<double> pos_noise(0.0, sigma_pos);

        Vector3 r{ t.getPosition() };

        return Vector3{
            r.getX() + pos_noise(m_generator),
            r.getY() + pos_noise(m_generator),
            r.getZ() + pos_noise(m_generator)
        };
    }
};