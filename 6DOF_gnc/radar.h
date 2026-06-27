#pragma once
#include <iostream>
#include <cmath>
#include <random>
#include "constants.h"
#include "motion.h"
#include "gnc.h"

struct SeekerMeasurement {
    Vector3 position;
    double  sigma;
};

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
    double m_beamwidth{};    // antenna beamwidth (radians), derived from gain
    RCSModel m_rcs;
    mutable std::default_random_engine m_generator;

public:
    Seeker(double power, double gain, double wavelength, double bandwidth,
        double noise_fig, double noise_temp, const RCSModel& rcs)
        : m_power{ power }, m_gain{ gain }, m_wavelength{ wavelength },
        m_bandwidth{ bandwidth }, m_noise_fig{ noise_fig },
        m_noise_temp{ noise_temp }, m_rcs{ rcs }, m_beamwidth{ Constants::pi / std::sqrt(m_gain / 0.6) } {
    }

    double computeSNR(double rcs, double range) const
    {
        if (range < 1e-10) return 1e10; // avoid division by zero

        double P_n{ Constants::k * m_noise_temp * m_bandwidth * m_noise_fig };
        double P_r{ m_power * m_gain * m_gain * m_wavelength * m_wavelength * rcs /
                       (std::pow(4.0 * Constants::pi, 3.0) * std::pow(range, 4.0)) };
        return P_r / P_n;
    }

    // Range sigma: derived from bandwidth and SNR (SI internally, output in ft)
    double sigmaRange(double snr) const
    {
        double snr_clamped{ std::max(snr, 0.01) };
        double c_ft{ 9.836e8 };  // speed of light in ft/s
        return c_ft / (2.0 * m_bandwidth * std::sqrt(2.0 * snr_clamped));
    }

    // Cross-range sigma: derived from beamwidth, range, SNR (monopulse)
    double sigmaAngle(double snr) const
    {
        double snr_clamped{ std::max(snr, 0.01) };
        return m_beamwidth / (1.6 * std::sqrt(snr_clamped));
    }

    SeekerMeasurement measure(const ProjectileMotion& p, const Target& t) const
    {
        double rcs{ m_rcs.compute(p.getPosition(), t.getPosition(), t.getOrientation()) };

        Vector3 r_TM{ t.getPosition() - p.getPosition() };
        double  range{ r_TM.magnitude() };

        // Convert RCS ft² -> m² for SNR (radar eq uses SI)
        double rcs_m2{ rcs / (3.28084 * 3.28084) };
        double range_m{ range / 3.28084 };
        double P_n{ Constants::k * m_noise_temp * m_bandwidth * m_noise_fig };
        double lam_m{ m_wavelength / 3.28084 };
        double P_r{ m_power * m_gain * m_gain * lam_m * lam_m * rcs_m2 /
                    (std::pow(4.0 * Constants::pi, 3.0) * std::pow(range_m, 4.0)) };
        double snr{ P_r / P_n };

        // Physically correct sigmas
        double sr{ std::max(sigmaRange(snr), 0.01) };   // along LOS (range)
        double sc{ std::max(range * sigmaAngle(snr), 0.01) }; // cross-range

        // LOS unit vector and two perpendicular axes for noise decomposition
        Vector3 los{ r_TM.normalize() };

        // Build perpendicular axes to LOS
        Vector3 ref{ std::abs(los.getX()) < 0.9 ? Vector3{1,0,0} : Vector3{0,1,0} };
        Vector3 perp1{ los.crossP(ref).normalize() };
        Vector3 perp2{ los.crossP(perp1).normalize() };

        std::normal_distribution<double> range_noise(0.0, sr);
        std::normal_distribution<double> cross_noise(0.0, sc);

        // Noise along LOS (range error) and perpendicular (angle errors)
        double nr{ range_noise(m_generator) };
        double nc1{ cross_noise(m_generator) };
        double nc2{ cross_noise(m_generator) };

        Vector3 pos_true{ t.getPosition() };
        Vector3 pos_noisy{
            pos_true.getX() + los.getX() * nr + perp1.getX() * nc1 + perp2.getX() * nc2,
            pos_true.getY() + los.getY() * nr + perp1.getY() * nc1 + perp2.getY() * nc2,
            pos_true.getZ() + los.getZ() * nr + perp1.getZ() * nc1 + perp2.getZ() * nc2
        };

        // Representative scalar sigma for R (geometric mean of range and cross-range)
        double sigma_rep{ std::sqrt((sr * sr + sc * sc + sc * sc) / 3.0) };

        return SeekerMeasurement{ pos_noisy, sigma_rep };
    }
};

