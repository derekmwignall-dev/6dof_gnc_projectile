#pragma once
#include <iostream>
#include <fstream>
#include "motion.h"
#include "linalg.h"
#include "integrator.h"
#include "gnc.h"
#include "filters.h"
#include "settings.h"
#include "runsimulation.h"

class MonteCarlo {
private:
	int m_runs{};
	double m_sigma_pos{};
	double m_sigma_vel{};
	std::default_random_engine m_generator{ std::random_device{}() };
	std::vector<SimResult> m_results{}; // stores all run results

public:
	MonteCarlo(int runs, double sigma_pos, double sigma_vel)
		: m_runs{ runs }, m_sigma_pos{ sigma_pos }, m_sigma_vel{ sigma_vel } {}

	void run(const SimConfig& baseCfg)
	{
		m_results.clear();
		m_results.reserve(m_runs);

		// define distributions once outside the loop
		std::normal_distribution<double> pos_dist(0.0, m_sigma_pos);
		std::normal_distribution<double> vel_dist(0.0, m_sigma_vel);


		for (int runNum{ 0 }; runNum < m_runs; ++runNum)
		{
			// copy base config
			SimConfig cfg{ baseCfg };

			// perturb initial conditions
			cfg.missilePos = cfg.missilePos + Vector3{
				pos_dist(m_generator),
				pos_dist(m_generator),
				pos_dist(m_generator)
			};
			cfg.missileVel = cfg.missileVel + Vector3{
				vel_dist(m_generator),
				vel_dist(m_generator),
				vel_dist(m_generator)
			};

			// recompute Kalman X0 from perturbed positions
			cfg.derive();

			// run simulation and store result
			SimResult result{ runSimulation(cfg, false) };
			result.initialPos = cfg.missilePos;  // store perturbed values
			result.initialVel = cfg.missileVel;
			m_results.push_back(result);

			// progress output
			std::cout << "Run " << runNum + 1 << "/" << m_runs
				<< "  miss=" << result.missDistance << " ft"
				<< "  hit=" << (result.hit ? "YES" : "NO") << "\n";
		}

		// print summary
		std::cout << "\n--- Monte Carlo Summary ---\n";
		std::cout << "Runs: " << m_runs << "\n";
		std::cout << "Hits: " << countHits() << "/" << m_runs << "\n";
		std::cout << "CEP:  " << computeCEP() << " ft\n";
		std::cout << "Max miss: " << maxMiss() << " ft\n";
	}

	// 50th percentile of miss distances
	double computeCEP() const
	{
		std::vector<double> distances{};
		for (const auto& r : m_results)
			distances.push_back(r.missDistance);
		std::sort(distances.begin(), distances.end());
		return distances[distances.size() / 2];
	}

	int countHits() const
	{
		int hits{ 0 };
		for (const auto& r : m_results)
			if (r.hit) ++hits;
		return hits;
	}

	double maxMiss() const
	{
		double max{ 0.0 };
		for (const auto& r : m_results)
			if (r.missDistance > max) max = r.missDistance;
		return max;
	}

	void writeResults(const std::string& filename) const
	{
		std::ofstream outFile{ filename };
		outFile << "run,dx,dy,dz,dvx,dvy,dvz,missDistance,interceptTime,hit\n";
		for (int i{ 0 }; i < static_cast<int>(m_results.size()); ++i)
		{
			outFile << i + 1 << ","
				<< m_results[i].initialPos.getX() << ","
				<< m_results[i].initialPos.getY() << ","
				<< m_results[i].initialPos.getZ() << ","
				<< m_results[i].initialVel.getX() << ","
				<< m_results[i].initialVel.getY() << ","
				<< m_results[i].initialVel.getZ() << ","
				<< m_results[i].missDistance << ","
				<< m_results[i].interceptTime << ","
				<< m_results[i].hit << "\n";
		}
	}
};