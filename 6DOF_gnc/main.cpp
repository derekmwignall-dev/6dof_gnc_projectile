#include <iostream>
#include <fstream>
#include "settings.h"
#include "runsimulation.h"
#include "montecarlo.h"

int main()
{
    SimConfig cfg{ SimConfig::load("settings.cfg") };

	MonteCarlo mc{ cfg.mc_runs, cfg.mc_sigma_pos, cfg.mc_sigma_vel };
	mc.run(cfg);
	mc.writeResults("montecarlo.csv");
    // Single logged run
    //SimResult result{ runSimulation(cfg, true) };
    //std::cout << "Miss distance: " << result.missDistance << " ft\n";
    //std::cout << "Intercept time: " << result.interceptTime << "s\n";

    return 0;
}