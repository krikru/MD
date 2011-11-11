// Test2.cpp : Defines the entry point for the console application.
//

// Standard includes
#include <tchar.h>

// Own includes
#include "mdsystem.h"

int _tmain(int argc, _TCHAR* argv[])
{
    int nrparticles_in = 32;
    float sigma_in = 1;
    float epsilon_in = 2;
    float inner_cutoff_in = 1;
    float outer_coutoff_in = 6;
    float mass_in = 2;
    float dt_in = 3;
    int nrinst_in = 3;
    float temperature_in = 300;
    int nrtimesteps_in = 500;
    float latticeconstant_in = 5;
    mdsystem simulation(nrparticles_in, sigma_in, epsilon_in, inner_cutoff_in, outer_coutoff_in, mass_in, dt_in, nrinst_in, temperature_in, nrtimesteps_in, latticeconstant_in, LT_FCC);
    simulation.run_simulation();
    return 0;
}

