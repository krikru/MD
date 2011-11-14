// main.cpp: Defines the entry point for the console application.
//

// Standard includes
#include <tchar.h>

// Own includes
#include "definitions.h"
#include "mdsystem.h"

int _tmain(int argc, _TCHAR* argv[])
{
    //srand((unsigned int)time(NULL)); // Pick a random seed based on the current time

#if 1
    //Let's use the Xenon (Xe) atom for the fcc lattice (Melting point 161.4 K)
    // Element constants
    float sigma_in = 3.98f * P_ANGSTROM;
    float epsilon_in = 320e-16f * P_ERG; //1 erg = 10^-7 J
    float mass_in = 131.293f * P_U;
    float latticeconstant_in = float((pow(2.0, 1.0/6.0)*sigma_in) * M_SQRT2);

    //Simulation constants
    float dt_in = .001f * P_PS; // [s]
    float temperature_in = 100; // [K]
#endif

    int nrparticles_in = 1000;    // The number of particles
    int nrinst_in = 10;          // Number of timesteps between measurements of properties
    int nrtimesteps_in = 500; // Desired (or minimum) total number of timesteps

    float inner_cutoff_in = 2.5f * sigma_in;
    float outer_cutoff_in = 1.2f * inner_cutoff_in; //Decreased the skin thickness -> fewer neighbors -> faster, but too thin skin not good.
    bool diff_c_on_in = true;
    bool Cv_on_in = true;
    bool pressure_on_in = true;
    bool msd_on_in = true;
    bool Ep_on_in = true;
    bool Ek_on_in = true;
    mdsystem simulation(nrparticles_in, sigma_in, epsilon_in, inner_cutoff_in, outer_cutoff_in, mass_in, dt_in, nrinst_in, temperature_in, nrtimesteps_in, latticeconstant_in, LT_FCC, diff_c_on_in, Cv_on_in, pressure_on_in, msd_on_in, Ep_on_in, Ek_on_in);
    simulation.run_simulation();
    return 0;
}

