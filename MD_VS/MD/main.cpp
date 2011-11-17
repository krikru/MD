// main.cpp: Defines the entry point for the console application.
//

// Standard includes
#include <iostream>
#include <tchar.h>

// Own includes
#include "definitions.h"
#include "mdsystem.h"

int _tmain(int argc, _TCHAR* argv[])
{
    // Randomize the simulation with a random seed based on the current time
   // uint random_seed = (unsigned int)time(NULL);
    //srand(random_seed);

    // Init element specific constants
#if 1
    //Let's use the Xenon (Xe) atom in an fcc lattice (Melting point 161.4 K)
    // Element constants
    uint  lattice_type_in = LT_FCC; // (enum_lattice_types)
    ftype sigma_in = 3.98f * P_ANGSTROM;
    ftype epsilon_in = 320e-16f * P_ERG; //1 erg = 10^-7 J
    ftype mass_in = 131.293f * P_U;
    ftype latticeconstant_in = ftype((pow(2.0, 1.0/6.0)*sigma_in) * M_SQRT2);

    // Simulation constants
    ftype dt_in = .1f * P_PS; // [s]
    ftype temperature_in = 100; // [K]
#endif

    // Init simulation specific constants
    uint nrparticles_in = 100; // The number of particles
    uint nrinst_in = 100;       // Number of timesteps between measurements of properties
    uint nrtimesteps_in = 5000; // Desired (or minimum) total number of timesteps
    ftype inner_cutoff_in = 2.0f * sigma_in;
    ftype outer_cutoff_in = 1.01f * inner_cutoff_in; //Fewer neighbors -> faster, but too thin skin is not good either.

    // Init flags
    bool diff_c_on_in = true;
    bool Cv_on_in = true;
    bool pressure_on_in = true;
    bool msd_on_in = true;
    bool Ep_on_in = true;
    bool Ek_on_in = true;

    // Init system and run simulation
    mdsystem simulation(nrparticles_in, sigma_in, epsilon_in, inner_cutoff_in, outer_cutoff_in, mass_in, dt_in, nrinst_in, temperature_in, nrtimesteps_in, latticeconstant_in, lattice_type_in, diff_c_on_in, Cv_on_in, pressure_on_in, msd_on_in, Ep_on_in, Ek_on_in);
    simulation.run_simulation();
    //std::cout << "Random seed " << random_seed << std::endl;
	system("PAUSE");
    return 0;
}

