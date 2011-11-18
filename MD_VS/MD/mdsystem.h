#ifndef  MDSYSTEM_H
#define  MDSYSTEM_H

//Standard includes
#include <vector>
#include <time.h>
using namespace std;

// Own includes
#include "definitions.h"
#include "base_float_vec3.h"
#include "particle.h"

enum enum_lattice_types
{
    LT_NO_LATTICE,
    LT_FCC,
    NUM_LATTICE_TYPES
};

class mdsystem
{
 public:
    // Constructor 
    mdsystem(uint nrparticles_in, ftype sigma_in, ftype epsilon_in, ftype inner_cutoff_in, ftype outer_cutoff_in, ftype mass_in, ftype dt_in, uint nrinst_in, ftype temperature_in, uint nrtimesteps_in, ftype latticeconstant_in, uint lattice_type_in, ftype desiredtemp_in, ftype thermostattime_in, bool thermostat_on_in, bool diff_c_on_in, bool Cv_on_in, bool pressure_on_in, bool msd_on_in, bool Ep_on_in, bool Ek_on_in);

    // Public functions
    void run_simulation();

private:
    //Private variables
    bool         cells_used;            //Flag to tell is the cell list is used or not
    vector<uint> cell_linklist;         //Contains the particle index of the next particle (with decreasing order of the particles) that is in the same cell as the particle the list entry corresponds to. If these is no more particle in the cell, the entry will be 0.
    vector<uint> cell_list;             //Contains the largest particle index each cell contains. The list is coded as if each cell would contain particle zero (although it is probably not located there!)
    vector<particle> particles;         //The elements in the vector particles are particle objects
    vector<uint> verlet_particles_list; //List of integernumber, each index points to an element in the verlet_neighbors_list which is the first neighbor to corresponding particle.
    vector<uint> verlet_neighbors_list; //List with index numbers to neighbors.
    vector<ftype> temp;     // Temperature
    vector<ftype> insttemp; // Instant temperature
    vector<ftype> Cv;       // Heat capacity
    vector<ftype> pressure; // Pressure
    vector<ftype> msd;      // Mean square distance
    vector<ftype> Ek;       // Kinetic energy
    vector<ftype> Ep;       // Potential energy
    vector<ftype> instEk;   // Instat kinetic energy
    vector<ftype> instEp;   // Instat potential energy
    vector<ftype> diffusion_coefficient;
    uint lattice_type; // (enum_lattice_types)
    ftype dt; //length of each timestep
    ftype sqr_inner_cutoff; // Square of the inner cut-off radius in the Verlet list
    ftype sqr_outer_cutoff; // Square of the outer cut-off radius in the Verlet list
    uint loop_num; //gives the current iteration
    uint nrparticles;
    ftype mass;
    ftype sqr_sigma; // Square of sigma in the Lennard Jones potential
    ftype four_epsilon; // Four times epsilon in the Lennard Jones potential
    uint nrinst; //nr of instantaneously measured values before taking the average...
    uint nrcells; //given in one dimension TODO: Change name?
    ftype cellsize;//Could be the same as outer_cutoff but perhaps we should think about that...
    ftype init_temp;
    uint nrtimesteps;
    ftype distanceforcesum;
    ftype a;        // The lattice constant
    uint n;         // Length of one side of the box in conventional unit cells
    ftype box_size; // Length of one side of the box in length units
    ftype p_half_box_size; // Half box side
    ftype n_half_box_size; // Negated half box side
    ftype thermostat;
    ftype desiredtemp;
    ftype thermostattime;
    bool thermostat_on;
    bool diff_c_on;
    bool Cv_on;
    bool pressure_on;
    bool msd_on;
    bool Ep_on;
    bool Ek_on;
    ftype largest_sqr_displacement; // == 2 * maximal displacement for any single atom.

    // Private functions
    void init();
    void init_particles();
    void leapfrog();
    void create_verlet_list();
    void create_linked_cells();
    void create_verlet_list_using_linked_cell_list();
    void force_calculation();
    void calculate_temperature();
    void calculate_Ep();
    void calculate_Ek();
    void calculate_properties();
    void calculate_specific_heat();
    void calculate_pressure();
    void calculate_mean_square_displacement();
    void calculate_diffusion_coefficient();
    vec3 modulos_distance(vec3 pos1, vec3 pos2) const;
    void calculate_largest_sqr_displacement();
    void update_non_modulated_particle_positions();
    inline void update_single_non_modulated_particle_position(uint i);
};

#endif  /* MDSYSTEM_H */
