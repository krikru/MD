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
    mdsystem();

    /********************
     * Public functions *
     ********************/
    void init(void (*event_handler_in)(void), uint nrparticles_in, ftype sigma_in, ftype epsilon_in, ftype inner_cutoff_in, ftype outer_cutoff_in, ftype mass_in, ftype dt_in, uint nrinst_in, ftype temperature_in, uint nrtimesteps_in, ftype latticeconstant_in, uint lattice_type_in, ftype desiredtemp_in, ftype thermostat_time_in, bool thermostat_on_in, bool diff_c_on_in, bool Cv_on_in, bool pressure_on_in, bool msd_on_in, bool Ep_on_in, bool Ek_on_in);
    void run_simulation();
    void abort_activities();
    bool is_operating() const;

private:
    /*********************
     * Private variables *
     *********************/
    // Thread safety
    bool operating;

    // Comunication with the application
    void (*event_handler)(void);
    bool abort_activities_requested;

    // The time
    ftype            dt;          // The length of each timestep
    uint             loop_num;    // How many timesteps that has been taken in the simulation
    uint             nrtimesteps; // How many timesteps the simulation will take in total
    // The particles
    uint             nrparticles;  // The number of particles in the system
    uint             lattice_type; // (enum_lattice_types)
    ftype            mass;         // The mass of one atom
    vector<particle> particles;    // The elements in the vector particles are particle objects
    // Initialization
    ftype init_temp; // The temperature the system has when it is initialized
    ftype a;         // The lattice constant
    uint  n;         // Length of one side of the box in conventional unit cells //TODO: Move away this variable
    // The box
    ftype box_size;        // Length of one side of the box in length units
    ftype p_half_box_size; // Half box side
    ftype n_half_box_size; // Negated half box side
    // Verlet list
    bool         cells_used;            // Flag to tell is the cell list is used or not
    uint         nrcells;               // Given in one dimension TODO: Change name?
    ftype        cellsize;              // Could be the same as outer_cutoff but perhaps we should think about that...
    vector<uint> cell_linklist;         // Contains the particle index of the next particle (with decreasing order of the particles) that is in the same cell as the particle the list entry corresponds to. If these is no more particle in the cell, the entry will be 0.
    vector<uint> cell_list;             // Contains the largest particle index each cell contains. The list is coded as if each cell would contain particle zero (although it is probably not located there!)
    vector<uint> verlet_particles_list; // List of integernumber, each index points to an element in the verlet_neighbors_list which is the first neighbor to corresponding particle.
    vector<uint> verlet_neighbors_list; // List with index numbers to neighbors.
    ftype        sqr_inner_cutoff;      // Square of the inner cut-off radius in the Verlet list
    ftype        sqr_outer_cutoff;      // Square of the outer cut-off radius in the Verlet list
    // Measurements
    uint          nrinst;   // Number of timesteps between each measurement
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
    // Constrol
    ftype         thermostat;      // Varying parameter telling how the velocities should change to adjust the temperature
    ftype         desiredtemp;     // The temperature the system strives to obtain
    ftype         thermostat_time; // The half time for the existing temperature deviation
    vector<ftype> therm;           // To store the values of the thermostat
    // Lennard Jones potential
    ftype sqr_sigma;    // Square of sigma in the Lennard Jones potential
    ftype four_epsilon; // Four times epsilon in the Lennard Jones potential
    //
    ftype distanceforcesum;
    // Flags
    bool thermostat_on;
    bool diff_c_on;
    bool Cv_on;
    bool pressure_on;
    bool msd_on;
    bool Ep_on;
    bool Ek_on;

    /*********************
     * Private functions *
     *********************/
    // Initialization
    void init_particles();
    // Verlet list
    void update_verlet_list_if_necessary();
    void create_verlet_list();
    void create_linked_cells();
    void create_verlet_list_using_linked_cell_list();
    void update_non_modulated_particle_positions();
    inline void update_single_non_modulated_particle_position(uint i);
    // Simulation
    void leapfrog();
    void force_calculation();
    // Measurements
    void calculate_temperature();
    void calculate_Ep();
    void calculate_Ek();
    void calculate_properties();
    void calculate_specific_heat();
    void calculate_pressure();
    void calculate_mean_square_displacement();
    void calculate_diffusion_coefficient();

    // Arithmetic operations
    vec3 modulos_distance(vec3 pos1, vec3 pos2) const;

    // Communication with the application
    void process_events();

    // Thread safety
    void start_operation();
    void finish_operation();
};

#endif  /* MDSYSTEM_H */
