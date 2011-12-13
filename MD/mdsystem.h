#ifndef  MDSYSTEM_H
#define  MDSYSTEM_H

//Standard includes
#include <vector>
#include <time.h>
#include <sstream>
using namespace std;

// Own includes
#include "definitions.h"
#include "callback.h"
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
    // Functions that affect the system
    void set_event_callback (callback<void (*)(void*        )> event_callback_in );
    void set_output_callback(callback<void (*)(void*, string)> output_callback_in);
    void init(uint nrparticles_in, ftype sigma_in, ftype epsilon_in, ftype inner_cutoff_in, ftype outer_cutoff_in, ftype mass_in, ftype dt_in, uint nrinst_in, ftype temperature_in, uint nrtimesteps_in, ftype latticeconstant_in, uint lattice_type_in, ftype desiredtemp_in, ftype thermostat_time_in, ftype deltaEp_in, bool thermostat_on_in, bool diff_c_on_in, bool Cv_on_in, bool pressure_on_in, bool msd_on_in, bool Ep_on_in, bool Ek_on_in);
    void run_simulation();
    void abort_activities();
    // Functins that not affect the system
    bool is_operating() const;
    uint get_loop_num() const;
    uint get_max_loops_num() const;

private:
    /*********************
     * Private variables *
     *********************/
    // Thread safety
    bool operating;

    // Comunication with the application
    callback<void (*)(void*        )> event_callback ;
    callback<void (*)(void*, string)> output_callback;
    bool abort_activities_requested;
    stringstream output;
    // The time
    ftype            dt;          // The length of each timestep
    uint             loop_num;    // How many timesteps that has been taken in the simulation
    uint             num_time_steps; // How many timesteps the simulation will take in total
    // The particles
    uint             num_particles; // The number of particles in the system
    uint             lattice_type;  // (enum_lattice_types)
    ftype            particle_mass; // The mass of one atom
    vector<particle> particles;     // The elements in the vector particles are particle objects
    // Initialization
    ftype init_temp;                     // The temperature the system has when it is initialized
    ftype lattice_constant;              // The lattice constant
    uint  box_size_in_lattice_constants; // Length of one side of the box in conventional unit cells //TODO: Move away this variable
    // The box
    ftype box_size;        // Length of one side of the box in length units
    ftype pos_half_box_size; // Half box side
    ftype neg_half_box_size; // Negated half box side
    // Verlet list
    bool         cells_used;            // Flag to tell is the cell list is used or not
    uint         box_size_in_cells;     // Given in one dimension TODO: Change name?
    ftype        cell_size;             // Could be the same as outer_cutoff but perhaps we should think about that...
    vector<uint> cell_linklist;         // Contains the particle index of the next particle (with decreasing order of the particles) that is in the same cell as the particle the list entry corresponds to. If these is no more particle in the cell, the entry will be 0.
    vector<uint> cell_list;             // Contains the largest particle index each cell contains. The list is coded as if each cell would contain particle zero (although it is probably not located there!)
    vector<uint> verlet_particles_list; // List of integernumber, each index points to an element in the verlet_neighbors_list which is the first neighbor to corresponding particle.
    vector<uint> verlet_neighbors_list; // List with index numbers to neighbors.
    ftype        sqr_inner_cutoff;      // Square of the inner cut-off radius in the Verlet list
    ftype        sqr_outer_cutoff;      // Square of the outer cut-off radius in the Verlet list
    // Graphs & measurements
    uint          sample_period;       // Number of timesteps between each measurement
    uint          num_sampling_points; // The number of samples taken for each property
    vector<ftype> temperature;         // Temperature
    vector<ftype> insttemp;            // Instant temperature
    vector<ftype> Cv;                  // Heat capacity
    vector<ftype> pressure;            // Pressure
    vector<ftype> msd;                 // Mean square distance
    vector<ftype> Ek;                  // Kinetic energy
    vector<ftype> Ep;                  // Potential energy
    vector<ftype> cohesive_energy;     // Negative potential energy per atom
    vector<ftype> instEk;              // Instat kinetic energy
    vector<ftype> instEp;              // Instat potential energy
    vector<ftype> diffusion_coefficient;
    // Constrol
    ftype         thermostat_value;  // Varying parameter telling how the velocities should change to adjust the temperature
    ftype         desired_temp;      // The temperature the system strives to obtain
    ftype         thermostat_time;   // The half time for the existing temperature deviation
    vector<ftype> thermostat_values; // To store the values of the thermostat
    // Lennard Jones potential
    ftype sqr_sigma;    // Square of sigma in the Lennard Jones potential
    // 
    ftype distance_force_sum; // For pressure calculation
    ftype dEp_tolerance;      //equilibrium is reached when abs((Ep(current)-Ep(previous))/Ep(current)) is below this value
    bool  equilibrium;
    // Flags
    bool thermostat_on;
    bool diff_c_on;
    bool Cv_on;
    bool pressure_on;
    bool msd_on;
    bool Ep_on;
    bool Ek_on;
    //Reduced unit
    ftype epsilon;
    ftype sigma;
    ftype outer_cutoff;
    ftype inner_cutoff;
    ftype Ep_shift;
    ftype E_cutoff;
    ftype sqr_distance;
    ftype sqr_distance_inv;
    ftype distance_inv;

    /*********************
     * Private functions *
     *********************/
    // Initialization
    void init_particles();
    void potential_energy_shift();
    void potential_energy_cutoff();
    // Verlet list
    void update_verlet_list_if_necessary();
    void create_verlet_list();
    void create_linked_cells();
    void create_verlet_list_using_linked_cell_list();
    void reset_non_modulated_relative_particle_positions();
    inline void reset_single_non_modulated_relative_particle_positions(uint i);
    void update_non_modulated_relative_particle_positions();
    inline void update_single_non_modulated_relative_particle_position(uint i);
    // Simulation
    void leapfrog();
    void force_calculation();
    // Measurements
    void calculate_temperature();
    void calculate_Ep();
    void calculate_cohesive_energy();
    void calculate_Ek();
    void calculate_properties();
    void calculate_specific_heat();
    void calculate_pressure();
    void calculate_mean_square_displacement();
    void calculate_diffusion_coefficient();
    // Output
    ofstream* open_ofstream_file(ofstream &o, const char* path) const;

    // Arithmetic operations
    vec3 modulus_position_minus(vec3 pos1, vec3 pos2) const;

    // Communication with the application
    void print_output_and_process_events();
    void process_events();
    void print_output();

    // Thread safety
    void start_operation();
    void finish_operation();
};

#endif  /* MDSYSTEM_H */
