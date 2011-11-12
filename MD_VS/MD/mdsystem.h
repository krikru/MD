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
	LT_NUM_LATTICE_TYPES
};

class mdsystem
{
 public:
    // Constructor 
    mdsystem(int nrparticles_in, float sigma_in, float epsilon_in, float inner_cutoff_in, float outer_cutoff_in, float mass_in, float dt_in, int nrinst_in, float temperature_in, int nrtimesteps_in, float latticeconstant_in, enum_lattice_types lattice_type_in, bool diff_c_on_in, bool Cv_on_in, bool pressure_on_in, bool msd_on_in, bool Ep_on_in, bool Ek_on_in);

    // Public functions
    void init();
    void leapfrog();
    void create_linked_cells();//see .cpp-file
    void create_verlet_list_using_linked_cell_list();//see .cpp-file
    void force_calculation();//not done
    void run_simulation();
    void calculate_temperature();
    void calculate_Ep();
    void calculate_Ek();
    void calculate_properties();
    void calculate_specific_heat();
    void calculate_pressure();
    void calculate_mean_square_displacement();

private:
    //Private variables
    vector<uint> cell_linklist;         //List with indexnumbers reached through the cell_list, each index number corresponds to a particleindex but also next element in the cell_linklist which are in the same cell. a index equal to zero means end off particle-chain in one cell.
    vector<uint> cell_list;             //List of indexnumbers, each index corresponds to an element in the cell_linklist, it is also the "head"-particle in the particle-chain in corresponding cell.
    vector<particle> particles;         //The elements in the vector particles are particle objects
    vector<uint> verlet_particles_list; //List of integernumber, each index points to an element in the verlet_neighbors_list which is the first neighbor to corresponding particle.
    vector<uint> verlet_neighbors_list; //List with index numbers to neighbors.
    vector<float> temp;     // Temperature
    vector<float> insttemp; // Instant temperature
    vector<float> Cv;       // Heat capacity
    vector<float> pressure; // Pressure
    vector<float> msd;      // Mean square distance
    vector<float> Ek;       // Kinetic energy
    vector<float> Ep;       // Potential energy
    vector<float> instEk;   // Instat kinetic energy
    vector<float> instEp;   // Instat potential energy
    enum_lattice_types lattice_type;
    float dt; //length of each timestep
    float sqr_inner_cutoff; // Square of the inner cut-off radius in the Verlet list
    float sqr_outer_cutoff; // Square of the outer cut-off radius in the Verlet list
    uint loop_num; //gives the current iteration
    uint nrparticles;
    float mass;
    float sqr_sigma; // Square of sigma in the Lennard Jones potential
    float four_epsilon; // Four times epsilon in the Lennard Jones potential
    uint nrinst; //nr of instantaneously measured values before taking the average...
    uint nrcells; //given in one dimension TODO: Change name?
    float cellsize;//Could be the same as outer_cutoff but perhaps we should think about that...
    float init_temp;
    uint nrtimesteps;
    float distanceforcesum;
    float kB;
    float a;        // The lattice constant
    uint n;         // Length of one side of the box in conventional unit cells
    float box_size; // Length of one side of the box in length units
    float p_half_box_size; // Half box side
    float n_half_box_size; // Negated half box side
    float diffusion_coefficient;
    bool diff_c_on;
    bool Cv_on;
    bool pressure_on;
    bool msd_on;
    bool Ep_on;
    bool Ek_on;

    // Private functions
    void init_particles();
    fvec3 modulos_distance(fvec3 pos1, fvec3 pos2) const;
};

#endif  /* MDSYSTEM_H */
