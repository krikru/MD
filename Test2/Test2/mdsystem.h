//Standard includes
#include <vector>
#include <time.h>
using namespace std;

// Own includes
#include "definitions.h"
#include "base_float_vec3.h"

enum enum_lattice_types
{
	LT_NO_LATTICE,
	LT_FCC,
	LT_NUM_LATTICE_TYPES
};

class mdsystem
{
 public:
    void leapfrog();
    void create_linked_cells();//see .cpp-file
    void create_verlet_list_using_linked_cell_list();//see .cpp-file
    void force_calculation();//not done
    void md();
    void init();
    void calculate_diffusion_coefficient();
    void calculate_temperature();
    void calculate_Ep();
    void calculate_Ek();
    void initpos();
    void calculate_properties();
    void calculate_specific_heat();
    void calculate_pressure();
    void calculate_mean_square_displacement();
    
    mdsystem(int nrparticles_in, float sigma_in, float epsilon_in, float inner_cutoff_in, float outer_coutoff_in, float mass_in, float dt_in, int nrinst_in, float temperature_in, int nrtimesteps_in, float latticeconstant_in, enum_lattice_types lattice_type_in);

 private:
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
    float inner_cutoff; // Parameter for the Verlet list
    float outer_cutoff; // Parameter for the Verlet list
    uint timestep; //gives the current iteration
    uint n;            //length of lattice in conventional unit cells
    uint nrparticles;
    float mass;
    float sigma;
    float epsilon;
    uint nrinst; //nr of instantaneously measured values before taking the average...
    uint nrcells; //given in one dimension TODO: Change name?
    float cellsize;//Could be the same as outer_cutoff but perhaps we should think about that...
    float init_temp;
    uint nrtimesteps;
    float distanceforcesum;
    float kB;
    float a;
};
