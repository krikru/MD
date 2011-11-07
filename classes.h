#include <vector>
using namespace std;

class particle {//Maybe we should have an additional prev_position which we use when we calculate how far the particles has moved when deciding if we should update the verlet lists
3d_vector pos;  // and if we update the verlet list we should also update the prev_position...
3d_vector prev_pos;
3d_vector vel;
3d_vector prev_vel;
3d_vector acc;
int part_nr; // Might not be needed (part_nr = i is the i-1 element of the vector particles (if part_nr goes from 1 to total number of particles) ) 
};

class system {
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
	
	system(int nrparticles_in, int nrsteps_in, float sigma_in, float epsilon_in, float inner_cutoff_in, float outer_coutoff_in, float mass_in, float dt_in, int nrcells_in, int nrinst_in);
	~system();

private:
	vector<int> cell_linklist;//List with indexnumbers reached through the cell_list, each index number corresponds to a particleindex but also next element in the cell_linklist which are in the same cell. a index equal to zero means end off particle-chain in one cell.
	vector<int> cell_list;//List of integernumber, each index corresponds to an element in the cell_linklist, it is also the "head"-particle in the particle-chain in corresponding cell.
	vector<particle> particles; //The elements in the vector particles are particle objects
	vector<int> verlet_particles_list;//List of integernumber, each index points to an element in the verlet_neighbors_list which is the first neighbor to corresponding particle.
	vector<int> verlet_neighbors_list;//List with index numbers to neighbors.
	vector<float> temp;
	vector<float> insttemp;
	vector<float> Cv;
	vector<float> pressure;
	vector<float> msd;
	vector<float> Ek;
	vector<float> Ep;
	vector<float> instEk;
	vector<float> instEp;
	float dt;//length of each timestep
	float inner_cutoff;
	float outer_cutoff;
	int nrsteps;
	int timestep;//gives the current iteration
	int nrparticles;
	float mass;
	float sigma;
	float epsilon;
	int nrinst;//nr of instantaneously measured values before taking the average...
	int nrcells;//given in one dimension
	float cellsize;//Could be the same as outer_cutoff but perhaps we should think about that...
};