#include "classes.h"

void system::leapfrog() {
3d_vector sumv = 0;
float sumvsq = 0;

for (int i = 0; i < nrparticles; i++) {
	(*particles[i]).vel = (*particles[i]).vel + dt*(*particles[i]).acc;
	(*particles[i]).pos = (*particles[i]).pos + dt*(*particles[i]).vel;
	sumv = sumv + (*particles[i]).vel;
	sumvsq = sumvsq + ((*particles[i]).vel)*((*particles[i]).vel);
	}
insttemp[timestep % nrinst] = mass*sumvsq/(3*nrparticles*epsilon);
instEk[timestep % nrinst] = mass*sumvsq/(2*epsilon);
}

void system::create_linked_cells() {//Assuming origo in the corner of the bulk, and positions given according to boundaryconditions i.e. between zero and lenght of the bulk.
	int cellindex = 0;
	for (int i = 0; i < nrcells; i++) {
		cell_list[i] = 0;
	}
	for (int i = 0; i < nrparticles; i++) {
		cellindex = 1 + (int) particles[i].pos[0]/cellsize
						+ ((int) particles[i].pos[1]/cellsize) * nrcells
						+ ((int) particles[i].pos[2]/cellsize) * nrcells * nrcells;
		cell_linklist[i] = cell_list[cellindex];
		cell_list[cellindex] = i;
	}
}
void system::create_verlet_list_using_linked_cell_list() { // This function ctreates the verlet_lists (verlet_vectors) using the linked cell lists
	int cellindex = 0;
	int particle_index = 0;
	for (int i = 0; i < nrparticles; i++) {
		verlet_particles_list[i] = 0;
	}
	for (int i = 0; i < nrparticles; i++) {
		if (i==0)
			verlet_neighbors_list[0] = 0;
		else
			verlet_neighbors_list[verlet_particles_list[i]] = 0;
		if (i =! (nrparticles-1))
			verlet_particles_list[i+1] = verlet_particles_list[i];
		for (float x =(particles[i].pos[0]-cellsize); x <= (particles[i].pos[0]-cellsize); x += cellsize) {
			for (float y =(particles[i].pos[1]-cellsize); y <= (particles[i].pos[1]-cellsize); y += cellsize) {
				for (float z =(particles[i].pos[2]-cellsize); z <= (particles[i].pos[2]-cellsize); z += cellsize) {
					if (x < 0)
						x = x + nrcells*cellsize;
					if (x > nrcells*cellsize)
						x = x - nrcells*cellsize;
					if (y < 0)
						y = y + nrcells*cellsize;
					if (y > nrcells*cellsize)
						y = y - nrcells*cellsize;
					if (z < 0)
						z = z + nrcells*cellsize;
					if (z > nrcells*cellsize)
						z = z - nrcells*cellsize;
					cellindex = 1 + (int) x/cellsize
								  + ((int) y/cellsize) * nrcells
								  + ((int) z/cellsize) * nrcells * nrcells;
					particle_index = cell_list[cellindex];
					int j = 0;
					while (particle_index != 0) {
						float distance = (particles[i].pos-particles[particle_index].pos).lenght; //Asuming 3d_vector has function lenght that calculates the lenght of a vector.
						if(distance < outer_cutoff) {
							j += 1;
							if (i =! (nrparticles-1))
								verlet_particles_list[i+1]=verlet_particles_list[i+1]+1;
							verlet_neighbors_list[i] += 1;
							verlet_neighbors_list[i+j] = particle_index;
						}
						particle_index = cell_linklist[particle_index];
					}
				}
			}
		}
	}
}

void system::force_calculation() {//Not done

}

system::system(int nrparticles_in, int nrsteps_in, float sigma_in, float epsilon_in, float outer_cutoff_in, float inner_cutoff_in, float mass_in, float dt_in, int nrcells_in, int nrinst_in):
cell_linklist(nrparticles_in),
cell_list(nrcells_in), 
particles(nrparticles_in), 
insttemp(nrinst_in), 
instEk(nrinst_in), 
instEp(nrinst_in), 
temp((int) nrsteps_in/nrinst_in), 
Ek((int) nrsteps_in/nrinst_in), 
Ep((int) nrsteps_in/nrinst_in),
Cv((int) nrsteps_in/nrinst_in),
pressure((int) nrsteps_in/nrinst_in),
msd((int) nrsteps_in/nrinst_in),
verlet_particles_list(), 
verlet_neighbors_list()
{
dt = dt_in;
outer_cutoff = outer_cutoff_in;
inner_cutoff = inner_cutoff_in;
nrsteps = nrsteps_in;
timestep = 1;
nrparticles = nrparticles_in;
mass = mass_in;
sigma = sigma_in;
epsilon = epsilon_in;
nrinst = nrinst_in;
}

system::~system() {

}