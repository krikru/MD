
// Own includes
#include "mdsystem.h"

mdsystem::mdsystem(int nrparticles_in, float sigma_in, float epsilon_in, float inner_cutoff_in, float outer_cutoff_in, float mass_in, float dt_in, int nrinst_in, float temperature_in, int nrtimesteps_in, float latticeconstant_in, enum_lattice_types lattice_type_in):
    cell_linklist(1),
    cell_list(1),
    particles(1), //TODO: we will resize it later (remove rwo?)
    insttemp(nrinst_in), 
    instEk(nrinst_in), 
    instEp(nrinst_in), 
    temp((int) nrtimesteps_in/nrinst_in), 
    Ek((int) nrtimesteps_in/nrinst_in), 
    Ep((int) nrtimesteps_in/nrinst_in),
    Cv((int) nrtimesteps_in/nrinst_in),
    pressure((int) nrtimesteps_in/nrinst_in),
    msd((int) nrtimesteps_in/nrinst_in),
    verlet_particles_list(), 
    verlet_neighbors_list()
{
	lattice_type = lattice_type_in;
    dt = dt_in;
    outer_cutoff = outer_cutoff_in;
    inner_cutoff = inner_cutoff_in;
    timestep = 1;
	if (lattice_type == LT_FCC) {
		n = int(std::pow(float(nrparticles_in / 4 ), float( 1.0 / 3.0 )));
		nrparticles = 4*n*n*n;   // Calculate the new number of atoms; all can't fit in the box since n is an integer
	}
    mass = mass_in;
    sigma = sigma_in;
    epsilon = epsilon_in;
    nrinst = nrinst_in;
    init_temp = temperature_in;
    nrtimesteps = nrtimesteps_in;
    distanceforcesum = 0;
    kB = 1.381e-23f;
    a = latticeconstant_in;
    nrcells = int(n*a/outer_cutoff);
    if (!nrcells) {
        nrcells = 1;
    }
    cellsize = n*a/nrcells;
}
/*
  mdsystem::mdsystem(int nrparticles_in, float sigma_in, float epsilon_in, float inner_cutoff_in, float outer_cutoff_in, float mass_in, float dt_in, int nrinst_in, float temperature_in, int nrtimesteps_in, float latticeconstant_in):
  cell_linklist(4*((int) pow(float (nrparticles_in/4), 1/3))*((int) pow(float (nrparticles_in/4), 1/3))*((int) pow(float (nrparticles_in/4), 1/3))),
  cell_list((int((int (pow(float (nrparticles_in/4),float (1/3))))*a/outer_cutoff))*(int((int (pow(float (nrparticles_in/4),float (1/3))))*a/outer_cutoff))*(int((int (pow(float (nrparticles_in/4),float (1/3))))*a/outer_cutoff))),
  //particles(4*((int) pow(float (nrparticles_in/4), 1/3))*((int) pow(float (nrparticles_in/4), 1/3))*((int) pow(float (nrparticles_in/4), 1/3))), 
  verlet_neighbors_list(),
  insttemp(nrinst_in)
  {
  dt = dt_in;
  outer_cutoff = outer_cutoff_in;
  inner_cutoff = inner_cutoff_in;
  timestep = 1;
  n = int (pow(float (nrparticles_in/4),float (1/3))); 
  nrparticles = 4*n*n*n;
  mass = mass_in;
  sigma = sigma_in;
  epsilon = epsilon_in;
  nrinst = nrinst_in;
  init_temp = temperature_in;
  nrtimesteps = nrtimesteps_in;
  distanceforcesum = 0;
  kB = float (1.381e-23);
  a = latticeconstant_in;
  nrcells = int(n*a/outer_cutoff);
  cellsize =float( n*a/nrcells);
  }
*/

void mdsystem::init() {
    fvec3 sumv = fvec3(0, 0, 0);
    float sumvsq = 0;
    initpos();
    srand((unsigned int)time(NULL));
    for (uint i = 0; i < nrparticles; i++) {
        for (uint j = 0; j < 3; j++) {
            particles[i].vel[j] = ((float) rand())/((float) RAND_MAX) - 0.5f;
        }
        sumv += particles[i].vel;
        sumvsq += particles[i].vel*particles[i].vel;
    }
    sumv = sumv/(float(nrparticles));
    sumvsq = sumvsq/nrparticles;
    float s = sqrt(3*init_temp/sumvsq);
    for (uint i = 0; i < nrparticles; i++) {
        particles[i].vel = (particles[i].vel - sumv)*s;
    }
}

void mdsystem::run_simulation() {
    init();
    create_linked_cells();
    create_verlet_list_using_linked_cell_list();
    while (timestep <= nrtimesteps) {
        force_calculation();
        leapfrog();
        calculate_properties();
        //if (neighbors should be updated) {
        create_linked_cells();
        create_verlet_list_using_linked_cell_list();
        // }
        timestep++;
    }
}
void mdsystem::leapfrog()
{
    fvec3 zero_vector = fvec3(0, 0, 0);
    float sumvsq = 0;
	float box_size = a*n; //TODO

    for (uint i = 0; i < nrparticles; i++) {
		// Update velocities
        particles[i].vel += dt * particles[i].acc;

        diffusion_coefficient += dt*particles[i].vel*particles[i].init_vel/(3*nrparticles);

		// Update positions
        particles[i].pos += dt * particles[i].vel;

		// Check boundaries in x-dir
        if (particles[i].pos[0] >= box_size) {
            particles[i].pos[0] -= box_size;
			while (particles[i].pos[0] >= box_size) {
				particles[i].pos[0] -= box_size;
			}
		}
        else if (particles[i].pos[0] < 0) {
            particles[i].pos[0] += box_size;
			while (particles[i].pos[0] < 0) {
				particles[i].pos[0] += box_size;
			}
		}
		// Check boundaries in y-dir
        if (particles[i].pos[1] >= box_size) {
            particles[i].pos[1] -= box_size;
			while (particles[i].pos[1] >= box_size) {
				particles[i].pos[1] -= box_size;
			}
		}
        else if (particles[i].pos[1] < 0) {
            particles[i].pos[1] += box_size;
			while (particles[i].pos[1] < 0) {
				particles[i].pos[1] += box_size;
			}
		}
		// Check boundaries in z-dir
        if (particles[i].pos[2] >= box_size) {
            particles[i].pos[2] -= box_size;
			while (particles[i].pos[2] >= box_size) {
				particles[i].pos[2] -= box_size;
			}
		}
        else if (particles[i].pos[2] < 0) {
            particles[i].pos[2] += box_size;
			while (particles[i].pos[2] < 0) {
				particles[i].pos[2] += box_size;
			}
		}

        sumvsq = sumvsq + particles[i].vel.sqr_length();
	}
    insttemp[timestep % nrinst] = mass*sumvsq/(3*nrparticles*epsilon);
    instEk[timestep % nrinst] = mass*sumvsq/(2*epsilon);
}

void mdsystem::create_linked_cells() {//Assuming origo in the corner of the bulk, and positions given according to boundaryconditions i.e. between zero and lenght of the bulk.
    int cellindex = 0;
    cell_list.resize(nrcells*nrcells*nrcells);
    cell_linklist.resize(nrparticles);
    for (uint i = 0; i < cell_list.size() ; i++) {
        cell_list[i] = 0;
    }
    for (uint i = 0; i < nrparticles; i++) {  //stops here
        int help_x = int(particles[i].pos[0] / cellsize);
        int help_y = int(particles[i].pos[1] / cellsize);
        int help_z = int(particles[i].pos[2] / cellsize);
        cellindex = help_x + help_y * nrcells + help_z * nrcells * nrcells;
        cell_linklist[i] = cell_list[cellindex];
        cell_list[cellindex] = i;
    }
}

void mdsystem::create_verlet_list_using_linked_cell_list() { // This function ctreates the verlet_lists (verlet_vectors) using the linked cell lists
    int cellindex = 0;
    int particle_index = 0;
    for (uint i = 0; i < nrparticles; i++) {
        verlet_particles_list[i] = 0;
    }
    for (uint i = 0; i < nrparticles; i++) {
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
                    cellindex = int(x/cellsize)
                        + (int(y/cellsize)) * nrcells
                        + (int(z/cellsize)) * nrcells * nrcells;
                    particle_index = cell_list[cellindex];
                    int j = 0;
                    while (particle_index != 0) {
                        float distance = (particles[i].pos-particles[particle_index].pos).length(); //Asuming 3d_vector has function lenght that calculates the lenght of a vector.
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

void mdsystem::force_calculation() { //using reduced unit
    vector<float> force_x;
    vector<float> force_y;
    vector<float> force_z;
    float distance = inner_cutoff ;
    float distance_inv = 1/distance ;
    float distance6_inv = pow(distance_inv,6) ;
    float E_cutoff = 4 * distance6_inv * (distance6_inv - 1);
    fvec3 x_hat = fvec3(1, 0, 0);
    fvec3 y_hat = fvec3(0, 1, 0);
    fvec3 z_hat = fvec3(0, 0, 1);                
    for (uint i=0; i < nrparticles ; i++) { 
        for (uint j = verlet_particles_list[i] + 1; j < verlet_particles_list[i] + verlet_neighbors_list[verlet_particles_list[i]] + 1 ; j++) { 
            fvec3 dr = particles[i].pos-particles[j].pos ;
            distance = dr.length() ;
            distance_inv = 1/distance ;
            distance6_inv = pow(distance_inv,6) ;
            float force = 48 * distance_inv * distance6_inv * (distance6_inv - 0.5f) ;
            dr.normalize();
            force_x[i] +=  force * (dr * x_hat) ;
            force_x[j] -=  force * (dr * x_hat) ;
            force_y[i] +=  force * (dr * y_hat) ;
            force_y[j] -=  force * (dr * y_hat) ;
            force_z[i] +=  force * (dr * z_hat) ;
            force_z[j] -=  force * (dr * z_hat) ;
            Ep[i] += 4 * distance6_inv * (distance6_inv - 1) - E_cutoff ; 
            Ep[j] += 4 * distance6_inv * (distance6_inv - 1) - E_cutoff ;
            distanceforcesum += force * distance ; // for pressure
        }
    }
}

    diffusion_coefficient = 0;
void mdsystem::calculate_temperature() {
    float sum = 0;
    for (uint i = 0; i < nrinst; i++) {
        sum += insttemp[i];
    }
    temp[timestep/nrinst] = sum/nrinst; 
}

void mdsystem::calculate_Ep() {
    float sum = 0;
    for (uint i = 0; i < nrinst; i++) {
        sum += instEp[i];
    }
    Ep[timestep/nrinst] = sum/nrinst; 
}

void mdsystem::calculate_Ek() {
    float sum = 0;
    for (uint i = 0; i < nrinst; i++) {
        sum += instEk[i];
    }
    Ek[timestep/nrinst] = sum/nrinst;
}

void mdsystem::initpos() {
    particles.resize(nrparticles);
	if (lattice_type == LT_FCC) {
		for (uint i = 0; i < n; i++) {
			for (uint j = 0; j < n; j++) {
				for (uint k = 0; k < n; k++) {
					int help_index = 4*(i*n*n + j*n + k);
					(particles[help_index + 0]).start_pos[0] = i*a;
					(particles[help_index + 0]).start_pos[1] = j*a;
					(particles[help_index + 0]).start_pos[2] = k*a;

					(particles[help_index + 1]).start_pos[0] = i*a;
					(particles[help_index + 1]).start_pos[1] = (j + 0.5f)*a;
					(particles[help_index + 1]).start_pos[2] = (k + 0.5f)*a;

					(particles[help_index + 2]).start_pos[0] = (i + 0.5f)*a;
					(particles[help_index + 2]).start_pos[1] = j*a;
					(particles[help_index + 2]).start_pos[2] = (k + 0.5f)*a;

					(particles[help_index + 3]).start_pos[0] = (i + 0.5f)*a;
					(particles[help_index + 3]).start_pos[1] = (j + 0.5f)*a;
					(particles[help_index + 3]).start_pos[2] = k*a;
				}
			}
		}
	}
    for (uint i = 0; i < nrparticles; i++) {
        particles[i].pos = particles[i].start_pos;
    }
}


void mdsystem::calculate_properties() {
    if ((timestep % nrinst) == 0) {
        calculate_specific_heat();            
        calculate_pressure();
        calculate_mean_square_displacement();
        calculate_temperature();
        calculate_Ep();
        calculate_Ek();
    }
}

void mdsystem::calculate_specific_heat() {
    float T2 = 0;
    for (uint i = 0; i < nrinst; i++){
        T2 += insttemp[i];
    }
    T2 = T2/nrinst;
    Cv[timestep/nrinst] = 9*kB/(6/nrparticles+4-4*T2/temp[timestep/nrinst]);
}

void mdsystem::calculate_pressure() {
    float V = n*a*n*a*n*a;
    pressure[timestep/nrinst] = nrparticles*kB*temp[timestep/nrinst]/V + distanceforcesum/(6*V*nrinst);
    distanceforcesum = 0;
}

void mdsystem::calculate_mean_square_displacement() {
    float sum = 0;
    for (uint i = 0; i < nrparticles;i++) {
        sum += (particles[i].pos - particles[i].start_pos)*(particles[i].pos - particles[i].start_pos);
    }
    sum = sum/nrparticles;
    msd[timestep/nrinst] = sum;
}
