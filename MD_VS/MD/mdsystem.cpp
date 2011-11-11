
// Standard includes
#include <iostream>

// Own includes
#include "mdsystem.h"

////////////////////////////////////////////////////////////////
// PUBLIC FUNCTIONS
////////////////////////////////////////////////////////////////

mdsystem::mdsystem(int nrparticles_in, float sigma_in, float epsilon_in, float inner_cutoff_in, float outer_cutoff_in, float mass_in, float dt_in, int nrinst_in, float temperature_in, int nrtimesteps_in, float latticeconstant_in, enum_lattice_types lattice_type_in, bool diff_c_on_in, bool Cv_on_in, bool pressure_on_in, bool msd_on_in, bool Ep_on_in, bool Ek_on_in):
    cell_linklist(),
    cell_list(),
    particles(), //TODO: we will resize it later (remove rwo?)
    insttemp(nrinst_in), 
    instEk  (nrinst_in), 
    instEp  (nrinst_in), 
    // Measure values
    temp    (0),
    Ek      (0),
    Ep      (0),
    Cv      (0),
    pressure(0),
    msd     (0),
    verlet_particles_list(), 
    verlet_neighbors_list()
{
	lattice_type = lattice_type_in; // One of the supported lattice types listed in enum_lattice_types
    dt = dt_in;                     // Delta time, the time step to be taken when solving the diff.eq.
    outer_cutoff = outer_cutoff_in; // Parameter for the Verlet list
    inner_cutoff = inner_cutoff_in; // Parameter for the Verlet list

    loop_num = 0;
    nrtimesteps = ((nrtimesteps_in - 1) / nrinst_in + 1) * nrinst_in; // Make the smallest multiple of nrinst_in that has at least the specified size

    temp    .resize(nrtimesteps/nrinst_in + 1); 
    Ek      .resize(nrtimesteps/nrinst_in + 1); 
    Ep      .resize(nrtimesteps/nrinst_in + 1);
    Cv      .resize(nrtimesteps/nrinst_in + 1);
    pressure.resize(nrtimesteps/nrinst_in + 1);
    msd     .resize(nrtimesteps/nrinst_in + 1);
    if (lattice_type == LT_FCC) {
        n = int(std::pow(float(nrparticles_in / 4 ), float( 1.0 / 3.0 )));
        nrparticles = 4*n*n*n;   // Calculate the new number of atoms; all can't fit in the box since n is an integer
    }
    mass = mass_in;
    sigma = sigma_in;
    epsilon = epsilon_in;
    nrinst = nrinst_in;
    init_temp = temperature_in;
    distanceforcesum = 0;
    kB = 1.381e-23f;
    a = latticeconstant_in;
    nrcells = int(n*a/outer_cutoff);
    if (!nrcells) {
        nrcells = 1;
    }
    cellsize = n*a/nrcells;
    diffusion_coefficient = 0;
    diff_c_on = diff_c_on_in;
    Cv_on = Cv_on_in;
    pressure_on = pressure_on_in;
    msd_on = msd_on_in;
    Ep_on = Ep_on_in;
    Ek_on = Ek_on_in;
}

void mdsystem::init() {
    init_particles();
    create_linked_cells();
    create_verlet_list_using_linked_cell_list();
}

void mdsystem::run_simulation() {
    init();
    while (loop_num <= nrtimesteps) {
        force_calculation();
        leapfrog();
        calculate_properties();
        //if (neighbors should be updated) {
        create_linked_cells();
        create_verlet_list_using_linked_cell_list();
        // }
        std::cout << loop_num++ << std::endl;
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

        if (diff_c_on) diffusion_coefficient += dt*particles[i].vel*particles[i].start_vel/(3*nrparticles);

		// Update positions
        particles[i].pos += dt * particles[i].vel;
        if (msd_on) particles[i].no_bound_pos += dt * particles[i].vel;

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
    insttemp[loop_num % nrinst] = mass*sumvsq/(3*nrparticles*epsilon);
    if (Ek_on) instEk[loop_num % nrinst] = mass*sumvsq/(2*epsilon);
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
    uint particle_index = 0;
	verlet_particles_list.resize(nrparticles);
	verlet_neighbors_list.resize(nrparticles*nrparticles);
    for (uint i = 0; i < nrparticles; i++) {
        verlet_particles_list[i] = 0;
    }
	for (uint i = 0; i < verlet_neighbors_list.size(); i++) {
        verlet_neighbors_list[i] = 0;
    }
    for (uint i = 0; i < nrparticles; i++) {
        if (i==0)
            verlet_neighbors_list[0] = 0;
        else
            verlet_neighbors_list[verlet_particles_list[i]] = 0;
        if (i < (nrparticles-1))
            verlet_particles_list[i+1] = verlet_particles_list[i]+1;
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
                        if((distance < outer_cutoff)&&(particle_index > i)) {
                            j += 1;
                            if (i < (nrparticles-1))
                                verlet_particles_list[i+1]=verlet_particles_list[i+1]+1;
                            verlet_neighbors_list[verlet_particles_list[i]] += 1;
                            verlet_neighbors_list[verlet_particles_list[i]+j] = particle_index;
                        }
                        particle_index = cell_linklist[particle_index];
                    }
                }
            }
        }
    }
}

void mdsystem::force_calculation() { //using reduced unit
	for (uint k = 0; k < nrparticles; k++)
	{
		particles[k].acc = fvec3(0, 0, 0);
	}
	float distance = inner_cutoff ;
    float distance_inv = 1/distance ;
    float distance6_inv = pow(distance_inv,6) ;
    float E_cutoff = 4 * distance6_inv * (distance6_inv - 1);
    float mass_inv=1/mass;             
    instEp[loop_num % nrinst] = 0;
    for (uint i=0; i < nrparticles ; i++) { 
        for (uint j = verlet_particles_list[i] + 1; j < verlet_particles_list[i] + verlet_neighbors_list[verlet_particles_list[i]] + 1 ; j++) { 
            fvec3 dr = particles[i].pos-particles[verlet_neighbors_list[j]].pos;
            distance = dr.length();
            if (distance >= inner_cutoff) {
                continue; // Skip this interaction and continue with the next one
            }
            distance_inv = 1 / distance;
            distance6_inv = pow(distance_inv, 6);
			float acceleration = 48 * distance_inv * distance6_inv * (distance6_inv - 0.5f) * mass_inv;
            dr.normalize();
            particles[i].acc +=  acceleration * dr;
			particles[verlet_neighbors_list[j]].acc -=  acceleration * dr;

            if (Ep_on) instEp[loop_num % nrinst] += 4 * distance6_inv * (distance6_inv - 1) - E_cutoff;
			 
            if (pressure_on) distanceforcesum += mass * acceleration * distance;
        }
    }
}

    
void mdsystem::calculate_temperature() {
    float sum = 0;
    for (uint i = 0; i < nrinst; i++) {
        sum += insttemp[i];
    }
    temp[loop_num/nrinst] = sum/nrinst; 
}

void mdsystem::calculate_Ep() {
    float sum = 0;
    for (uint i = 0; i < nrinst; i++) {
        sum += instEp[i];
    }
    Ep[loop_num/nrinst] = sum/nrinst; 
}

void mdsystem::calculate_Ek() {
    float sum = 0;
    for (uint i = 0; i < nrinst; i++) {
        sum += instEk[i];
    }
    Ek[loop_num/nrinst] = sum/nrinst;
}

void mdsystem::calculate_properties() {
    if ((loop_num % nrinst) == 0) {
        if (Cv_on) calculate_specific_heat();            
        if (pressure_on) calculate_pressure();
        if (msd_on) calculate_mean_square_displacement();
        calculate_temperature();
        if (Ep_on) calculate_Ep();
        if (Ek_on) calculate_Ek();
    }
}

void mdsystem::calculate_specific_heat() {
    float T2 = 0;
    for (uint i = 0; i < nrinst; i++){
        T2 += insttemp[i]*insttemp[i];
    }
    T2 = T2/nrinst;
    Cv[loop_num/nrinst] = 9*kB/(6/nrparticles+4-4*T2/(temp[loop_num/nrinst]*temp[loop_num/nrinst]));
}

void mdsystem::calculate_pressure() {
    float V = n*a*n*a*n*a;
    pressure[loop_num/nrinst] = nrparticles*kB*temp[loop_num/nrinst]/V + distanceforcesum/(6*V*nrinst);
    distanceforcesum = 0;
}

void mdsystem::calculate_mean_square_displacement() {
    float sum = 0;
    for (uint i = 0; i < nrparticles;i++) {
        sum += (particles[i].no_bound_pos - particles[i].start_pos)*(particles[i].no_bound_pos - particles[i].start_pos);
    }
    sum = sum/nrparticles;
    msd[loop_num/nrinst] = sum;
}

////////////////////////////////////////////////////////////////
// PRIVATE FUNCTIONS
////////////////////////////////////////////////////////////////

void mdsystem::init_particles() {
    // Allocate space for particles
    particles.resize(nrparticles);

    //Place out particles according to the lattice pattern
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
    
    //Randomixe the velocities
    fvec3 sum_vel = fvec3(0, 0, 0);
    float sum_sqr_vel = 0;
    for (uint i = 0; i < nrparticles; i++) {
        for (uint j = 0; j < 3; j++) {
            particles[i].start_vel[j] = ((float) rand())/((float) RAND_MAX) - 0.5f;
        }
        sum_vel    += particles[i].start_vel;
        sum_sqr_vel += particles[i].start_vel.sqr_length();
    }

    // Compensate for incorrect start temperature and total velocities and finalize the initialization values
    fvec3 average_vel = sum_vel/float(nrparticles);
    float vel_variance = sum_sqr_vel/nrparticles - average_vel.sqr_length();
    float scale_factor = sqrt(3 * init_temp / vel_variance);
    for (uint i = 0; i < nrparticles; i++) {
        particles[i].start_vel = (particles[i].start_vel - sum_vel)*scale_factor;
        particles[i].vel = particles[i].start_vel;
        particles[i].pos = particles[i].start_pos;
    }
}
