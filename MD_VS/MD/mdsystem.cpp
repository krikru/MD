
// Standard includes
#include <iostream>
using std::cout;
using std::endl;

// Own includes
#include "mdsystem.h"

////////////////////////////////////////////////////////////////
// PUBLIC FUNCTIONS
////////////////////////////////////////////////////////////////

mdsystem::mdsystem(int nrparticles_in, float sigma_in, float epsilon_in, float inner_cutoff_in, float outer_cutoff_in, float mass_in, float dt_in, int nrinst_in, float temperature_in, int nrtimesteps_in, float latticeconstant_in, enum_lattice_types lattice_type_in, bool diff_c_on_in, bool Cv_on_in, bool pressure_on_in, bool msd_on_in, bool Ep_on_in, bool Ek_on_in):
    cell_linklist(),
    cell_list    (),
    particles    (),
    insttemp(nrinst_in), 
    instEk  (nrinst_in), 
    instEp  (nrinst_in), 
    // Measure values
    temp    (),
    Ek      (),
    Ep      (),
    Cv      (),
    pressure(),
    msd     (),
    // For the verlet list
    verlet_particles_list(), 
    verlet_neighbors_list()
{
    lattice_type = lattice_type_in; // One of the supported lattice types listed in enum_lattice_types
    dt = dt_in;                     // Delta time, the time step to be taken when solving the diff.eq.
    sqr_outer_cutoff = outer_cutoff_in*outer_cutoff_in; // Parameter for the Verlet list
    sqr_inner_cutoff = inner_cutoff_in*inner_cutoff_in; // Parameter for the Verlet list

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
    sqr_sigma = sigma_in*sigma_in;
    four_epsilon = 4*epsilon_in;
    nrinst = nrinst_in;
    init_temp = temperature_in;
    distanceforcesum = 0;
    kB = 1.381e-23f;
    a = latticeconstant_in;
    box_size = a*n;
    p_half_box_size = 0.5f * box_size;
    n_half_box_size = -p_half_box_size;
    nrcells = int(n*a/outer_cutoff_in);
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
#if 1 //TODO
        cout << "loop number = " << loop_num << endl;
        

        if (loop_num == 9) {
            loop_num = loop_num;
        /*cout << "Ep = "        << Ep[loop_num/nrinst] << endl;
        cout << "Ek = "            << Ek[loop_num/nrinst] << endl;
        cout << "pressure = "    << pressure[loop_num/nrinst] << endl;
        cout << "MSD = "        << msd[loop_num/nrinst] << endl;
        cout << "T = "            << temp[loop_num/nrinst] << endl;
        cout << "Cv = "            << Cv[loop_num/nrinst] << endl;
        */
        }
#endif
        force_calculation();
        //cout << "Box size: " << box_size << endl;
        //cout << "particles[4] x-pos: " << particles[4].pos[0] << endl;
        //cout << "particles[4] y-pos: " << particles[4].pos[1] << endl;
        //cout << "particles[4] z-pos: " << particles[4].pos[2] << endl;
        //cout << endl;
        leapfrog();
        calculate_properties();
        if (1) {
        //if (neighbors should be updated) {
            create_linked_cells();
            create_verlet_list_using_linked_cell_list();
        }
        loop_num++;
    }
}

void mdsystem::leapfrog()
{
    fvec3 zero_vector = fvec3(0, 0, 0);
    float sum_sqr_vel = 0;
    for (uint i = 0; i < nrparticles; i++) {
        //cout << "\ti = " << i << endl;
        if (loop_num == 2) {
            //cout << "i = " << i << endl;
            if (i == 22) {
                i = i; //TODO
            }
        }

        //TODO: Check if vel and pos are stored for the same time or not, in that case, compensate for that

        // Update velocities
        particles[i].vel += dt * particles[i].acc;

        // Update positions
        particles[i].pos += dt * particles[i].vel;

        // Check boundaries in x-dir
        if (particles[i].pos[0] >= box_size) {
            particles[i].pos[0] -= box_size;
            while (particles[i].pos[0] >= box_size) {
                particles[i].pos[0] -= box_size;
            }
        }
        else {
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
        else {
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
        else {
            while (particles[i].pos[2] < 0) {
                particles[i].pos[2] += box_size;
            }
        }

        sum_sqr_vel = sum_sqr_vel + particles[i].vel.sqr_length();
    }
    insttemp[loop_num % nrinst] = mass * sum_sqr_vel / (.75f * nrparticles * four_epsilon);
    if (Ek_on) instEk[loop_num % nrinst] = mass * sum_sqr_vel / (.5f * four_epsilon);
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
                        float sqr_distance = modulos_distance(particles[particle_index].pos, particles[i].pos).sqr_length();
                        if((sqr_distance < sqr_outer_cutoff) && (particle_index > i)) {
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
    // Reset accelrations for all particles
    for (uint k = 0; k < nrparticles; k++) {
        particles[k].acc = fvec3(0, 0, 0);
    }
    float distance;
    float sqr_distance;
    float distance_inv;
    float p = pow(sqr_sigma / sqr_inner_cutoff, 3); // For calculating the cutoff energy
    float E_cutoff = four_epsilon * p * (p - 1);
    float mass_inv = 1/mass;             
    instEp[loop_num % nrinst] = 0;
    for (uint i1 = 0; i1 < nrparticles ; i1++) { // Loop through all particles
        for (uint j = verlet_particles_list[i1] + 1; j < verlet_particles_list[i1] + verlet_neighbors_list[verlet_particles_list[i1]] + 1 ; j++) { 
            uint i2 = verlet_neighbors_list[j]; // Get index of the second (possibly) interacting particle 
            fvec3 r = modulos_distance(particles[i2].pos, particles[i1].pos); // Calculate the closest distance
            sqr_distance = r.sqr_length();
            if (sqr_distance >= sqr_inner_cutoff) {
                continue; // Skip this interaction and continue with the next one
            }
            distance = sqrt(sqr_distance);
            if (loop_num == 2 && (i1 == 22 || i2 == 22)) {
                i1 = i1; //TODO
            }
#if 0 //Emil's code
            distance_inv = sigma / distance;
            distance6_inv = pow(distance_inv, 6);
            float acceleration = 48 * epsilon * distance_inv * distance6_inv * (distance6_inv - 0.5f) * mass_inv; // Emil's formula is incorrect!!! ;P
#endif //Kristofer's code
            distance_inv = 1 / distance;
            p = pow(sqr_sigma * distance_inv * distance_inv, 3);
            float acceleration = 12 * four_epsilon * distance_inv * p * (p - 0.5f) * mass_inv;
            if (acceleration > 1e16f) {
                acceleration = acceleration; //TODO
            }

            // Update accelerations of interacting particles
            fvec3 r_hat = r * distance_inv;
            particles[i1].acc +=  acceleration * r_hat;
            particles[i2].acc -=  acceleration * r_hat;

            // Update properties
            if (Ep_on) instEp[loop_num % nrinst] += four_epsilon * p * (p - 1) - E_cutoff;
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
        sum += (particles[i].pos - particles[i].start_pos).sqr_length();
    }
    sum = sum/nrparticles;
    msd[loop_num/nrinst] = sum;
}

void mdsystem::calculate_diffusion_coefficient() {
    for (int i = 0; i < nrparticles; i++) {
        diffusion_coefficient += (particles[i].pos - particles[i].start_pos)*particles[i].start_vel;
    }
    diffusion_coefficient /= (3*nrparticles);
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
    float scale_factor = sqrt(1.5f * P_KB * init_temp / (0.5f * vel_variance * mass)); // Termal energy = 1.5 * P_KB * init_temp
    for (uint i = 0; i < nrparticles; i++) {
        particles[i].start_vel = (particles[i].start_vel - sum_vel)*scale_factor;
        particles[i].vel = particles[i].start_vel;
        particles[i].pos = particles[i].start_pos;
    }
}

fvec3 mdsystem::modulos_distance(fvec3 pos1, fvec3 pos2) const
{
    fvec3 d = pos2 - pos1;

    // Check boundaries in x-direction
    if (d[0] >= p_half_box_size) {
        d[0] -= box_size;
        while (d[0] >= p_half_box_size) {
            d[0] -= box_size;
        }
    }
    else {
        while (d[0] < n_half_box_size) {
            d[0] += box_size;
        }
    }

    // Check boundaries in y-direction
    if (d[1] >= p_half_box_size) {
        d[1] -= box_size;
        while (d[1] >= p_half_box_size) {
            d[1] -= box_size;
        }
    }
    else {
        while (d[1] < n_half_box_size) {
            d[1] += box_size;
        }
    }

    // Check boundaries in z-direction
    if (d[2] >= p_half_box_size) {
        d[2] -= box_size;
        while (d[2] >= p_half_box_size) {
            d[2] -= box_size;
        }
    }
    else {
        while (d[2] < n_half_box_size) {
            d[2] += box_size;
        }
    }

    return d;
}
