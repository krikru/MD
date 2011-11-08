#include "classes.h"

void system::leapfrog()
{
    3d_vector sumv = 0;
    float sumvsq = 0;

    for (int i = 0; i < nrparticles; i++) {
        particles[i].vel = particles[i].vel + dt*particles[i].acc;
        base_float_vec3 newpos = particles[i].pos + dt*particles[i].vel;
        if (newpos[0] < a*n)
            particles[i].pos[0] = newpos[0];
        else
            particles[i].pos[0] = newpos[0] - a*n;
        if (newpos[1] < a*n)
            particles[i].pos[1] = newpos[1];
        else
            particles[i].pos[1] = newpos[1] - a*n;
        if (newpos[2] < a*n)
            particles[i].pos[2] = newpos[2];
        else
            particles[i].pos[2] = newpos[2] - a*n;
        sumv = sumv + particles[i].vel;
        sumvsq = sumvsq + (particles[i].vel)*(particles[i].vel);
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


void system::force_calculation() { //using reduced unit
    vector<float> force_x ;
    vector<float> force_y ;
    vector<float> force_z ;
    float distance = inner_cutoff ;
    float distance_inv = 1/distance ;
    float distance6_inv = pow(distance_i,6) ;
    float E_c = 4 * distance6_inv * (distance6_inv - 1);
                    
    for (i=0, i<nrparticles , i++) { 
        for (j = verlet_particles_list[i] + 1, j < verlet_particles_list[i+1] , j++) { 
            base_float_vec3 dr = particles[i].pos-particles[j].pos ;
            distance = (dr).lenght ;
            distance_inv = 1/distance ;
            distance6_inv = pow(distance_i,6) ;
            float force = 48 * distance_inv * distance6_inv * (distance6_inv - 0.5) ;
            dr.normalize();
            fvec3 x_hat = fvec3(1, 0, 0);
            fvec3 y_hat = fvec3(0, 1, 0);
            fvec3 z_hat = fvec3(0, 0, 1);
            force_x[i]+ =    force * (dr *  x_hat) ;
            force_x[j]- =    force * (dr *  x_hat) ;
            force_y[i]+ =    force * (dr *  y_hat) ;
            force_y[j]- =    force * (dr *  y_hat) ;
            force_z[i]+ =    force * (dr *  z_hat) ;
            force_z[j]- =    force * (dr *  z_hat) ;
            Ep[i] += 4 * distance6_i * (distance6_i - 1) - E_c ; 
            Ep[j] += 4 * distance6_i * (distance6_i - 1) - E_c ;
            float sigma_r_dot_f += force * distance ; // for pressure
        }
    }
}

system::system(int nrparticles_in, int nrsteps_in, float sigma_in, float epsilon_in, float outer_cutoff_in, float inner_cutoff_in, float mass_in, float dt_in, int nrinst_in, float temperature_in, int nrtimesteps_in, float latticeconstant_in):
    cell_linklist(4*((int) pow(nrparticles_in/4, 1/3))*((int) pow(nrparticles_in/4, 1/3))*((int) pow(nrparticles_in/4, 1/3))),
    cell_list((floor(n*a/outer_cutoff))*(floor(n*a/outer_cutoff))*(floor(n*a/outer_cutoff))), 
    particles(4*((int) pow(nrparticles_in/4, 1/3))*((int) pow(nrparticles_in/4, 1/3))*((int) pow(nrparticles_in/4, 1/3))), 
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
    n = (int) pow(nrparticles_in/4, 1/3); 
    nrparticles = 4*n*n*n;
    mass = mass_in;
    sigma = sigma_in;
    epsilon = epsilon_in;
    nrinst = nrinst_in;
    init_temp = temperature_in;
    nrtimesteps = nrtimesteps_in;
    distanceforcesum = 0;
    kB = 1.381e-23;
    a = latticeconstant_in;
    nrcells = floor(n*a/outer_cutoff);
    cellsize = n*a/nrcells;
}

void system::md() {
    init();
    create_linked_cells();
    create_verlet_list_using_linked_cell_list();
    while (timestep <= nrtimesteps) {
        force_calculation();
        leapfrog();
        calculate_properties();
        if (neighbors should be updated) {
            create_linked_cells();
            create_verlet_list_using_linked_cell_list();
        }
        timestep++;
    }
}

void system::init() {
    3d_vector sumv = 0;
    float sumvsq = 0;
    initpos();
    srand(time(NULL));
    for (i = 0; i < nrparticles; i++) {
        for (j = 0; j < 3; j++) {
            particles[i].vel[j] = ((float) rand())/((float) RAND_MAX) - 0.5;
        }
        sumv += particles[i].vel;
        sumvsq += particles[i].vel*particles[i].vel;
    }
    sumv = sumv/nrparticles;
    sumvsq = sumvsq/nrparticles;
    float s = sqrt(3*init_temp/sumvsq);
    for (i = 0; i < nrparticles; i++) {
        particles[i].vel = (particles[i].vel - sumv)*s;
    }
}

void system::calculate_diffusion_coefficient() {                //not finished

}

void system::calculate_temperature() {
    float sum = 0;
    for (i = 0; i < nrinst; i++) {
        sum += insttemp[i];
    }
    temp[timestep/nrinst] = sum/nrinst; 
}

void system::calculate_Ep() {
    float sum = 0;
    for (i = 0; i < nrinst; i++) {
        sum += instEp[i];
    }
    Ep[timestep/nrinst] = sum/nrinst; 
}

void system::calculate_Ek() {
    float sum = 0;
    for (i = 0; i < nrinst; i++) {
        sum += instEk[i];
    }
    Ek[timestep/nrinst] = sum/nrinst; 
}

void system::initpos() {
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                (*particles[4*(i*n*n + j*n + k)]).fcc_pos[0] = i*a;
                (*particles[4*(i*n*n + j*n + k)]).fcc_pos[1] = j*a;
                (*particles[4*(i*n*n + j*n + k)]).fcc_pos[2] = k*a;
                (*particles[4*(i*n*n + j*n + k) + 1]).fcc_pos[0] = (i + 0.5)*a;
                (*particles[4*(i*n*n + j*n + k) + 1]).fcc_pos[1] = (j + 0.5)*a;
                (*particles[4*(i*n*n + j*n + k) + 1]).fcc_pos[2] = k*a;
                (*particles[4*(i*n*n + j*n + k) + 2]).fcc_pos[0] = (i + 0.5)*a;
                (*particles[4*(i*n*n + j*n + k) + 2]).fcc_pos[1] = j*a;
                (*particles[4*(i*n*n + j*n + k) + 2]).fcc_pos[2] = (k + 0.5)*a;
                (*particles[4*(i*n*n + j*n + k) + 3]).fcc_pos[0] = i*a;
                (*particles[4*(i*n*n + j*n + k) + 3]).fcc_pos[1] = (j + 0.5)*a;
                (*particles[4*(i*n*n + j*n + k) + 3]).fcc_pos[2] = (k + 0.5)*a;
            }
        }
    }
    for (i = 0; i < nrparticles; i++) {
        (*particles[i]).pos = (*particles[i]).fcc_pos;
    }
}


void system::calculate_properties() {
    if ((timestep % nrinst) == 0) {
        calculate_specific_heat();            
        calculate_pressure();
        calculate_mean_square_displacement();
        calculate_diffusion_coefficient();
        calculate_temperature();
        calculate_Ep();
        calculate_Ek();
    }
}

void system::calculate_specific_heat() {
    float T2 = 0;
    for (i = 0; i < nrinst; i++){
        T2 += insttemp[i];
    }
    T2 = T2/nrinst;
    Cv[timestep/nrinst] = 9*kB/(6/nrparticles+4-4*T2/temp[timestep/nrinst]);
}

void system::calculate_pressure() {
    float V = nrcells*nrcells*nrcells*outer_cutoff*outer_cutoff*outer_cutoff;
    pressure[timestep/nrinst] = nrparticles*kB*temp[timestep/nrinst]/V + distanceforcesum/(6*V*nrinst);
    distanceforcesum = 0;
}

void system::calculate_mean_square_displacement() {
    float sum = 0;
    for (i = 0; i < nrparticles;i++) {
        sum += ((*particles[i]).pos - (*particles[i]).fcc_pos)*((*particles[i]).pos - (*particles[i]).fcc_pos);
    }
    sum = sum/nrparticles;
    msd[timestep/nrinst] = sum;
}
