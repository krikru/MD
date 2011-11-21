
// Standard includes
#include <cstdlib>
#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
#include <fstream>
using std::ofstream;
// Own includes
#include "mdsystem.h"

////////////////////////////////////////////////////////////////
// CONSTRUCTOR
////////////////////////////////////////////////////////////////

mdsystem::mdsystem(uint nrparticles_in, ftype sigma_in, ftype epsilon_in, ftype inner_cutoff_in, ftype outer_cutoff_in, ftype mass_in, ftype dt_in, uint nrinst_in, ftype temperature_in, uint nrtimesteps_in, ftype latticeconstant_in, uint lattice_type_in, ftype desiredtemp_in, ftype thermostat_time_in, bool thermostat_on_in, bool diff_c_on_in, bool Cv_on_in, bool pressure_on_in, bool msd_on_in, bool Ep_on_in, bool Ek_on_in)
{
    lattice_type = lattice_type_in; // One of the supported lattice types listed in enum_lattice_types
    dt = dt_in;                     // Delta time, the time step to be taken when solving the diff.eq.
    sqr_outer_cutoff = outer_cutoff_in*outer_cutoff_in; // Parameter for the Verlet list
    sqr_inner_cutoff = inner_cutoff_in*inner_cutoff_in; // Parameter for the Verlet list

    loop_num = 0;
    nrtimesteps = ((nrtimesteps_in - 1) / nrinst_in + 1) * nrinst_in; // Make the smallest multiple of nrinst_in that has at least the specified size

    insttemp.resize(nrinst_in);
    instEk  .resize(nrinst_in);
    instEp  .resize(nrinst_in);
    temp                 .resize(nrtimesteps/nrinst_in + 1);
    therm                .resize(nrtimesteps/nrinst_in + 1);
    Ek                   .resize(nrtimesteps/nrinst_in + 1);
    Ep                   .resize(nrtimesteps/nrinst_in + 1);
    Cv                   .resize(nrtimesteps/nrinst_in + 1);
    pressure             .resize(nrtimesteps/nrinst_in + 1);
    msd                  .resize(nrtimesteps/nrinst_in + 1);
    diffusion_coefficient.resize(nrtimesteps/nrinst_in + 1);
    if (lattice_type == LT_FCC) {
        n = int(pow(ftype(nrparticles_in / 4 ), ftype( 1.0 / 3.0 )));
        nrparticles = 4*n*n*n;   // Calculate the new number of atoms; all can't fit in the box since n is an integer
    }
    mass = mass_in;
    sqr_sigma = sigma_in*sigma_in;
    four_epsilon = 4*epsilon_in;
    nrinst = nrinst_in;
    init_temp = temperature_in;
    distanceforcesum = 0;
    a = latticeconstant_in;
    box_size = a*n;
    p_half_box_size = 0.5f * box_size;
    n_half_box_size = -p_half_box_size;
    nrcells = int(n*a/outer_cutoff_in);
    if (nrcells > 3) {
        cells_used = true;
    }
    else {
        cells_used = false;
    }
    cellsize = n*a/nrcells;
    diff_c_on = diff_c_on_in;
    Cv_on = Cv_on_in;
    pressure_on = pressure_on_in;
    msd_on = msd_on_in;
    Ep_on = Ep_on_in;
    Ek_on = Ek_on_in;
    desiredtemp = desiredtemp_in;
    thermostat_time = thermostat_time_in;
    thermostat_on = thermostat_on_in;
}

////////////////////////////////////////////////////////////////
// PUBLIC FUNCTIONS
////////////////////////////////////////////////////////////////

void mdsystem::run_simulation(void (*event_handler_in)(void)) {
    // Set the event handler
    event_handler = event_handler_in;

    // Initialize the system
    init();

    // Start simulating
    for (loop_num = 0; loop_num <= nrtimesteps; loop_num++) {
        cout << "loop number = " << loop_num << endl; //TODO: Remove

        // Evolve the system in time
        force_calculation();
        leapfrog(); // TODO: Compensate for half time steps

        // Calculate properties each nrinst loops
        if (loop_num % nrinst == 0 && loop_num != 0) {
            calculate_properties();
        }

        // Update Verlet list if necessary
        calculate_largest_sqr_displacement();
        if (4 * largest_sqr_displacement > (sqr_outer_cutoff + sqr_inner_cutoff - 2*sqrt(sqr_outer_cutoff*sqr_inner_cutoff))) {
            cout<<int(100*loop_num/nrtimesteps)<<" % done"<<endl;
            create_verlet_list();
            cout<<int(100*loop_num/nrtimesteps)<<" % done"<<endl;
        }

        // Process events
        event_handler();
    }
    cout << "Simulation completed" << endl;

    // The out files are like cin
    ofstream out_etot_data ;
    ofstream out_temp_data ;
    ofstream out_therm_data;

    // Open the output files
    out_etot_data .open("TotalEnergy.dat");
    out_temp_data .open("Temperature.dat");
    out_therm_data.open("Thermostat.dat" );
    if( !out_etot_data || !out_temp_data ) { // file couldn't be opened
        cerr << "Error: Output files could not be opened" << endl;
    }
    else {
        for (uint i = 1; i < temp.size(); i++) {
            out_etot_data  << setprecision(9) << Ek   [i] + Ep[i] << endl;
            out_temp_data  << setprecision(9) << temp [i]         << endl;
            out_therm_data << setprecision(9) << therm[i]         << endl;

            // Process events
            event_handler();
        }
        out_etot_data .close();
        out_temp_data .close();
        out_therm_data.close();
    }

    for (uint i = 1; i < temp.size();i++)
    {
        cout << "Temp    = " <<setprecision(9) << temp[i]         << endl;
        cout << "Ek + Ep = " <<setprecision(9) << Ek  [i] + Ep[i] << endl;
        cout << "Ek      = " <<setprecision(9) << Ek  [i]         << endl;
        cout << "Ep      = " <<setprecision(9) << Ep  [i]         << endl;
        cout << "Cv      = " <<setprecision(9) << Cv  [i]         << endl;
        
        // Process events
        event_handler();
    }
}

////////////////////////////////////////////////////////////////
// PRIVATE FUNCTIONS
////////////////////////////////////////////////////////////////

void mdsystem::init() {
    init_particles();
    create_verlet_list();
}

void mdsystem::init_particles() {
    // Allocate space for particles
    particles.resize(nrparticles);

    //Place out particles according to the lattice pattern
    if (lattice_type == LT_FCC) {
        for (uint z = 0; z < n; z++) {
            for (uint y = 0; y < n; y++) {
                for (uint x = 0; x < n; x++) {
                    int help_index = 4*(x + n*(y + n*z));

                    (particles[help_index + 0]).start_pos[0] = x*a;
                    (particles[help_index + 0]).start_pos[1] = y*a;
                    (particles[help_index + 0]).start_pos[2] = z*a;

                    (particles[help_index + 1]).start_pos[0] = x*a;
                    (particles[help_index + 1]).start_pos[1] = (y + 0.5f)*a;
                    (particles[help_index + 1]).start_pos[2] = (z + 0.5f)*a;

                    (particles[help_index + 2]).start_pos[0] = (x + 0.5f)*a;
                    (particles[help_index + 2]).start_pos[1] = y*a;
                    (particles[help_index + 2]).start_pos[2] = (z + 0.5f)*a;

                    (particles[help_index + 3]).start_pos[0] = (x + 0.5f)*a;
                    (particles[help_index + 3]).start_pos[1] = (y + 0.5f)*a;
                    (particles[help_index + 3]).start_pos[2] = z*a;
                } // X
            } // Y
        } // Z
    }
    
    //Randomixe the velocities
    vec3 sum_vel = vec3(0, 0, 0);
    ftype sum_sqr_vel = 0;
    for (uint i = 0; i < nrparticles; i++) {
        for (uint j = 0; j < 3; j++) {
            particles[i].start_vel[j] = 0;
            for (uint terms = 0; terms < 5; terms++) { //This will effectivelly create a distribution very similar to normal distribution. (If you want to see what the distribution looks like, go to www.wolframalpha.com/input/?i=fourier((sinc(x))^n) and replace n by the number of terms)
                particles[i].start_vel[j] += ftype(rand());
            }
        }
        sum_vel    += particles[i].start_vel;
        sum_sqr_vel += particles[i].start_vel.sqr_length();
    }

    // Compensate for incorrect start temperature and total velocities and finalize the initialization values
    vec3 average_vel = sum_vel/ftype(nrparticles);
    ftype vel_variance = sum_sqr_vel/nrparticles - average_vel.sqr_length();
    ftype scale_factor = sqrt(1.5f * P_KB * init_temp / (0.5f * vel_variance * mass)); // Termal energy = 1.5 * P_KB * init_temp
    for (uint i = 0; i < nrparticles; i++) {
        particles[i].start_vel = (particles[i].start_vel - average_vel)*scale_factor;
        particles[i].vel = particles[i].start_vel;
        particles[i].pos = particles[i].start_pos;
        particles[i].non_modulated_relative_pos = vec3(0, 0, 0);
        particles[i].pos_when_non_modulated_relative_pos_was_calculated = particles[i].start_pos;
    }
}

void mdsystem::create_verlet_list()
{
    //Updating pos_when_verlet_list_created and non_modulated_relative_pos for all particles
    for (uint i = 0; i < nrparticles; i++) {
        update_single_non_modulated_particle_position(i);
        particles[i].pos_when_verlet_list_created = particles[i].pos;
    }

    if (cells_used) {
        create_linked_cells();
    }
    create_verlet_list_using_linked_cell_list();
}

void mdsystem::create_linked_cells() {//Assuming origo in the corner of the bulk, and positions given according to boundaryconditions i.e. between zero and lenght of the bulk.
    int cellindex = 0;
    cell_list.resize(nrcells*nrcells*nrcells);
    cell_linklist.resize(nrparticles);
    for (uint i = 0; i < cell_list.size() ; i++) {
        cell_list[i] = 0; // Beware! Particle zero is a member of all cells!
    }
    for (uint i = 0; i < nrparticles; i++) {
        uint help_x = int(particles[i].pos[0] / cellsize);
        uint help_y = int(particles[i].pos[1] / cellsize);
        uint help_z = int(particles[i].pos[2] / cellsize);
        if (help_x == nrcells || help_y == nrcells || help_z == nrcells) { // This actually occationally happens
            help_x -= help_x == nrcells;
            help_y -= help_y == nrcells;
            help_z -= help_z == nrcells;
        }
        cellindex = help_x + nrcells * (help_y + nrcells * help_z);
        cell_linklist[i] = cell_list[cellindex];
        cell_list[cellindex] = i;
    }
}

void mdsystem::create_verlet_list_using_linked_cell_list() { // This function ctreates the verlet_lists (verlet_vectors) using the linked cell lists
    uint cellindex = 0;
    uint neighbour_particle_index = 0;
    verlet_particles_list.resize(nrparticles);
    verlet_neighbors_list.resize(nrparticles * 100); //This might be unnecessarily large //TODO: CHANGE THIS AS SOON AS POSSIBLE!!!

    //Creating new verlet_list
    verlet_particles_list[0] = 0;
    for (uint i = 0; i < nrparticles;) { // Loop through all particles
        // Init this neighbour list and point to the next list
        verlet_neighbors_list[verlet_particles_list[i]] = 0; // Reset number of neighbours
        int next_particle_list = verlet_particles_list[i] + 1; // Link to the next particle list

        if (cells_used) { //Loop through all neighbour cells
            // Calculate cell indexes
            uint cellindex_x = int(particles[i].pos[0]/cellsize);
            uint cellindex_y = int(particles[i].pos[1]/cellsize);
            uint cellindex_z = int(particles[i].pos[2]/cellsize);
            if (cellindex_x == nrcells || cellindex_y == nrcells || cellindex_z == nrcells) { // This actually occationally happens
                cellindex_x -= cellindex_x == nrcells;
                cellindex_y -= cellindex_y == nrcells;
                cellindex_z -= cellindex_z == nrcells;
            }
            for (int index_z = int(cellindex_z) - 1; index_z <= int(cellindex_z) + 1; index_z++) {
                for (int index_y = int(cellindex_y) - 1; index_y <= int(cellindex_y) + 1; index_y++) {
                    for (int index_x = int(cellindex_x) - 1; index_x <= int(cellindex_x) + 1; index_x++) {
                        int modulated_x = index_x;
                        int modulated_y = index_y;
                        int modulated_z = index_z;
                        // Control boundaries
                        if (modulated_x == -1) {
                            modulated_x = int(nrcells) - 1;
                        }
                        else if (modulated_x == int(nrcells)) {
                            modulated_x = 0;
                        }
                        if (modulated_y == -1) {
                            modulated_y = int(nrcells) - 1;
                        }
                        else if (modulated_y == int(nrcells)) {
                            modulated_y = 0;
                        }
                        if (modulated_z == -1) {
                            modulated_z = int(nrcells) - 1;
                        }
                        else if (modulated_z == int(nrcells)) {
                            modulated_z = 0;
                        }
                        cellindex = uint(modulated_x + nrcells * (modulated_y + nrcells * modulated_z)); // Calculate neighbouring cell index
                        neighbour_particle_index = cell_list[cellindex]; // Get the largest particle index of the particles in this cell
                        while (neighbour_particle_index > i) { // Loop though all particles in the cell with greater index
                            ftype sqr_distance = modulos_distance(particles[neighbour_particle_index].pos, particles[i].pos).sqr_length();
                            if(sqr_distance < sqr_outer_cutoff) {
                                verlet_neighbors_list[verlet_particles_list[i]] += 1;
                                verlet_neighbors_list[next_particle_list] = neighbour_particle_index;
                                next_particle_list++;
                            }
                            neighbour_particle_index = cell_linklist[neighbour_particle_index]; // Get the next particle in the cell
                        }
                    } // X
                } // Y
            } // Z
        } // if (cells_used)
        else {
            for (neighbour_particle_index = i+1; neighbour_particle_index < nrparticles; neighbour_particle_index++) { // Loop though all particles with greater index
                ftype sqr_distance = modulos_distance(particles[neighbour_particle_index].pos, particles[i].pos).sqr_length();
                if(sqr_distance < sqr_outer_cutoff) {
                    verlet_neighbors_list[verlet_particles_list[i]] += 1;
                    verlet_neighbors_list[next_particle_list] = neighbour_particle_index;
                    next_particle_list++;
                }
            }
        }
        i++; // Continue with the next particle (if there exists any)
        if (i < nrparticles) { // Point to the next particle list
            verlet_particles_list[i] = next_particle_list;
        }
    }
}

//Updating largest square displacement.
void mdsystem::calculate_largest_sqr_displacement()
{
    largest_sqr_displacement = 0;
    for (uint i = 0; i < nrparticles; i++) {
        ftype sqr_displacement = modulos_distance(particles[i].pos, particles[i].pos_when_verlet_list_created).sqr_length();
        if (sqr_displacement > largest_sqr_displacement) {
            largest_sqr_displacement = sqr_displacement;
        }
    }
}

void mdsystem::update_non_modulated_particle_positions()
{
    for (uint i = 0; i < nrparticles; i++) {
        update_single_non_modulated_particle_position(i);
    }
}

inline void mdsystem::update_single_non_modulated_particle_position(uint i)
{
    particles[i].non_modulated_relative_pos += modulos_distance(particles[i].pos_when_non_modulated_relative_pos_was_calculated, particles[i].pos);
    particles[i].pos_when_non_modulated_relative_pos_was_calculated = particles[i].pos;
}

void mdsystem::leapfrog()
{
    ftype sum_sqr_vel = 0;
    
    //TODO: What if the temperature (insttemp) is zero? Has to randomize new velocities in that case.
#if THERMOSTAT == LASSES_THERMOSTAT
	thermostat = loop_num > 0 ? (1 - desiredtemp/insttemp[(loop_num-1) % nrinst]) / (2*thermostat_time) : 0;
#elif THERMOSTAT == CHING_CHIS_THERMOSTAT
	/////Using Smooth scaling Thermostat (Berendsen et. al, 1984)/////
	thermostat = loop_num > 0 ? sqrt(1 +  dt / thermostat_time * ((desiredtemp) / insttemp[(loop_num-1) % nrinst] - 1)) : 1;
#endif

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
#if THERMOSTAT == CHING_CHIS_THERMOSTAT
        particles[i].vel = particles[i].vel * thermostat;
#endif
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
    insttemp[loop_num % nrinst] = mass * sum_sqr_vel / (3 * nrparticles * P_KB);
    
    if (Ek_on) instEk[loop_num % nrinst] = 0.5f * mass * sum_sqr_vel;
}

void mdsystem::force_calculation() { //Using si-units
    // Reset accelrations for all particles
    for (uint k = 0; k < nrparticles; k++) {
        particles[k].acc = vec3(0, 0, 0);
    }
    ftype sqr_distance;
    ftype sqr_distance_inv;
    ftype distance_inv;
    ftype p = sqr_sigma / sqr_inner_cutoff; // For calculating the cutoff energy
    p = p*p*p;
    ftype E_cutoff = four_epsilon * p * (p - 1);
    ftype mass_inv = 1/mass;             
    instEp[loop_num % nrinst] = 0;
    for (uint i1 = 0; i1 < nrparticles ; i1++) { // Loop through all particles
        for (uint j = verlet_particles_list[i1] + 1; j < verlet_particles_list[i1] + verlet_neighbors_list[verlet_particles_list[i1]] + 1 ; j++) { 
            uint i2 = verlet_neighbors_list[j]; // Get index of the second (possibly) interacting particle 
            vec3 r = modulos_distance(particles[i2].pos, particles[i1].pos); // Calculate the closest distance
            sqr_distance = r.sqr_length();
            if (sqr_distance >= sqr_inner_cutoff) {
                continue; // Skip this interaction and continue with the next one
            }
            sqr_distance_inv = 1/sqr_distance;

            //Calculating acceleration
            distance_inv = sqrt(sqr_distance_inv);
            p = sqr_sigma * sqr_distance_inv;
            p = p*p*p;
            ftype acceleration = 12 * four_epsilon * distance_inv * p * (p - 0.5f) * mass_inv;

            // Update accelerations of interacting particles
            vec3 r_hat = r * distance_inv;
            particles[i1].acc +=  acceleration * r_hat;
            particles[i2].acc -=  acceleration * r_hat;

            // Update properties
            //TODO: Remove these two from force calculation and place them somewhere else
            if (Ep_on) instEp[loop_num % nrinst] += four_epsilon * p * (p - 1) - E_cutoff;
            if (pressure_on) distanceforcesum += mass * acceleration / distance_inv;
        }
    }
#if THERMOSTAT == LASSES_THERMOSTAT
    if (thermostat_on) {
        for (uint i = 0; i < nrparticles; i++) {
            particles[i].acc -= thermostat * particles[i].vel;
        }
    }
#endif
}

void mdsystem::calculate_properties() {
    update_non_modulated_particle_positions();
    calculate_temperature();
    if (Cv_on) calculate_specific_heat();            
    if (pressure_on) calculate_pressure();
    if (msd_on) calculate_mean_square_displacement();
        
    if (Ep_on) calculate_Ep();
    if (Ek_on) calculate_Ek();
    if (diff_c_on) calculate_diffusion_coefficient();
}

void mdsystem::calculate_temperature() {
    ftype sum = 0;
    for (uint i = 0; i < nrinst; i++) {
        sum += insttemp[i];
    }
    temp[loop_num/nrinst] = sum/nrinst;
    therm[loop_num/nrinst] = thermostat;
}

void mdsystem::calculate_Ep() {
    ftype sum = 0;
    for (uint i = 0; i < nrinst; i++) {
        sum += instEp[i];
    }
    Ep[loop_num/nrinst] = sum/nrinst; 
}

void mdsystem::calculate_Ek() {
    ftype sum = 0;
    for (uint i = 0; i < nrinst; i++) {
        sum += instEk[i];
    }
    Ek[loop_num/nrinst] = sum/nrinst;
}

void mdsystem::calculate_specific_heat() {
    ftype T2 = 0;
    for (uint i = 0; i < nrinst; i++){
        T2 += insttemp[i]*insttemp[i];
    }
    T2 = T2/nrinst;
    Cv[loop_num/nrinst] = 9*P_KB/(6/nrparticles+4-4*T2/(temp[loop_num/nrinst]*temp[loop_num/nrinst]));
}

void mdsystem::calculate_pressure() {
    ftype V = n*a*n*a*n*a;
    pressure[loop_num/nrinst] = nrparticles*P_KB*temp[loop_num/nrinst]/V + distanceforcesum/(6*V*nrinst);
    distanceforcesum = 0;
}

void mdsystem::calculate_mean_square_displacement() {
    ftype sum = 0;
    for (uint i = 0; i < nrparticles;i++) {
        sum += (particles[i].non_modulated_relative_pos - particles[i].start_pos).sqr_length();
    }
    sum = sum/nrparticles;
    msd[loop_num/nrinst] = sum;
}

void mdsystem::calculate_diffusion_coefficient() {
    diffusion_coefficient[loop_num/nrinst] = msd[loop_num/nrinst]/(6*dt*loop_num);
}

vec3 mdsystem::modulos_distance(vec3 pos1, vec3 pos2) const
{
    vec3 d = pos2 - pos1;

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
