////////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////////

// Standard includes
#include <stdexcept>
using std::runtime_error;
#include <cstdlib>
#include <iostream>
#include <iomanip>
using std::endl;
#include <fstream>
using std::ofstream;

// Own includes
#include "mdsystem.h"

////////////////////////////////////////////////////////////////
// CONSTRUCTOR
////////////////////////////////////////////////////////////////

mdsystem::mdsystem()
{
    operating = false;
}

////////////////////////////////////////////////////////////////
// PUBLIC FUNCTIONS
////////////////////////////////////////////////////////////////

void mdsystem::set_event_callback(callback<void (*)(void*)> event_callback_in)
{
    start_operation();
    event_callback = event_callback_in;
    finish_operation();
}

void mdsystem::set_output_callback(callback<void (*)(void*, string)> output_callback_in)
{
    start_operation();
    output_callback = output_callback_in;
    finish_operation();
}

void mdsystem::init(uint nrparticles_in, ftype sigma_in, ftype epsilon_in, ftype inner_cutoff_in, ftype outer_cutoff_in, ftype mass_in, ftype dt_in, uint nrinst_in, ftype temperature_in, uint nrtimesteps_in, ftype latticeconstant_in, uint lattice_type_in, ftype desiredtemp_in, ftype thermostat_time_in, ftype deltaEp_in, bool thermostat_on_in, bool diff_c_on_in, bool Cv_on_in, bool pressure_on_in, bool msd_on_in, bool Ep_on_in, bool Ek_on_in)
{
#if RU_ON == 1
        // The system is *always* operating when running non-const functions
        start_operation();

        mass = mass_in;
        sqr_sigma = sigma_in * sigma_in;
        four_epsilon = 4 * epsilon_in;
        sigma = sigma_in;
        epsilon = epsilon_in;

        //reduced unit
        dt = dt_in / sqrt(mass * sqr_sigma / epsilon);
        init_temp = temperature_in * P_KB/ epsilon;
        desiredtemp = desiredtemp_in * P_KB/ epsilon;
        thermostat_time = thermostat_time_in / sqrt(mass * sqr_sigma / epsilon);
        a = latticeconstant_in / sigma;
        inner_cutoff = inner_cutoff_in / sigma;
        outer_cutoff = outer_cutoff_in / sigma;
        sqr_outer_cutoff =outer_cutoff*outer_cutoff ; // Parameter for the Verlet list
        sqr_inner_cutoff =inner_cutoff*inner_cutoff ; // Parameter for the Verlet list

        //
        lattice_type = lattice_type_in; // One of the supported lattice types listed in enum_lattice_types

        loop_num = 0;
        nrtimesteps = ((nrtimesteps_in - 1) / nrinst_in + 1) * nrinst_in; // Make the smallest multiple of nrinst_in that has at least the specified size
        insttemp.resize(nrinst_in);
        instEk  .resize(nrinst_in);
        instEp  .resize(nrinst_in);
        temp                 .resize(nrtimesteps/nrinst_in + 1);
        therm                .resize(nrtimesteps/nrinst_in + 1);
        Ek                   .resize(nrtimesteps/nrinst_in + 1);
        Ep                   .resize(nrtimesteps/nrinst_in + 1);
        cohesive_energy      .resize(nrtimesteps/nrinst_in + 1);
        Cv                   .resize(nrtimesteps/nrinst_in + 1);
        pressure             .resize(nrtimesteps/nrinst_in + 1);
        msd                  .resize(nrtimesteps/nrinst_in + 1);
        diffusion_coefficient.resize(nrtimesteps/nrinst_in + 1);
        if (lattice_type == LT_FCC) {
            n = int(pow(ftype(nrparticles_in / 4 ), ftype( 1.0 / 3.0 )));
            nrparticles = 4*n*n*n;   // Calculate the new number of atoms; all can't fit in the box since n is an integer
        }

        nrinst = nrinst_in;
        distanceforcesum = 0;

        box_size = a*n;

        p_half_box_size = 0.5f * box_size;
        n_half_box_size = -p_half_box_size;

        nrcells = int(box_size/outer_cutoff);
        if (nrcells > 3) {
            cells_used = true;
        }
        else {
            cells_used = false;
        }
        cellsize = box_size/nrcells;

        // turned off some functions for the time being
        deltaEp = deltaEp_in;
        diff_c_on = diff_c_on_in;
        Cv_on = Cv_on_in;
        pressure_on = pressure_on_in;
        msd_on = msd_on_in;
        Ep_on = Ep_on_in;
        Ek_on = Ek_on_in;
        desiredtemp = desiredtemp_in;
        thermostat_time = thermostat_time_in;
        thermostat_on = 0;
        equilibrium = false;

        abort_activities_requested = false;

        init_particles();
        create_verlet_list();

        // Finish the operation
        finish_operation();

#else
        // The system is *always* operating when running non-const functions
        start_operation();

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
        cohesive_energy      .resize(nrtimesteps/nrinst_in + 1);
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
        deltaEp = deltaEp_in;
        diff_c_on = diff_c_on_in;
        Cv_on = Cv_on_in;
        pressure_on = pressure_on_in;
        msd_on = msd_on_in;
        Ep_on = Ep_on_in;
        Ek_on = Ek_on_in;
        desiredtemp = desiredtemp_in;
        thermostat_time = thermostat_time_in;
        thermostat_on = thermostat_on_in;
        equilibrium = false;

        abort_activities_requested = false;

        init_particles();
        create_verlet_list();

        // Finish the operation
        finish_operation();

#endif
}

void mdsystem::run_simulation()
{
    // The system is *always* operating when running non-const functions
    start_operation();

    // Open the output files. They work like cin
    ofstream out_etot_data ;
    ofstream out_ep_data ;
    ofstream out_ek_data;
    ofstream out_cv_data;
    ofstream out_temp_data ;
    ofstream out_therm_data;
    ofstream out_msd_data  ;
    ofstream out_cohe_data ;
    ofstream out_posx      ;
    ofstream out_posy      ;
    ofstream out_posz      ;

    vector<float> posx(nrtimesteps + 1); // TODO: This code shouldn't be here
    vector<float> posy(nrtimesteps + 1); // TODO: This code shouldn't be here
    vector<float> posz(nrtimesteps + 1); // TODO: This code shouldn't be here

    // Start simulating
    for (loop_num = 0; loop_num <= nrtimesteps; loop_num++) {

        // Check if the simulation has been requested to abort
        if (abort_activities_requested) {
            goto operation_finished;
        }
        //output <<"loop number = " << loop_num << endl;

        // Evolve the system in time
        //cout << loop_num << endl;
        force_calculation();
        leapfrog(); // TODO: Compensate for half time steps
        /*
        posx[loop_num]=particles[nrparticles/2].pos[0]/a;
        posy[loop_num]=particles[nrparticles/2].pos[1]/a;
        posz[loop_num]=particles[nrparticles/2].pos[2]/a;
        */
        // Calculate properties each nrinst loops
        if (loop_num % nrinst == 0 && loop_num != 0) {
            calculate_properties();
        }

        // Update Verlet list if necessary
        update_verlet_list_if_necessary();

        // Process events
        print_output_and_process_events();
    }
    output << "Simulation completed." << endl;

    /*
     * TODO: The following code should be moved into another public function.
     * This function should *just* run the simulation since that is what it
     * says it does.
     */
    output << "Opening output files..." << endl;
    if (!(open_ofstream_file(out_etot_data , "TotalEnergy.dat") &&
          open_ofstream_file(out_ep_data   , "Potential.dat") &&
          open_ofstream_file(out_ek_data   , "Kinetic.dat") &&
          open_ofstream_file(out_cv_data   , "Cv.dat") &&
          open_ofstream_file(out_temp_data , "Temperature.dat") &&
          open_ofstream_file(out_therm_data, "Thermostat.dat" ) &&
          open_ofstream_file(out_msd_data  , "MSD.dat"        ) &&
          open_ofstream_file(out_cohe_data , "cohesive.dat"       ) &&
          open_ofstream_file(out_posx      , "posx.dat"       ) &&
          open_ofstream_file(out_posy      , "posy.dat"       ) &&
          open_ofstream_file(out_posz      , "posz.dat"       )
          )) {
        cerr << "Error: Output files could not be opened" << endl;
    }
    else {
        output << "Writing to output files..." << endl;
        print_output_and_process_events();
        /*for (uint i = 1; i < posx.size(); i++)
        {
            if (abort_activities_requested) {
                break;
            }
            out_posx << setprecision(9) << posx[i] << endl;
            out_posy << setprecision(9) << posy[i] << endl;
            out_posz << setprecision(9) << posz[i] << endl;
        }*/
        for (uint i = 1; i < temp.size(); i++) {
            if (abort_activities_requested) {
                break;
            }
            out_etot_data  << setprecision(9) << Ep[i] + Ek[i]    << endl;
            out_ep_data    << setprecision(9) << Ep[i]            << endl;
            out_ek_data    << setprecision(9) << Ek[i]            << endl;
            out_cv_data    << setprecision(9) << Cv[i]            << endl;
            out_temp_data  << setprecision(9) << temp [i]         << endl;
            out_therm_data << setprecision(9) << therm[i]         << endl;
            out_msd_data   << setprecision(9) << msd  [i]         << endl;
            out_cohe_data  << setprecision(9) << cohesive_energy  [i]<< endl;

            // Process events
            print_output_and_process_events();
        }
        out_etot_data .close();
        out_ep_data   .close();
        out_ek_data   .close();
        out_cv_data   .close();
        out_temp_data .close();
        out_therm_data.close();
        out_msd_data  .close();
        out_cohe_data  .close();
        out_posx.close();
        out_posy.close();
        out_posz.close();
    }
    output << "Writing to output files done." << endl;

    for (uint i = 1; i < temp.size();i++)
    {
        if (abort_activities_requested) {
            goto operation_finished;
        }
        output << "Temp            = " <<setprecision(9) << temp[i]               << endl;
        output << "Ek + Ep         = " <<setprecision(9) << Ek  [i] + Ep[i]       << endl;
        output << "Ek              = " <<setprecision(9) << Ek  [i]               << endl;
        output << "Ep              = " <<setprecision(9) << Ep  [i]               << endl;
        output << "Cohesive energy = " <<setprecision(9) << (cohesive_energy [i])/P_EV   << endl;
        output << "Cv              = " <<setprecision(9) << Cv  [i]               << endl;
        output << "msd             = " <<setprecision(9) << msd [i]               << endl;

        // Process events
        print_output_and_process_events();
    }
    output<< "a=" << a<<endl;
    output<< "boxsize=" << box_size<<endl;
    output<<"dt="<< dt << endl;
    output<<"init_temp= "<<init_temp<<endl;
    output << "Complete" << endl;

operation_finished:
    // Finish the operation
    finish_operation();
}

void mdsystem::abort_activities()
{
    abort_activities_requested = true;
}

bool mdsystem::is_operating() const
{
    return operating;
}

uint mdsystem::get_loop_num() const
{
    return loop_num;
}

uint mdsystem::get_max_loops_num() const
{
    return nrtimesteps;
}

////////////////////////////////////////////////////////////////
// PRIVATE FUNCTIONS
////////////////////////////////////////////////////////////////

void mdsystem::init_particles() {
    // Allocate space for particles
    particles.resize(nrparticles);

    //Place out particles according to the lattice pattern
    if (lattice_type == LT_FCC) {
        for (uint z = 0; z < n; z++) {
            for (uint y = 0; y < n; y++) {
                for (uint x = 0; x < n; x++) {
                    int help_index = 4*(x + n*(y + n*z));

                    (particles[help_index + 0]).pos[0] = x*a;
                    (particles[help_index + 0]).pos[1] = y*a;
                    (particles[help_index + 0]).pos[2] = z*a;

                    (particles[help_index + 1]).pos[0] = x*a;
                    (particles[help_index + 1]).pos[1] = (y + 0.5f)*a;
                    (particles[help_index + 1]).pos[2] = (z + 0.5f)*a;

                    (particles[help_index + 2]).pos[0] = (x + 0.5f)*a;
                    (particles[help_index + 2]).pos[1] = y*a;
                    (particles[help_index + 2]).pos[2] = (z + 0.5f)*a;

                    (particles[help_index + 3]).pos[0] = (x + 0.5f)*a;
                    (particles[help_index + 3]).pos[1] = (y + 0.5f)*a;
                    (particles[help_index + 3]).pos[2] = z*a;
                } // X
            } // Y
        } // Z
    }
    
    //Randomize the velocities
    vec3 sum_vel = vec3(0, 0, 0);
    ftype sum_sqr_vel = 0;
    for (uint i = 0; i < nrparticles; i++) {
        for (uint j = 0; j < 3; j++) {
            particles[i].vel[j] = 0;
            for (uint terms = 0; terms < 5; terms++) { //This will effectivelly create a distribution very similar to normal distribution. (If you want to see what the distribution looks like, go to www.wolframalpha.com/input/?i=fourier((sinc(x))^n) and replace n by the number of terms)
                particles[i].vel[j] += ftype(rand());
            }
        }
        sum_vel    += particles[i].vel;
        sum_sqr_vel += particles[i].vel.sqr_length();
    }

    // Compensate for incorrect start temperature and total velocities and finalize the initialization values
    vec3 average_vel = sum_vel/ftype(nrparticles);
    ftype vel_variance = sum_sqr_vel/nrparticles - average_vel.sqr_length();
    ftype scale_factor;
#if RU_ON == 1
        scale_factor = sqrt(3.0f  * init_temp  / (vel_variance)); // Termal energy = 1.5 * P_KB * init_temp = 0.5 m v*v
#else
        scale_factor = sqrt(3.0f * P_KB * init_temp / (vel_variance * mass)); // Termal energy = 1.5 * P_KB * init_temp = 0.5 m v*v
#endif
    for (uint i = 0; i < nrparticles; i++) {
        particles[i].vel = (particles[i].vel - average_vel)* scale_factor;
    }

    reset_non_modulated_relative_particle_positions();
}

void mdsystem::update_verlet_list_if_necessary()
{
    // Check if largest displacement too large for not updating the Verlet list
    ftype sqr_limit = (sqr_outer_cutoff + sqr_inner_cutoff - 2*sqrt(sqr_outer_cutoff*sqr_inner_cutoff));
    uint i;
    for (i = 0; i < nrparticles; i++) {
        ftype sqr_displacement = modulus_position_minus(particles[i].pos, particles[i].pos_when_verlet_list_created).sqr_length();
        if (sqr_displacement > sqr_limit) {
            break;
        }
    }
    if (i < nrparticles) {
        // Displacement that is to large was found
        output << int(100*loop_num/nrtimesteps) << " % done" <<endl;
        create_verlet_list();
    }
}

void mdsystem::create_verlet_list()
{
    //Updating pos_when_verlet_list_created and non_modulated_relative_pos for all particles
    for (uint i = 0; i < nrparticles; i++) {
        update_single_non_modulated_relative_particle_position(i);
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
                            // TODO: The modolus can be removed if
                            ftype sqr_distance = modulus_position_minus(particles[i].pos, particles[neighbour_particle_index].pos).sqr_length();
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
                ftype sqr_distance = modulus_position_minus(particles[i].pos, particles[neighbour_particle_index].pos).sqr_length();
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

void mdsystem::reset_non_modulated_relative_particle_positions()
{
    for (uint i = 0; i < nrparticles; i++) {
        reset_single_non_modulated_relative_particle_positions(i);
    }
}

inline void mdsystem::reset_single_non_modulated_relative_particle_positions(uint i)
{
    particles[i].non_modulated_relative_pos = vec3(0, 0, 0);
    particles[i].pos_when_non_modulated_relative_pos_was_calculated = particles[i].pos;
}

void mdsystem::update_non_modulated_relative_particle_positions()
{
    for (uint i = 0; i < nrparticles; i++) {
        update_single_non_modulated_relative_particle_position(i);
    }
}

inline void mdsystem::update_single_non_modulated_relative_particle_position(uint i)
{
    particles[i].non_modulated_relative_pos += modulus_position_minus(particles[i].pos, particles[i].pos_when_non_modulated_relative_pos_was_calculated);
    particles[i].pos_when_non_modulated_relative_pos_was_calculated = particles[i].pos;
}

void mdsystem::leapfrog()
{
    ftype sum_sqr_vel = 0;
    
    //TODO: What if the temperature (insttemp) is zero? Has to randomize new velocities in that case.
#if THERMOSTAT == LASSES_THERMOSTAT
    thermostat = thermostat_on && loop_num > 0 ? (1 - desiredtemp/insttemp[(loop_num-1) % nrinst]) / (2*thermostat_time) : 0;
#elif THERMOSTAT == CHING_CHIS_THERMOSTAT
    /////Using Smooth scaling Thermostat (Berendsen et. al, 1984)/////
    thermostat = thermostat_on && loop_num > 0 ? sqrt(1 +  dt / thermostat_time * ((desiredtemp) / insttemp[(loop_num-1) % nrinst] - 1)) : 1;
#endif

    for (uint i = 0; i < nrparticles; i++) {
        //output << "\ti = " << i << endl;
        if (loop_num == 2) {
            //output << "i = " << i << endl;
            if (i == 22) {
                i = i; //TODO
            }
        }

        //TODO: Check if vel and pos are stored for the same time or not, in that case, compensate for that

        // Update velocities
#if THERMOSTAT == CHING_CHIS_THERMOSTAT
        if (thermostat_on) {
            particles[i].vel = particles[i].vel * thermostat;
        }
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

        // TODO: Remove this from leapfrog and place it somewhere else
        sum_sqr_vel = sum_sqr_vel + particles[i].vel.sqr_length();



    }
#if RU_ON ==1
    insttemp[loop_num % nrinst] =  sum_sqr_vel / (3 * nrparticles );
    //cout <<"insttemp= "<<insttemp[loop_num % nrinst]<<endl;
    if (Ek_on) instEk[loop_num % nrinst] = 0.5f * sum_sqr_vel;
#else
    insttemp[loop_num % nrinst] = mass * sum_sqr_vel / (3 * nrparticles * P_KB);
    if (Ek_on) instEk[loop_num % nrinst] = 0.5f * mass * sum_sqr_vel;
#endif

}

void mdsystem::force_calculation() { //Using si-units
    // Reset accelrations for all particles
    for (uint k = 0; k < nrparticles; k++) {
        particles[k].acc = vec3(0, 0, 0);
    }
    ftype sqr_distance;
    ftype sqr_distance_inv;
    ftype distance_inv;
    ftype p;
#if RU_ON == 1
        p = 1/sqr_inner_cutoff;
        p = p * p * p;
        ftype E_cutoff = 4 * p * (p - 1);
#else
        p = sqr_sigma / sqr_inner_cutoff; // For calculating the cutoff energy
        ftype mass_inv = 1/mass;
        p = p * p * p;
        ftype E_cutoff = four_epsilon * p * (p - 1);
#endif

    instEp[loop_num % nrinst] = 0;
    for (uint i1 = 0; i1 < nrparticles ; i1++) { // Loop through all particles
        for (uint j = verlet_particles_list[i1] + 1; j < verlet_particles_list[i1] + verlet_neighbors_list[verlet_particles_list[i1]] + 1 ; j++) { 
            // TODO: automatically detect if a boundary is crossed and compensate for that in this function
            // Calculate the closest distance to the second (possibly) interacting particle
            uint i2 = verlet_neighbors_list[j];
            vec3 r = modulus_position_minus(particles[i1].pos, particles[i2].pos);
            sqr_distance = r.sqr_length();
            if (sqr_distance >= sqr_inner_cutoff) {
                continue; // Skip this interaction and continue with the next one
            }
            sqr_distance_inv = 1/sqr_distance;

            //Calculating acceleration
            distance_inv = sqrt(sqr_distance_inv);
            ftype acceleration;
#if RU_ON == 1
            p = sqr_distance_inv;
            p = p*p*p;
            acceleration = 48  * distance_inv * p * (p - 0.5f);
#else
            p = sqr_sigma * sqr_distance_inv;; // For calculating the cutoff energy
            p = p*p*p;
            acceleration = 12 * four_epsilon * distance_inv * p * (p - 0.5f) * mass_inv;
#endif


            // Update accelerations of interacting particles
            vec3 r_hat = r * distance_inv;
            particles[i1].acc +=  acceleration * r_hat;
            particles[i2].acc -=  acceleration * r_hat;

            // Update properties
            //TODO: Remove these two from force calculation and place them somewhere else
#if RU_ON ==1
            if (Ep_on) instEp[loop_num % nrinst] += 4 * p * (p - 1) - E_cutoff;
            if (pressure_on) distanceforcesum += acceleration / distance_inv;
#else
            if (Ep_on) instEp[loop_num % nrinst] += four_epsilon * p * (p - 1) - E_cutoff;
            if (pressure_on) distanceforcesum += mass * acceleration / distance_inv;
#endif
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
    update_non_modulated_relative_particle_positions();
    calculate_temperature();
    if (Cv_on) calculate_specific_heat();            
    if (pressure_on) calculate_pressure();
    if (Ep_on) {
        calculate_Ep();
        calculate_cohesive_energy();
    }
    if (Ek_on) calculate_Ek();
    if (msd_on) calculate_mean_square_displacement();
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

void mdsystem::calculate_cohesive_energy() {
    cohesive_energy[loop_num/nrinst] = -Ep[loop_num/nrinst]/nrparticles;
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
#if RU_ON == 1
    ftype sqr_avgT = temp[loop_num/nrinst]*temp[loop_num/nrinst];
    //cout<< "<T>2=" << sqr_avgT << endl;
    //cout<< "<T2>=" << T2 << endl;
    ftype Cv_inv = (2/3.0/nrparticles - 4/9*((T2/sqr_avgT)-1)) ;
    //cout<< "Cv_inv=" << Cv_inv << endl;
    Cv[loop_num/nrinst] = 1/Cv_inv;
#else
    Cv[loop_num/nrinst] = P_KB/(ftype(2)/3 + nrparticles*(1 - T2/(temp[loop_num/nrinst]*temp[loop_num/nrinst]))) * P_AVOGADRO;
    //Cv[loop_num/nrinst] = 9*P_KB/(6.0f/nrparticles+4.0f-4*T2/(temp[loop_num/nrinst]*temp[loop_num/nrinst])) * P_AVOGADRO;
#endif

}

void mdsystem::calculate_pressure() {
    ftype V = n*a*n*a*n*a;
#if RU_ON == 1
    pressure[loop_num/nrinst] = nrparticles*temp[loop_num/nrinst]/V + distanceforcesum/(6*V*nrinst);
#else
    pressure[loop_num/nrinst] = nrparticles*P_KB*temp[loop_num/nrinst]/V + distanceforcesum/(6*V*nrinst);
#endif
    distanceforcesum = 0;
}

void mdsystem::calculate_mean_square_displacement() {
    ftype sum = 0;
    if (equilibrium == false && loop_num/nrinst) {
        // Check if equilibrium has been reached
        ftype variation = (Ep[loop_num/nrinst] - Ep[loop_num/nrinst - 1]) / Ep[loop_num/nrinst];
        variation = variation >= 0 ? variation : -variation;
        if (variation < deltaEp) {
            // The requirements for equilibrium has been reached
            equilibrium = true;
            // Consider the particles to "start" now
            reset_non_modulated_relative_particle_positions();
        }
    }
    if (equilibrium) {
        // Calculate mean square displacement
        for (uint i = 0; i < nrparticles;i++) {
            sum += particles[i].non_modulated_relative_pos.sqr_length();
        }
        sum = sum/nrparticles;
        msd[loop_num/nrinst] = sum;
    }
    else {
        // Equilibrium not reached; don't calculate this property.
        msd[loop_num/nrinst] = 0;
    }
}

void mdsystem::calculate_diffusion_coefficient()
{
    diffusion_coefficient[loop_num/nrinst] = msd[loop_num/nrinst]/(6*dt*loop_num);
}

ofstream* mdsystem::open_ofstream_file(ofstream &o, const char* path) const
{
    o.open(path);
    return &o;
}

vec3 mdsystem::modulus_position_minus(vec3 pos1, vec3 pos2) const
{
    vec3 d = pos1 - pos2;

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

void mdsystem::print_output_and_process_events()
{
    print_output();
    process_events();
}

void mdsystem::process_events()
{
    // Let the application process its events
    if (event_callback.func) {
        event_callback.func(event_callback.param);
    }
}

void mdsystem::print_output()
{
    if (output.str().empty()) {
        // Nothing to write
        return;
    }

    // Print the contents of the output buffer and then empty it
    if (output_callback.func) {
        output_callback.func(output_callback.param, output.str());
    }
    output.str("");
}

void mdsystem::start_operation()
{
    while (operating) {
        // Wait for the other operation to finish
        process_events();
    }
    operating = true;
}

void mdsystem::finish_operation()
{
    print_output();
    if (!operating) {
        throw runtime_error("Tried to finish operation that was never started");
    }
    operating = false;
}
