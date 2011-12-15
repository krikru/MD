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
    abort_activities_requested = false;
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

void mdsystem::init(uint nrparticles_in, ftype sigma_in, ftype epsilon_in, ftype inner_cutoff_in, ftype outer_cutoff_in, ftype particle_mass_in, ftype dt_in, uint ensemblesize_in, uint sample_period_in, ftype temperature_in, uint nrtimesteps_in, ftype lattice_constant_in, uint lattice_type_in, ftype desired_temp_in, ftype thermostat_time_in, ftype dEp_tolerance_in, ftype impulse_response_decay_time_in, bool thermostat_on_in, bool diff_c_on_in, bool Cv_on_in, bool pressure_on_in, bool msd_on_in, bool Ep_on_in, bool Ek_on_in)
{
    // The system is *always* operating when running non-const functions
    start_operation();

#if THERMOSTAT == LASSES_THERMOSTAT
    thermostat_value = 0;
#endif
    // Copy in parameters to member variables
    particle_mass    = particle_mass_in;
    sigma            = sigma_in;
    epsilon          = epsilon_in;
    ensemblesize     = ensemblesize_in;
    sampling_period  = sample_period_in;
    init_temp        = temperature_in;
    lattice_constant = lattice_constant_in;
    lattice_type     = lattice_type_in; // One of the supported lattice types listed in enum_lattice_types
    dt               = dt_in;           // Delta time, the time step to be taken when solving the diff.eq.
    outer_cutoff     = outer_cutoff_in;
    inner_cutoff     = inner_cutoff_in;
    dEp_tolerance    = dEp_tolerance_in;
    diff_c_on        = diff_c_on_in;
    Cv_on            = Cv_on_in;
    pressure_on      = pressure_on_in;
    msd_on           = msd_on_in;
    Ep_on            = Ep_on_in;
    Ek_on            = Ek_on_in;
    desired_temp     = desired_temp_in;
    thermostat_time  = thermostat_time_in;
    impulse_response_decay_time = impulse_response_decay_time_in;

    /*
     * Reduced units
     *
     * Length unit: sigma
     * Energy unit: epsilon
     * Mass unit: particle mass
     * Temperature unit: epsilon/KB
     *
     * Time unit: sigma * (particle mass / epsilon)^.5
     */
    // Lengths
    lattice_constant /= sigma;
    inner_cutoff     /= sigma;
    outer_cutoff     /= sigma;
    // Temperatures
    init_temp        *= P_KB/ epsilon;
    desired_temp     *= P_KB/ epsilon;
    // Times
    dt               /= sqrt(particle_mass * sigma * sigma / epsilon);
    thermostat_time  /= sqrt(particle_mass * sigma * sigma / epsilon);

    sqr_outer_cutoff = outer_cutoff*outer_cutoff; // Parameter for the Verlet list
    sqr_inner_cutoff = inner_cutoff*inner_cutoff; // Parameter for the Verlet list

    loop_num = 0;
#if FILTER == 0
    num_time_steps = ((nrtimesteps_in - 1) / sampling_period + 1) * sampling_period; // Make the smallest multiple of sample_period that has at least the specified size
    num_sampling_points = num_time_steps/sampling_period + 1;
#elif FILTER == 1
    num_sampling_points = ((nrtimesteps_in - 1) / ensemblesize + 1) * ensemblesize ; // Make the smallest multiple of ensemblesize  that has at least the specified size
    num_time_steps = num_sampling_points + 1;
#endif

    insttemp             .resize(num_sampling_points);
    instEk               .resize(num_sampling_points);
    instEp               .resize(num_sampling_points);
    temperature          .resize(num_sampling_points);
    thermostat_values    .resize(num_sampling_points);
    Ek                   .resize(num_sampling_points);
    Ep                   .resize(num_sampling_points);
    cohesive_energy      .resize(num_sampling_points);
    Cv                   .resize(num_sampling_points);
    pressure             .resize(num_sampling_points);
    msd                  .resize(num_sampling_points);
    diffusion_coefficient.resize(num_sampling_points);
    distanceforcesum     .resize(num_sampling_points);

    if (lattice_type == LT_FCC) {
        box_size_in_lattice_constants = int(pow(ftype(nrparticles_in / 4 ), ftype( 1.0 / 3.0 )));
        num_particles = 4*box_size_in_lattice_constants*box_size_in_lattice_constants*box_size_in_lattice_constants;   // Calculate the new number of atoms; all can't fit in the box since n is an integer
    }
    else {
        cerr << "Lattice type unknown" << endl;
        return;
    }

    // Verlet list and box cells
    box_size = lattice_constant*box_size_in_lattice_constants;
    pos_half_box_size = 0.5f * box_size;
    neg_half_box_size = -pos_half_box_size;
    box_size_in_cells = int(box_size/outer_cutoff);
    if (box_size_in_cells > 3) {
        cells_used = true;
        cell_size = box_size/box_size_in_cells;
    }
    else {
        cells_used = false;
    }

    // Thermostat
    thermostat_on = thermostat_on_in;
    equilibrium = false;

    //
    init_particles();
    create_verlet_list();
    calculate_potential_energy_cutoff();
    calculate_potential_energy_shift();

    // Finish the operation
    finish_operation();
}

void mdsystem::run_simulation()
{
    // The system is *always* operating when running non-const functions
    start_operation();
    /*
     * All variables define in this function has to defined here since we use
     * goto's.
     */
    // Open the output files. They work like cin
    ofstream out_etot_data ;
    ofstream out_ep_data ;
    ofstream out_ek_data;
    ofstream out_cv_data;
    ofstream out_temp_data ;
    ofstream out_therm_data;
    ofstream out_msd_data  ;
    ofstream out_cohe_data ;
    ofstream out_pressure_data;
    // For calculating the average specific heat
    ftype Cv_sum;
    ftype Cv_num;

    ftype sum_sqr_vel;
    // Start simulating
    for (loop_num = 0; loop_num <= num_time_steps; loop_num++) {
        // Check if the simulation has been requested to abort
        if (abort_activities_requested) {
            goto operation_finished;
        }
        //output <<"loop number = " << loop_num << endl;

        // Evolve the system in time
        //cout << loop_num << endl;
        force_calculation();
        leapfrog(); // TODO: Compensate for half time steps

        if (loop_num % sampling_period == 0) {
            update_non_modulated_relative_particle_positions();
            sum_sqr_vel = 0;
            for (uint i = 0; i < num_particles; i++) {
                sum_sqr_vel = sum_sqr_vel + particles[i].vel.sqr_length();
            }
            insttemp[loop_num / sampling_period] =  sum_sqr_vel / (3 * num_particles );
            if (Ek_on) instEk[loop_num/sampling_period] = 0.5f * sum_sqr_vel;
            thermostat_values[loop_num/sampling_period] = thermostat_value;
            if (msd_on) calculate_mean_square_displacement();
            if (diff_c_on) calculate_diffusion_coefficient();
        }

        // Update Verlet list if necessary
        update_verlet_list_if_necessary();

        // Process events
        print_output_and_process_events();
    }
    calculate_properties();
    output << "Simulation completed." << endl;

    /*
     * TODO: The following code should be moved into another public function.
     * This function should *just* run the simulation since that is what it
     * says it does.
     */
    output << "Opening output files..." << endl;
    if (!(open_ofstream_file(out_etot_data , "TotalEnergy.dat") &&
          open_ofstream_file(out_ep_data   , "Potential.dat"  ) &&
          open_ofstream_file(out_ek_data   , "Kinetic.dat"    ) &&
          open_ofstream_file(out_cv_data   , "Cv.dat"         ) &&
          open_ofstream_file(out_temp_data , "Temperature.dat") &&
          open_ofstream_file(out_therm_data, "Thermostat.dat" ) &&
          open_ofstream_file(out_msd_data  , "MSD.dat"        ) &&
          open_ofstream_file(out_pressure_data,"Pressure.dat"     ) &&
          open_ofstream_file(out_cohe_data , "cohesive.dat"       )
          )) {
        cerr << "Error: Output files could not be opened" << endl;
    }
    else {
        output << "Writing to output files..." << endl;
        print_output_and_process_events();

/////////////////Start writing files////////////////////////////////////////////////////////

        for (uint i = 1; i < temperature.size(); i++) {
            if (abort_activities_requested) {
                break;
            }
        out_temp_data  << setprecision(9) << temperature[i] *epsilon/P_KB          << endl;
            // Process events
            print_output_and_process_events();
        }
        for (uint i = 1; i < Ek.size(); i++) {
            if (abort_activities_requested) {
                break;
            }
        out_etot_data  << setprecision(9) << (Ek[i] + (Ep[i]-Ep_shift))*epsilon/P_EV          << endl;
            // Process events
            print_output_and_process_events();
        }
        for (uint i = 1; i < Ek.size(); i++) {
            if (abort_activities_requested) {
                break;
            }
        out_ek_data    << setprecision(9) << Ek[i]*epsilon/P_EV                    << endl;
            // Process events
            print_output_and_process_events();
        }
        for (uint i = 1; i < Ep.size(); i++) {
            if (abort_activities_requested) {
                break;
            }
        out_ep_data    << setprecision(9) << (Ep[i]-Ep_shift)*epsilon/P_EV                    << endl;
            // Process events
            print_output_and_process_events();
        }
        for (uint i = 1; i < cohesive_energy.size(); i++) {
            if (abort_activities_requested) {
                break;
            }
        out_cohe_data  << setprecision(9) << (cohesive_energy[i])/P_EV*epsilon     << endl;
            // Process events
            print_output_and_process_events();
        }
        for (uint i = 1; i < Cv.size(); i++) {
            if (abort_activities_requested) {
                break;
            }
        out_cv_data    << setprecision(9) << Cv[i]*P_KB/(1000 * particle_mass)     << endl;
            // Process events
            print_output_and_process_events();
        }
        for (uint i = 1; i < msd.size(); i++) {
            if (abort_activities_requested) {
                break;
            }
        out_msd_data   << setprecision(9) << msd[i]*sigma*sigma                    << endl;
            // Process events
            print_output_and_process_events();
        }
        for (uint i = 1; i < thermostat_values.size(); i++) {
            if (abort_activities_requested) {
                break;
            }
        out_therm_data << setprecision(9) << thermostat_values[i]                  << endl;
            // Process events
            print_output_and_process_events();
        }
        for (uint i = 1; i < pressure.size(); i++) {
            if (abort_activities_requested) {
                break;
            }
        out_pressure_data<<setprecision(9)<< pressure[i]*epsilon/(sigma*sigma*sigma)<< endl;
            // Process events
            print_output_and_process_events();
        }

/////////////////Finish writing files///////////////////////////////////////////////////////

        out_etot_data .close();
        out_ep_data   .close();
        out_ek_data   .close();
        out_cv_data   .close();
        out_temp_data .close();
        out_therm_data.close();
        out_msd_data  .close();
        out_cohe_data .close();
        out_pressure_data.close();
    }
    output << "Writing to output files done." << endl;

    for (uint i = 1; i < temperature.size();i++)
    {
        if (abort_activities_requested) {
            goto operation_finished;
        }

        output << "Temp            (K)   = " <<setprecision(9) << temperature[i] *epsilon/P_KB       << endl;
        output << "Ek + Ep         (eV)  = " <<setprecision(9) << (Ek[i] + (Ep[i]-Ep_shift))*epsilon/P_EV      << endl;
        output << "Ek              (eV)  = " <<setprecision(9) << Ek[i]*epsilon/P_EV                << endl;
        output << "Ep              (eV)  = " <<setprecision(9) << (Ep[i]-Ep_shift)*epsilon/P_EV                << endl;
        output << "Cohesive energy (eV)  = " <<setprecision(9) << (cohesive_energy[i])/P_EV*epsilon  << endl;
        output << "Cv              (J/K) = " <<setprecision(9) << Cv[i]*P_KB/(1000 * particle_mass)  << endl;
        output << "msd             (m^2) = " <<setprecision(9) << msd[i]*sigma*sigma        << endl;
        output << "Pressure        (Pa)  = " <<setprecision(9) << pressure[i]*epsilon/(sigma*sigma*sigma)     << endl;

        // Process events
        print_output_and_process_events();
    }

        if (abort_activities_requested) {
            goto operation_finished;
        }
        Cv_sum = 0;
        Cv_num = 0;
        for(uint i = 0; i < Cv.size();i++)//I know it's an ugly filtering method but it actually gives a very nice result sometimes...,
        {
            if (abort_activities_requested) {
                goto operation_finished;
            }
            if((Cv[i]*P_KB/(1000 * particle_mass)< 1.0f)&& (Cv[i]*P_KB/(1000 * particle_mass) > 0.0f))
            {
                Cv_sum += Cv[i]*P_KB/(1000 * particle_mass);
                Cv_num++;
            }
        }
    output <<"*******************"<<endl;
    output<<"Cv = "<<Cv_sum/Cv_num<<endl;
    output<< "a=" << lattice_constant<<endl;
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
    /*
     * This is not an *operation* in that sence since the variable being
     * changed cannot be locked for writing to a single thread
     */
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
    return num_time_steps;
}

////////////////////////////////////////////////////////////////
// PRIVATE FUNCTIONS
////////////////////////////////////////////////////////////////

void mdsystem::init_particles() {
    // Allocate space for particles
    particles.resize(num_particles);

    //Place out particles according to the lattice pattern
    if (lattice_type == LT_FCC) {
        for (uint z = 0; z < box_size_in_lattice_constants; z++) {
            for (uint y = 0; y < box_size_in_lattice_constants; y++) {
                for (uint x = 0; x < box_size_in_lattice_constants; x++) {
                    int help_index = 4*(x + box_size_in_lattice_constants*(y + box_size_in_lattice_constants*z));

                    (particles[help_index + 0]).pos[0] = x*lattice_constant;
                    (particles[help_index + 0]).pos[1] = y*lattice_constant;
                    (particles[help_index + 0]).pos[2] = z*lattice_constant;

                    (particles[help_index + 1]).pos[0] = x*lattice_constant;
                    (particles[help_index + 1]).pos[1] = (y + 0.5f)*lattice_constant;
                    (particles[help_index + 1]).pos[2] = (z + 0.5f)*lattice_constant;

                    (particles[help_index + 2]).pos[0] = (x + 0.5f)*lattice_constant;
                    (particles[help_index + 2]).pos[1] = y*lattice_constant;
                    (particles[help_index + 2]).pos[2] = (z + 0.5f)*lattice_constant;

                    (particles[help_index + 3]).pos[0] = (x + 0.5f)*lattice_constant;
                    (particles[help_index + 3]).pos[1] = (y + 0.5f)*lattice_constant;
                    (particles[help_index + 3]).pos[2] = z*lattice_constant;
                } // X
            } // Y
        } // Z
    }
    
    //Randomize the velocities
    vec3 sum_vel = vec3(0, 0, 0);
    ftype sum_sqr_vel = 0;
    for (uint i = 0; i < num_particles; i++) {
        for (uint j = 0; j < 3; j++) {
            particles[i].vel[j] = 0;
            for (uint terms = 0; terms < 5; terms++) { //This will effectivelly create a distribution very similar to normal distribution. (If you want to see what the distribution looks like, go to www.wolframalpha.com/input/?i=fourier((sinc(x))^n) and replace n by the number of terms)
                particles[i].vel[j] += ftype(rand());
            }
        }
        sum_vel     += particles[i].vel;
        sum_sqr_vel += particles[i].vel.sqr_length();
    }

    // Compensate for incorrect start temperature and total velocities and finalize the initialization values
    vec3 average_vel = sum_vel/ftype(num_particles);
    ftype vel_variance = sum_sqr_vel/num_particles - average_vel.sqr_length();
    ftype scale_factor;

        scale_factor = sqrt(3.0f  * init_temp  / (vel_variance)); // Termal energy = 1.5 * P_KB * init_temp = 0.5 m v*v

    for (uint i = 0; i < num_particles; i++) {
        particles[i].vel = (particles[i].vel - average_vel)* scale_factor;
    }

    reset_non_modulated_relative_particle_positions();
}

void mdsystem::calculate_potential_energy_shift(){
    Ep_shift = 0;
    if (!SHIFT_EP) return;
    for (uint i1 = 0; i1 < num_particles ; i1++) { // Loop through all particles
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
            ftype p;

            p = sqr_distance_inv;
            p = p*p*p;

            if (Ep_on) Ep_shift += 4 * p * (p - 1) - E_cutoff;

        }
    }
}
void mdsystem::calculate_potential_energy_cutoff()
{
    ftype q;
    q = 1/sqr_inner_cutoff;
    q = q * q * q;
    E_cutoff = 4 * q * (q - 1);
}

void mdsystem::update_verlet_list_if_necessary()
{
    // Check if largest displacement too large for not updating the Verlet list
    ftype sqr_limit = (sqr_outer_cutoff + sqr_inner_cutoff - 2*sqrt(sqr_outer_cutoff*sqr_inner_cutoff));
    uint i;
    for (i = 0; i < num_particles; i++) {
        ftype sqr_displacement = modulus_position_minus(particles[i].pos, particles[i].pos_when_verlet_list_created).sqr_length();
        if (sqr_displacement > sqr_limit) {
            break;
        }
    }
    if (i < num_particles) {
        // Displacement that is to large was found
        output << "Verlet list updated. " << int(100*loop_num/num_time_steps) << " % done" <<endl;
        create_verlet_list();
    }
}

void mdsystem::create_verlet_list()
{
    //Updating pos_when_verlet_list_created and non_modulated_relative_pos for all particles
    for (uint i = 0; i < num_particles; i++) {
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
    cell_list.resize(box_size_in_cells*box_size_in_cells*box_size_in_cells);
    cell_linklist.resize(num_particles);
    for (uint i = 0; i < cell_list.size() ; i++) {
        cell_list[i] = 0; // Beware! Particle zero is a member of all cells!
    }
    for (uint i = 0; i < num_particles; i++) {
        uint help_x = int(particles[i].pos[0] / cell_size);
        uint help_y = int(particles[i].pos[1] / cell_size);
        uint help_z = int(particles[i].pos[2] / cell_size);
        if (help_x == box_size_in_cells || help_y == box_size_in_cells || help_z == box_size_in_cells) { // This actually occationally happens
            help_x -= help_x == box_size_in_cells;
            help_y -= help_y == box_size_in_cells;
            help_z -= help_z == box_size_in_cells;
        }
        cellindex = help_x + box_size_in_cells * (help_y + box_size_in_cells * help_z);
        cell_linklist[i] = cell_list[cellindex];
        cell_list[cellindex] = i;
    }
}

void mdsystem::create_verlet_list_using_linked_cell_list() { // This function ctreates the verlet_lists (verlet_vectors) using the linked cell lists
    uint cellindex = 0;
    uint neighbour_particle_index = 0;
    verlet_particles_list.resize(num_particles);
    verlet_neighbors_list.resize(0); //The elements will be push_back'ed to the Verlet list

    //Creating new verlet_list
    verlet_particles_list[0] = 0;
    for (uint i = 0; i < num_particles;) { // Loop through all particles
        // Init this neighbour list and point to the next list
        verlet_neighbors_list.push_back(0); // Reset number of neighbours
        int next_particle_list = verlet_particles_list[i] + 1; // Link to the next particle list

        if (cells_used) { //Loop through all neighbour cells
            // Calculate cell indexes
            uint cellindex_x = int(particles[i].pos[0]/cell_size);
            uint cellindex_y = int(particles[i].pos[1]/cell_size);
            uint cellindex_z = int(particles[i].pos[2]/cell_size);
            if (cellindex_x == box_size_in_cells || cellindex_y == box_size_in_cells || cellindex_z == box_size_in_cells) { // This actually occationally happens
                cellindex_x -= cellindex_x == box_size_in_cells;
                cellindex_y -= cellindex_y == box_size_in_cells;
                cellindex_z -= cellindex_z == box_size_in_cells;
            }
            for (int index_z = int(cellindex_z) - 1; index_z <= int(cellindex_z) + 1; index_z++) {
                for (int index_y = int(cellindex_y) - 1; index_y <= int(cellindex_y) + 1; index_y++) {
                    for (int index_x = int(cellindex_x) - 1; index_x <= int(cellindex_x) + 1; index_x++) {
                        int modulated_x = index_x;
                        int modulated_y = index_y;
                        int modulated_z = index_z;
                        // Control boundaries
                        if (modulated_x == -1) {
                            modulated_x = int(box_size_in_cells) - 1;
                        }
                        else if (modulated_x == int(box_size_in_cells)) {
                            modulated_x = 0;
                        }
                        if (modulated_y == -1) {
                            modulated_y = int(box_size_in_cells) - 1;
                        }
                        else if (modulated_y == int(box_size_in_cells)) {
                            modulated_y = 0;
                        }
                        if (modulated_z == -1) {
                            modulated_z = int(box_size_in_cells) - 1;
                        }
                        else if (modulated_z == int(box_size_in_cells)) {
                            modulated_z = 0;
                        }
                        cellindex = uint(modulated_x + box_size_in_cells * (modulated_y + box_size_in_cells * modulated_z)); // Calculate neighbouring cell index
                        neighbour_particle_index = cell_list[cellindex]; // Get the largest particle index of the particles in this cell
                        while (neighbour_particle_index > i) { // Loop though all particles in the cell with greater index
                            // TODO: The modolus can be removed if
                            ftype sqr_distance = modulus_position_minus(particles[i].pos, particles[neighbour_particle_index].pos).sqr_length();
                            if(sqr_distance < sqr_outer_cutoff) {
                                verlet_neighbors_list[verlet_particles_list[i]] += 1;
                                verlet_neighbors_list.push_back(neighbour_particle_index);
                                next_particle_list++;
                            }
                            neighbour_particle_index = cell_linklist[neighbour_particle_index]; // Get the next particle in the cell
                        }
                    } // X
                } // Y
            } // Z
        } // if (cells_used)
        else {
            for (neighbour_particle_index = i+1; neighbour_particle_index < num_particles; neighbour_particle_index++) { // Loop though all particles with greater index
                ftype sqr_distance = modulus_position_minus(particles[i].pos, particles[neighbour_particle_index].pos).sqr_length();
                if(sqr_distance < sqr_outer_cutoff) {
                    verlet_neighbors_list[verlet_particles_list[i]] += 1;
                    verlet_neighbors_list.push_back(neighbour_particle_index);
                    next_particle_list++;
                }
            }
        }
        i++; // Continue with the next particle (if there exists any)
        if (i < num_particles) { // Point to the next particle list
            verlet_particles_list[i] = next_particle_list;
        }
    }
}

void mdsystem::reset_non_modulated_relative_particle_positions()
{
    for (uint i = 0; i < num_particles; i++) {
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
    for (uint i = 0; i < num_particles; i++) {
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
    //TODO: What if the temperature (insttemp) is zero? Has to randomize new velocities in that case.
#if THERMOSTAT == LASSES_THERMOSTAT
    thermostat_value = thermostat_on && loop_num > 0 && loop_num % sampling_period == 0 ? (1 - desired_temp/insttemp[loop_num / sampling_period - 1]) / (2*thermostat_time) : 0;

#elif THERMOSTAT == CHING_CHIS_THERMOSTAT

    /////Using Smooth scaling Thermostat (Berendsen et. al, 1984)/////
    thermostat_value = thermostat_on && loop_num > 0 ? sqrt(1 +  dt / thermostat_time * ((desired_temp) / insttemp[(loop_num-1) / sampling_period] - 1)) : 1;
#endif

    for (uint i = 0; i < num_particles; i++) {

        //TODO: Check if vel and pos are stored for the same time or not, in that case, compensate for that

        // Update velocities
#if THERMOSTAT == CHING_CHIS_THERMOSTAT
        if (thermostat_on) {
            particles[i].vel = particles[i].vel * thermostat_value;
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
    }
}
void mdsystem::force_calculation() {


    // Reset accelrations for all particles
    for (uint k = 0; k < num_particles; k++) {
        particles[k].acc = vec3(0, 0, 0);
    }
    if (loop_num % sampling_period == 0) {
        instEp[loop_num/sampling_period] = 0;
        distanceforcesum[loop_num/sampling_period] = 0;
    }


    for (uint i1 = 0; i1 < num_particles ; i1++) { // Loop through all particles
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
            ftype p;

            p = sqr_distance_inv;
            p = p*p*p;
            acceleration = 48  * distance_inv * p * (p - ftype(0.5));

            // Update accelerations of interacting particles
            vec3 r_hat = r * distance_inv;
            particles[i1].acc +=  acceleration * r_hat;
            particles[i2].acc -=  acceleration * r_hat;

            // Update properties
            //TODO: Remove these two from force calculation and place them somewhere else

            if (Ep_on && loop_num % sampling_period == 0) {
                instEp[loop_num / sampling_period] += 4 * p * (p - 1) - E_cutoff;
            }
            if (pressure_on && loop_num % sampling_period == 0) {
                distanceforcesum[loop_num/sampling_period] += epsilon*acceleration / distance_inv;
            }

        }
    }
#if THERMOSTAT == LASSES_THERMOSTAT
    if (thermostat_on) {
        for (uint i = 0; i < num_particles; i++) {
            particles[i].acc -= thermostat_value * particles[i].vel/sampling_period;
        }
    }

#endif
}

void mdsystem::calculate_properties() {
    calculate_temperature();
    if (Cv_on) calculate_specific_heat();            
    if (pressure_on) calculate_pressure();
    if (Ep_on) {
        calculate_Ep();
        calculate_cohesive_energy();
    }
    if (Ek_on) calculate_Ek();
}

void mdsystem::calculate_temperature() {
    filter(insttemp, temperature, impulse_response_decay_time);
}

void mdsystem::calculate_Ep() {
    filter(instEp, Ep, impulse_response_decay_time);
}

void mdsystem::calculate_cohesive_energy() {
    for (uint i = 0; i < cohesive_energy.size(); i++) {
        cohesive_energy[i] = -Ep[i]/num_particles;
    }
}

void mdsystem::calculate_Ek() {
    filter(instEk, Ek, impulse_response_decay_time);
}

void mdsystem::calculate_specific_heat() {
    vector<ftype> instT2(insttemp.size());
    vector<ftype> T2(temperature.size());
    for (uint i = 0; i < insttemp.size(); i++){
        instT2[i] = insttemp[i]*insttemp[i];
    }
    filter(instT2, T2, impulse_response_decay_time);

    Cv.resize(T2.size());
    for (uint i = 0; i < Cv.size(); i++) {
        Cv[i] = 1/(ftype(2)/3 + num_particles*(1 - T2[i]/(temperature[i]*temperature[i])));
    }
}

void mdsystem::calculate_pressure() {
    ftype V = box_size*box_size*box_size;
    vector<ftype> filtereddistanceforcesum(temperature.size());
    filter(distanceforcesum, filtereddistanceforcesum, impulse_response_decay_time);


    for (uint i = 0; i < pressure.size(); i++) {
        pressure[i] = epsilon*num_particles*temperature[i]/V+ filtereddistanceforcesum[i]/(3*V);
    }
}

void mdsystem::calculate_mean_square_displacement() {
    ftype sum = 0;
    if (equilibrium == false) {
        // Check if equilibrium has been reached
        ftype variation = (instEp[loop_num/sampling_period] - instEp[loop_num/sampling_period - 1]) / instEp[loop_num/sampling_period];
        variation = variation >= 0 ? variation : -variation;
        if (variation < dEp_tolerance) { //TODO: Is this a sufficient check? Probably not
            // The requirements for equilibrium has been reached
            equilibrium = true;
            // Consider the particles to "start" now
            reset_non_modulated_relative_particle_positions();
        }
    }
    if (equilibrium) {
        // Calculate mean square displacement
        for (uint i = 0; i < num_particles;i++) {
            sum += particles[i].non_modulated_relative_pos.sqr_length();
        }
        sum = sum/num_particles;
        msd[loop_num/sampling_period] = sum;
    }
    else {
        // Equilibrium not reached; don't calculate this property.
        msd[loop_num/sampling_period] = 0;
    }
}

void mdsystem::calculate_diffusion_coefficient()
{
    diffusion_coefficient[loop_num/sampling_period] = msd[loop_num/sampling_period]/(6*dt*loop_num);
}

void mdsystem::filter(vector<ftype> &unfiltered, vector<ftype> &filtered, ftype impulse_response_decay_time) {

#if FILTER == 0
    ftype f = exp(-dt*sampling_period/impulse_response_decay_time);
    ftype k = 1 - f;
    ftype a, w;
    int vector_size = unfiltered.size();
    vector<ftype> total_weight(vector_size);
    filtered.resize(vector_size);

    // Left side exponential decay
    a = w = 0;
    for (int i = 0; i < vector_size; i++) {
        a = f*a + k*unfiltered[i];
        w = f*w + k              ;
        filtered    [i] = a;
        total_weight[i] = w;
    }

    // Right side exponential decay
    a = w = 0;
    for (int i = vector_size - 1; i >= 0; i--) {
        a = f*a + k*unfiltered[i];
        w = f*w + k              ;
        filtered    [i] += a;
        total_weight[i] += w;

        // Compensate for weights at the same time
        filtered[i] /= total_weight[i];
    }
#endif
#if FILTER == 1
    filtered.resize(unfiltered.size()/ensemblesize);
    for(uint i = 0; i < filtered.size(); i++){
        ftype sum = 0;
        for(uint j = 0; j < ensemblesize; j++){
            sum += unfiltered[i*ensemblesize+j];
        }
        filtered[i] = sum / ensemblesize;
    }

#endif
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
    if (d[0] >= pos_half_box_size) {
        d[0] -= box_size;
        while (d[0] >= pos_half_box_size) {
            d[0] -= box_size;
        }
    }
    else {
        while (d[0] < neg_half_box_size) {
            d[0] += box_size;
        }
    }

    // Check boundaries in y-direction
    if (d[1] >= pos_half_box_size) {
        d[1] -= box_size;
        while (d[1] >= pos_half_box_size) {
            d[1] -= box_size;
        }
    }
    else {
        while (d[1] < neg_half_box_size) {
            d[1] += box_size;
        }
    }

    // Check boundaries in z-direction
    if (d[2] >= pos_half_box_size) {
        d[2] -= box_size;
        while (d[2] >= pos_half_box_size) {
            d[2] -= box_size;
        }
    }
    else {
        while (d[2] < neg_half_box_size) {
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
