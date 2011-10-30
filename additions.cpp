void system::calculate_specific_heat() {							//obs, no reduced units
float T2 = 0;
for (i = 0; i < nrinst; i++){
	T2 += insttemp[i];
	}
T2 = T2/nrinst;
Cv[timestep/nrinst] = 9*kB/(6/nrparticles+4-4*T2/temp[timestep/nrinst]);
}

void system::calculate_pressure() {												//not finished
float V = nrcells*nrcells*nrcells*outer_cutoff*outer_cutoff*outer_cutoff;
pressure[timestep/nrinst] = nrparticles*kB*temp[timestep/nrinst]/V
}

void system::calculate_mean_square_displacement() {
float sum = 0;
for (i = 0; i < nrparticles;i++) {
	sum += ((*particles[i]).pos - (*particles[i]).fcc_pos)*((*particles[i]).pos - (*particles[i]).fcc_pos);
	}
sum = sum/nrparticles;
msd[timestep/nrinst] = sum;
}


//particle should have a member called fcc_pos that has the value that pos has at the start of the simulation 
//cellsize in create_verlet_list_using_linked_cell_list should be outer_cutoff