#include <vector>
#include "base_float_vec3.h"
using namespace std;

class particle {//Maybe we should have an additional prev_position which we use when we calculate how far the particles has moved when deciding if we should update the verlet lists
public:
	fvec3 pos;  // and if we update the verlet list we should also update the prev_position...
    fvec3 prev_pos;
    fvec3 vel;
    fvec3 prev_vel;
    fvec3 acc;
	fvec3 start_pos;
    int part_nr; // Might not be needed (part_nr = i is the i-1 element of the vector particles (if part_nr goes from 1 to total number of particles) ) 
};