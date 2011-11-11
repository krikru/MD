#ifndef  PARTICLE_H
#define  PARTICLE_H

#include <vector>
#include "base_float_vec3.h"
using namespace std;

class particle {//Maybe we should have an additional prev_position which we use when we calculate how far the particles has moved when deciding if we should update the verlet lists
public:
	fvec3 pos;  // and if we update the verlet list we should also update the prev_position...
    fvec3 no_bound_pos;
    fvec3 prev_pos;
    fvec3 vel;
    fvec3 acc;
	fvec3 start_pos;
    fvec3 start_vel;
};

#endif  /* PARTICLE_H */
