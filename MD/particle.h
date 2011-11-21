#ifndef  PARTICLE_H
#define  PARTICLE_H

#include <vector>
#include "base_float_vec3.h"
using namespace std;

class particle {//Maybe we should have an additional prev_position which we use when we calculate how far the particles has moved when deciding if we should update the verlet lists
public:
    vec3 pos;  // and if we update the verlet list we should also update the prev_position...
    vec3 non_modulated_relative_pos;
    vec3 pos_when_non_modulated_relative_pos_was_calculated;
    vec3 pos_when_verlet_list_created;
    vec3 vel;
    vec3 acc;
    vec3 start_pos;
    vec3 start_vel;
};

#endif  /* PARTICLE_H */
