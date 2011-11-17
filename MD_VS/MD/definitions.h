#ifndef  DEFINITIONS_H
#define  DEFINITIONS_H

////////////////////////////////////////////////////////////////
// INCLUDES
////////////////////////////////////////////////////////////////

#include "base_float_vec3.h"

////////////////////////////////////////////////////////////////
// TYPEDEFS
////////////////////////////////////////////////////////////////

typedef  unsigned int            uint ;
typedef  float                   ftype;
typedef  base_float_vec3<ftype>  vec3 ;

//Mathematical constants
const ftype M_SQRT2    = 1.41421356237309504880168872420969807856967187537695f; // Square root of 2
const ftype M_PI       = 3.14159265358979323846264338327950288419716939937510f; // pi
const ftype M_E        = 2.71828182845904523536028747135266249775724709369995f; // Euler's number

//Physical constants
const ftype P_ERG      = 1e-7f             ; // [J  ] Ergon
const ftype P_ANGSTROM = 1e-10f            ; // [m  ] Ångström
const ftype P_U        = 1.66053892173e-27f; // [kg ] The "unified atomic mass unit" (1.660538921(73)*10^-27)
const ftype P_PS       = 1e-12f            ; // [s  ] pikoseconds
const ftype P_KB       = 1.380648813e-23f  ; // [J/K] Boltzmann constant (1.3806488(13)*10^-23)

#endif  /* DEFINITIONS_H */
