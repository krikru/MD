#ifndef  DEFINITIONS_H
#define  DEFINITIONS_H

////////////////////////////////////////////////////////////////
// COMPILER DEFINITIONS
////////////////////////////////////////////////////////////////

//#define WIN32_LEAN_AND_MEAN

////////////////////////////////////////////////////////////////
// INCLUDES
////////////////////////////////////////////////////////////////

#include "base_float_vec3.h"

////////////////////////////////////////////////////////////////
// COMPILE TIME OPTIONS
////////////////////////////////////////////////////////////////

#define  NO_THERMOSTAT           0
#define  BERENDSEN_THERMOSTAT    1
#define  NOSE_HOOVER_THERMOSTAT  2

#define  CHING_CHIS_THERMOSTAT   BERENDSEN_THERMOSTAT
#define  LASSES_THERMOSTAT       NOSE_HOOVER_THERMOSTAT

//#define  THERMOSTAT              CHING_CHIS_THERMOSTAT
#define  THERMOSTAT              LASSES_THERMOSTAT


////////////////////////////////////////////////////////////////
// MISCELANEOUS DEFINITIONS
////////////////////////////////////////////////////////////////

#define  RU_ON                 0
#define  USE_DOUBLE_PRECITION  1

////////////////////////////////////////////////////////////////
// TYPEDEFS
////////////////////////////////////////////////////////////////

typedef  unsigned int            uint ;
#if USE_DOUBLE_PRECITION
typedef  double                  ftype;
#else
typedef  float                   ftype;
#endif
typedef  base_float_vec3<ftype>  vec3 ;

//Mathematical constants
#ifndef _MATH_H_
const ftype M_SQRT2  = ftype(1.41421356237309504880168872420969807856967187537695); // Square root of 2
const ftype M_PI     = ftype(3.14159265358979323846264338327950288419716939937510); // pi
const ftype M_E      = ftype(2.71828182845904523536028747135266249775724709369995); // Euler's number
#elif 0
#define  M_SQRT2  ftype(1.41421356237309504880168872420969807856967187537695) // Square root of 2
#define  M_PI     ftype(3.14159265358979323846264338327950288419716939937510) // pi
#define  M_E      ftype(2.71828182845904523536028747135266249775724709369995) // Euler's number
#endif

//Physical constants
const ftype P_ERG      = ftype(1e-7             ); // [J  ] Ergon
const ftype P_EV       = ftype(1.60217648740e-19); // [J  ] Electron volt
const ftype P_ANGSTROM = ftype(1e-10            ); // [m  ] A
const ftype P_U        = ftype(1.66053892173e-27); // [kg ] The "unified atomic mass unit" (1.660538921(73)*10^-27)
//const ftype P_PS       = ftype(1e-12            ); // [s  ] pikoseconds
const ftype P_FS       = ftype(1e-15            ); // [s  ] femtoseconds
const ftype P_KB       = ftype(1.380648813e-23  ); // [J/K] Boltzmann constant (1.3806488(13)*10^-23)
const ftype P_AVOGADRO = ftype(6.0221420e23     ); // [1/mol] Avogadro constant
#endif  /* DEFINITIONS_H */
