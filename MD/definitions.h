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

/***************
 * Thermostats *
 ***************/
/* Enumeration of thermostats */
#define  NO_THERMOSTAT           0
#define  BERENDSEN_THERMOSTAT    1
#define  NOSE_HOOVER_THERMOSTAT  2

#define  CHING_CHIS_THERMOSTAT   BERENDSEN_THERMOSTAT
#define  LASSES_THERMOSTAT       NOSE_HOOVER_THERMOSTAT

#define  THERMOSTAT              CHING_CHIS_THERMOSTAT /////Using Smooth scaling Thermostat (Berendsen et. al, 1984)/////
//#define  THERMOSTAT              LASSES_THERMOSTAT

/***********
 * Filters *
 ***********/
/* Enumeration of filters */
#define  NO_FILTER                           0
#define  TWO_SIDED_EXPONENTIAL_DECAY_FILTER  1
#define  ENSEMBLE_AVERAGE_FILTER             2

#define  KRISTOFERS_FILTER  TWO_SIDED_EXPONENTIAL_DECAY_FILTER
#define  EMILS_FILTER       ENSEMBLE_AVERAGE_FILTER

#define  FILTER  KRISTOFERS_FILTER
//#define  FILTER  EMILS_FILTER

/*********
 * Check *
 *********/
#if !THERMOSTAT
#error No thermostat chosen
#endif
#if !FILTER
#error No filter chosen
#endif

////////////////////////////////////////////////////////////////
// MISCELANEOUS DEFINITIONS
////////////////////////////////////////////////////////////////

#define  USE_DOUBLE_PRECISION  0
#define  SHIFT_EP              1
#define  PRINT_OUTPUT          0

////////////////////////////////////////////////////////////////
// TYPEDEFS
////////////////////////////////////////////////////////////////

typedef  unsigned int            uint ;
#if USE_DOUBLE_PRECISION
typedef  double                  ftype;
#else
typedef  float                   ftype;
#endif
typedef  base_float_vec3<ftype>  vec3 ;

/* Mathematical constants */
#ifndef _MATH_H_
#define M_E		2.7182818284590452354  // Euler's number
#define M_LOG2E		1.4426950408889634074  // log_2(e)
#define M_LOG10E	0.43429448190325182765 // log_10(e)
#define M_LN2		0.69314718055994530942 // ln(2)
#define M_LN10		2.30258509299404568402 // ln(10)
#define M_PI		3.14159265358979323846 // pi
#define M_PI_2		1.57079632679489661923 // pi/2
#define M_PI_4		0.78539816339744830962 // pi/4
#define M_1_PI		0.31830988618379067154 // 1/pi
#define M_2_PI		0.63661977236758134308 // 2/pi
#define M_2_SQRTPI	1.12837916709551257390 // 2/sqrt(pi)
#define M_SQRT2		1.41421356237309504880 // sqrt(2)
#define M_SQRT1_2	0.70710678118654752440 // sqrt(1/2)
#endif

/*
 * SI units
 */

/* Physical constants */
const ftype P_SI_ERG      = ftype(1e-7             ); // [J  ] Ergon
const ftype P_SI_EV       = ftype(1.60217648740e-19); // [J  ] Electron volt
const ftype P_SI_ANGSTROM = ftype(1.00000000000e-10); // [m  ] A
const ftype P_SI_U        = ftype(1.66053892173e-27); // [kg ] The "unified atomic mass unit" (1.660538921(73)*10^-27)
//const ftype P_SI_PS       = ftype(1e-12            ); // [s  ] pikoseconds
const ftype P_SI_FS       = ftype(1e-15            ); // [s  ] femtoseconds
const ftype P_SI_KB       = ftype(1.380648813e-23  ); // [J/K] Boltzmann constant (1.3806488(13)*10^-23)
const ftype P_SI_AVOGADRO = ftype(6.0221420e23     ); // [1/mol] Avogadro constant


/*
 * Reduced units
 */

/* Physical constants */
// Masses / particle_mass_in_kg
#define  P_RU_U         (P_SI_U       /particle_mass_in_kg)
// Lengths / sigma_in_m
#define  P_RU_ANGSTROM  (P_SI_ANGSTROM/sigma_in_m)
// Energies / epsilon_in_j
#define  P_RU_ERG       (P_SI_ERG     /epsilon_in_j)
#define  P_RU_EV        (P_SI_EV      /epsilon_in_j)
// Temperatures * P_KB / epsilon_in_j
// Times / sqrt(particle_mass_in_kg * sigma_in_m * sigma_in_m / epsilon_in_j)
#define  P_RU_PS        (P_SI_PS      /sqrt(particle_mass_in_kg*sigma_in_m*sigma_in_m/epsilon_in_j))
#define  P_RU_FS        (P_SI_FS      /sqrt(particle_mass_in_kg*sigma_in_m*sigma_in_m/epsilon_in_j))
// Pressures * sigma_in_m * sigma_in_m * sigma_in_m / epsilon_in_j
// Unitless
#define  P_RU_AVOGADRO  (P_SI_AVOGADRO)
// Other
#define  P_RU_KB        (1)

////////////////////////////////////////////////////////////////
// MACROS
////////////////////////////////////////////////////////////////

#define  TEMP_SWAP(x, y, temp) { \
    (temp) = (x);                \
    (x) = (y);                   \
    (y) = (temp);                \
    }

#endif  /* DEFINITIONS_H */
