#ifndef __INCLUDED_CONSTANT_HPP__
#define __INCLUDED_CONSTANT_HPP__

#include <math.h>

namespace constant{
    const double NA = 6.022140857 * pow(10,23); // Avogadoro constant [mol^-1]
    const double kB = 8.6173324 * pow(10,-5); // Boltzman constant [eV K^-1]
    const double h = 4.135667662 * pow(10,-15); // Planck constant [eV s]
    const double joule = 1.6021766208 * pow(10,-19); // [eV J^-1]
    const double angstrom = pow(10,10);          // [Å m^-1]
}

#endif
