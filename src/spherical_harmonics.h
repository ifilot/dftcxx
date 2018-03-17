/*************************************************************************
 *
 *  This file is part of DFTCXX.
 *
 *  Author: Ivo Filot <i.a.w.filot@tue.nl>
 *
 *  DFTCXX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  DFTCXX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with DFTCXX.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************/

#ifndef _SPHERICAL_HARMONIC_H
#define _SPHERICAL_HARMONIC_H

#include <cmath>
#include <boost/math/special_functions/factorials.hpp>

namespace SH {
    double spherical_harmonic(int l, int m, double pole, double azimuth);

    double prefactor_spherical_harmonic(int l, int m);

    double polar_function(int l, int m, double theta);

    double azimuthal_function(int m, double phi);

    double legendre (int n, double x);

    double legendre_p (int n, int m, double x);
};

#endif // _SPHERICAL_HARMONIC_H
