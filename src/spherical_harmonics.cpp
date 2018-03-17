/**************************************************************************
 *   This file is part of DFTCXX.                                         *
 *                                                                        *
 *   Author: Ivo Filot <ivo@ivofilot.nl>                                  *
 *                                                                        *
 *   DFTCXX is free software:                                             *
 *   you can redistribute it and/or modify it under the terms of the      *
 *   GNU General Public License as published by the Free Software         *
 *   Foundation, either version 3 of the License, or (at your option)     *
 *   any later version.                                                   *
 *                                                                        *
 *   DFTCXX is distributed in the hope that it will be useful,            *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 *   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 *   See the GNU General Public License for more details.                 *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 *                                                                        *
 **************************************************************************/

#include "spherical_harmonics.h"

double SH::spherical_harmonic(int l, int m, double pole, double azimuth) {
    return SH::polar_function(l, m, pole) * SH::azimuthal_function(m, azimuth);
}

double SH::prefactor_spherical_harmonic(int l, int m) {
    static const double pre = 1.0 / std::sqrt(4 * M_PI);
    return pre * (m == 0 ? 1 : std::sqrt(2.0)) *
            std::sqrt((double)(2 * l + 1) * boost::math::factorial<double>(l - std::abs(m)) /
                 boost::math::factorial<double>(l + std::abs(m)) );
}

double SH::polar_function(int l, int m, double theta) {
    return SH::legendre_p(l, std::abs(m), std::cos(theta));
}

double SH::azimuthal_function(int m, double phi) {
    if(m == 0) return 1.0;

    if(m > 0) {
        return std::cos((double)m * phi);
    } else {
        return std::sin(-(double)m * phi);
    }
}

/* legendre function */

double SH::legendre (int n, double x) {
    int i;

    if(n < 0) {
        return -1;
    }

    double v[n];
        v[0] = 1.0;

    if(n < 1) {
        return 1.0;
    }

    v[1] = x;

    for ( i = 2; i <= n; i++ ) {
          v[i] =     ( ( double ) ( 2 * i - 1 ) * x    * v[i-1]
                     - ( double ) (     i - 1 ) *        v[i-2] )
                     / ( double ) (     i     );
    }

    return v[n];
}

/* associated legendre function
 *
 * note that x should lie between -1 and 1 for this to work, else
 * a NAN will be returned
 */

double SH::legendre_p (int n, int m, double x) {
    double fact;
    int i;
    int j;
    int k;
    double v[n+1];

    for ( i = 0; i < n + 1; i++ ) {
        v[i] = 0.0;
    }

    //
    //  J = M is the first nonzero function.
    //
    if ( m <= n ) {
        v[m] = 1.0;

        fact = 1.0;
        for ( k = 0; k < m; k++ ) {
            v[m] *= - fact * std::sqrt ( 1.0 - x * x);
            fact += 2.0;
        }
    }

    //
    //  J = M + 1 is the second nonzero function.
    //
    if ( m + 1 <= n ) {
        v[m+1] = x * ( double ) ( 2 * m + 1 ) * v[m];
    }
    //
    //  Now we use a three term recurrence.
    //
    for ( j = m + 2; j <= n; j++ ) {
          v[j] = ( ( double ) ( 2 * j     - 1 ) * x * v[(j-1)]
              + ( double ) (   - j - m + 1 ) *        v[(j-2)] )
              / ( double ) (     j - m     );
    }

    return v[n];
}
