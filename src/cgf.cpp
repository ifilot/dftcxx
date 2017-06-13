/**************************************************************************
 *   cgf.cpp  --  This file is part of DFTCXX.                            *
 *                                                                        *
 *   Copyright (C) 2016, Ivo Filot                                        *
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

#include "cgf.h"

GTO::GTO(double _c,
         const vec3& _position,     // position (unit = Bohr)
         double _alpha,
         unsigned int _l,
         unsigned int _m,
         unsigned int _n):
    c(_c),
    alpha(_alpha),
    l(_l),
    m(_m),
    n(_n),
    position(_position) {

    // calculate the normalization constant
    this->calculate_normalization_constant();
}

/*
 * @fn get_amp
 * @brief Gets the amplitude of the GTO
 *
 * @param vec3 r    coordinates
 *
 * @return const double amplitude
 */
const double GTO::get_amp(const vec3& r) const {
    double r2 = (r - position).squaredNorm();

    return this->norm *
           std::pow(r[0]-position[0], l) *
           std::pow(r[1]-position[1], m) *
           std::pow(r[2]-position[2], n) *
           std::exp(- alpha * r2);
}

/*
 * @fn calculate_normalization_constant
 * @brief Calculates the normalization constant so that <GTO|GTO>=1
 *
 * @return void
 */
void GTO::calculate_normalization_constant() {
    static const double pi = 3.14159265359;

    double nom =   std::pow(2.0, 2.0 * (l + m + n) + 3.0 / 2.0) *
                   std::pow(alpha, (l + m + n) + 3.0 / 2.0);

    double denom = (l < 1 ? 1 : boost::math::double_factorial<double>(2 * l - 1) )*
                   (m < 1 ? 1 : boost::math::double_factorial<double>(2 * m - 1) )*
                   (n < 1 ? 1 : boost::math::double_factorial<double>(2 * n - 1) )*
                   std::pow(pi, 3.0 / 2.0);

    this->norm = std::sqrt(nom / denom);
}

/*
 * @fn CGF
 * @brief Constructor
 *
 * @return CGF
 */
CGF::CGF():
    r(vec3(0,0,0)) {
        // do nothing
}

/*
 * @fn CGF
 * @brief Default constructor
 *
 * @return CGF
 */
CGF::CGF(const vec3& _r):
    r(_r) {
        // do nothing
}

/*
 * @fn get_amp
 * @brief Gets the amplitude of the CGF
 *
 * @param vec3 r    coordinates
 *
 * @return const double amplitude
 */
const double CGF::get_amp(const vec3& r) const {
    double sum = 0.0;

    for(unsigned int i=0; i<this->gtos.size(); i++) {
        sum += gtos[i].get_coefficient() * gtos[i].get_amp(r);
    }

    return sum;
}

/*
 * @fn add_GTO
 * @brief Add a GTO to the CGF
 *
 * @param unsigned int type     type of the orbital (see above for the list)
 * @param double alpha          alpha value
 * @param double c              coefficient
 * @param const vec3& vec3      position
 *
 * @return void
 */
void CGF::add_gto(unsigned int type,  // type of the orbital (see above for the list)
                  double alpha,       // alpha value
                  double c,           // coefficient
                  const vec3& vec3) { // position

    switch(type) {
        // S ORBITAL
        case GTO_S:
            gtos.push_back(GTO(c, r, alpha, 0,0,0));
        break;

        // P ORBITALS
        case GTO_PX:
            gtos.push_back(GTO(c, r, alpha, 1,0,0));
        break;
        case GTO_PY:
            gtos.push_back(GTO(c, r, alpha, 0,1,0));
        break;
        case GTO_PZ:
            gtos.push_back(GTO(c, r, alpha, 0,0,1));
        break;

        // D ORBITALS
        case GTO_DX2:
            gtos.push_back(GTO(c, r, alpha, 2,0,0));
        break;
        case GTO_DXY:
            gtos.push_back(GTO(c, r, alpha, 1,1,0));
        break;
        case GTO_DXZ:
            gtos.push_back(GTO(c, r, alpha, 1,0,1));
        break;
        case GTO_DY2:
            gtos.push_back(GTO(c, r, alpha, 0,2,0));
        break;
        case GTO_DYZ:
            gtos.push_back(GTO(c, r, alpha, 0,1,1));
        break;
        case GTO_DZ2:
            gtos.push_back(GTO(c, r, alpha, 0,0,2));
        break;

        default:
            std::cerr << "Undefined orbital type. Exiting..." << std::endl;
            exit(-1);
        break;
    }
}

/*
 * @fn set_position
 * @brief Set a (new) center for the CGF
 *
 * @param pos   center of the CGF
 *
 * @return void
 */
void CGF::set_position(const vec3 &pos) {
    this->r = pos;

    for(unsigned int i=0; i<this->gtos.size(); i++) {
        this->gtos[i].set_position(pos);
    }
}
