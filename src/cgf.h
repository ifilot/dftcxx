/**************************************************************************
 *   cgf.h  --  This file is part of DFTCXX.                              *
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

#ifndef _CGF_H
#define _CGF_H

#include <iostream>
#include <Eigen/Dense>
#include <boost/math/special_functions/factorials.hpp>

typedef Eigen::Vector3d vec3;

/*
 * Gaussian Type Orbital
 *
 * N * (x-X)^l * (y-Y)^m * (z-Z)^n * exp(-alpha * r^2)
 *
 * where r = sqrt(x^2 + y^2 + z^2)
 * and N a normalization constant such that <GTO|GTO> = 1
 *
 */

class GTO { // Gaussian Type Orbital
private:
    double c;               // coefficient
    double alpha;           // alpha value in the exponent
    unsigned int l,m,n;     // powers of the polynomial
    vec3 position;          // position vector (unit = Bohr)
    double norm;            // normalization constant

public:
    /*
     * @fn GTO
     * @brief Construct Gaussian Type Orbital
     *
     * @param _c        coefficient
     * @param _position position of the Gaussian
     * @param _alpha    alpha value in the exponent
     * @param _l        power of x in the polynomial
     * @param _m        power of y in the polynomial
     * @param _n        power of z in the polynomial
     *
     * returns double value of the incomplete Gamma Function
     */
    GTO(double _c,
        const vec3& _position,
        double _alpha,
        unsigned int _l,
        unsigned int _m,
        unsigned int _n);

    /*
     * INLINE GETTERS
     */

    /*
     * @fn get_coefficient
     * @brief get the coefficient
     *
     * @return const double coefficient
     */
    inline const double get_coefficient() const {
        return this->c;
    }

    /*
     * @fn get_coefficient
     * @brief Gets the alpha value in the polynomial
     *
     * @return const double coefficient
     */
    inline const double get_alpha() const {
        return this->alpha;
    }

    /*
     * @fn get_l
     * @brief Gets power of the x component
     *
     * @return const double power
     */
    inline const unsigned int get_l() const {
        return this->l;
    }

    /*
     * @fn get_m
     * @brief Gets power of the y component
     *
     * @return const double power
     */
    inline const unsigned int get_m() const {
        return this->m;
    }

    /*
     * @fn get_n
     * @brief Gets power of the z component
     *
     * @return const double power
     */
    inline const unsigned int get_n() const {
        return this->n;
    }

    /*
     * @fn get normalization constant
     * @brief Returns the constant that ensures that the GTO is normalized
     *
     * @return const double norm constant
     */
    inline const double get_norm() const {
        return this->norm;
    }

    /*
     * @fn get_position
     * @brief Get the center of the Gaussian
     *
     * @return const double vec3
     */
    inline const vec3& get_position() const {
        return this->position;
    }

    /*
     * @fn get_amp
     * @brief Gets the amplitude of the GTO
     *
     * @param vec3 r    coordinates
     *
     * @return const double amplitude
     */
    const double get_amp(const vec3& r) const;

    /*
     * @fn set_position
     * @brief Set (new) position of GTO
     *
     * @return void
     */
    inline void set_position(const vec3& _position) {
        this->position = _position;
    }

private:
    /*
     * @fn calculate_normalization_constant
     * @brief Calculates the normalization constant so that <GTO|GTO>=1
     *
     * @return void
     */
    void calculate_normalization_constant();
};


class CGF { // Contracted Gaussian Function
private:
    std::vector<GTO> gtos;  // vector holding all gtos
    vec3 r;                 // position of the CGF

public:
    /*
     * @fn CGF
     * @brief Constructor
     *
     * @return CGF
     */
    CGF();

    /*
     * @fn CGF
     * @brief Default constructor
     *
     * @return CGF
     */
    CGF(const vec3& _r);

    // type of GTOs to add
    enum{
        GTO_S,
        GTO_PX,
        GTO_PY,
        GTO_PZ,
        GTO_DX2,
        GTO_DXY,
        GTO_DXZ,
        GTO_DY2,
        GTO_DYZ,
        GTO_DZ2,

        NUM_GTO
    };

    /*
     * @fn size
     * @brief Returns the length of the contraction
     *
     * @return unsigned int length
     */
    inline const unsigned int size() const {
        return this->gtos.size();
    }

    /*
     * @fn size
     * @brief Returns the normalization constant of a GTO
     *
     * @param unsigned int i        ith GTO in the CGF
     *
     * @return double normalization constant
     */
    inline const double get_norm_gto(const unsigned int i) const {
        return this->gtos[i].get_norm();
    }

    /*
     * @fn size
     * @brief Returns the coefficient of the GTO
     *
     * @param unsigned int i        ith GTO in the CGF
     *
     * @return double GTO coefficient
     */
    inline const double get_coefficient_gto(const unsigned int i) const {
        return this->gtos[i].get_coefficient();
    }

    /*
     * @fn size
     * @brief Returns the GTO
     *
     * @param unsigned int i        ith GTO in the CGF
     *
     * @return GTO
     */
    inline const GTO& get_gto(const unsigned int i) const {
        return this->gtos[i];
    }

    /*
     * @fn get_amp
     * @brief Gets the amplitude of the CGF
     *
     * @param vec3 r    coordinates
     *
     * @return const double amplitude
     */
    const double get_amp(const vec3& r) const;

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
    void add_gto(unsigned int type,
                 double alpha,
                 double c,
                 const vec3& vec3);

    /*
     * @fn set_position
     * @brief Set a (new) center for the CGF
     *
     * @param pos   center of the CGF
     *
     * @return void
     */
    void set_position(const vec3 &pos);
};

#endif //_CGF_H
