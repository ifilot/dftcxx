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

#ifndef _GRIDPOINT_H
#define _GRIDPOINT_H

#include <Eigen/Dense>

#include "molecule.h"

typedef Eigen::Vector3d vec3;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXXd;

/*
 * @class GridPoint
 * @brief Represents a point in the molecular grid
 *
 * Grid point at position r; stores local information such as the amplitude
 * of all the basis functions in the basis set and the local density.
 *
 * Numerical integrations of any kind of property is conducted by
 * evaluating the local value of the functional at the grid point and
 * multiplying this by the weight of the grid point. Integration is then
 * performed by summing over all the grid points.
 *
 * The weights take into account:
 * - the Jacobian for integration into spherical coordinates
 * - the weight for the Lebedev integration (angular part)
 * - the weight for the Gauss-Chebyshev integration (radial part)
 * - the weight for the Becke grid (see: http://dx.doi.org/10.1063/1.454033)
 *
 */
class GridPoint {
public:
    /**
     * @fn DFT
     * @brief DFT routine constructor
     *
     * @return DFT instance
     */
    GridPoint(const vec3& _r);

    /*
     * SETTERS
     */

     /**
     * @fn set_weight
     * @brief set the weight at the grid point
     *
     * @param _w      value of the weight
     *
     * @return void
     */
    inline void set_weight(double _w) {
        this->w = _w;
    }

    /**
     * @fn multiply_weight
     * @brief multiply weight by factor at the grid point
     *
     * @param _w      weight multiplication factor
     *
     * @return void
     */
    inline void multiply_weight(double _w) {
        this->w *= _w;
    }

    /**
     * @fn set_atom
     * @brief defines the atom to which this grid point is 'linked'
     *
     * @param at    pointer to Atom class
     *
     * @return void
     */
    inline void set_atom(const Atom* at) {
        this->atom = at;
    }

    /**
     * @fn set_basis_func_amp
     * @brief calculates the amplitudes at the grid point of all basis functions
     *
     * @param _mol      pointer to the molecule object
     *
     * @return void
     */
    void set_basis_func_amp(const std::shared_ptr<Molecule>& _mol);

    /**
     * @fn set_density
     * @brief calculates the density at the grid point using the density matrix
     *
     * @param _mol      reference to density matrix
     *
     * @return void
     */
    void set_density(const MatrixXXd& D);

    /*
     * GETTERS
     */

     /**
     * @fn get_position
     * @brief get the position of the grid point
     *
     * @return vec3 position of the grid point
     */
    inline const vec3& get_position() const {
        return this->r;
    }

    /**
     * @fn get_atom_position
     * @brief get the position of the atom to which this grid point is 'linked'
     *
     * @return vec3 position of the atom
     */
    inline const vec3& get_atom_position() const {
        return this->atom->get_position();
    }

    /**
     * @fn get_weight
     * @brief get the weight of the grid point in the numerical integration
     *
     * @return double weight
     */
    inline const double get_weight() const {
        return this->w;
    }

    /**
     * @fn get_density
     * @brief get the electron density at the grid point
     *
     * @return double density
     */
    inline const double get_density() const {
        return this->density;
    }

    /**
     * @fn get_basis_func_amp
     * @brief get the amplitude of the basis functions
     *
     * @return (Eigen3) vector containing basis function amplitudes
     */
    inline const VectorXd& get_basis_func_amp() const {
        return this->basis_func_amp;
    }

private:
    const vec3 r;               // position in 3D space
    double w;                   // weight
    const Atom* atom;           // atom this gridpoint adheres to
    VectorXd basis_func_amp;    // amplitude of basis functions at gridpoint
    double density;             // current density at grid point
};

#endif // _GRIDPOINT_H
