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

#ifndef _MOLECULAR_GRID_H
#define _MOLECULAR_GRID_H

#include <boost/math/special_functions/factorials.hpp>
#include <cmath>
#include <chrono>

#include "molecule.h"
#include "gridpoint.h"
#include "atomicgrid.h"

/*
 * @class MolecularGrid
 * @brief Set of grid points for the numerical integration
 *
 * For evaluating non-local properties of the system, a numerical integration
 * has to be performed. The numerical integration proceeds over a set of grid points.
 * The MolecularGrid routine constructs these grid points and sets appropriate
 * weights for these grid points according to the procedure by Becke:
 *
 * A multicenter numerical integration scheme for polyatomic molecules
 * A. D. Becke
 * The Journal of Chemical Physics 88, 2547 (1988); doi: 10.1063/1.454033
 * http://dx.doi.org/10.1063/1.454033
 *
 */
class MolecularGrid {
private:
    static constexpr double pi = 3.14159265358979323846;

    std::shared_ptr<Molecule> mol;                                // pointer to molecule object

    std::vector<std::unique_ptr<AtomicGrid>> atomic_grids;        // vector of pointers to atomic grids

    std::vector<double> atomic_densities;

    unsigned int angular_points;
    unsigned int radial_points;
    unsigned int lebedev_order;
    unsigned int lebedev_offset;
    int lmax;

    // density coefficients (MXN) matrix where M are the radial points and N the combined lm index
    MatrixXXd rho_lm;

    // hartree potential coefficient (MXN) matrix where M are the radial points and N the combined lm index
    MatrixXXd U_lm;

    // vector holding radial distances
    VectorXd r_n;

    size_t grid_size;       //!< total grid size

public:
    /**
     * @fn MolecularGrid
     * @brief MolecularGrid constructor
     *
     * @param _mol      pointer to Molecule class
     *
     * @return MolecularGrid instance
     */
    MolecularGrid(const std::shared_ptr<Molecule>& _mol);

    /**
     * @fn calculate_density
     * @brief calculate the total electron density (number of electrons)
     *
     * Calculates the total number of electrons by summing the weights
     * multiplied by the local value of the electron density.
     *
     * @return number of electrons
     */
    double calculate_density() const;

    /**
     * @brief      Sets the grid parameters (fine-tuning)
     *
     * @param[in]  _radial_points  number of radial points
     * @param[in]  _lebedev_order  The lebedev order
     * @param[in]  _lmax           maximum angular momentum for spherical harmonics
     */
    void set_grid_parameters(unsigned int _radial_points, unsigned int _lebedev_order, unsigned int _lmax);

    /**
     * @fn create_grid
     * @brief creates the molecular grid
     *
     * Creates a molecular grid given a set of atoms and a set of
     * basis functions. For each atom, an atomic grid is created. The
     * numerical integration is carried out over all the atomic grids. The
     * contribution of the atomic grid to the overall integration over the
     * whole molecule is controlled via a weight factor as defined by Becke.
     *
     * @return void
     */
    void create_grid();

    /**
     * @fn set_density
     * @brief set the density at the grid point given a density matrix
     *
     * @param P     reference to density matrix
     *
     * Loops over all the grid points and calculates the local density
     * using the amplitudes of the basis functions and the density matrix
     *
     * @return void
     */
    void set_density(const MatrixXXd& P);

    /**
     * @brief      Calculates the hartree potential.
     *
     * @return     The hartree potential.
     */
    MatrixXXd calculate_hartree_potential();

    /**
     * @fn get_weights
     * @brief get the weights of all the grid points as a vector
     *
     * @return vector containing all the weights
     */
    VectorXd get_weights() const;

    /**
     * @fn get_densities
     * @brief get the densities of all the grid points as a vector
     *
     * @return vector containing all the densities
     */
    VectorXd get_densities() const;

    /**
     * @fn get_amplitudes
     * @brief get the amplitudes of all the grid points and all the basis functions as a matrix
     *
     * @return matrix (basis functions x grid points)
     */
    MatrixXXd get_amplitudes() const;

    /**
     * @brief      correct the densities to match the total number of electrons
     */
    void correct_densities();

private:
    /**
     * @fn get_becke_weight_pn
     * @brief calculate Pn(r) (i.e. for a single gridpoint) for atom n (eq. 22 in A.D. Becke J.Chem.Phys. 88, 2547)
     *
     * Reference: A.D. Becke, J.Chem.Phys. 88, 2547, (1988)
     * Link: http://dx.doi.org/10.1063/1.454033
     *
     * @param atnr atomic number
     * @param p0   grid point position
     *
     * @return double Becke weight value
     */
    double get_becke_weight_pn(unsigned int i, const vec3& p0);

    /*
     * auxiliary functions for the Becke grid
     */

    double cutoff(double mu);
    double fk(unsigned int k, double mu);
};

#endif //_MOLECULAR_GRID_H
