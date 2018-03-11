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

#ifndef _ATOMIC_GRID_H
#define _ATOMIC_GRID_H

#include "gridpoint.h"
#include "quadrature.h"
#include "molecule.h"

class AtomicGrid {
private:
    static constexpr double pi = 3.14159265358979323846;

    // pointer to Atom object
    std::shared_ptr<Atom> atom;

    // pointer to molecule object
    std::shared_ptr<Molecule> mol;

    // density coefficients (MXN) matrix where M are the radial points and N the combined lm index
    MatrixXXd rho_lm;

    // hartree potential coefficient (MXN) matrix where M are the radial points and N the combined lm index
    MatrixXXd U_lm;

    // vector holding radial distances
    VectorXd r_n;

    std::vector<GridPoint> grid;                    // set of all gridpoints

    unsigned int angular_points;
    unsigned int radial_points;
    unsigned int lebedev_offset;
    unsigned int lmax;

    MatrixXXd Jj;                                   // numerical J matrix for this atom

public:
    AtomicGrid(const std::shared_ptr<Atom>& _atom);

    void create_atomic_grid(unsigned int _radial_points,
                            unsigned int _angular_points,
                            unsigned int _lebedev_offset,
                            unsigned int _lmax,
                            const std::shared_ptr<Molecule>& _mol);

    inline size_t get_grid_size() const {
        return this->grid.size();
    }

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
     * @fn get_positions
     * @brief get the positions of all the grid points as a vector
     *
     * @return vector containing all the positions
     */
    MatrixXXd get_positions() const;

    /**
     * @fn correct_weights
     * @brief correct the weights
     *
     */
    void correct_weights(const VectorXd& corr);

    /**
     * @fn get_amplitudes
     * @brief get the amplitudes of all the grid points and all the basis functions as a matrix
     *
     * @return matrix (basis functions x grid points)
     */
    MatrixXXd get_amplitudes() const;

    /**
     * @brief      get coulomb matrix
     *
     * @return     Coulomb matrix
     */
    inline const MatrixXXd& get_J() const {
        return this->Jj;
    }

    void calculate_rho_lm();
    void calculate_U_lm();

private:
    double spherical_harmonic(int l, int m, double pole, double azimuth) const;
    double prefactor_spherical_harmonic(int l, int m) const;
    double polar_function(int l, int m, double theta) const;
    double azimuthal_function(int m, double phi) const;
    double legendre (int n, double x) const;
    double legendre_p (int n, int m, double x) const;

    double d2zdr2(double r, double m);
    double dzdrsq(double r, double m);
};

#endif //_ATOMIC_GRID_H
