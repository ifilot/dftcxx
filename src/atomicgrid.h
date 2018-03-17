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
#include "cspline.h"
#include "spherical_harmonics.h"

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
    unsigned int lm;

    // hartree potential
    VectorXd V_fuzzy_cell;
    VectorXd V;

    // vector holding radial distances
    VectorXd r_n;

    std::vector<GridPoint> grid;                    // set of all gridpoints

    std::vector<Cspline> splines;

    unsigned int angular_points;
    unsigned int radial_points;
    unsigned int lebedev_offset;
    unsigned int lmax;

    MatrixXXd Jj;                                   // numerical J matrix for this atom
    double density;                                 // total electron density for this fuzzy cell
    bool density_cached;                            // whether density has been evaluated

public:

    /**
     * @brief      create atomic grid
     *
     * @param[in]  _atom  pointer to atom
     */
    AtomicGrid(const std::shared_ptr<Atom>& _atom);

    /**
     * @brief      Creates an atomic grid.
     *
     * @param[in]  _radial_points   number of radial points
     * @param[in]  _angular_points  number of angular points
     * @param[in]  _lebedev_offset  lebedev offset
     * @param[in]  _lmax            maximum angular quantum number l
     * @param[in]  _mol             pointer to molecule
     */
    void create_atomic_grid(unsigned int _radial_points,
                            unsigned int _angular_points,
                            unsigned int _lebedev_offset,
                            unsigned int _lmax,
                            const std::shared_ptr<Molecule>& _mol);

    /**
     * @brief      Gets the grid size.
     *
     * @return     The grid size.
     */
    inline size_t get_grid_size() const {
        return this->grid.size();
    }

    /**
     * @brief      get the center position of the fuzzy cell
     *
     * @return     fuzzy cell center position
     */
    inline const vec3& get_position() const {
        return this->atom->get_position();
    }

    /**
     * @brief      get the position of the grid point
     *
     * @param[in]  idx   index of grid point
     *
     * @return     grid point position
     */
    inline const vec3& get_position_grid_point(unsigned int idx) const {
        return this->grid[idx].get_position();
    }

    /**
     * @fn get_density
     * @brief get the total electron density
     *
     * @return total electron density
     */
    double get_density();

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
     * @brief      set fuzzy cell weights
     *
     * @param[in]  vw  fuzzy cell weights
     */
    void set_becke_weights(const VectorXd& vw);

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

    /**
     * @brief      calculate density at distance r
     */
    void calculate_rho_lm();

    /**
     * @brief      calculate spherical harmonic coefficients
     */
    void calculate_U_lm();

    /**
     * @brief      Calculates the j matrix.
     */
    void calculate_J_matrix();

    /**
     * @brief      Gets the sh value.
     *
     * @param[in]  r     radial distance from fuzzy cell center
     * @param[in]  lm    lm index
     *
     * @return     The sh value.
     */
    double get_sh_value(double r, unsigned int lm) const;

    /**
     * @brief      get the local hartree potential for this fuzzy cell
     *
     * @param[in]  idx   grid point index
     *
     * @return     local cell hartree potential Vn
     */
    inline double get_V(unsigned int idx) const {
        return this->V_fuzzy_cell(idx);
    }

    /**
     * @brief      Set total Hartree potential for gridpoint
     *
     * @param[in]  idx   grid point index
     * @param[in]  v     total hartree potential
     */
    inline void set_V_gp(unsigned int idx, double v) {
        this->V(idx) = v;
    }

private:
    double d2zdr2(double r, double m);

    double dzdrsq(double r, double m);

    /**
     * @fn calculate_density
     * @brief calculate the total electron density (number of electrons)
     *
     * Calculates the total number of electrons by summing the weights
     * multiplied by the local value of the electron density.
     *
     * @return number of electrons
     */
    void calculate_density();

    /**
     * @brief      perform interpolation on the spherical harmonic coefficients
     */
    void interpolate_sh_coeff();
};

#endif //_ATOMIC_GRID_H
