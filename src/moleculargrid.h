/**************************************************************************
 *   moleculargrid.h  --  This file is part of DFTCXX.                    *
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

#ifndef _MOLECULAR_GRID_H
#define _MOLECULAR_GRID_H

#ifdef HAS_OPENMP
#include <mutex>
#endif

#include <Eigen/Dense>
#include <boost/math/special_functions/factorials.hpp>

#include "molecule.h"
#include "quadrature.h"

#include "bragg.h"

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


    /**
     * @fn scale_density
     *
     * @brief multiplies density at gridpoint with factor
     *
     * @param factor multiplication factor
     *
     * @return void
     */
    void scale_density(double factor);

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

    enum {  // defines fineness of the numerical integration
        GRID_COARSE,
        GRID_MEDIUM,
        GRID_FINE,
        GRID_ULTRAFINE,

        NR_GRID_RESOLUTIONS
    };

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

    void calculate_hartree_potential();

    /*
     *  vector & matrix getters
     */

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
     * @fn renormalize_density
     * @brief normalize density at gridpoints so that sum equals all electrons
     *
     * @return void
     */
    void renormalize_density(unsigned int nr_elec);

private:
    static constexpr double pi = 3.14159265358979323846;

    std::shared_ptr<Molecule> mol;  // pointer to molecule object

    std::vector<GridPoint> grid;    // set of all gridpoints

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

    /**
     * @fn create_grid
     * @brief creates the molecular grid
     *
     * @param fineness      ENUM (unsigned int) defining the resolution of the grid
     *
     * Creates a molecular grid given a set of atoms and a set of
     * basis functions. For each atom, an atomic grid is created. The
     * numerical integration is carried out over all the atomic grids. The
     * contribution of the atomic grid to the overall integration over the
     * whole molecule is controlled via a weight factor as defined by Becke.
     *
     * @return void
     */
    void create_grid(unsigned int fineness = GRID_MEDIUM);

    void calculate_rho_lm();
    void calculate_U_lm();

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
    double get_becke_weight_pn(unsigned int atnr, const vec3& p0);

    /**
     * @brief      Gets a correction factor for the Becke mu
     *
     * @param[in]  at1   atomic number of atom 1
     * @param[in]  at2   atomic number of atom 2
     * @param[in]  mu    Becke weight value
     *
     * @return     mu correction factor
     */
    double get_becke_mu_correction_hetero_atoms(unsigned int at1, unsigned int at2, double mu);

    /*
     * auxiliary functions for the Becke grid
     */

    /**
     * @fn cutoff
     * @brief calculates the Becke weight
     *
     * @param fineness      mu (elliptical coordinate)
     *
     * @return double Becke weight value
     */
    double cutoff(double mu);

    /**
     * @fn fk
     * @brief recursive function to have smooth overlapping atomic grids
     *
     * @param k     number of iterations
     * @param mu    elliptical coordinate
     *
     * @return value
     */
    double fk(unsigned int k, double mu);

    double spherical_harmonic(int l, int m, double pole, double azimuth) const;
    double prefactor_spherical_harmonic(int l, int m) const;
    double polar_function(int l, int m, double theta) const;
    double azimuthal_function(int m, double phi) const;
    double legendre (int n, double x) const;
    double legendre_p (int n, int m, double x) const;

    double dz2dr2(double r, double m);
    double dzdrsq(double r, double m);
};

#endif //_MOLECULAR_GRID_H
