/**************************************************************************
 *   moleculargrid.cpp  --  This file is part of DFTCXX.                  *
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

#include "moleculargrid.h"

/*
 *
 * GRID POINTS
 *
 */

/**
 * @fn GridPoint
 * @brief GridPoint Constructor
 *
 * @param _r    vec3 position of the grid point
 *
 * @return GridPoint instance
 */
GridPoint::GridPoint(const vec3& _r):
    r(_r) {

}

/**
 * @fn set_basis_func_amp
 * @brief calculates the amplitudes at the grid point of all basis functions
 *
 * @param _mol      pointer to the molecule object
 *
 * @return void
 */
void GridPoint::set_basis_func_amp(const std::shared_ptr<Molecule>& _mol) {
    unsigned int size = _mol->get_nr_bfs();
    this->basis_func_amp = VectorXd(size);

    for(unsigned int i=0; i<size; i++) {
        this->basis_func_amp(i) = _mol->get_cgf(i).get_amp(this->r);
    }
}

/**
 * @fn set_density
 * @brief calculates the density at the grid point using the density matrix
 *
 * @param _mol      reference to density matrix
 *
 * @return void
 */
void GridPoint::set_density(const MatrixXXd& D) {
    this->density = this->basis_func_amp.dot(D * this->basis_func_amp);
}

/**
 * @fn scale_density
 *
 * @brief multiplies density at gridpoint with factor
 *
 * @param factor multiplication factor
 *
 * @return void
 */
void GridPoint::scale_density(double factor) {
    this->density *= factor;
}

/*
 *
 * MOLECULAR GRID
 *
 */

/**
 * @fn MolecularGrid
 * @brief MolecularGrid constructor
 *
 * @param _mol      pointer to Molecule class
 *
 * @return MolecularGrid instance
 */
MolecularGrid::MolecularGrid(const std::shared_ptr<Molecule>& _mol) {
    this->mol = _mol;
    this->create_grid(GRID_ULTRAFINE);
}

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
void MolecularGrid::set_density(const MatrixXXd& P) {
    for(unsigned int i=0; i<this->grid.size(); i++) {
        this->grid[i].set_density(P);
    }
}

/**
 * @fn get_weights
 * @brief get the weights of all the grid points as a vector
 *
 * @return vector containing all the weights
 */
VectorXd MolecularGrid::get_weights() const {
    VectorXd weights = VectorXd(this->grid.size());
    for(unsigned int i=0; i<this->grid.size(); i++) {
        weights(i) = this->grid[i].get_weight();
    }

    return weights;
}

/**
 * @fn get_densities
 * @brief get the densities of all the grid points as a vector
 *
 * @return vector containing all the densities
 */
VectorXd MolecularGrid::get_densities() const {
    VectorXd densities = VectorXd(this->grid.size());
    for(unsigned int i=0; i<this->grid.size(); i++) {
        densities(i) = this->grid[i].get_density();
    }

    return densities;
}

/**
 * @fn get_amplitudes
 * @brief get the amplitudes of all the grid points and all the basis functions as a matrix
 *
 * @return matrix (basis functions x grid points)
 */
MatrixXXd MolecularGrid::get_amplitudes() const {
    MatrixXXd amplitudes = MatrixXXd(this->mol->get_nr_bfs(), this->grid.size());

    #ifdef HAS_OPENMP
    #pragma omp parallel for
    #endif
    for(unsigned int i=0; i<this->grid.size(); i++) {
        const VectorXd bfs = this->grid[i].get_basis_func_amp();
        for(unsigned int j=0; j<this->mol->get_nr_bfs(); j++) {
            amplitudes(j,i) = bfs(j);
        }
    }

    return amplitudes;
}

/**
 * @fn renormalize_density
 * @brief normalize density at gridpoints so that sum equals all electrons
 *
 * @return void
 */
void MolecularGrid::renormalize_density(unsigned int nr_elec) {
    double factor = (double)nr_elec / this->calculate_density();

    #ifdef HAS_OPENMP
    #pragma omp parallel for
    #endif
    for(unsigned int i=0; i<this->grid.size(); i++) {
        this->grid[i].scale_density(factor);
    }
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
double MolecularGrid::calculate_density() const {
    double sum = 0.0;
    #ifdef HAS_OPENMP
    #pragma omp parallel for reduction ( + : sum)
    #endif
    for(unsigned int i=0; i<this->grid.size(); i++) {
        sum += this->grid[i].get_weight() * this->grid[i].get_density();
    }

    return sum;
}

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
 * Reference: A.D. Becke, J.Chem.Phys. 88, 2547, (1988)
 * Link: http://dx.doi.org/10.1063/1.454033
 *
 * @return void
 */
void MolecularGrid::create_grid(unsigned int fineness) {

    // variables controlling the resolution of the radial and angular grid
    unsigned int radial_points;
    unsigned int lebedev_order = 0;

    // set the resolution of the grid
    switch(fineness) {
        case GRID_COARSE:
            radial_points = 10;
            lebedev_order = Quadrature::LEBEDEV_50;
        break;
        case GRID_MEDIUM:
            radial_points = 15;
            lebedev_order = Quadrature::LEBEDEV_110;
        break;
        case GRID_FINE:
            radial_points = 20;
            lebedev_order = Quadrature::LEBEDEV_146;
        break;
        case GRID_ULTRAFINE:
            radial_points = 30;
            lebedev_order = Quadrature::LEBEDEV_194;
        break;
        default: // medium settings
            radial_points = 15;
            lebedev_order = Quadrature::LEBEDEV_110;
        break;
    }

    // calculate number of angular points
    const unsigned int angular_points = Quadrature::num_lebedev_points[lebedev_order];

    // The Lebedev coefficients and radial points are stored in a large matrix.
    // Given a particular order for the integration, we have to start reading
    // this matrix from a specific position. Here, we calculate that position.
    unsigned int start = 0;
    for(unsigned int i=0; i<lebedev_order; i++) {
        start += Quadrature::num_lebedev_points[i];
    }

    // generate for each atom an atomic grid
    for(unsigned int i=0; i<this->mol->get_nr_atoms(); i++) {
        // set the position of atom 1
        const vec3 p1 = this->mol->get_atom(i).get_position();

        // construct grid points
        unsigned int grid_start = this->grid.size();
        double f = pi / (double)(radial_points + 1);

        // loop over all the radial points
        for(unsigned int p=1; p<=radial_points; p++) {
            // calculate weight function of the Gauss-Chebyshev integration
            double w = f * std::pow(std::sin(f * (double)p), 2.0);

            // calculate x point in the interval [-1,1] for Gauss-Chebyshev integration
            const double x = std::cos(f * (double)p);

            // convert x to r so that the interval is converted from [-1,1] to [0, infinity]
            const double r = (x +1.0) / (x - 1.0);

            // calculate weight function of the Gauss-Chebyshev integration
            w =  w / std::sqrt(1.0 - std::pow(x, 2.0)) *
                 2.0 / std::pow(1.0 - x, 2.0);

            // loop over all the angular points
            for(unsigned int a=start; a<angular_points + start; a++) {
                // create a new grid point given the coordinates on the unit sphere; these
                // coordinates are multiplied by the radius as calculated above
                this->grid.push_back(GridPoint(p1 + vec3(Quadrature::lebedev_coefficients[a][0],
                                                         Quadrature::lebedev_coefficients[a][1],
                                                         Quadrature::lebedev_coefficients[a][2]) * r));

                // set the atom to which the grid point adheres to (for defining the atomic grid)
                this->grid.back().set_atom(&this->mol->get_atom(i));

                // set the total integration weight by multiplying the Gauss-Chebyshev weight by the Lebedev weight
                // note that the 'weight' from the Jacobian due to spherical coordinates and the weight
                // due to Becke grid are calculated later
                this->grid.back().set_weight(w * Quadrature::lebedev_coefficients[a][3]);

                // set the amplitudes of the basis functions at this grid point
                this->grid.back().set_basis_func_amp(this->mol);
            }
        }
        unsigned int grid_stop = this->grid.size();

        // calculate the Becke weights for the atomic grids
        #ifdef HAS_OPENMP
        #pragma omp parallel for
        #endif
        for(unsigned int g=grid_start; g<grid_stop; g++) {

            double denom = 0.0;
            double nom = 1.0;

            // loop over all atoms to get Pn(r)
            for(unsigned int j=0; j<this->mol->get_nr_atoms(); j++) {
                double term = this->get_becke_weight_pn(j, this->grid[g].get_position()); // obtain single Pn(r) term
                denom += term;
                if(i == j) {
                    nom = term;
                }
            }

            // set weight from cell function: wn(r) = Pn(r) / SUM_m Pm(r) (eq. 22)
            if(denom != 0.0) {
                this->grid[g].multiply_weight(nom / denom);
            }
        }
    }

    // finally loop over all gridpoints again and
    // correct the weights for the Jacobian r**2
    #ifdef HAS_OPENMP
    #pragma omp parallel for
    #endif
    for(unsigned int i=0; i<grid.size(); i++) {
        this->grid[i].multiply_weight((this->grid[i].get_position() - this->grid[i].get_atom_position()).squaredNorm() * 4.0 * pi);
    }
}

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
double MolecularGrid::get_becke_weight_pn(unsigned int atnr, const vec3& p0) {
    double wprod = 1.0;
    const vec3 p1 = this->mol->get_atom(atnr).get_position();

    for(unsigned int j=0; j<this->mol->get_nr_atoms(); j++) {
        if(atnr == j) {
            continue;
        }

        const vec3 p2 = this->mol->get_atom(j).get_position();

        double mu =  ((p0 - p1).norm()
                    - (p0 - p2).norm() )/
                      (p2-p1).norm();

        wprod *= this->cutoff(mu);
    }

    return wprod;
}

/**
 * @fn cutoff
 * @brief calculates the Becke weight
 *
 * @param fineness      mu (elliptical coordinate)
 *
 * @return double Becke weight value
 */
double MolecularGrid::cutoff(double mu) {
    return 0.5 * (1.0 - this->fk(3.0, mu));
}

/**
 * @fn fk
 * @brief recursive function to have smooth overlapping atomic grids
 *
 * @param k     number of iterations
 * @param mu    elliptical coordinate
 *
 * @return value
 */
double MolecularGrid::fk(unsigned int k, double mu) {
    for(unsigned int i=0; i<k; i++) {
        mu = 3.0/2.0 * mu - 0.5 * std::pow(mu, 3.0);
    }

    return mu;
}
