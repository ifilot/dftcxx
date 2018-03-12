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

#include "moleculargrid.h"

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
    this->grid_size = 0;
    this->create_grid(GRID_COARSE);
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
    for(unsigned int i=0; i<this->atomic_grids.size(); i++) {
        this->atomic_grids[i]->set_density(P);
    }
}

/**
 * @fn get_weights
 * @brief get the weights of all the grid points as a vector
 *
 * @return vector containing all the weights
 */
VectorXd MolecularGrid::get_weights() const {
    VectorXd weights = VectorXd::Zero(this->grid_size);
    unsigned int idx = 0;

    for(unsigned int i=0; i<this->atomic_grids.size(); i++) {
        VectorXd vec_weight = this->atomic_grids[i]->get_weights();
        for(unsigned int j=0; j<vec_weight.size(); j++) {
            weights(idx) = vec_weight(j);
            idx++;
        }
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
    VectorXd densities = VectorXd::Zero(this->grid_size);
    unsigned int idx = 0;

    for(unsigned int i=0; i<this->atomic_grids.size(); i++) {
        VectorXd vec_densities = this->atomic_grids[i]->get_densities();
        for(unsigned int j=0; j<vec_densities.size(); j++) {
            densities(idx) = vec_densities(j);
            idx++;
        }
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
    MatrixXXd amplitudes = MatrixXXd::Zero(this->mol->get_nr_bfs(), this->grid_size);
    unsigned int idx = 0;

    for(unsigned int i=0; i<this->atomic_grids.size(); i++) {
        MatrixXXd agrid_amplitudes = this->atomic_grids[i]->get_amplitudes();

        for(unsigned int j = 0; j<this->atomic_grids[i]->get_grid_size(); j++) {
            for(unsigned int k=0; k<this->mol->get_nr_bfs(); k++) {
                amplitudes(k,j) = agrid_amplitudes(k,j);
            }
            idx++;
        }
    }

    return amplitudes;
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
    for(unsigned int i=0; i<this->atomic_grids.size(); i++) {
        sum += this->atomic_grids[i]->calculate_density();
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
 * @return void
 */
void MolecularGrid::create_grid(unsigned int fineness) {

    // set the resolution of the grid
    switch(fineness) {
        case GRID_COARSE:
            this->radial_points = 10;
            this->lebedev_order = Quadrature::LEBEDEV_50;
            this->lmax = 5;
        break;
        case GRID_MEDIUM:
            this->radial_points = 15;
            this->lebedev_order = Quadrature::LEBEDEV_110;
            this->lmax = 8;
        break;
        case GRID_FINE:
            this->radial_points = 20;
            this->lebedev_order = Quadrature::LEBEDEV_146;
            this->lmax = 10;
        break;
        case GRID_ULTRAFINE:
            this->radial_points = 30;
            this->lebedev_order = Quadrature::LEBEDEV_194;
            this->lmax = 11;
        break;
        default: // medium settings
            this->radial_points = 15;
            this->lebedev_order = Quadrature::LEBEDEV_110;
            this->lmax = 8;
        break;
    }

    // calculate number of angular points
    this->angular_points = Quadrature::num_lebedev_points[this->lebedev_order];

    // The Lebedev coefficients and radial points are stored in a large matrix.
    // Given a particular order for the integration, we have to start reading
    // this matrix from a specific position. Here, we calculate that position.
    lebedev_offset = 0;
    for(unsigned int i=0; i<this->lebedev_order; i++) {
        lebedev_offset += Quadrature::num_lebedev_points[i];
    }

    // generate for each atom an atomic grid
    for(unsigned int i=0; i<this->mol->get_nr_atoms(); i++) {
        // add atomic density
        this->atomic_densities.push_back(0);

        this->atomic_grids.emplace_back(new AtomicGrid(this->mol->get_atom(i)));
        this->atomic_grids.back()->create_atomic_grid(radial_points,
                                                      angular_points,
                                                      lebedev_offset,
                                                      lmax,
                                                      mol);
        this->grid_size += this->atomic_grids.back()->get_grid_size();
    }

    // calculate the Becke weights for the atomic grids
    for(unsigned int i=0; i<this->mol->get_nr_atoms(); i++) {

        MatrixXXd positions = this->atomic_grids[i]->get_positions();
        VectorXd becke_coeff = VectorXd::Zero(positions.rows());

        for(unsigned int j=0; j<positions.rows(); j++) {

            double denom = 0.0;
            double nom = 1.0;

            // loop over all atoms to get Pn(r)
            for(unsigned int k=0; k<this->mol->get_nr_atoms(); k++) {
                double term = this->get_becke_weight_pn(k, positions.row(j)); // obtain single Pn(r) term
                denom += term;
                if(i == k) {
                    nom = term;
                }
            }

            // set weight from cell function: wn(r) = Pn(r) / SUM_m Pm(r) (eq. 22)
            if(denom != 0.0) {
                becke_coeff(j) = nom / denom;
            } else {
                becke_coeff(j) = 1.0;
            }
        }

        this->atomic_grids[i]->correct_weights(becke_coeff);
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

    // get position of central atom in fuzzy cell
    const vec3 p1 = this->mol->get_atom(atnr)->get_position();

    for(unsigned int j=0; j<this->mol->get_nr_atoms(); j++) {
        if(atnr == j) {
            continue;
        }

        // get position of other atom
        const vec3 p2 = this->mol->get_atom(j)->get_position();

        // calculate mu (eq. 11 in A.D. Becke J.Chem.Phys. 88, 2547)
        double mu =  ((p0 - p1).norm()
                    - (p0 - p2).norm() )/
                      (p2-p1).norm();

        // (eq. 13 in A.D. Becke J.Chem.Phys. 88, 2547)
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
    return 0.5 * (1.0 - this->fk(3, mu));
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

void MolecularGrid::calculate_hartree_potential() {
    MatrixXXd J = MatrixXXd::Zero(this->mol->get_nr_bfs(), this->mol->get_nr_bfs());

    for(unsigned int i=0; i<this->atomic_grids.size(); i++) {
        std::cout << "Atomic density on " << (i+1) << ": " << this->atomic_grids[i]->calculate_density() << std::endl;

        this->atomic_grids[i]->calculate_rho_lm();
        this->atomic_grids[i]->calculate_U_lm();

        J += this->atomic_grids[i]->get_J();
    }

    std::cout << J << std::endl;
}
