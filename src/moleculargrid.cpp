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
    #pragma omp parallel for
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
    unsigned int idx =0;

    for(unsigned int i=0; i<this->atomic_grids.size(); i++) {
        VectorXd vec_weight = this->atomic_grids[i]->get_weights();

        #pragma omp parallel for
        for(unsigned int j=0; j<vec_weight.size(); j++) {
            weights(idx + j) = vec_weight(j);
        }

        idx += vec_weight.size();
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

        #pragma omp parallel for
        for(unsigned int j=0; j<vec_densities.size(); j++) {
            densities(idx + j) = vec_densities(j);
        }

        idx += vec_densities.size();
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

        #pragma omp parallel for
        for(unsigned int j = 0; j<this->atomic_grids[i]->get_grid_size(); j++) {
            for(unsigned int k=0; k<this->mol->get_nr_bfs(); k++) {
                amplitudes(k, idx + j) = agrid_amplitudes(k,j);
            }
        }

        idx += this->atomic_grids[i]->get_grid_size();
    }

    return amplitudes;
}

/**
 * @brief      correct the densities to match the total number of electrons
 */
void MolecularGrid::correct_densities() {
    double sum = 0.0;
    double charge = 0.0;

    for(unsigned int i=0; i<this->atomic_grids.size(); i++) {
        sum += this->atomic_grids[i]->get_density();
        charge += this->mol->get_atomic_charge(i);
    }

    double correction_factor = charge / sum;

    for(unsigned int i=0; i<this->atomic_grids.size(); i++) {
        this->atomic_grids[i]->correct_density(correction_factor);
    }
}

/**
 * @brief      calculate the Hellmann-Feynman forces
 *
 * @return     atomic forces
 */
MatrixXXd MolecularGrid::get_forces_atoms() const {
    MatrixXXd ans = MatrixXXd::Zero(this->mol->get_nr_atoms(), 3);

    for(unsigned int i=0; i<this->mol->get_nr_atoms(); i++) {

        for(unsigned int j=0; j<this->mol->get_nr_atoms(); j++) {

            // calculate attraction of electron densities
            const vec3 f = this->atomic_grids[j]->get_force(this->mol->get_atom(i)->get_position(), this->mol->get_atom(i)->get_charge());
            ans(i,0) += f[0];
            ans(i,1) += f[1];
            ans(i,2) += f[2];

            // add repulsion of other atoms
            if(i != j) {
                double strength = this->mol->get_atom(i)->get_charge() * this->mol->get_atom(j)->get_charge();
                vec3 direction = this->mol->get_atom(i)->get_position() - this->mol->get_atom(j)->get_position();
                strength /= std::pow(direction.squaredNorm(), 1.5);

                ans(i,0) -= strength * direction[0];
                ans(i,1) -= strength * direction[1];
                ans(i,2) -= strength * direction[2];
            }
        }
    }

    return ans;
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

    #pragma omp parallel for reduction(+:sum)
    for(unsigned int i=0; i<this->atomic_grids.size(); i++) {
        sum += this->atomic_grids[i]->get_density();
    }

    return sum;
}

/**
 * @brief      Sets the grid parameters (fine-tuning)
 *
 * @param[in]  _radial_points  number of radial points
 * @param[in]  _lebedev_order  The lebedev order
 * @param[in]  _lmax           maximum angular momentum for spherical harmonics
 */
void MolecularGrid::set_grid_parameters(unsigned int _radial_points, unsigned int _lebedev_order, unsigned int _lmax) {
    this->radial_points = _radial_points;
    this->lebedev_order = _lebedev_order;
    this->lmax = _lmax;
}

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
void MolecularGrid::create_grid() {

    auto start = std::chrono::system_clock::now(); //tic

    std::cout << "      Constructing molecular grid       " << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Number of radial points: " << this->radial_points << std::endl;
    std::cout << "Lebedev order: " << this->lebedev_order << std::endl;
    std::cout << "Lmax value: " << this->lmax << std::endl;

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
        VectorXd becke_weight_coeff = VectorXd::Zero(positions.rows());

        // loop over positions
        for(unsigned int j=0; j<positions.rows(); j++) {

            double denom = 0.0;
            double nom = 1.0;

            // loop over all atoms to get Pn(r)
            for(unsigned int k=0; k<this->mol->get_nr_atoms(); k++) {
                const double term = this->get_becke_weight_pn(k, positions.row(j)); // obtain single Pn(r) term
                denom += term;
                if(i == k) {
                    nom = term;
                }
            }

            // set weight from cell function: wn(r) = Pn(r) / SUM_m Pm(r) (eq. 22)
            becke_weight_coeff(j) = nom / denom;
        }

        this->atomic_grids[i]->set_becke_weights(becke_weight_coeff);
    }

    auto end = std::chrono::system_clock::now(); //toc
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << boost::format("Total time: %f ms\n") % elapsed.count();
    std::cout << "========================================" << std::endl;
    std::cout << std::endl;
}

/**
 * @fn get_becke_weight_pn
 * @brief calculate Pi(r) (i.e. for a single gridpoint) for atom n (eq. 13 in A.D. Becke J.Chem.Phys. 88, 2547)
 *
 * Reference: A.D. Becke, J.Chem.Phys. 88, 2547, (1988)
 * Link: http://dx.doi.org/10.1063/1.454033
 *
 * @param i    atom id i
 * @param p0   grid point position
 *
 * @return double Becke weight value
 */
double MolecularGrid::get_becke_weight_pn(unsigned int i, const vec3& p0) {
    double wprod = 1.0;

    // get position of central atom in fuzzy cell
    const vec3 p1 = this->mol->get_atom(i)->get_position();

    // loop over all other atoms
    for(unsigned int j=0; j<this->mol->get_nr_atoms(); j++) {
        if(i == j) {
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

/**
 * @brief      Calculates the hartree potential.
 *
 * @return     The hartree potential.
 */
MatrixXXd MolecularGrid::calculate_hartree_potential() {
    for(unsigned int i=0; i<this->atomic_grids.size(); i++) {
        this->atomic_grids[i]->calculate_rho_lm();
        this->atomic_grids[i]->calculate_U_lm();
    }

    for(unsigned int i=0; i<this->atomic_grids.size(); i++) {

        #pragma omp parallel for
        for(unsigned int j=0; j<this->atomic_grids[i]->get_grid_size(); j++) {

            double V = 0.0;

            for(unsigned int k=0; k<this->atomic_grids.size(); k++) {
                if(i == k) {    // obtain Vn for fuzzy cell
                    V += this->atomic_grids[k]->get_V(j);
                } else {        // interpolatie Vn

                    // set lm index
                    unsigned int lm = 0;

                    // loop over l and m
                    for(int l=0; l<=lmax; l++) {
                        for(int m=-l; m<=l; m++) {

                            const double pre = SH::prefactor_spherical_harmonic(l, m);
                            const vec3 pos = this->atomic_grids[i]->get_position_grid_point(j) - this->atomic_grids[k]->get_position();

                            // convert to spherical coordinates
                            const double r = pos.norm(); // due to unit sphere
                            const double azimuth = std::atan2(pos(1),pos(0));
                            const double pole = std::acos(pos(2) / r); // due to unit sphere else z/r

                            const double y_lm = pre * SH::spherical_harmonic(l, m, pole, azimuth);
                            V += 1.0 / r * y_lm * this->atomic_grids[k]->get_sh_value(r, lm);

                            lm++;
                        }
                    }
                }
            }

            this->atomic_grids[i]->set_V_gp(j, V);
        }
    }

    MatrixXXd J = MatrixXXd::Zero(this->mol->get_nr_bfs(), this->mol->get_nr_bfs());
    for(unsigned int i=0; i<this->atomic_grids.size(); i++) {
        this->atomic_grids[i]->calculate_J_matrix();
        J += this->atomic_grids[i]->get_J();
    }

    return J;
}
