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

#include "dft.h"

/**
 * @brief      default constructor
 *
 * @param[in]  _mol       molecule
 * @param[in]  _settings  settings
 */
DFT::DFT(const std::shared_ptr<Molecule>& _mol, const std::shared_ptr<Settings>& _settings) :
settings(_settings),
flag_has_energy(false),
flag_has_forces(false) {
    is_first = true;
    this->add_molecule(_mol);
    this->integrator = std::make_unique<Integrator>();
    this->functional = std::make_unique<Functional>();
}

/**
 * @fn add_molecule
 * @brief add molecule to the DFT routine
 *
 * @param _mol  pointer to Molecule class
 *
 * Links molecule to the DFT routine by a pointer. The
 * Contracted Gaussian Functions of the molecule are copied
 * to the DFT routine class and the routine matrices are
 * constructed.
 *
 * @return void
 */
void DFT::add_molecule(const std::shared_ptr<Molecule>& _mol) {
    // set molecule
    this->mol = _mol;
    this->nelec = this->mol->get_nr_elec();

    // construct molecular grid based on molecule and basis functions
    this->molgrid = std::make_unique<MolecularGrid>(this->mol);
    this->molgrid->set_grid_parameters(this->settings->get_radial_points(),
                                       this->settings->get_lebedev_order(),
                                       this->settings->get_lmax());
    this->molgrid->create_grid();

    // copy the contracted gaussian functions to the DFT class
    this->copy_cgfs_from_molecule();

    std::cout << "Loading molecule and constructing matrices." << std::endl;
    auto start = std::chrono::system_clock::now();

    // construct all matrices
    this->construct_matrices();

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << boost::format("Total time: %f ms\n") % elapsed.count();
    std::cout << std::endl;
}

/**
 * @fn scf
 * @brief perform the self-consistent field iterations
 *
 * @return void
 */
void DFT::scf(bool verbose) {
    if(this->flag_has_energy) {
        return;
    }

    static const unsigned int max_iterations = 100;

    // print coordinates
    this->mol->print_geometry();

    double old_energy = this->et;
    double difference = 1.0;
    unsigned int iteration = 0;

    // print header
    std::cout << " Performing electronic structure calculation " << std::endl;
    std::cout << "=============================================" << std::endl;

    while(difference > 1e-5 || iteration < 3) {
        iteration++;

        auto start = std::chrono::system_clock::now(); //tic

        this->calculate_density_matrix();
        this->calculate_electronic_repulsion_matrix();
        this->calculate_exchange_correlation_matrix();
        this->calculate_energy();

        auto end = std::chrono::system_clock::now(); //toc
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        difference = std::abs(this->et - old_energy);
        old_energy = this->et;

        // detailed output if verbose is enabled
        if(verbose) {
            std::cout << boost::format("%3i    %9.7f    %4.2f (%3i) \n")
                                % iteration
                                % this->et
                                % this->molgrid->calculate_density()
                                % this->nelec;

            std::cout << boost::format("\tE_XC \t= %9.7f\n\tE_NUC \t= %9.7f\n\tE_ONE \t= %9.7f\n\tE_J \t= %9.7f\n\tt \t=%9.7f ms\n")
                                % this->exc
                                % this->enuc
                                % this->single_electron_energy
                                % this->electronic_repulsion
                                % (elapsed.count());

            std::cout << "----------------------------------------" << std::endl;
        } else {    // regular output
            std::cout << boost::format("%4i | %12.8f | %12.8f | %4.1f ms")
                                % iteration
                                % this->et
                                % difference
                                % (elapsed.count());
            std::cout << std::endl;
        }

        if(iteration >= max_iterations) {
            std::cout << "=============================================" << std::endl;
            std::cout << "Stopping because maximum number of iterations has been reached." << std::endl;
            std::cout << std::endl;
            break;
        }
    }

    if(iteration < max_iterations) {
        std::cout << "=============================================" << std::endl;
        std::cout << "Stopping because energy criterion is reached." << std::endl;
        std::cout << std::endl;
    }

    this->flag_has_energy = true;
}

/**
 * @brief      get the total energy of the molecule
 *
 * @return     The energy.
 */
double DFT::get_energy() {
    if(!this->flag_has_energy) {
        this->scf();
    }

    return this->et;
}

/**
 * @brief      Pre-initialize a wave function
 *
 * @param[in]  _P    density matrix P
 */
void DFT::set_wavefunction(const MatrixXXd& _P) {
    // set the new density
    this->P = _P;

    // calculate properties based on the density
    this->molgrid->set_density(this->P);
    this->molgrid->correct_densities();
    this->calculate_electronic_repulsion_matrix();
    this->calculate_exchange_correlation_matrix();
}

/**
 * @brief      get the force vector
 */
VectorXd DFT::get_force_vector() {
    if(!this->flag_has_energy) {
        this->scf();
    }

    if(!this->flag_has_forces) {
        this->calculate_hellmann_feynman_forces();
    }

    auto Ft = F;
    Ft.transposeInPlace();
    return Eigen::Map<VectorXd>(Ft.data(), Ft.cols() * Ft.rows());
}

/**
 * @fn copy_cgfs_from_molecule
 * @brief copy cgfs from the molecule class to the dft class
 *
 * @return void
 */
void DFT::copy_cgfs_from_molecule() {
    this->cgfs = this->mol->get_cgfs();
}

/**
 * @fn construct_matrices
 * @brief construct the matrices and two-electron integrals for the SCF procedure
 *
 * @return void
 */
void DFT::construct_matrices() {
    const unsigned int size = this->mol->get_nr_bfs();

    // overlap matrix
    this->S = MatrixXXd::Zero(size, size);

    // kinetic energy matrix
    this->T = MatrixXXd::Zero(size, size);

    // nuclear attraction matrix
    this->V = MatrixXXd::Zero(size, size);

    // single-electron matrix
    this->H = MatrixXXd::Zero(size, size);

    // coulombic repulsion matrix
    this->J = MatrixXXd::Zero(size, size);

    // exchange-correlation matrix
    this->XC = MatrixXXd::Zero(size, size);

    // density matrix
    this->P = MatrixXXd::Zero(size, size);

    // calculate values of matrix elements
    #pragma omp parallel for schedule(dynamic)
    for(unsigned int i=0; i<this->cgfs->size(); i++) {
        for(unsigned int j=i; j<this->cgfs->size(); j++) {
            S(i,j) = S(j,i) = this->integrator->overlap(cgfs->at(i), cgfs->at(j));
            T(i,j) = T(j,i) = this->integrator->kinetic(cgfs->at(i), cgfs->at(j));
            double v_sum = 0.0;
            for(unsigned int k=0; k<this->mol->get_nr_atoms(); k++) {
                 v_sum += this->integrator->nuclear(cgfs->at(i), cgfs->at(j),
                                    this->mol->get_atomic_position(k),
                                    this->mol->get_atomic_charge(k));
            }
            V(i,j) = V(j,i) = v_sum;
        }
    }

    // calculate single-electron matrix by summing T and V matrices
    this->H = this->T + this->V;

    // calculate nuclear repulsion energy
    this->calculate_nuclear_repulsion();

    // calculate all two-electron integrals
    if(this->settings->get_hartree_evaluation_method() == Settings::TWO_ELECTRON_INTEGRALS) {
        auto start = std::chrono::system_clock::now(); //tic
        std::cout << "Explicitly calculating two electron integrals... " << std::flush;
        this->calculate_two_electron_integrals();
        auto end = std::chrono::system_clock::now(); //toc
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << boost::format("(%f ms)\n") % elapsed.count();
    }

    // calculate transformation matrix from the overlap matrix
    this->calculate_transformation_matrix();

    // calculate the density matrix from the coefficient matrices
    this->calculate_density_matrix();

    // calculate the electronic repulsion matrix from the density
    this->calculate_electronic_repulsion_matrix();

    // calculate the total energy of the molecule
    this->calculate_energy();
}

/**
 * @fn calculate_nuclear_repulsion
 * @brief calculate the nuclear repulsion energy
 *
 * @return void
 */
void DFT::calculate_nuclear_repulsion() {
    this->enuc = 0.0;

    for(unsigned int i=0; i<this->mol->get_nr_atoms(); i++) {
        for(unsigned int j=i+1; j<this->mol->get_nr_atoms(); j++) {
            this->enuc += this->mol->get_atomic_charge(i) *
                          this->mol->get_atomic_charge(j) /
                          (this->mol->get_atomic_position(i) -
                          this->mol->get_atomic_position(j)).norm();
        }
    }
}

/**
 * @fn calculate_two_electron_integrals
 * @brief calculate the two-electron integrals
 *
 * @return void
 */
void DFT::calculate_two_electron_integrals() {
    const unsigned int size = this->mol->get_nr_bfs();

    unsigned int nrints = this->integrator->teindex(size,size,size,size);
    this->ints = VectorXd::Ones(nrints);
    this->ints *= -1.0;

    #ifdef HAS_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for(unsigned int i=0; i<size; i++) {
        for(unsigned int j=0; j<size; j++) {
            unsigned int ij = i*(i+1)/2 + j;
            for(unsigned int k=0; k<size; k++) {
                for(unsigned int l=0; l<size; l++) {
                    unsigned int kl = k * (k+1)/2 + l;
                    if(ij <= kl) {
                        unsigned int idx = this->integrator->teindex(i,j,k,l);

                        // this avoids recalculating an integral which has
                        // already been evaluated
                        if(this->ints(idx) != -1.0) {
                            continue;
                        }

                        this->ints(idx) = this->integrator->repulsion(this->cgfs->at(i),
                                                                      this->cgfs->at(j),
                                                                      this->cgfs->at(k),
                                                                      this->cgfs->at(l));
                    }
                }
            }
        }
    }
}

/**
 * @fn calculate_transformation_matrix
 * @brief calculate the transformation matrix X from the overlap matrix S
 *
 * The transformation matrix X is calculated from the overlap matrix S by
 * canonical orthogonalization. This transformation matrix X ensures that
 * the Slater determinant consists of orthogonal spin-orbitals.
 *
 * @return void
 */
void DFT::calculate_transformation_matrix() {
    const unsigned int size = this->mol->get_nr_bfs(); // nr of cgfs

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
    Eigen::MatrixXd D = es.eigenvalues().real().asDiagonal();
    for(unsigned int i=0; i<size; i++) {
        D(i,i) = 1.0 / std::sqrt(D(i,i));
    }
    Eigen::MatrixXd U = es.eigenvectors().real();

    // Calculate the transformation matrix
    X = MatrixXXd::Zero(size, size);
    Xp = MatrixXXd::Zero(size, size);

    this->X = U * D;
    this->Xp = X.transpose();
}

/**
 * @fn calculate_density_matrix
 * @brief calculate the density matrix from the coefficient matrix
 *
 * Performs eigenvector decomposition on the Fock matrix F. From these
 * the coefficients in the transformed basis are calculated. The basis
 * function coefficients in the original basis can be obtained by
 * using the transformation matrix X. From these coefficients in the
 * original representation, the density matrix P is calculated.
 *
 * @return void
 */
void DFT::calculate_density_matrix() {
    static const double alpha = 0.50; // mixing parameter alpha (NOTE: obtain this value from input file...)
    const unsigned int size = this->mol->get_nr_bfs(); // nr of cgfs

    MatrixXXd F = this->H + 2.0 * this->J + this->XC;

    MatrixXXd Fp = this->Xp * F * this->X;

    // compute eigenvectors and eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Fp);
    this->Cc = es.eigenvectors().real();

    // obtain true coefficient matrix using the transformation matrix
    this->C = this->X * this->Cc;

    MatrixXXd Pnew  = MatrixXXd::Zero(size, size);
    for(unsigned int i=0; i<size; i++) {
        for(unsigned int j=0; j<size; j++) {
            for(unsigned int k=0; k < this->nelec / 2; k++) {
                Pnew(i,j) += C(i,k) * C(j,k);
            }
        }
    }

    if(is_first) {
        this->P = Pnew;
        is_first = false;
    } else {
        this->P = (1.0 - alpha) * Pnew + alpha * P;
    }

    // set the new density
    this->molgrid->set_density(P);

    // correct density
    this->molgrid->correct_densities();
}

/**
 * @fn calculate_electronic_repulsion_matrix
 * @brief calculates electronic repulsion matrix
 *
 * Calculate the electronic repulsion matrix using the
 * density matrix and the two-electron integrals.
 *
 * @return void
 */
void DFT::calculate_electronic_repulsion_matrix() {
    switch(this->settings->get_hartree_evaluation_method()) {
        case Settings::BECKE_GRID:
            this->calculate_hartree_potential_becke_poisson();
        break;
        case Settings::TWO_ELECTRON_INTEGRALS:
            this->calculate_hartree_potential_te_int();
        break;
    }
}

/**
 * @fn calculate_exchange_correlation_matrix
 * @brief calculates the exchange-correlation matrix
 *
 * @return void
 */
void DFT::calculate_exchange_correlation_matrix() {
    VectorXd  densities   = this->molgrid->get_densities();
    VectorXd  weights     = this->molgrid->get_weights();
    MatrixXXd amplitudes  = this->molgrid->get_amplitudes();

    VectorXd ex;    // exchange energy
    VectorXd vxa;   // exchange potential (spin-up)
    VectorXd vxb;   // exchange potential (spin-down)

    VectorXd ec;    // correlation energy
    VectorXd vca;   // correlation potential (spin-up)
    VectorXd vcb;   // correlation potential (spin-down)

    VectorXd densitiesa = densities * 0.5;
    VectorXd densitiesb = densities * 0.5;

    // calculate exchange energies and potential
    this->functional->xalpha_x_functional(densitiesa, densitiesb, ex, vxa, vxb);

    // calculate correlation energy and potential
    this->functional->vwm_c_functional(densitiesa, densitiesb, ec, vca, vcb);

    // calculate total exchange-correlation energy
    this->exc = weights.dot(ex + ec);

    // calculate gridpoint-wise xc potential

    // I still don't understand why I need to put a 0.5 here... (see also the function above)
    // Does it have anything to do how I construct the density matrix P?
    // Can someone who knows the answer send me an e-mail?
    VectorXd wva = weights.cwiseProduct((vxa + vxb + vca + vcb) * 0.5);

    for(unsigned int i=0; i<this->cgfs->size(); i++) {
        VectorXd row = amplitudes.row(i);
        VectorXd wva_i = wva.cwiseProduct(row);
        for(unsigned int j=0; j<this->cgfs->size(); j++) {
            this->XC(i,j) = wva_i.dot(amplitudes.row(j));
        }
    }
}

/**
 * @fn calculate_energy
 * @brief calculate the energy of the molecule
 *
 * @return void
 */
void DFT::calculate_energy() {
    this->single_electron_energy = 2.0 * (this->P * this->H).trace();
    this->electronic_repulsion = 2.0 * (this->P * this->J).trace();

    // sum all terms
    this->et = this->single_electron_energy + this->electronic_repulsion + this->enuc + this->exc;
}

/**
 * @brief      Calculates the hartree potential from two electrons integrals
 */
void DFT::calculate_hartree_potential_te_int() {
    const unsigned int size = this->mol->get_nr_bfs();

    // Populate two-electron repulsion matrix
    #pragma omp parallel for collapse(2)
    for(unsigned int i=0; i<size; i++) {
        for(unsigned int j=0; j<size; j++) {
            this->J(i,j) = 0.; /* reset J matrix */
            for(unsigned int k=0; k<size; k++) {
                for(unsigned int l=0; l<size; l++) {
                    const unsigned int index = this->integrator->teindex(i,j,k,l);
                    this->J(i,j) += P(k,l) * ints(index);
                }
            }
        }
    }
}

/**
 * @brief      Calculates the hartree potential over Becke grid using Poisson equation
 */
void DFT::calculate_hartree_potential_becke_poisson() {
    this->J = this->molgrid->calculate_hartree_potential();
}

/**
 * @brief      Calculate the hellmann feynman forces.
 */
void DFT::calculate_hellmann_feynman_forces() {
    if(this->flag_has_forces) {
        return;
    }

    this->F = this->molgrid->get_forces_atoms();

    this->flag_has_forces = true;
}
