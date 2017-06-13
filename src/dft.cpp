/**************************************************************************
 *   dft.cpp  --  This file is part of DFTCXX.                            *
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

#include "dft.h"

/**
 * @fn DFT
 * @brief DFT routine constructor
 *
 * @return DFT instance
 */
DFT::DFT() {
    is_first = true;
};

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
void DFT::add_molecule(const Molecule* _mol) {
    // set molecule
    this->mol = _mol;
    this->nelec = this->mol->get_nr_elec();

    // construct molecular grid based on molecule and basis functions
    this->molgrid = new MolecularGrid(this->mol);

    // copy the contracted gaussian functions to the DFT class
    this->copy_cgfs_from_molecule();

    // construct all matrices
    this->construct_matrices();
}

/**
 * @fn scf
 * @brief perform the self-consistent field iterations
 *
 * @return void
 */
void DFT::scf() {
    std::cout << "          Starting calculation          " << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "  #        energy    elec" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    double old_energy = this->et;
    double difference = 1.0;
    unsigned int iteration = 0;

    while(difference > 1e-4) {
        iteration++;

        this->calculate_density_matrix();
        this->calculate_electronic_repulsion_matrix();
        this->calculate_exchange_correlation_matrix();
        this->calculate_energy();

        std::cout << boost::format("%3i    %9.7f    %4.2f (%3i) \n")
                            % iteration
                            % et
                            % this->molgrid->calculate_density()
                            % this->nelec;

        difference = std::abs(this->et - old_energy);
        old_energy = this->et;
    }

    std::cout << "========================================" << std::endl;
    std::cout << "Stopping because energy criterion is reached." << std::endl;
    std::cout << std::endl;
}

/**
 * @fn copy_cgfs_from_molecule
 * @brief copy cgfs from the molecule class to the dft class
 *
 * @return void
 */
void DFT::copy_cgfs_from_molecule() {
    for(unsigned int i=0; i<this->mol->get_nr_bfs(); i++) {
        this->cgfs.push_back(&this->mol->get_cgf(i));
    }
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
    this->S = MatrixXXd(size, size);

    // kinetic energy matrix
    this->T = MatrixXXd(size, size);

    // nuclear attraction matrix
    this->V = MatrixXXd(size, size);

    // single-electron matrix
    this->H = MatrixXXd(size, size);

    // coulombic repulsion matrix
    this->J = MatrixXXd(size, size);

    // exchange-correlation matrix
    this->XC = MatrixXXd::Zero(size, size);

    // density matrix
    this->P = MatrixXXd::Zero(size, size);

    // calculate values of matrix elements
    for(unsigned int i=0; i<this->cgfs.size(); i++) {
        for(unsigned int j=i; j<this->cgfs.size(); j++) {
            S(i,j) = S(j,i) = this->integrator->overlap(*cgfs[i], *cgfs[j]);
            T(i,j) = T(j,i) = this->integrator->kinetic(*cgfs[i], *cgfs[j]);
            double v_sum = 0.0;
            for(unsigned int k=0; k<this->mol->get_nr_atoms(); k++) {
                 v_sum += this->integrator->nuclear(*cgfs[i], *cgfs[j],
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
    this->calculate_two_electron_integrals();

    // calculate transformation matrix from the overlap matrix
    this->calculate_transformation_matrix();

    // calculate the density matrix from the coefficient matrices
    this->calculate_density_matrix();

    // calculate the electronic repulsion matrix from the density
    // matrix and the two-electron integrals
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

    unsigned int nrints = integrator->teindex(size,size,size,size);
    this->ints = VectorXd::Ones(nrints);
    this->ints *= -1.0;
    for(unsigned int i=0; i<size; i++) {
        for(unsigned int j=0; j<size; j++) {
            unsigned int ij = i*(i+1)/2 + j;
            for(unsigned int k=0; k<size; k++) {
                for(unsigned int l=0; l<size; l++) {
                    unsigned int kl = k * (k+1)/2 + l;
                    if(ij <= kl) {
                        unsigned int idx = integrator->teindex(i,j,k,l);

                        // this avoids recalculating an integral which has
                        // already been evaluated
                        if(this->ints(idx) != -1.0) {
                            continue;
                        }

                        this->ints(idx) = integrator->repulsion(*this->cgfs[i],
                                                                *this->cgfs[j],
                                                                *this->cgfs[k],
                                                                *this->cgfs[l]);
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

    Eigen::EigenSolver<MatrixXXd> es(this->S, true);
    MatrixXXd D = MatrixXXd::Zero(size, size);
    for(unsigned int i=0; i<size; i++) {
        D(i,i) = 1.0 / std::sqrt(es.eigenvalues()(i, 0).real() );
    }

    MatrixXXd U = es.eigenvectors().real();

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
    static const double alpha = 0.5; // mixing parameter alpha
    const unsigned int size = this->mol->get_nr_bfs(); // nr of cgfs

    MatrixXXd F = this->T + this->V + 2.0 * this->J + this->XC;

    MatrixXXd Fp = this->Xp * F * this->X;

    // compute eigenvectors and eigenvalues
    Eigen::EigenSolver<MatrixXXd> es(Fp, true);
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

    this->molgrid->set_density(P);
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
    const unsigned int size = this->mol->get_nr_bfs();

    // Populate two-electron repulsion matrix
    for(unsigned int i=0; i<size; i++) {
        for(unsigned int j=0; j<size; j++) {
            this->J(i,j) = 0.; /* reset G matrix */
            for(unsigned int k=0; k<size; k++) {
                for(unsigned int l=0; l<size; l++) {
                    unsigned int index = integrator->teindex(i,j,k,l);
                    this->J(i,j) += P(k,l) * ints(index);
                }
            }
        }
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

    VectorXd ex;
    VectorXd vx;

    VectorXd ec;
    VectorXd vca;
    VectorXd vcb;

    VectorXd densitiesa = densities / 2.0;
    VectorXd densitiesb = densities / 2.0;

    // calculate exchange energies and potential
    this->functional->xalpha_x_functional(densitiesa, densitiesb, ex, vx);

    // calculate correlation energy and potential
    this->functional->vwm_c_functional(densitiesa, densitiesb, ec, vca, vcb);

    // calculate total exchange-correlation energy
    this->exc = weights.dot(ex + ec);

    // calculate gridpoint-wise xc potential
    VectorXd wva = weights.cwiseProduct(vx + vca);

    for(unsigned int i=0; i<this->cgfs.size(); i++) {
        VectorXd row = amplitudes.row(i);
        VectorXd wva_i = wva.cwiseProduct(row);
        for(unsigned int j=i; j<this->cgfs.size(); j++) {
            XC(i,j) = XC(j,i) = wva_i.dot(amplitudes.row(j));
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
    double single_electron_energy = 2.0 * (this->P * this->H).trace();
    double electronic_repulsion = 2.0 * (this->P * this->J).trace();

    // sum all terms
    this->et = single_electron_energy + electronic_repulsion + this->enuc + this->exc;
}
