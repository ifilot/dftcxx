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
void GridPoint::set_basis_func_amp(std::shared_ptr<Molecule> _mol) {
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
    this->density = 2.0 * this->basis_func_amp.dot(
                        D * this->basis_func_amp);
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

    for(unsigned int i=0; i<this->grid.size(); i++) {
        const VectorXd bfs = this->grid[i].get_basis_func_amp();
        for(unsigned int j=0; j<this->mol->get_nr_bfs(); j++) {
            amplitudes(j,i) = bfs(j);
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

        // set the position of atom 1
        const vec3 p1 = this->mol->get_atom(i).get_position();

        // construct grid points
        unsigned int grid_start = this->grid.size();
        const double f = pi / (double)(this->radial_points + 1);

        // loop over all the radial points
        for(unsigned int p=1; p<=this->radial_points; p++) {
            // calculate weight function of the Gauss-Chebyshev integration
            double w = f * std::pow(std::sin(f * (double)p), 2.0);

            // calculate x point in the interval [-1,1] for Gauss-Chebyshev integration
            double x = std::cos(f * (double)p);

            // convert x to r so that the interval is converted from [-1,1] to [0, infinity]
            const double r = (1.0 + x) / (1.0 - x);

            // calculate weight function of the Gauss-Chebyshev integration
            w =  w / std::sqrt(1.0 - std::pow(x, 2.0)) *
                 2.0 / std::pow(1.0 - x, 2.0);

            // loop over all the angular points
            for(unsigned int a=lebedev_offset; a<this->angular_points + lebedev_offset; a++) {
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

        // if there is only one atom in the system, just ignore the
        // the Becke grid part below and continue
        if(this->mol->get_nr_atoms() == 1) {
            continue;
        }

        // calculate the Becke weights for the atomic grids
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

    // correct the weights for the Jacobian r**2 and the 4PI from the Lebedev quadrature
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

void MolecularGrid::output_density(const MatrixXXd& D) {
    const double distx = 10.0;
    const double disty = 10.0;
    const double distz = 10.0;

    const unsigned nrpoints = 100;

    const double dx = distx / (double)(nrpoints - 1);
    const double dy = disty / (double)(nrpoints - 1);
    const double dz = distz / (double)(nrpoints - 1);

    unsigned int counter = 0;

    std::cout << "electron density" << std::endl;
    std::cout << "    1.00000" << std::endl;
    std::cout << distx << "    0    0" << std::endl;
    std::cout << "0    " << disty << "    0" << std::endl;
    std::cout << "0    0    " << distz << std::endl;
    std::cout << "0" << std::endl;
    std::cout << "Direct" << std::endl;
    std::cout << std::endl;
    std::cout << nrpoints << "  " << nrpoints << "  " << nrpoints << std::endl;

    for(unsigned int k=0; k<nrpoints; k++) {
        for(unsigned int l=0; l<nrpoints; l++) {
            for(unsigned int m=0; m<nrpoints; m++) {
                double x = m * dx - distx / 2.0;
                double y = l * dy - disty / 2.0;
                double z = k * dz - distz / 2.0;

                double value = 0.0;
                for(unsigned int i=0; i<this->mol->get_nr_bfs(); i++) {
                    for(unsigned int j=0; j<this->mol->get_nr_bfs(); j++) {
                        value = 2.0 * D(i,j) *
                                this->mol->get_cgf(i).get_amp(vec3(x,y,z)) *
                                this->mol->get_cgf(j).get_amp(vec3(x,y,z));
                    }
                }

                std::cout << value << "\t";

                if(counter == 6) {
                    counter = 0;
                    std::cout << std::endl;
                    continue;
                }

                counter++;
            }
        }
    }
}

void MolecularGrid::calculate_hartree_potential() {
    this->calculate_rho_lm();
    this->calculate_U_lm();
}

/**
 * @fn calculate_rho_lm
 * @brief Calculate rho_lm from rho_n
 *
 * For each radial point rho_r(r), the values rho(r, theta, phi) can be calculated
 * by summing over rho_lm(r) * Y_lm(theta, phi). Here, the coefficients rho_lm
 * are calculated for the spherical harmonic expansion of the electron density.
 * From this expansion, the Hartree potential can be evaluated.
 *
 * @return void
 */
void MolecularGrid::calculate_rho_lm() {
    // calculate the size of the vector given a maximum angular quantum number
    unsigned int n = 0;
    for(int l=0; l <= lmax; l++) {
        n += 2 * l + 1;
    }

    // expand vector have appropriate size
    rho_lm = MatrixXXd::Zero(this->grid.size() / this->angular_points, n);

    // create vector holding the radial electron density
    VectorXd rho_r = VectorXd::Zero(this->grid.size() / this->angular_points);

    // create vector holding r values (i.e. distance to atomic center)
    this->r_n = VectorXd(this->grid.size() / this->angular_points);

    VectorXd weights = this->get_weights();
    VectorXd densities = this->get_densities();
    VectorXd rho_n = weights.cwiseProduct(densities);

    for(unsigned int i=0; i<this->grid.size() / this->angular_points; i++) {
        n = 0;

        // calculate radial distance for each radial grid point
        this->r_n(i) = this->grid[i * this->angular_points].get_position().norm();

        // calculate total radial density
        for(unsigned int j=0; j<this->angular_points; j++) {
            unsigned int idx = i * angular_points + j;
            rho_r(i) += rho_n(idx);
        }

        // calculate rho_lm from the radial density and the spherical harmonics
        for(int l=0; l<=lmax; l++) {
            for(int m=-l; m<=l; m++) {
                double pre = this->prefactor_spherical_harmonic(l, m);
                for(unsigned int j=0; j<this->angular_points; j++) {
                    // get index positions
                    unsigned int idx = i * angular_points + j;

                    // get cartesian coordinates
                    vec3 pos = this->grid[idx].get_position();
                    const double w = Quadrature::lebedev_coefficients[j + lebedev_offset][3];

                    // convert to spherical coordinates
                    double r = pos.norm(); // due to unit sphere
                    double azimuth = std::atan2(pos(1),pos(0));
                    double pole = std::acos(pos(2) / r); // due to unit sphere else z/r

                    double y_lm = pre * spherical_harmonic(l, m, pole, azimuth);

                    this->rho_lm(i, n) += densities(idx) * y_lm * w;
                }
                this->rho_lm(i, n) *= 4.0 * M_PI;

                n++;
            }
        }
    }
}

void MolecularGrid::calculate_U_lm() {
    unsigned int N = this->grid.size() / this->angular_points;

    static const double sqrt4pi = std::sqrt(4.0 * M_PI);

    std::vector<double> atomic_density;
    atomic_density.resize(this->mol->get_nr_atoms(), 0.0);

    VectorXd weights = this->get_weights();
    VectorXd densities = this->get_densities();

    // calculate total electron density per atom
    for(unsigned int i=0; i<this->mol->get_nr_atoms(); i++) {
        for(unsigned int j=0; j<this->radial_points; j++) {
            for(unsigned int k=0; k<this->angular_points; k++) {
                unsigned int idx = j * angular_points + k;
                atomic_density[i] += weights(idx) * densities(idx);
            }
        }
        std::cout << "Atomic density on " << i+1 << ": " << atomic_density[i] << std::endl;
    }

    double q_n = atomic_density[0];

    double c1 = 0.0;    // constant for d2U/dz2
    double c2 = 0.0;    // constant for (dU/dz)^2

    const double h = 1.0 / (double)(N + 1);

    // calculate the size of the vector given a maximum angular quantum number
    unsigned int lm = 0;
    for(int l=0; l <= lmax; l++) {
        lm += 2 * l + 1;
    }

    // construct the finite difference matrix
    MatrixXXd A = MatrixXXd::Zero(N+2,N+2);

    // start constructing matrix
    for(unsigned int i=0; i<N+2; i++) {

        if(i > 0 && i < N+1) {
            c1 = dzdrsq(r_n(i-1), 1.0);
            c2 = d2zdr2(r_n(i-1), 1.0);
        }

        if(i == 0) {
            A(0,0) = 1.0;
            continue;
        }
        if(i == 1) {
            c1 /= 12.0 * h * h;
            c2 /= 12.0 * h;
            A(i,0) = 11.0 * c1 -3.0 * c2;
            A(i,1) = -20.0 * c1 - 10.0 * c2;
            A(i,2) = 6.0 * c1 + 18.0 * c2;
            A(i,3) = 4.0 * c1 - 6.0 * c2;
            A(i,4) = -1.0 * c1 + 1.0 * c2;
            continue;
        }
        if(i == 2) {
            c1 /= 12.0 * h*h;
            c2 /= 60.0 * h;
            A(i,0) = -1.0 * c1 + 3.0 * c2;
            A(i,1) = 16.0 * c1 - 30.0 * c2;
            A(i,2) = -30.0 * c1 - 20.0 * c2;
            A(i,3) = 16.0 * c1 + 60.0 * c2;
            A(i,4) = -1.0 * c1 - 15.0 * c2;
            A(i,5) = 0.0 * c1 + 2.0 * c2;
            continue;
        }
        if(i == N-1) {
            c1 /= 12.0 * h*h;
            c2 /= 60.0 * h;
            A(i,N-4) = 0.0 * c1 - 2.0 * c2;
            A(i,N-3) = -1.0 * c1 + 15.0 * c2;
            A(i,N-2) = 16.0 * c1 - 60.0 * c2;
            A(i,N-1) = -30.0 * c1 + 20.0 * c2;
            A(i,N ) = 16.0 * c1 + 30.0 * c2;
            A(i,N+1) = -1.0 * c1 - 3.0 * c2;
            continue;
        }
        if(i == N) {
            c1 /= 12.0 * h * h;
            c2 /= 12.0 * h;
            A(i,N-3) = -1.0 * c1 - 1.0 * c2;
            A(i,N-2) = 4.0 * c1 + 6.0 * c2;
            A(i,N-1) = 6.0 * c1 - 18.0 * c2;
            A(i,N  ) = -20.0 * c1 + 10.0 * c2;
            A(i,N+1) = 11.0 * c1 + 3.0 * c2;
            continue;
        }
        if(i == N+1) {
            A(i,i) = 1.0;
            continue;
        }
        c1 /= 180.0 * h*h;
        c2 /= 60.0 * h;
        A(i,i-3) = 2.0 * c1 - 1.0 * c2;
        A(i,i-2) = -27.0 * c1 + 9.0 * c2;
        A(i,i-1) = 270.0 * c1 - 45.0 * c2;
        A(i,i)   = -490.0 * c1;
        A(i,i+1) = 270.0 * c1 + 45.0 * c2;
        A(i,i+2) = -27.0 * c1 - 9.0 * c2;
        A(i,i+3) = 2.0 * c1 + 1.0 * c2;
    } // end constructing matrix

    // expand vector have appropriate size
    U_lm = MatrixXXd::Zero(N, lm);

    unsigned int n=0;
    // for each lm value, calculate the Ulm value
    for(int l=0; l<=lmax; l++) {
        for(int m=-l; m<=l; m++) {
            // construct the divergence vector
            VectorXd g = VectorXd::Zero(N+2);
            MatrixXXd M(A); // make a copy of matrix A

            // construct divergence vector
            for(unsigned int i=0; i<N+2; i++) {
                if(i == 0) {
                    if(n == 0) {
                        g(i) = sqrt4pi * q_n;
                    } else {
                        g(i) = 0.0;
                    }
                    continue;
                }
                if(i == N+1) {
                    g(i) = 0.0;
                    continue;
                }
                if(i > 0 && i < N+1) {
                    M(i,i) -= (double)l * (double)(l+1) / (r_n(i-1) * r_n(i-1));
                    g(i) = -4.0 * M_PI * r_n(i-1) * this->rho_lm(i-1,n);
                }
            } // end constructing matrix

            Eigen::PartialPivLU<MatrixXXd> dec(M);
            VectorXd b = dec.solve(g);

            // place b vector back into U_lm
            for(unsigned int i=1; i<N+1; i++) {
                this->U_lm(i-1,n) = b(i);
            }
            n++;
        }
    }

    VectorXd V = VectorXd::Zero(this->grid.size());

    // calculate V(r, theta, phi) from the radial density and the spherical harmonics
    // loop over radial points
    for(unsigned int i=0; i<this->grid.size() / this->angular_points; i++) {
        n = 0;
        for(int l=0; l<=lmax; l++) {
            for(int m=-l; m<=l; m++) {
                const double pre = this->prefactor_spherical_harmonic(l, m);
                for(unsigned int j=0; j<this->angular_points; j++) {
                    // get grid index positions
                    unsigned int idx = i * angular_points + j;

                    // get cartesian coordinates
                    vec3 pos = this->grid[idx].get_position();

                    // convert to spherical coordinates
                    double r = pos.norm(); // due to unit sphere
                    double azimuth = std::atan2(pos(1),pos(0));
                    double pole = std::acos(pos(2) / r); // due to unit sphere else z/r

                    double y_lm = pre * spherical_harmonic(l, m, pole, azimuth);
                    V(idx) += 1.0 / r * y_lm * U_lm(i,n);
                }
                n++;
            }
        }
    }

    MatrixXXd amplitudes  = this->get_amplitudes();
    VectorXd Vp = weights.cwiseProduct(V);
    VectorXd row = amplitudes.row(0);
    const unsigned int size = this->mol->get_nr_bfs();
    MatrixXXd Jj = MatrixXXd::Zero(size, size);

    for(unsigned int i=0; i<size; i++) {
        VectorXd row = amplitudes.row(i);
        VectorXd row_i = Vp.cwiseProduct(row);
        for(unsigned int j=i; j<size; j++) {
            Jj(i,j) = Jj(j,i) = row_i.dot(amplitudes.row(j)) / 2.0;
        }
    }

    std::cout << "Numerical approximation of Jj: " << Jj << std::endl;
}

double MolecularGrid::spherical_harmonic(int l, int m, double pole, double azimuth) const {
    return this->polar_function(l, m, pole) * this->azimuthal_function(m, azimuth);
}

double MolecularGrid::prefactor_spherical_harmonic(int l, int m) const {
    static const double pre = 1.0 / sqrt(4 * M_PI);
    return pre * (m == 0 ? 1 : sqrt(2.0)) *
            sqrt((double)(2 * l + 1) * boost::math::factorial<double>(l - std::abs(m)) /
                 boost::math::factorial<double>(l + std::abs(m)) );
}

double MolecularGrid::polar_function(int l, int m, double theta) const {
    return this->legendre_p(l, std::abs(m), cos(theta));
}

double MolecularGrid::azimuthal_function(int m, double phi) const {
    if(m == 0) return 1.0;

    if(m > 0) {
        return cos((double)m * phi);
    } else {
        return sin(-(double)m * phi);
    }
}

/* legendre function */

double MolecularGrid::legendre (int n, double x) const {
    int i;

    if(n < 0) {
        return -1;
    }

    double v[n];
        v[0] = 1.0;

    if(n < 1) {
        return 1.0;
    }

    v[1] = x;

    for ( i = 2; i <= n; i++ ) {
          v[i] =     ( ( double ) ( 2 * i - 1 ) * x    * v[i-1]
                     - ( double ) (     i - 1 ) *        v[i-2] )
                     / ( double ) (     i     );
    }

    return v[n];
}

/* associated legendre function
 *
 * note that x should lie between -1 and 1 for this to work, else
 * a NAN will be returned
 */

double MolecularGrid::legendre_p (int n, int m, double x) const {
    double fact;
    int i;
    int j;
    int k;
    double v[n+1];

    for ( i = 0; i < n + 1; i++ ) {
        v[i] = 0.0;
    }

    //
    //  J = M is the first nonzero function.
    //
    if ( m <= n ) {
        v[m] = 1.0;

        fact = 1.0;
        for ( k = 0; k < m; k++ ) {
            v[m] *= - fact * sqrt ( 1.0 - x * x);
            fact += 2.0;
        }
    }

    //
    //  J = M + 1 is the second nonzero function.
    //
    if ( m + 1 <= n ) {
        v[m+1] = x * ( double ) ( 2 * m + 1 ) * v[m];
    }
    //
    //  Now we use a three term recurrence.
    //
    for ( j = m + 2; j <= n; j++ ) {
          v[j] = ( ( double ) ( 2 * j     - 1 ) * x * v[(j-1)]
              + ( double ) (   - j - m + 1 ) *        v[(j-2)] )
              / ( double ) (     j - m     );
    }

    return v[n];
}

double MolecularGrid::d2zdr2(double r, double m) {
    double nom = m*m * (m + 3.0 * r);
    double denom = 2.0 * M_PI * std::pow((m * r) / ((m+r)*(m+r)),1.5) * std::pow(m+r,5.0);
    return nom/denom;
}

double MolecularGrid::dzdrsq(double r, double m) {
    return m / (M_PI * M_PI * r * (m + r)*(m + r));
}
