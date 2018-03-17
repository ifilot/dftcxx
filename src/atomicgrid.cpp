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

#include "atomicgrid.h"

AtomicGrid::AtomicGrid(const std::shared_ptr<Atom>& _atom) :
atom(_atom)
{

}

void AtomicGrid::create_atomic_grid(unsigned int _radial_points,
                                    unsigned int _angular_points,
                                    unsigned int _lebedev_offset,
                                    unsigned int _lmax,
                                    const std::shared_ptr<Molecule>& _mol) {
    this->radial_points = _radial_points;
    this->angular_points = _angular_points;
    this->lebedev_offset = _lebedev_offset;
    this->lmax = _lmax;
    this->mol = _mol;
    this->density_cached = false;

    // construct grid points
    const vec3 p1 = this->atom->get_position();
    const double f = pi / (double)(radial_points + 1);
    this->r_n = VectorXd::Zero(this->radial_points);

    // loop over all the radial points
    for(unsigned int p=1; p<=radial_points; p++) {
        // calculate weight function of the Gauss-Chebyshev integration
        double w = f * std::pow(std::sin(f * (double)p), 2.0);

        // calculate x point in the interval [-1,1] for Gauss-Chebyshev integration
        double x = std::cos(f * (double)p);

        // convert x to r so that the interval is converted from [-1,1] to [0, infinity]
        const double r = (1.0 + x) / (1.0 - x);
        this->r_n(p-1) = r;

        // calculate weight function of the Gauss-Chebyshev integration
        w =  w / std::sqrt(1.0 - std::pow(x, 2.0)) *
             2.0 / std::pow(1.0 - x, 2.0);

        // loop over all the angular points
        for(unsigned int a=lebedev_offset; a<angular_points + lebedev_offset; a++) {
            // create a new grid point given the coordinates on the unit sphere; these
            // coordinates are multiplied by the radius as calculated above
            const vec3 p2 = vec3(Quadrature::lebedev_coefficients[a][0],
                                 Quadrature::lebedev_coefficients[a][1],
                                 Quadrature::lebedev_coefficients[a][2]) * r;
            this->grid.push_back(GridPoint(p1 + p2, p2));

            // set the total integration weight by multiplying the Gauss-Chebyshev weight by the Lebedev weight
            // note that the 'weight' from the Jacobian due to spherical coordinates and the weight
            // due to Becke grid are calculated later
            this->grid.back().set_weight(w * Quadrature::lebedev_coefficients[a][3]);

            // set the amplitudes of the basis functions at this grid point
            this->grid.back().set_basis_func_amp(mol);
        }
    }

    // correct the weights for the Jacobian r**2 and the 4PI from the Lebedev quadrature
    for(unsigned int i=0; i<grid.size(); i++) {
        this->grid[i].multiply_weight((this->grid[i].get_position_fuzzy_cell()).squaredNorm() * 4.0 * pi);
    }
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
void AtomicGrid::set_density(const MatrixXXd& P) {
    for(unsigned int i=0; i<this->grid.size(); i++) {
        this->grid[i].set_density(P);
    }
    this->density_cached = false;
}

/**
 * @fn get_density
 * @brief get the total electron density
 *
 * @return total electron density
 */
double AtomicGrid::get_density() {
    if(!this->density_cached) {
        this->calculate_density();
    }

    return this->density;
}

/**
 * @fn get_weights
 * @brief get the weights of all the grid points as a vector
 *
 * @return vector containing all the weights
 */
VectorXd AtomicGrid::get_weights() const {
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
VectorXd AtomicGrid::get_densities() const {
    VectorXd densities = VectorXd(this->grid.size());
    for(unsigned int i=0; i<this->grid.size(); i++) {
        densities(i) = this->grid[i].get_density();
    }

    return densities;
}

/**
 * @fn get_positions
 * @brief get the positions of all the grid points as a vector
 *
 * @return vector containing all the positions
 */
MatrixXXd AtomicGrid::get_positions() const {
    MatrixXXd positions = MatrixXXd(this->grid.size(), 3);
    for(unsigned int i=0; i<this->grid.size(); i++) {
        for(unsigned int j=0; j<3; j++) {
            positions(i,j) = this->grid[i].get_position()[j];
        }
    }

    return positions;
}

/**
 * @brief      set fuzzy cell weights
 *
 * @param[in]  vw  fuzzy cell weights
 */
void AtomicGrid::set_becke_weights(const VectorXd& vw) {
    if(vw.size() != this->grid.size()) {
        throw std::runtime_error("Dimension of correction vector does not match dimension of grid");
    }

    for(unsigned int i=0; i<vw.size(); i++) {
        this->grid[i].multiply_weight(vw(i));
        this->grid[i].set_becke_weight(vw(i));
    }
}

/**
 * @fn get_amplitudes
 * @brief get the amplitudes of all the grid points and all the basis functions as a matrix
 *
 * @return matrix (basis functions x grid points)
 */
MatrixXXd AtomicGrid::get_amplitudes() const {
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
void AtomicGrid::calculate_rho_lm() {
    // calculate the size of the vector given a maximum angular quantum number
    unsigned int n = 0;
    for(int l=0; l <= lmax; l++) {
        n += 2 * l + 1;
    }

    // expand vector have appropriate size
    rho_lm = MatrixXXd::Zero(this->grid.size() / this->angular_points, n);

    // create vector holding the radial electron density
    VectorXd rho_r = VectorXd::Zero(this->grid.size() / this->angular_points);

    VectorXd weights = this->get_weights();
    VectorXd densities = this->get_densities();
    VectorXd rho_n = weights.cwiseProduct(densities);

    #pragma omp parallel for
    for(unsigned int i=0; i<this->grid.size() / this->angular_points; i++) {
        unsigned int cnt = 0;

        // calculate total radial density
        for(unsigned int j=0; j<this->angular_points; j++) {
            unsigned int idx = i * angular_points + j;
            rho_r(i) += rho_n(idx);
        }

        // calculate rho_lm from the radial density and the spherical harmonics
        for(int l=0; l<=lmax; l++) {
            for(int m=-l; m<=l; m++) {
                const double pre = SH::prefactor_spherical_harmonic(l, m);
                for(unsigned int j=0; j<this->angular_points; j++) {
                    // get index positions
                    unsigned int idx = i * angular_points + j;

                    // get cartesian coordinates
                    const vec3& pos = this->grid[idx].get_position_fuzzy_cell();
                    const double w = Quadrature::lebedev_coefficients[j + lebedev_offset][3];
                    const double wb = this->grid[idx].get_becke_weight();

                    // convert to spherical coordinates
                    const double r = pos.norm(); // due to unit sphere
                    const double azimuth = std::atan2(pos(1),pos(0));
                    const double pole = std::acos(pos(2) / r); // due to unit sphere else z/r

                    const double y_lm = pre * SH::spherical_harmonic(l, m, pole, azimuth);

                    this->rho_lm(i, cnt) += densities(idx) * y_lm * w * wb;
                }
                this->rho_lm(i, cnt) *= 4.0 * M_PI;   // correct for unit sphere

                cnt++;
            }
        }
    }
}

/**
 * @brief      calculate spherical harmonic coefficients
 */
void AtomicGrid::calculate_U_lm() {
    const size_t N = this->grid.size() / this->angular_points;

    static const double sqrt4pi = std::sqrt(4.0 * M_PI);
    VectorXd weights = this->get_weights();
    const double q_n = this->get_density();
    const double h = 1.0 / (double)(N + 1);

    // calculate the size of the vector given a maximum angular quantum number
    this->lm = 0;
    for(int l=0; l <= lmax; l++) {
        this->lm += 2 * l + 1;
    }

    // construct the finite difference matrix
    MatrixXXd A = MatrixXXd::Zero(N+2,N+2);

    // start constructing matrix
    #pragma omp parallel for
    for(unsigned int i=0; i<N+2; i++) {

        double c1 = 0.0;    // constant for d2U/dz2
        double c2 = 0.0;    // constant for (dU/dz)^2

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
    U_lm = MatrixXXd::Zero(N, this->lm);

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

    this->V_fuzzy_cell = VectorXd::Zero(this->grid.size());
    this->V = VectorXd::Zero(this->grid.size());

    // calculate V(r, theta, phi) from the radial density and the spherical harmonics
    // loop over radial points
    #pragma omp parallel for
    for(unsigned int i=0; i<this->grid.size() / this->angular_points; i++) {
        for(int l=0; l<=lmax; l++) {
            for(int m=-l; m<=l; m++) {
                unsigned int n = this->calculate_lm(l,m);
                const double pre = SH::prefactor_spherical_harmonic(l, m);
                for(unsigned int j=0; j<this->angular_points; j++) {
                    // get grid index positions
                    unsigned int idx = i * angular_points + j;

                    // get cartesian coordinates
                    vec3 pos = this->grid[idx].get_position_fuzzy_cell();

                    // convert to spherical coordinates
                    const double r = pos.norm(); // due to unit sphere
                    const double azimuth = std::atan2(pos(1),pos(0));
                    const double pole = std::acos(pos(2) / r); // due to unit sphere else z/r

                    const double y_lm = pre * SH::spherical_harmonic(l, m, pole, azimuth);
                    this->V_fuzzy_cell(idx) += 1.0 / r * y_lm * U_lm(i,n);
                }
            }
        }
    }

    // construct splines for interpolation
    this->interpolate_sh_coeff();
}

/**
 * @brief      Calculates the j matrix.
 */
void AtomicGrid::calculate_J_matrix() {
    VectorXd weights = this->get_weights();
    MatrixXXd amplitudes = this->get_amplitudes();     // get NxM amplitude matrix where N are the number of BF and M the gridpoints

    VectorXd Vp = weights.cwiseProduct(this->V);

    const unsigned int size = this->mol->get_nr_bfs();
    this->Jj = MatrixXXd::Zero(size, size);

    for(unsigned int i=0; i<size; i++) {
        VectorXd row = amplitudes.row(i);
        VectorXd row_i = Vp.cwiseProduct(row);
        for(unsigned int j=i; j<size; j++) {
            this->Jj(i,j) = this->Jj(j,i) = row_i.dot(amplitudes.row(j)) * 0.5;
        }
    }
}

/**
 * @brief      Gets the sh value.
 *
 * @param[in]  r     radial distance from fuzzy cell center
 * @param[in]  lm    lm index
 *
 * @return     The sh value.
 */
double AtomicGrid::get_sh_value(double r, unsigned int lm) const {
    return this->splines[lm].eval(r);
}

double AtomicGrid::d2zdr2(double r, double m) {
    double nom = m*m * (m + 3.0 * r);
    double denom = 2.0 * M_PI * std::pow((m * r) / ((m+r)*(m+r)),1.5) * std::pow(m+r,5.0);
    return nom/denom;
}

double AtomicGrid::dzdrsq(double r, double m) {
    return m / (M_PI * M_PI * r * (m + r)*(m + r));
}

/**
 * @fn calculate_density
 * @brief calculate the total electron density (number of electrons)
 *
 * Calculates the total number of electrons by summing the weights
 * multiplied by the local value of the electron density.
 *
 */
void AtomicGrid::calculate_density() {
    this->density = 0.0;
    for(unsigned int i=0; i<this->grid.size(); i++) {
        this->density += this->grid[i].get_weight() * this->grid[i].get_density();
    }

    this->density_cached = true;
}

/**
 * @brief      perform interpolation on the spherical harmonic coefficients
 */
void AtomicGrid::interpolate_sh_coeff() {
    // clear any previously constructed splines
    this->splines.clear();

    // construct x values
    std::vector<double> x;
    for(int i=this->radial_points-1; i>=0; i--) {
        x.push_back(this->r_n(i));
    }

    // calculate for each lm value the y values
    for(unsigned int i=0; i<this->lm; i++) {
        std::vector<double> y;
        for(int j=this->radial_points-1; j>=0; j--) {
            y.push_back(this->U_lm(j,i));
        }

        this->splines.push_back(Cspline());
        this->splines.back().set_values(x,y);
        this->splines.back().generate_spline();
    }
}

/**
 * @brief      calculate lm index
 *
 * @param[in]  l     l quantum number
 * @param[in]  m     m quantum number
 *
 * @return     the lm index
 */
unsigned int AtomicGrid::calculate_lm(unsigned int l, unsigned int m) {
    unsigned int lm_idx = 0;
    for(int i=0; i < l; i++) {
        lm_idx += 2 * i + 1;
    }

    return lm_idx + l + m;
}
