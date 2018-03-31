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

#include "gridpoint.h"

/**
 * @fn GridPoint
 * @brief GridPoint Constructor
 *
 * @param _r    vec3 position of the grid point
 *
 * @return GridPoint instance
 */
GridPoint::GridPoint(const vec3& _r, const vec3& _r_at):
    r(_r),
    r_at(_r_at)
    {}

/**
 * @fn set_basis_func_amp
 * @brief calculates the amplitudes at the grid point of all basis functions
 *
 * @param _mol      pointer to the molecule object
 *
 * @return void
 */
void GridPoint::set_basis_func_amp(const std::shared_ptr<Molecule>& _mol) {
    const unsigned int size = _mol->get_nr_bfs();
    this->basis_func_amp = VectorXd(size);

    for(unsigned int i=0; i<size; i++) {
        this->basis_func_amp(i) = _mol->get_cgf(i).get_amp(this->r);
    }
}

/**
 * @fn set_basis_func_grad
 * @brief calculates the gradient at the grid point of all basis functions
 *
 * @param _mol      pointer to the molecule object
 *
 * @return void
 */
void GridPoint::set_basis_func_grad(const std::shared_ptr<Molecule>& _mol) {
    const unsigned int size = _mol->get_nr_bfs();
    this->basis_func_grad = MatrixXXd(size, 3);

    for(unsigned int i=0; i<size; i++) {
        vec3 grad = _mol->get_cgf(i).get_grad(this->r);
        for(unsigned int j=0; j<3; j++) {
            this->basis_func_grad(i,j) = grad(j);
        }
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
    this->density = 2.0 * this->basis_func_amp.dot(D * this->basis_func_amp);
}

/**
 * @fn set_density
 * @brief calculates the gradient at the grid point using the density matrix
 *
 * @param _mol      reference to density matrix
 *
 * @return void
 */
void GridPoint::set_gradient(const MatrixXXd& D) {
    VectorXd x = this->basis_func_grad.col(0);
    VectorXd y = this->basis_func_grad.col(1);
    VectorXd z = this->basis_func_grad.col(2);

    // apply product rule for the gradient
    double gx = 2.0 * this->basis_func_amp.transpose().dot(D * x) +
                2.0 * x.transpose().dot(D * this->basis_func_amp);
    double gy = 2.0 * this->basis_func_amp.transpose().dot(D * y) +
                2.0 * y.transpose().dot(D * this->basis_func_amp);
    double gz = 2.0 * this->basis_func_amp.transpose().dot(D * z) +
                2.0 * z.transpose().dot(D * this->basis_func_amp);

    this->grad = vec3(gx, gy, gz);
}
