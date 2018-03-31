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

#include "rectangulargrid.h"

RectangularGrid::RectangularGrid(const std::shared_ptr<Molecule>& _mol) {
    this->mol = _mol;
}

/**
 * @brief      build a grid
 *
 * @param[in]  size  size of the grid in angstrom
 * @param[in]  dp    number of grid points in each Cartesian direction
 */
void RectangularGrid::build_grid(double size, unsigned int dp) {
    const double gdsize = size / (double)(dp - 1);
    const vec3 ctr(0,0,0);
    this->gridsize = dp;

    for(unsigned int i=0; i<dp; i++) {
        for(unsigned int j=0; j<dp; j++) {
            for(unsigned int k=0; k<dp; k++) {
                vec3 pos((double)k * gdsize - size / 2.0,
                         (double)j * gdsize - size / 2.0,
                         (double)i * gdsize - size / 2.0);

                // add new grid point
                this->grid.emplace_back(pos, ctr);

                // set basis function amplitude
                this->grid.back().set_basis_func_amp(mol);

                // set basis function gradient
                this->grid.back().set_basis_func_grad(mol);
            }
        }
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
void RectangularGrid::set_density(const MatrixXXd& P) {
    for(unsigned int i=0; i<this->grid.size(); i++) {
        this->grid[i].set_density(P);
        this->grid[i].set_gradient(P);
    }
}

/**
 * @brief      write the gradient to a data file
 *
 * @param[in]  filename  path to data file
 */
void RectangularGrid::write_gradient(const std::string& filename) {
    std::ofstream out(filename);

    for(unsigned int i=0; i<this->grid.size(); i++) {

        const vec3& pos = this->grid[i].get_position();
        const vec3& grad = this->grid[i].get_gradient();

        out << boost::format("%12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f\n")
                % pos[0] % pos[1] % pos[2] % grad[0] % grad[1] % grad[2];
    }

    out.close();
}
