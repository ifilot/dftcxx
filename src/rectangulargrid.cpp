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
                this->grid.emplace_back(pos, ctr);
                this->grid.back().set_basis_func_amp(mol);
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

void RectangularGrid::write_density() {
    std::vector<uint8_t> buffer(this->gridsize * this->gridsize * 4, 0);

    for(unsigned int i=0; i<this->gridsize; i++) {
        for(unsigned int j=0; j<this->gridsize; j++) {
            const unsigned int idx = i * this->gridsize * this->gridsize + (gridsize/2) * gridsize + j;

            vec3 grad = this->grid[idx].get_gradient();
            double valx = (std::log10(std::abs(grad[0])) + 5) * 255.0;
            double valy = (std::log10(std::abs(grad[1])) + 5) * 255.0;
            double valz = (std::log10(std::abs(grad[2])) + 5) * 255.0;

            const unsigned int canvas_idx = i * this->gridsize + j;
            valx = std::max(0.0, valx);
            valx = std::min(255.0, valx);
            valy = std::max(0.0, valy);
            valy = std::min(255.0, valy);
            valz = std::max(0.0, valz);
            valz = std::min(255.0, valz);

            buffer[canvas_idx * 4] = (uint8_t)valx;
            buffer[canvas_idx * 4 + 1] = (uint8_t)valy;
            buffer[canvas_idx * 4 + 2] = (uint8_t)valz;
            buffer[canvas_idx * 4 + 3] = 255;
        }
    }

    PNG::write_image_buffer_to_png("test.png", buffer, this->gridsize, this->gridsize, PNG_COLOR_TYPE_RGBA);
}
