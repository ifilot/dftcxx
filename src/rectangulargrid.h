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

#ifndef _RECTANGULAR_GRID
#define _RECTANGULAR_GRID

#include "gridpoint.h"
#include "pngfuncs.h"

class RectangularGrid {
private:
    std::shared_ptr<Molecule> mol;  //!< pointer to molecule class
    std::vector<GridPoint> grid;    //!< set of all gridpoints
    unsigned int gridsize;

public:
    RectangularGrid(const std::shared_ptr<Molecule>& mol);

    /**
     * @brief      build a grid
     *
     * @param[in]  size  size of the grid in angstrom
     * @param[in]  dp    number of grid points in each Cartesian direction
     */
    void build_grid(double size, unsigned int dp);

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
    void set_density(const MatrixXXd& P);

    void write_density();

private:
};

#endif // _RECTANGULAR_GRID
