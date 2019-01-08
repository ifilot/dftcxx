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

#ifndef _CONJUGATE_GRADIENT_H
#define _CONJUGATE_GRADIENT_H

#include "dft.h"

class ConjugateGradient {
private:
    std::shared_ptr<Molecule> molecule;
    std::shared_ptr<Settings> settings;
    std::unique_ptr<DFT> dft;
    MatrixXXd P;                        // copy of the density matrix
    std::vector<double> energies;       // store energies
    unsigned int iter;

public:
    /**
     * @brief      Constructs the object.
     *
     * @param[in]  _dft  pointer to dft object
     * @param[in]  _dft  pointer to settings object
     */
    ConjugateGradient(const std::shared_ptr<Molecule> _mol,
                      const std::shared_ptr<Settings> _settings);

    /**
     * @brief      optimize DFT structure using conjugate gradient
     */
    void optimize();

    /**
     * @brief      line search backtrack
     */
    double line_search_backtrack(const VectorXd& g, const VectorXd& h, double rho, double max);

    double line_search_interpolation(const VectorXd& g, const VectorXd& h, double rho, double max);

    /**
     * @brief      perform scf and obtain energy
     *
     * @return     The energy.
     */
    double get_energy();

    /**
     * @brief      perform scf and get gradient
     *
     * @return     The gradient
     */
    VectorXd get_gradient();

    /**
     * @brief      Gets the energy after perturbation, but restore configuation
     */
    double get_energy_perturbation(const VectorXd& h);

private:
    void output_iteration(unsigned int iter);

};

#endif // _CONJUGATE_GRADIENT_H
