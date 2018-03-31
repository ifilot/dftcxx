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

#include "conjugate_gradient.h"

/**
 * @brief      Constructs the object.
 *
 * @param[in]  _dft  pointer to dft object
 */
ConjugateGradient::ConjugateGradient(const std::shared_ptr<DFT>& _dft) :
dft(_dft) {

}

/**
 * @brief      optimize DFT structure using conjugate gradient
 */
void ConjugateGradient::optimize() {
    static const double tol = 1e-3;
    static const unsigned int reset_cg = 10;

    unsigned int iter = 1;              // keep track of iterations
    unsigned int reset_counter = 0;     // nr of cycles since last cg reset

    std::cout << "*** Calculate initial energy" << std::endl;
    double enew = this->get_energy();
    double eold = 1000; // some large number

    VectorXd g = this->dft->get_force_vector();
    VectorXd h = -g;
    h.normalize();

    while((std::abs(eold - enew) > tol || iter < 2) && iter < 30) {
        std::cout << "*** Performing line search" << std::endl;
        const double alpha = this->line_search_backtrack(g, h, .5, 0.5);
        dft->perturb_atoms(alpha * h);

        // get new energy and gradient
        eold = enew;
        std::cout << "*** Using new line search value" << std::endl;
        enew = this->get_energy();

        if(enew > eold) {
            std::cout << "*** Performing conjugate gradient step" << std::endl;
            this->dft->perturb_atoms(-alpha * h);
            this->dft->scf();
            eold = 1000;
            h = -g;
            continue;
        }

        const VectorXd oldg = g;
        g = this->dft->get_force_vector();

        // calculate new search direction
        double beta = g.dot(g - oldg) / oldg.dot(oldg);

        // reset algorithm every so many iterations
        if(reset_counter > reset_cg) {
            reset_counter = 0;
            beta = 0;
        } else {
            reset_counter++;
        }

        h = -g + beta * h;
        h.normalize();

        iter++;
    }
}

/**
 * @brief      line search backtrack
 */
double ConjugateGradient::line_search_backtrack(const VectorXd& g, const VectorXd& h, double rho, double max) {
    static const double t = 1e-4;
    static const double c_ideal = .5;
    unsigned int iter = 0;

    double alpha = max;
    double estart = this->get_energy();
    double enew = this->get_energy_perturbation(alpha * h);

    const double e_ideal = estart + c_ideal * t * g.dot(h);

    while(enew > e_ideal && iter < 20) {
        alpha *= rho;
        double eold = enew;
        enew = this->get_energy_perturbation(alpha * h);
        iter++;


    }

    return alpha;
}

/**
 * @brief      perform scf and obtain energy
 *
 * @return     The energy.
 */
double ConjugateGradient::get_energy() {
    std::cout << "---------------------------------------------------" << std::endl;
    this->dft->scf();
    std::cout << "---------------------------------------------------" << std::endl;
    return this->dft->get_energy();
}

/**
 * @brief      Gets the energy after perturbation, but restore configuration
 */
double ConjugateGradient::get_energy_perturbation(const VectorXd& h) {
    this->dft->perturb_atoms(h);
    std::cout << "---------------------------------------------------" << std::endl;
    this->dft->scf();
    std::cout << "---------------------------------------------------" << std::endl;
    double e = this->dft->get_energy();
    this->dft->perturb_atoms(-h);

    return e;
}
