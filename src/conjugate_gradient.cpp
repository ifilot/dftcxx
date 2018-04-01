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
ConjugateGradient::ConjugateGradient(const std::shared_ptr<Molecule> _mol, const std::shared_ptr<Settings> _settings) :
molecule(_mol),
settings(_settings) {

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
    // h.normalize();

    while((std::abs(eold - enew) > tol || iter < 2) && iter < 30) {
        std::cout << "*** Performing line search" << std::endl;
        const double alpha = this->line_search_backtrack(g, h, .5, 0.5);

        // get new energy and gradient
        eold = enew;
        std::cout << "*** Using new line search value" << std::endl;
        enew = this->get_energy();

        if(enew > eold) {
            std::cout << "*** Performing conjugate gradient step" << std::endl;
            this->molecule->perturb_atoms(-alpha * h);
            this->dft.release();
            eold = 1000;
            h = -g;
            continue;
        }

        const VectorXd oldg = g;
        g = this->get_gradient();

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
        // h.normalize();

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
 * @brief      line search backtrack
 */
double ConjugateGradient::line_search_interpolation(const VectorXd& g, const VectorXd& h, double rho, double max) {
    static const double t = 1e-4;

    // calculate energy and gradient at 0
    double e0 = this->get_energy();
    auto f0 = this->get_gradient();

    // calculate energy at a0
    double a0 = max;
    double e1 = this->get_energy_perturbation(a0 * h);
    double deriv = f0.dot(h);

    std::cout << a0 << std::endl;

    // check whether this step was sufficient
    double e_ideal = e0 + t * a0 * deriv;
    std::cout << "Expecting to get at least: " << e_ideal << ", got: " << e1 << std::endl;
    if(e1 <= e_ideal) {
        return a0;
    }

    // establish predictor step
    double a1 = -f0.dot(h) * a0 * a0 / (2.0 * (e1 - e0 - deriv * a0));

    std::cout << a1 << std::endl;

    double e2 = this->get_energy_perturbation((a1 - a0) * h);

    // check whether this step was sufficient
    e_ideal = e0 + t * a1 * deriv;
    std::cout << "Expecting to get at least: " << e_ideal << ", got: " << e2 << std::endl;
    if(e2 <= e_ideal) {
        return a1;
    }

    double a = 1.0 / (a0 * a0 * a1 * a1 * (a1 - a0)) * (a0 * a0 * (e2 - e0 - deriv * a1) - a1 * a1 * (e1 - e0 - deriv * a0));
    double b = 1.0 / (a0 * a0 * a1 * a1 * (a1 - a0)) * (-a0 * a0 * a0 * (e2 - e0 - deriv * a1) + a1 * a1 * a1 * (e1 - e0 - deriv * a0));
    double a2 = -b + std::sqrt(b * b - 3.0 * a * deriv) / (3.0 * a);

    if(a2 < 0) {
        a2 = a1 / 2.0;
    }

    std::cout << a2 << std::endl;

    double e3 = this->get_energy_perturbation((a2 - a1) * h);

    return a2;
}

/**
 * @brief      perform scf and obtain energy
 *
 * @return     The energy.
 */
double ConjugateGradient::get_energy() {
    if(!this->dft) {
        this->dft = std::make_unique<DFT>(this->molecule, this->settings);
    }
    return this->dft->get_energy();
}

/**
 * @brief      perform scf and get gradient
 *
 * @return     The gradient
 */
VectorXd ConjugateGradient::get_gradient() {
    if(!this->dft) {
        this->dft = std::make_unique<DFT>(this->molecule, this->settings);
    }
    return this->dft->get_force_vector();
}

/**
 * @brief      Gets the energy after perturbation, but restore configuration
 */
double ConjugateGradient::get_energy_perturbation(const VectorXd& h) {
    this->dft.release();
    this->molecule->perturb_atoms(h);
    return this->get_energy();
}
