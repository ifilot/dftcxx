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
ConjugateGradient::ConjugateGradient(const std::shared_ptr<Molecule> _mol,
                                     const std::shared_ptr<Settings> _settings) :
molecule(_mol),
settings(_settings) {

}

/**
 * @brief      optimize DFT structure using conjugate gradient
 */
void ConjugateGradient::optimize() {
    static const double tol = 1e-4;
    static const unsigned int reset_cg = 10;
    const unsigned int itermax = 50;

    this->iter = 1;                     // keep track of iterations
    unsigned int reset_counter = 0;     // nr of cycles since last cg reset

    this->output_iteration(iter);

    std::cout << "[CG-OPT] Calculate initial energy" << std::endl << std::endl;
    double enew = this->get_energy();
    this->energies.push_back(enew);
    double eold = 1000; // some large number

    VectorXd g = this->dft->get_force_vector();
    VectorXd h = -g;
    // h.normalize();

    // keep looping until energy criterion (tol) is met or maximum
    // number of iterations is reached
    while(std::abs(eold - enew) > tol && iter < itermax) {
        this->iter++;
        this->output_iteration(iter);

        std::cout << "[CG-OPT] Performing line search" << std::endl << std::endl;
        const double alpha = this->line_search_backtrack(g, h, 0.5, 0.5);

        // get new energy and gradient
        eold = enew;
        std::cout << "[CG-OPT] Using new line search value" << std::endl << std::endl;
        enew = this->get_energy();
        this->energies.push_back(enew);

        if(enew > eold) {
            std::cout << "[CG-OPT] Performing conjugate gradient step" << std::endl << std::endl;
            this->molecule->perturb_atoms(-alpha * h);
            this->dft.release();    // destroy DFT object
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
    }

    if(iter < itermax) {
        std::cout << "------------------------------------------------------" << std::endl;
        std::cout << " Stopping CG-opt: Energy criterion is reached" << std::endl;
        std::cout << "------------------------------------------------------" << std::endl;
    } else {
        std::cout << "------------------------------------------------------" << std::endl;
        std::cout << " Stopping CG-opt: Maximum number of iterations" << std::endl;
        std::cout << "------------------------------------------------------" << std::endl;
    }

    // print energy convergence
    for(unsigned int i=0; i<this->energies.size(); i++) {
        std::cout << boost::format("%3i  %12.6f") % (i+1) % this->energies[i] << std::endl;
    }
    std::cout << "------------------------------------------------------" << std::endl;

    // print final geometry
    //this->molecule->print_geometry();
}

/**
 * @brief      line search backtrack
 */
double ConjugateGradient::line_search_backtrack(const VectorXd& g, const VectorXd& h, double rho, double max) {
    static const double t = 1e-4;
    static const double c_ideal = .5;
    unsigned int maxsteps = 4;
    unsigned int steps = 0;

    double alpha = max;
    double estart = this->get_energy();
    double enew = this->get_energy_perturbation(alpha * h);

    const double e_ideal = estart + c_ideal * t * g.dot(h);

    while(enew > e_ideal && steps < maxsteps) {
        alpha *= rho;
        double eold = enew;
        this->iter++;
        this->output_iteration(this->iter);
        enew = this->get_energy_perturbation(alpha * h);
        steps++;
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
    // check if DFT object is destroyed
    if(!this->dft) {
        this->dft = std::make_unique<DFT>(this->molecule, this->settings);

        // if the density matrix is not empty, re-use the old wave function
        if(this->P.cols() > 0) {
            this->dft->set_wavefunction(this->P);
        }
    }

    double energy = this->dft->get_energy();
    this->P = this->dft->get_density_matrix();

    return energy;
}

/**
 * @brief      perform scf and get gradient
 *
 * @return     The gradient
 */
VectorXd ConjugateGradient::get_gradient() {
    if(!this->dft) {
        this->dft = std::make_unique<DFT>(this->molecule, this->settings);

        // if the density matrix is not empty, re-use the old wave function
        if(this->P.cols() > 0) {
            this->dft->set_wavefunction(this->P);
        }
    }

    auto forces = this->dft->get_force_vector();
    this->P = this->dft->get_density_matrix();

    return forces;
}

/**
 * @brief      Gets the energy after perturbation, but restore configuration
 */
double ConjugateGradient::get_energy_perturbation(const VectorXd& h) {
    this->dft.release();
    this->molecule->perturb_atoms(h);

    return this->get_energy();
}

void ConjugateGradient::output_iteration(unsigned int iter) {
    std::cout << boost::format("------------------------------------------------------") << std::endl << std::endl;
    std::cout << boost::format("                         %02i                         ") % iter << std::endl << std::endl;
    std::cout << boost::format("------------------------------------------------------") << std::endl << std::endl;
}
