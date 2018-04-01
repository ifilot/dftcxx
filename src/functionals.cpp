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

#include "functionals.h"

Functional::Functional() {

}

void Functional::xalpha_x_functional(const VectorXd& densitiesa,
                                     const VectorXd& densitiesb,
                                     VectorXd& ex,
                                     VectorXd& vxa,
                                     VectorXd& vxb) {

    static const double tol = 1e-10;
    static const double xalpha = 2.0 / 3.0;
    static const double fac = -2.25 * xalpha * std::pow(3.0 / 4.0 / pi, 1.0 / 3.0);

    ex  = VectorXd::Zero(densitiesa.size());
    vxa = VectorXd::Zero(densitiesa.size());
    vxb = VectorXd::Zero(densitiesb.size());

    #ifdef HAS_OPENMP
    #pragma omp parallel for
    #endif
    for(unsigned int i=0; i<densitiesa.size(); i++) {
        if(densitiesa(i) < tol) {
            continue;
        } else {
            double rho3 = std::pow(densitiesa(i), 1.0 / 3.0);
            ex(i) += fac * densitiesa(i) * rho3;
            vxa(i) += 4.0 / 3.0 * fac * rho3;
        }
    }

    #ifdef HAS_OPENMP
    #pragma omp parallel for
    #endif
    for(unsigned int i=0; i<densitiesb.size(); i++) {
        if(densitiesb(i) < tol) {
            continue;
        } else {
            double rho3 = std::pow(densitiesb(i), 1.0 / 3.0);
            ex(i) += fac * densitiesb(i) * rho3;
            vxb(i) += 4.0 / 3.0 * fac * rho3;
        }
    }
}

void Functional::vwm_c_functional(const VectorXd& densitiesa,
                                  const VectorXd& densitiesb,
                                  VectorXd& ec,
                                  VectorXd& vca,
                                  VectorXd& vcb) {

    static const double tol = 1e-10;

    ec.resize(densitiesa.size());
    vca.resize(densitiesa.size());
    vcb.resize(densitiesb.size());

    #ifdef HAS_OPENMP
    #pragma omp parallel for
    #endif
    for(unsigned int i=0; i<densitiesa.size(); i++) {
        double dens = densitiesa(i) + densitiesb(i);
        if(dens < tol) {
            ec(i) = 0.0;
            vca(i) = 0.0;
            vcb(i) = 0.0;
            continue;
        } else {
            double zeta = (densitiesa(i) - densitiesb(i)) / dens;
            double x = std::pow(3.0 / 4.0 / pi / dens, 1.0 / 6.0);
            double epsp = vwn_epsp(x);
            double depsp = vwn_depsp(x);
            double g = vwn_g(zeta);

            if(g < tol) {
                ec(i) = epsp * dens;
                vca(i) = epsp - (x / 6.0) * depsp;
                vcb(i) = epsp - (x / 6.0) * depsp;
                continue;
            }
            double epsf = vwn_epsf(x);
            double eps = epsp + g * (epsf - epsp);
            ec(i) = eps * dens;

            double depsf = vwn_depsf(x);
            double dg = vwn_dg(zeta);

            double deps_dx = depsp + g * (depsf - depsp);
            double deps_dg = (epsf - epsp) * dg;

            vca(i) = eps - (x / 6.0) * deps_dx + deps_dg * (1.0 - zeta);
            vcb(i) = eps - (x / 6.0) * deps_dx - deps_dg * (1.0 + zeta);
        }
    }
}

double Functional::vwn_xx(double x, double b, double c) {
    return x*x+b*x+c;
}

double Functional::vwn_epsp(double x) {
    return vwn_eps(x,0.0310907,-0.10498,3.72744,12.9352);
}

double Functional::vwn_epsf(double x) {
    return vwn_eps(x,0.01554535,-0.32500,7.06042,13.0045);
}

double Functional::vwn_eps(double x, double a, double x0, double b, double c) {
    double q = std::sqrt(4.0*c-b*b);
    double eps = a*(log(x*x/vwn_xx(x,b,c))
             - b*(x0/vwn_xx(x0,b,c))*log(pow(x-x0,2.0)/vwn_xx(x,b,c))
             + (2.0*b/q)*(1.0-(x0*(2.0*x0+b)/vwn_xx(x0,b,c))) * atan(q/(2.0*x+b)));
    return eps;
}

double Functional::vwn_depsp(double x) {
    return vwn_deps(x,0.0310907,-0.10498,3.72744,12.9352);
}

double Functional::vwn_depsf(double x) {
    return vwn_deps(x,0.01554535,-0.32500,7.06042,13.0045);
}

double Functional::vwn_deps(double x, double a, double x0, double b, double c) {
    double q = sqrt(4.0*c-b*b);
    double deps = a*(2.0/x - (2.0*x+b)/vwn_xx(x,b,c)
              - 4.0*b/(pow(2.0*x+b,2.0)+q*q) - (b*x0/vwn_xx(x0,b,c))
              * (2.0/(x-x0)-(2.0*x+b)/vwn_xx(x,b,c)-4.0*(2.0*x0+b)/(pow(2.0*x+b,2.0)+q*q)));
    return deps;
}

double Functional::vwn_g(double z) {
    return 1.125*(pow(1.0+z,4.0/3.0)+pow(1.0-z,4.0/3.0)-2.0);
}

double Functional::vwn_dg(double z) {
    return 1.5*(pow(1.0+z,1.0/3.0)-pow(1.0-z,1.0/3.0));
}
