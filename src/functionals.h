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

#ifndef _FUNCTIONALS_H
#define _FUNCTIONALS_H

#include <Eigen/Dense>
#include <boost/math/special_functions/factorials.hpp>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXXd;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;

class Functional {
public:
    Functional();

    // functionals
    void xalpha_x_functional(const VectorXd& densitiesa,
                             const VectorXd& densitiesb,
                             VectorXd& ex,
                             VectorXd& vxa,
                             VectorXd& vxb);

    void vwm_c_functional(const VectorXd& densitiesa,
                          const VectorXd& densitiesb,
                          VectorXd& ec,
                          VectorXd& vca,
                          VectorXd& vcb);

private:
    static constexpr double pi = 3.14159265358979323846;

    // auxiliary functions
    double vwn_xx(double x, double b, double c);
    double vwn_epsp(double x);
    double vwn_epsf(double x);
    double vwn_eps(double x, double a, double x0, double b, double c);
    double vwn_depsp(double x);
    double vwn_depsf(double x);
    double vwn_deps(double x, double a, double x0, double b, double c);
    double vwn_g(double z);
    double vwn_dg(double z);
};

#endif //_FUNCTIONALS_H
