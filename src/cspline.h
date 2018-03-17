/*************************************************************************
 *
 *  This file is part of DFTCXX.
 *
 *  Author: Ivo Filot <i.a.w.filot@tue.nl>
 *
 *  DFTCXX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  DFTCXX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with DFTCXX.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************/

#ifndef _CSPLINE_H
#define _CSPLINE_H

#include <vector>
#include <math.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <sstream>

/**
 * @brief      Cubic spline class.
 */
class Cspline {
private:
std::vector<double> x;  // x and y values with the same indices are datapoints
std::vector<double> y;
std::vector< std::vector<double> > spline;  //contains n-rows of polynomial coefficients
bool spline_generated;

public:
    /**
     * @brief      default constructor
     */
    Cspline();

    /**
     * @brief      Set the x and y values used for spline construction.
     *
     * @param[in]  xx    X-axis values.
     * @param[in]  yy    Y-axis values.
     */
    void set_values(const std::vector<double>& xx, const std::vector<double>& yy);

    /**
     * @brief      Generate a spline through the (x,y) points.
     *
     *             The spline generation is done by setting up a tridiagonal
     *             matrix problem. This matrix problem can be solved with the
     *             Thomas algorithm. The algorithm is generalized to work with
     *             points that are irregularly spaced over the x-axis.
     */
    void generate_spline();

    /**
     * @brief      Evaluate the spline at a point on the x-axis.
     *
     * @param[in]  xx    Point of evaluation.
     *
     * @return     Spline value at point of evaluation.
     */
    double eval(double x) const;

    /**
     * @brief      Gets the minimum.
     *
     * @return     The minimum.
     */
    inline double get_min() const {
        return x[0];
    }

    /**
     * @brief      Gets the maximum.
     *
     * @return     The maximum.
     */
    inline double get_max() const {
        return x.back();
    }

    /**
     * @brief      print the data points
     */
    void print() const;

};

#endif //_CSPLINE_H
