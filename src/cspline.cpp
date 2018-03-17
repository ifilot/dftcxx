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

#include "cspline.h"

/**
 * @brief      default constructor
 */
Cspline::Cspline() :
spline_generated(false) {}

/**
 * @brief      Set the x and y values used for spline construction.
 *
 * @param[in]  xx    X-axis values.
 * @param[in]  yy    Y-axis values.
 */
void Cspline::set_values(const std::vector<double>& xx, const std::vector<double>& yy) {
    x = xx;
    y = yy;
    if(x.size() != y.size()) {
        std::stringstream error;
        error << "Cspline x and y data ranges do not match: "<< x.size() << " != " << y.size();
        throw std::runtime_error(error.str());
    }
    if(x.size() < 3) {
        std::stringstream error;
        error << "Cspline data range has to contain 3 of more items. Range = "<< x.size();
        throw std::runtime_error(error.str());
    }
    for(unsigned int i=1; i<x.size(); ++i) {
        if(x[i]<= x[i-1]) {
            std::stringstream error;
            error << "Spline x-data should be continuously increasing. @i="<< i+1 << ": " << x[i-1] << "->" << x[i];
            throw std::runtime_error(error.str());
        }
    }
}

/**
 * @brief      Generate a spline through the (x,y) points.
 *
 *             The spline generation is done by setting up a tridiagonal matrix
 *             problem. This matrix problem can be solved with the Thomas
 *             algorithm. The algorithm is generalized to work with points that
 *             are irregularly spaced over the x-axis.
 */
void Cspline::generate_spline() {
    //spline contains n-rows of polynomial coefficients
    spline = std::vector< std::vector<double> >(x.size()-1,std::vector<double>(4,0.0));

    //tridiagonal matrix problem to be solved with the Thomas algorithm
    //|bc00||D0||Y0|
    //|abc0||D1||Y1|
    //|0abc||D2||Y2|
    //|00ab||D3||Y3|

    std::vector<double> A(x.size(),0.0);
    std::vector<double> B(x.size(),0.0);
    std::vector<double> C(x.size(),0.0);
    std::vector<double> D(x.size(),0.0);
    std::vector<double> Y(x.size(),0.0);
    double h0;
    double h1;
    double r0;
    double r1;
    // FILE_LOG(logDEBUG) << "generate b_0, c_0, and Y_0";
    h0 = x[1]-x[0];
    h1 = x[2]-x[1];
    r0 = (y[1]-y[0])/h0;
    r1 = (y[2]-y[1])/h1;
    B[0] = h1*(h0+h1);
    C[0] = (h0+h1)*(h0+h1);
    // C[0] = h0;
    Y[0] = r0*(3*h0*h1+2*h1*h1)+r1*h0*h0;
    //generate a_i, b_i, c_i, and Y_i
    for(unsigned int i=1; i<x.size()-1; ++i) {
        h0 = x[i]-x[i-1];
        h1 = x[i+1]-x[i];
        r0 = (y[i]-y[i-1])/h0;
        r1 = (y[i+1]-y[i])/h1;
        A[i] = h1;
        B[i] = 2*(h0+h1);
        C[i] = h0;
        Y[i] =3*(r0*h1+r1*h0);
    }
    //generate a_n, b_n, and Y_n
    A[x.size()-1] = (h0+h1)*(h0+h1);
    // A[x.size()-1] = h1;
    B[x.size()-1] = h0*(h0+h1);
    Y[x.size()-1] = r0*h1*h1+r1*(3*h0*h1+2*h0*h0);

    //get modified c's after swipe
    C[0] = C[0]/B[0];
    for(unsigned int i=1; i<x.size()-1; ++i) {
        C[i] = C[i] / (B[i] - A[i]*C[i-1]);
    }

    //get modified Y's after swipe
    Y[0] = Y[0]/B[0];
    for(unsigned int i=1; i<x.size(); ++i) {
        Y[i] = (Y[i]- A[i]*Y[i-1]) / (B[i] - A[i]*C[i-1]);
    }

    //back substitute to get D's (which are the derivatives at each data point)
    D[x.size()-1] = Y[x.size()-1];
    for(unsigned int i=x.size()-1; i>0; i--) {
        D[i-1] = Y[i-1] - C[i-1] * D[i];
    }
    //the matrix problem has been solved and can be used to get the polynomial coefficients
    //yi(x) = ai + bi*x + ci*x^2 + di*x^3
    double dx;
    double dy;
    for(unsigned int i=0; i<x.size()-1; ++i) {
        dx = 1.0 / (x[i+1]-x[i]);
        dy = (y[i+1]-y[i])*dx;
        spline[i][0] = y[i]; //ai
        spline[i][1] = D[i]; //bi
        spline[i][2] = dx*(3*dy - 2*D[i] - D[i+1]); //ci
        spline[i][3] = dx*dx*(-2*dy + D[i] + D[i+1]); //di
    }

    this->spline_generated = true;
}

/**
 * @brief      Evaluate the spline at a point on the x-axis.
 *
 * @param[in]  xx    Point of evaluation.
 *
 * @return     Spline value at point of evaluation.
 */
double Cspline::eval(double xx) const {
    if(!this->spline_generated) {
        throw std::runtime_error("Spline interpolation was not performed, cannot evaluate.");
    }

    if(xx<x[0]) {
        std::stringstream error;
        error << "Spline starts at x="<< x[0] << " but is evaluated at x=" << xx;
        throw std::runtime_error(error.str());
    }
    if(xx >= x.back()) {
        return y.back();
    }

    //find correct x-data interval to interpolate and return the interpolated value
    for(unsigned int i=1; i<x.size(); ++i) {
        if(xx<=x[i]) {
            double xpar = (xx - x[i-1]);
            return spline[i-1][0] + spline[i-1][1]*xpar + spline[i-1][2]*xpar*xpar + spline[i-1][3]*xpar*xpar*xpar;
        }
    }

    return 0.0;
}

/**
 * @brief      print the data points
 */
void Cspline::print() const {
    for(unsigned int i=0; i<x.size(); i++) {
        std::cout << x[i] << "\t" << y[i] << std::endl;
    }
}
