/**************************************************************************
 *   integrals.cpp  --  This file is part of DFTCXX.                      *
 *                                                                        *
 *   Copyright (C) 2016, Ivo Filot                                        *
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

#include "integrals.h"

/**
 * @fn Integrator
 * @brief Integrator constructor method
 *
 * @return Integrator class
 */
Integrator::Integrator(){}

/**
 * @fn overlap
 * @brief Calculates overlap integral of two CGF
 *
 * @param const CGF& cgf1   Contracted Gaussian Function
 * @param const CGF& cgf2   Contracted Gaussian Function
 *
 * Calculates the value of < cgf1 | cgf2 >
 *
 * @return double value of the overlap integral
 */
double Integrator::overlap(const CGF& cgf1, const CGF& cgf2) {
    double sum = 0.0;

    // loop over all GTOs inside the CGF, calculate the overlap integrals
    // and sum all the integral values
    for(unsigned int k = 0; k < cgf1.size(); k++) {
        for(unsigned int l = 0; l < cgf2.size(); l++) {
            sum += cgf1.get_norm_gto(k) *
                   cgf2.get_norm_gto(l) *
                   cgf1.get_coefficient_gto(k) *
                   cgf2.get_coefficient_gto(l) *
                   this->overlap(cgf1.get_gto(k), cgf2.get_gto(l) );
        }
    }

    return sum;
}

/**
 * @fn overlap
 * @brief Calculates overlap integral of two GTO
 *
 * @param const GTO& gto1   Gaussian Type Orbital
 * @param const GTO& gto2   Gaussian Type Orbital
 *
 * Calculates the value of < gto1 | gto2 >
 *
 * @return double value of the overlap integral
 */
double Integrator::overlap(const GTO& gto1, const GTO& gto2) {
    return this->overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_position(),
                         gto2.get_alpha(), gto2.get_l(), gto2.get_m(), gto2.get_n(), gto2.get_position());
}

/**
 * @fn kinetic
 * @brief Calculates kinetic integral of two CGF
 *
 * @param const CGF& cgf1   Contracted Gaussian Function
 * @param const CGF& cgf2   Contracted Gaussian Function
 *
 * Calculates the value of < cgf1 | T | cgf2 >
 *
 * @return double value of the kinetic integral
 */
double Integrator::kinetic(const CGF& cgf1, const CGF& cgf2) {
    double sum = 0.0;

    // loop over all GTOs inside the CGF, calculate the kinetic integrals
    // and sum all the integral values
    for(unsigned int k = 0; k < cgf1.size(); k++) {
        for(unsigned int l = 0; l < cgf2.size(); l++) {
            sum += cgf1.get_norm_gto(k) *
                   cgf2.get_norm_gto(l) *
                   cgf1.get_coefficient_gto(k) *
                   cgf2.get_coefficient_gto(l) *
                   this->kinetic(cgf1.get_gto(k), cgf2.get_gto(l) );
        }
    }

    return sum;
}

/**
 * @fn kinetic
 * @brief Calculates kinetic integral of two GTO
 *
 * @param const GTO& gto1   Gaussian Type Orbital
 * @param const GTO& gto2   Gaussian Type Orbital
 *
 * Calculates the value of < gto1 | H | gto2 >
 *
 * @return double value of the kinetic integral
 */
double Integrator::kinetic(const GTO& gto1, const GTO& gto2) {
    double term0 = gto2.get_alpha() *
                   (2 * ( gto2.get_l() + gto2.get_m() + gto2.get_n() ) + 3 ) *
                   this->overlap(gto1, gto2);

    double term1 = -2 * pow(gto2.get_alpha(), 2.0) * (
        this->overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_position(), gto2.get_alpha(), gto2.get_l()+2, gto2.get_m(), gto2.get_n(), gto2.get_position()) +
        this->overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_position(), gto2.get_alpha(), gto2.get_l(), gto2.get_m()+2, gto2.get_n(), gto2.get_position()) +
        this->overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_position(), gto2.get_alpha(), gto2.get_l(), gto2.get_m(), gto2.get_n()+2, gto2.get_position())
    );
    double term2 = -0.5 * (gto2.get_l() * (gto2.get_l() - 1) *
        this->overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_position(), gto2.get_alpha(), gto2.get_l()-2, gto2.get_m(), gto2.get_n(), gto2.get_position()) +
                                                 gto2.get_m() * (gto2.get_m() - 1) *
        this->overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_position(), gto2.get_alpha(), gto2.get_l(), gto2.get_m()-2, gto2.get_n(), gto2.get_position()) +
                                                 gto2.get_n() * (gto2.get_n() - 1) *
        this->overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_position(), gto2.get_alpha(), gto2.get_l(), gto2.get_m(), gto2.get_n()-2, gto2.get_position()) );

    return term0 + term1 + term2;
}

/**
 * @fn nuclear
 * @brief Calculates nuclear integral of two CGF
 *
 * @param const CGF& cgf1       Contracted Gaussian Function
 * @param const CGF& cgf2       Contracted Gaussian Function
 * @param unsigned int charge   charge of the nucleus in a.u.
 *
 * Calculates the value of < cgf1 | V | cgf2 >
 *
 * @return double value of the nuclear integral
 */
double Integrator::nuclear(const CGF& cgf1, const CGF& cgf2, const vec3 &nucleus, unsigned int charge) {
    double sum = 0.0;

    for(unsigned int k = 0; k < cgf1.size(); k++) {
        for(unsigned int l = 0; l < cgf2.size(); l++) {
            sum += cgf1.get_norm_gto(k) *
                   cgf2.get_norm_gto(l) *
                   cgf1.get_coefficient_gto(k) *
                   cgf2.get_coefficient_gto(l) *
                   this->nuclear(cgf1.get_gto(k), cgf2.get_gto(l), nucleus);
        }
    }

    return sum * (double)charge;
}

/**
 * @fn nuclear
 * @brief Calculates nuclear integral of two CGF
 *
 * @param const GTO& gto1       Contracted Gaussian Function
 * @param const GTO& gto2       Contracted Gaussian Function
 * @param unsigned int charge   charge of the nucleus in a.u.
 *
 * Calculates the value of < gto1 | V | gto2 >
 *
 * @return double value of the nuclear integral
 */
double Integrator::nuclear(const GTO& gto1, const GTO& gto2, const vec3 &nucleus) {
    return nuclear(gto1.get_position(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_alpha(),
                   gto2.get_position(), gto2.get_l(), gto2.get_m(), gto2.get_n(), gto2.get_alpha(), nucleus);
}

/**
 * @fn repulsion
 * @brief Calculates two-electron repulsion integral of four CGF
 *
 * @param const CGF& cgf1       Contracted Gaussian Function
 * @param const CGF& cgf2       Contracted Gaussian Function
 * @param const CGF& cgf3       Contracted Gaussian Function
 * @param const CGF& cgf4       Contracted Gaussian Function
 *
 * Calculates the value of < cgf1 | cgf2 | cgf3 | cgf4 >
 *
 * @return double value of the repulsion integral
 */
double Integrator::repulsion(const CGF &cgf1,const CGF &cgf2,const CGF &cgf3,const CGF &cgf4) {
    double sum = 0;

    for(unsigned int i=0; i< cgf1.size(); i++) {
        for(unsigned int j=0; j< cgf2.size(); j++) {
            for(unsigned int k=0; k < cgf3.size(); k++) {
                for(unsigned int l=0; l < cgf4.size(); l++) {
                    double pre = cgf1.get_coefficient_gto(i) * cgf2.get_coefficient_gto(j) * cgf3.get_coefficient_gto(k) * cgf4.get_coefficient_gto(l);
                    sum += pre * repulsion(cgf1.get_gto(i),cgf2.get_gto(j),cgf3.get_gto(k),cgf4.get_gto(l));
                }
            }
        }
    }

    return sum;
}

/**
 * @fn repulsion
 * @brief Calculates two-electron repulsion integral of four CGF
 *
 * @param const GTO& gto1       Contracted Gaussian Function
 * @param const GTO& gto2       Contracted Gaussian Function
 * @param const GTO& gto3       Contracted Gaussian Function
 * @param const GTO& gto4       Contracted Gaussian Function
 *
 * Calculates the value of < gto1 | gto2 | gto3 | gto4 >
 *
 * @return double value of the repulsion integral
 */
double Integrator::repulsion(const GTO &gto1, const GTO &gto2, const GTO &gto3, const GTO &gto4) {

    return repulsion(gto1.get_position(), gto1.get_norm(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_alpha(),
                     gto2.get_position(), gto2.get_norm(), gto2.get_l(), gto2.get_m(), gto2.get_n(), gto2.get_alpha(),
                     gto3.get_position(), gto3.get_norm(), gto3.get_l(), gto3.get_m(), gto3.get_n(), gto3.get_alpha(),
                     gto4.get_position(), gto4.get_norm(), gto4.get_l(), gto4.get_m(), gto4.get_n(), gto4.get_alpha());
}

/**
 * @fn overlap
 * @brief Performs overlap integral evaluation
 *
 * @param double alpha1     Gaussian exponent of the first GTO
 * @param unsigned int l1   Power of x component of the polynomial of the first GTO
 * @param unsigned int m1   Power of y component of the polynomial of the first GTO
 * @param unsigned int n1   Power of z component of the polynomial of the first GTO
 * @param vec3 a            Center of the Gaussian orbital of the first GTO
 * @param double alpha2     Gaussian exponent of the second GTO
 * @param unsigned int l2   Power of x component of the polynomial of the second GTO
 * @param unsigned int m2   Power of y component of the polynomial of the second GTO
 * @param unsigned int n2   Power of z component of the polynomial of the second GTO
 * @param vec3 b            Center of the Gaussian orbital of the second GTO
 *
 * Calculates the value of < gto1 | gto2 >
 *
 * @return double value of the overlap integral
 */
double Integrator::overlap(double alpha1, unsigned int l1, unsigned int m1, unsigned int n1, const vec3 &a,
                           double alpha2, unsigned int l2, unsigned int m2, unsigned int n2, const vec3 &b) {

    static const double pi = 3.14159265359;

    double rab2 = (a-b).squaredNorm();
    double gamma = alpha1 + alpha2;
    vec3 p = this->gaussian_product_center(alpha1, a, alpha2, b);

    double pre = std::pow(pi / gamma, 1.5) * std::exp(-alpha1 * alpha2 * rab2 / gamma);
    double wx = this->overlap_1D(l1, l2, p[0]-a[0], p[0]-b[0], gamma);
    double wy = this->overlap_1D(m1, m2, p[1]-a[1], p[1]-b[1], gamma);
    double wz = this->overlap_1D(n1, n2, p[2]-a[2], p[2]-b[2], gamma);

    return pre * wx * wy * wz;
}

/**
 * @fn overlap_1D
 * @brief Calculates one dimensional overlap integral
 *
 * @param int l1        Power of 'x' component of the polynomial of the first GTO
 * @param int l2        Power of 'x' component of the polynomial of the second GTO
 * @param double x1     'x' component of the position of the first GTO
 * @param double x2     'x' component of the position of the second GTO
 * @param double gamma  Sum of the two Gaussian exponents
 *
 * NOTE: in contrast to other places, here int has to be used rather than unsigned int
 *       because sometimes negative numbers are parsed
 *
 * @return double value of the one dimensional overlap integral
 */
double Integrator::overlap_1D(int l1, int l2, double x1, double x2, double gamma) {
    double sum = 0.0;

    for(int i=0; i < (1 + std::floor(0.5 * (l1 + l2))); i++) {
        sum += this->binomial_prefactor(2*i, l1, l2, x1, x2) *
                     (i == 0 ? 1 : boost::math::double_factorial<double>(2 * i - 1) ) /
                     std::pow(2 * gamma, i);
    }

    return sum;
}

/**
 * @fn gaussian_product_center
 * @brief Calculates the Gaussian product center of two GTOs
 *
 * @param double alpha1     Gaussian exponent of the first GTO
 * @param double alpha2     Gaussian exponent of the second GTO
 * @param const vec3 a      Center of the first GTO
 * @param const vec3 b      Center of the second GTO
 *
 *
 * @return new gaussian product center
 */
vec3 Integrator::gaussian_product_center(double alpha1, const vec3& a,
                                         double alpha2, const vec3& b) {
    return (alpha1 * a + alpha2 * b) / (alpha1 + alpha2);
}

double Integrator::binomial_prefactor(int s, int ia, int ib,
                                      double xpa, double xpb) {
    double sum = 0.0;

    for (int t=0; t < s+1; t++) {
        if ((s-ia <= t) && (t <= ib)) {
            sum += this->binomial(ia, s-t)   *
                   this->binomial(ib, t)     *
                   std::pow(xpa, ia - s + t) *
                   std::pow(xpb, ib - t);
        }
    }

    return sum;
}

double Integrator::binomial(int a, int b) {
    if( (a < 0) | (b < 0) | (a-b < 0) ) {
        return 1.0;
    }
    return boost::math::factorial<double>(a) / (boost::math::factorial<double>(b) * boost::math::factorial<double>(a-b));
}

double Integrator::nuclear(const vec3& a,
               int l1, int m1, int n1, double alpha1,
               const vec3& b,
               int l2, int m2, int n2,
               double alpha2, const vec3& c) {

    static const double pi = 3.14159265359;

    double gamma = alpha1 + alpha2;

    vec3 p = gaussian_product_center(alpha1, a, alpha2, b);
    double rab2 = (a-b).squaredNorm();
    double rcp2 = (c-p).squaredNorm();

    std::vector<double> ax = A_array(l1, l2, p[0]-a[0], p[0]-b[0], p[0]-c[0], gamma);
    std::vector<double> ay = A_array(m1, m2, p[1]-a[1], p[1]-b[1], p[1]-c[1], gamma);
    std::vector<double> az = A_array(n1, n2, p[2]-a[2], p[2]-b[2], p[2]-c[2], gamma);

    double sum = 0.0;

    for(int i=0; i<=l1+l2;i++) {
        for(int j=0; j<=m1+m2;j++) {
            for(int k=0; k<=n1+n2;k++) {
                sum += ax[i] * ay[j] * az[k] * this->gamma_inc.Fgamma(i+j+k,rcp2*gamma);
            }
        }
    }

    return -2.0 * pi / gamma * std::exp(-alpha1*alpha2*rab2/gamma) * sum;
}

std::vector<double> Integrator::A_array(const int l1, const int l2, const double pa, const double pb, const double cp, const double g) {
    int imax = l1 + l2 + 1;
    std::vector<double> arrA(imax, 0);

    for(int i=0; i<imax; i++) {
        for(int r=0; r<=i/2; r++) {
            for(int u=0; u<=(i-2*r)/2; u++) {
                int iI = i - 2 * r - u;
                arrA[iI] += A_term(i, r, u, l1, l2, pa, pb, cp, g);
            }
        }
    }

    return arrA;
}

double Integrator::A_term(const int i, const int r, const int u, const int l1, const int l2, const double pax, const double pbx, const double cpx, const double gamma) {
    return  std::pow(-1,i) * this->binomial_prefactor(i,l1,l2,pax,pbx)*
                    std::pow(-1,u) * boost::math::factorial<double>(i)*std::pow(cpx,i-2*r-2*u)*
                    std::pow(0.25/gamma,r+u)/boost::math::factorial<double>(r)/boost::math::factorial<double>(u)/boost::math::factorial<double>(i-2*r-2*u);
}

double Integrator::repulsion(const vec3 &a, const double norma, const int la, const int ma, const int na, const double alphaa,
                         const vec3 &b, const double normb, const int lb, const int mb, const int nb, const double alphab,
                         const vec3 &c, const double normc, const int lc, const int mc, const int nc, const double alphac,
                         const vec3 &d, const double normd, const int ld, const int md, const int nd, const double alphad) {

    static const double pi = 3.14159265359;
    double rab2 = (a-b).squaredNorm();
    double rcd2 = (c-d).squaredNorm();

    vec3 p = gaussian_product_center(alphaa, a, alphab, b);
    vec3 q = gaussian_product_center(alphac, c, alphad, d);

    double rpq2 = (p-q).squaredNorm();

    double gamma1 = alphaa + alphab;
    double gamma2 = alphac + alphad;
    double delta = 0.25 * (1.0 / gamma1 + 1.0 / gamma2);

    std::vector<double> bx = B_array(la, lb, lc, ld, p[0], a[0], b[0], q[0], c[0], d[0], gamma1, gamma2, delta);
    std::vector<double> by = B_array(ma, mb, mc, md, p[1], a[1], b[1], q[1], c[1], d[1], gamma1, gamma2, delta);
    std::vector<double> bz = B_array(na, nb, nc, nd, p[2], a[2], b[2], q[2], c[2], d[2], gamma1, gamma2, delta);

    double sum = 0.0;
    for(int i=0; i<=(la+lb+lc+ld); i++) {
        for(int j=0; j<=(ma+mb+mc+md); j++) {
            for(int k=0; k<=(na+nb+nc+nd); k++) {
                sum += bx[i]*by[j]*bz[k]*this->gamma_inc.Fgamma(i+j+k,0.25*rpq2/delta);
            }
        }
    }

    return 2.0 * std::pow(pi,2.5)/(gamma1*gamma2*std::sqrt(gamma1+gamma2))*
        std::exp(-alphaa*alphab*rab2/gamma1)*
        std::exp(-alphac*alphad*rcd2/gamma2)*sum*norma*normb*normc*normd;
}

std::vector<double> Integrator::B_array(const int l1,const int l2,const int l3,const int l4,
        const double p, const double a, const double b, const double q, const double c, const double d,
        const double g1, const double g2, const double delta) {

    int imax = l1 + l2 + l3 + l4 + 1;
    std::vector<double> arrB(imax,0);

    for(int i1=0; i1<l1+l2+1; i1++) {
        for(int i2=0; i2<l3+l4+1; i2++) {
            for(int r1=0; r1 < i1/2+1; r1++) {
                for(int r2=0; r2 < i2/2+1; r2++) {
                    for(int u=0; u<(i1+i2)/2-r1-r2+1; u++) {
                        int i = i1+i2-2*(r1+r2)-u;
                        arrB[i] += B_term(i1,i2,r1,r2,u,l1,l2,l3,l4,
                                                            p,a,b,q,c,d,g1,g2,delta);
                    }
                }
            }
        }
    }
    return arrB;
}

double Integrator::B_term(const int i1, const int i2, const int r1, const int r2, const int u, const int l1, const int l2, const int l3, const int l4,
        const double px, const double ax, const double bx, const double qx, const double cx, const double dx, const double gamma1,
        const double gamma2, const double delta) {
    return fB(i1,l1,l2,px,ax,bx,r1,gamma1)*
        pow(-1,i2) * fB(i2,l3,l4,qx,cx,dx,r2,gamma2)*
        pow(-1,u)*fact_ratio2(i1+i2-2*(r1+r2),u)*
        pow(qx-px,i1+i2-2*(r1+r2)-2*u)/
        pow(delta,i1+i2-2*(r1+r2)-u);
}

double Integrator::fB(const int i, const int l1, const int l2, const double p, const double a, const double b, const int r, const double g) {
    return binomial_prefactor(i, l1, l2, p-a, p-b) * B0(i, r, g);
}

double Integrator::B0(int i, int r, double g) {
    return fact_ratio2(i,r) * pow(4*g,r-i);
}

double Integrator::fact_ratio2(const int a, const int b) {
    return boost::math::factorial<double>(a) / boost::math::factorial<double>(b) / boost::math::factorial<double>(a - 2*b);
}

const unsigned int Integrator::teindex(unsigned int i, unsigned int j, unsigned int k, unsigned int l) {
    if(i < j) {
        std::swap(i,j);
    }
    if(k < l) {
        std::swap(k,l);
    }

    unsigned int ij = i * (i + 1) / 2 + j;
    unsigned int kl = k * (k + 1) / 2 + l;

    if(ij < kl) {
        std::swap(ij,kl);
    }

    return ij * (ij + 1) / 2 + kl;
}
