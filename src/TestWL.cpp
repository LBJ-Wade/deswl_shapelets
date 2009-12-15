
#include <valarray>
#include <fstream>
#include <iostream>
#include <cmath>

#include "BVec.h"
#include "Ellipse.h"
#include "EllipseSolver.h"
#include "TMV.h"
#include "TMV_Small.h"
#include "dbg.h"
#include "TestHelper.h"
#include "Log.h"
#include "InputCatalog.h"
#include "StarCatalog.h"
#include "PsfCatalog.h"
#include "FittedPsf.h"
#include "ShearCatalog.h"

#ifdef _OPENMP
#include "omp.h"
#endif

bool shouldShowTests = false;
bool shouldThrow = true;
std::string lastSuccess = "";

const double PI = 3.14159265359;
const double sqrtpi = sqrt(PI);

std::ostream* dbgout = 0;
bool XDEBUG = true;
//bool XDEBUG = false;

#if 0
#define TEST1
#define TEST2
#define TEST3
#define TEST4
#define TEST5
#endif
#define TEST6

#ifdef TEST1
#define TEST12
#define TEST12345
#endif
#ifdef TEST2
#define TEST12
#define TEST23
#define TEST234
#define TEST12345
#endif
#ifdef TEST3
#define TEST23
#define TEST234
#define TEST12345
#endif
#ifdef TEST4
#define TEST234
#define TEST45
#define TEST12345
#endif
#ifdef TEST5
#define TEST45
#define TEST12345
#endif

// The tests of the jacobians are pretty time consuming.
// They are worth testing, but once the code is working, I usually turn
// turn them off (by commenting out the next line) to save time.
#define TESTJ

inline void getFakePixList(
    PixelList& pix,
    double xcen, double ycen,
    const tmv::SmallMatrix<double,2,2>& D,
    double aperture, const BVec& b, const Ellipse& ell)
{
    const size_t border = b.getOrder();
    Assert(border <= 8); // This is as high as I have implemented.

    double det = std::abs(D.Det());
    //double pixScale = sqrt(det); // arcsec/pixel

    // xap,yap are the maximum deviation from the center in x,y
    // such that u^2+v^2 = aperture^2
    double xap = aperture / det * sqrt(D(0,0)*D(0,0) + D(0,1)*D(0,1));
    double yap = aperture / det * sqrt(D(1,0)*D(1,0) + D(1,1)*D(1,1));

    // xcen,ycen are given on a 1-based grid.
    // ie. where the lower left corner pixel is (1,1), rather than (0,0).
    // The easiest way to do this is to just decrease xcen,ycen by 1 each:
    --xcen; --ycen;

    int i1 = int(floor(xcen-xap));
    int i2 = int(ceil(xcen+xap));
    int j1 = int(floor(ycen-yap));
    int j2 = int(ceil(ycen+yap));

    double apsq = aperture*aperture;

    std::complex<double> zc = ell.getCen();
    std::complex<double> gamma = ell.getGamma();
    double mu = ell.getMu();

    double m = exp(-mu)/b.getSigma()/sqrt(1.-std::norm(gamma));

    const int NV = 1;
    double varlist[NV] = {0.1};
    // MJ: The non-unifiorm variance seems to lead to problems with 
    // things like applyZ and such.  Not sure why yet.
    // This could be a sign of a problem...
    //const int NV = 10;
    //double varlist[NV] = {1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9};
    //double varlist[NV] = {0.2,1.3,5.3,2.1,3.8,0.3,1.6,2.6,1.9,0.8};
    int ivar = 0;

    double chipx = i1-xcen;
    for(int i=i1;i<=i2;++i,chipx+=1.) {
        double chipy = j1-ycen;
        double u = D(0,0)*chipx+D(0,1)*chipy;
        double v = D(1,0)*chipx+D(1,1)*chipy;
        for(int j=j1;j<=j2;++j,u+=D(0,1),v+=D(1,1)) {
            std::complex<double> globz(u,v);
            // globz is in arcsec
            double rsq = std::norm(globz);
            if (rsq > apsq) continue;
            // Set up test intensity pattern:
            double var = varlist[ivar];
            if (++ivar == NV) ivar = 0;

            std::complex<double> zz = globz - zc;
            std::complex<double> z = m * (zz - gamma*std::conj(zz));
            rsq = std::norm(z);
            double expfactor = exp(-rsq/2.)/sqrtpi;

            // The factors of 2 and -2 for m!-0 terms below are due to the fact
            // that we want real I values, so we assume b(q,p) = b*(p,q).
            // So I_pq = b(p,q) psi(p,q) + b(q,p) psi(q,p) 
            //         = 2 * Re(b(p,q) psi(p,q))
            //         = 2 * Re(b(p,q)) Re(psi(p,q)) - 2 * Im(b(p,q)) Im(psi(p,q))

            tmv::Vector<double>::const_iterator bi = b.vec().begin();
            // 00
            double I = *bi;
            do {
                // This construct is just a means to have something to 
                // break out of from various places depending on the value of n.
                // It is a bit more readable than a bunch of if statements.

                if (border==0) break;

                // 10r
                I += *(++bi)*(2.*std::real(z));
                // 10i
                I += *(++bi)*(-2.*std::imag(z));
                if (border==1) break;

                // 20r
                std::complex<double> z2 = z*z;
                double temp = 2./sqrt(2.);
                I += *(++bi)*(temp*std::real(z2));
                // 20i
                I += *(++bi)*(-temp*std::imag(z2));
                // 11
                I += *(++bi)*(-(1.-rsq));
                if (border==2) break;

                // 30r
                std::complex<double> z3 = z2*z;
                temp = 2./sqrt(6.);
                I += *(++bi)*(temp*std::real(z3));
                // 30i
                I += *(++bi)*(-temp*std::imag(z3));
                // 21r
                temp = -2.*(2.-rsq)/sqrt(2.);
                I += *(++bi)*(temp*std::real(z));
                // 21i
                I += *(++bi)*(-temp*std::imag(z));
                if (border==3) break;

                // 40r
                std::complex<double> z4 = z3*z;
                double rsq2 = rsq*rsq;
                temp = 2./sqrt(24.);
                I += *(++bi)*(temp*std::real(z4));
                // 40i
                I += *(++bi)*(-temp*std::imag(z4));
                // 31r
                temp = -2.*(3.-rsq)/sqrt(6.);
                I += *(++bi)*(temp*std::real(z2));
                // 31i
                I += *(++bi)*(-temp*std::imag(z2));
                // 22
                I += *(++bi)*(rsq2-4*rsq+2.)/2.;
                if (border==4) break;

                // 50r
                std::complex<double> z5 = z4*z;
                temp = 2./sqrt(120.);
                I += *(++bi)*(temp*std::real(z5));
                // 50i
                I += *(++bi)*(-temp*std::imag(z5));
                // 41r
                temp = -2.*(4.-rsq)/sqrt(24.);
                I += *(++bi)*(temp*std::real(z3));
                // 41i
                I += *(++bi)*(-temp*std::imag(z3));
                // 32r
                temp = 2.*(3.-3.*rsq+0.5*rsq2)/sqrt(3.);
                I += *(++bi)*(temp*std::real(z));
                // 32i
                I += *(++bi)*(-temp*std::imag(z));
                if (border==5) break;

                // 60r
                std::complex<double> z6 = z5*z;
                double rsq3 = rsq2*rsq;
                temp = 2./sqrt(720.);
                I += *(++bi)*(temp*std::real(z6));
                // 60i
                I += *(++bi)*(-temp*std::imag(z6));
                // 51r
                temp = -2.*(5.-rsq)/sqrt(120.);
                I += *(++bi)*(temp*std::real(z4));
                // 51i
                I += *(++bi)*(-temp*std::imag(z4));
                // 42r
                temp = 2.*(6.-4.*rsq+0.5*rsq2)/sqrt(12.);
                I += *(++bi)*(temp*std::real(z2));
                // 42i
                I += *(++bi)*(-temp*std::imag(z2));
                // 33
                I += *(++bi)*(-(1.-3.*rsq+1.5*rsq2-rsq3/6.));
                if (border==6) break;

                // 70r
                std::complex<double> z7 = z6*z;
                temp = 2./sqrt(5040.);
                I += *(++bi)*(temp*std::real(z7));
                // 70i
                I += *(++bi)*(-temp*std::imag(z7));
                // 61r
                temp = -2.*(6.-rsq)/sqrt(720.);
                I += *(++bi)*(temp*std::real(z5));
                // 61i
                I += *(++bi)*(-temp*std::imag(z5));
                // 52r
                temp = 2.*(10.-5.*rsq+0.5*rsq2)/sqrt(60.);
                I += *(++bi)*(temp*std::real(z3));
                // 52i
                I += *(++bi)*(-temp*std::imag(z3));
                // 43r
                temp = -2.*(4.-6.*rsq+2.*rsq2-rsq3/6.)/2.;
                I += *(++bi)*(temp*std::real(z));
                // 43i
                I += *(++bi)*(-temp*std::imag(z));
                if (border==7) break;

                // 80r
                std::complex<double> z8 = z7*z;
                double rsq4 = rsq3*rsq;
                temp = 2./sqrt(40320.);
                I += *(++bi)*(temp*std::real(z8));
                // 80i
                I += *(++bi)*(-temp*std::imag(z8));
                // 71r
                temp = -2.*(7.-rsq)/sqrt(5040.);
                I += *(++bi)*(temp*std::real(z6));
                // 71i
                I += *(++bi)*(-temp*std::imag(z6));
                // 62r
                temp = 2.*(15.-6.*rsq+0.5*rsq2)/sqrt(360.);
                I += *(++bi)*(temp*std::real(z4));
                // 62i
                I += *(++bi)*(-temp*std::imag(z4));
                // 53r
                temp = -2.*(10.-10.*rsq+2.5*rsq2-rsq3/6.)/sqrt(20.);
                I += *(++bi)*(temp*std::real(z2));
                // 53i
                I += *(++bi)*(-temp*std::imag(z2));
                // 44
                I += *(++bi)*(1.-4.*rsq+3.*rsq2-2.*rsq3/3.+rsq4/24.);
                Assert(border==8);
            } while (false);
            // breaks go to here.

            I *= expfactor;

            pix.push_back(Pixel(std::real(globz),std::imag(globz),I,var));
        }
    }
}

inline void RtoC(const BVec& b, tmv::Vector<std::complex<double> >& bc)
{
    const size_t border = b.getOrder();
    Assert(bc.size() == (border+2)*(border+2)/4);

    size_t kr=0, kc=0;
    for(size_t ppq=0;ppq<=border;++ppq) {
        for(int p=ppq,q=0;p>=q;--p,++q,++kc,++kr) {
            if (p>q) {
                bc(kc) = std::complex<double>(b(kr),b(kr+1)); ++kr;
            } else {
                bc(kc) = b(kr);
            }
        }
    }
    Assert(kr == b.size());
    Assert(kc == bc.size());
}

inline void CtoR(const tmv::Vector<std::complex<double> >& bc, BVec& b)
{
    const size_t border = b.getOrder();
    Assert(bc.size() == (border+2)*(border+2)/4);

    size_t kr=0, kc=0;
    for(size_t ppq=0;ppq<=border;++ppq) {
        for(int p=ppq,q=0;p>=q;--p,++q,++kc,++kr) {
            b(kr) = std::real(bc[kc]);
            if (p>q) b(++kr) = std::imag(bc[kc]);
            else Assert(std::abs(std::imag(bc[kc])) < 1.e-5);
        }
    }
    Assert(kr == b.size());
    Assert(kc == bc.size());
}

inline void DirectConvolveB(const BVec& rbi, const BVec& rbp, BVec& rb)
    // This is based on direct convolutions of the various psi_pq o psi_st
    // combinations in Maple.  So it is independent of the formulae in 
    // BJ02.
{
    const double twosqrtpi = 2.*sqrtpi;
    Assert(rbi.getOrder() == 4);
    Assert(rbp.getOrder() == 4);
    Assert(rb.getOrder() == 8);

    tmv::Vector<std::complex<double> > bi(9);
    tmv::Vector<std::complex<double> > bp(9);
    tmv::Vector<std::complex<double> > b(25);
    RtoC(rbi,bi);
    RtoC(rbp,bp);

    const size_t K00 = 0;
    const size_t K10 = 1;
    const size_t K20 = 2, K11 = 3;
    const size_t K30 = 4,  K21 = 5;
    const size_t K40 = 6, K31 = 7, K22 = 8;
    const size_t K50 = 9, K41 = 10, K32 = 11;
    const size_t K60 = 12, K51 = 13, K42 = 14, K33 = 15;
    const size_t K70 = 16, K61 = 17, K52 = 18, K43 = 19;
    const size_t K80 = 20, K71 = 21, K62 = 22, K53 = 23, K44 = 24;

    double sigma = rb.getSigma();
    double sigma_i = rbi.getSigma();
    double sigmaPsf = rbp.getSigma();
    Assert(std::abs(sigma-sqrt(pow(sigma_i,2) + pow(sigmaPsf,2)))<1.e-5);

    double fact = 1.;
    for(int ppq=0,k=0;ppq<=4;++ppq) {
        if (ppq>0) fact *= sigma_i/sigma;
        for(int p=ppq,q=0;p>=q;--p,++q,++k)
            bi[k] *= fact;
    }
    fact = 1.;
    for(int spt=0,m=0;spt<=4;++spt) {
        if (spt>0) fact *= sigmaPsf/sigma;
        for(int s=spt,t=0;s>=t;--s,++t,++m)
            bp[m] *= fact;
    }

    double ri = sigma_i/sigmaPsf; ri *= ri;
    double ri2 = ri*ri;
    double rp = sigmaPsf/sigma_i; rp *= rp;
    double rp2 = rp*rp;

    b.Zero();

    b[K00] += bi[K00] * bp[K00];
    b[K10] += bi[K00] * bp[K10];
    b[K20] += bi[K00] * bp[K20];
    b[K30] += bi[K00] * bp[K30];
    b[K40] += bi[K00] * bp[K40];
    b[K11] += bi[K00] * bp[K11];
    b[K00] += ri * bi[K00] * bp[K11];
    b[K21] += bi[K00] * bp[K21];
    b[K10] += sqrt(2.)*ri * bi[K00] * bp[K21];
    b[K31] += bi[K00] * bp[K31];
    b[K20] += sqrt(3.)*ri * bi[K00] * bp[K31];
    b[K22] += bi[K00] * bp[K22];
    b[K11] += 2.*ri * bi[K00] * bp[K22];
    b[K00] += ri*ri * bi[K00] * bp[K22];

    b[K10] += bi[K10] * bp[K00];
    b[K20] += sqrt(2.) * bi[K10] * bp[K10];
    b[K30] += sqrt(3.) * bi[K10] * bp[K20];
    b[K40] += 2. * bi[K10] * bp[K30];
    b[K50] += sqrt(5.) * bi[K10] * bp[K40];
    b[K11] += bi[K10] * std::conj(bp[K10]);
    b[K00] -= bi[K10] * std::conj(bp[K10]);
    b[K21] += sqrt(2.) * bi[K10] * bp[K11];
    b[K10] -= (1.-ri) * bi[K10] * bp[K11];
    b[K31] += sqrt(3.) * bi[K10] * bp[K21];
    b[K20] -= (1.-2.*ri) * bi[K10] * bp[K21];
    b[K41] += 2. * bi[K10] * bp[K31];
    b[K30] -= (1.-3.*ri) * bi[K10] * bp[K31];
    b[K22] += sqrt(2.) * bi[K10] * std::conj(bp[K21]);
    b[K11] -= sqrt(2.)*(1.-ri) * bi[K10] * std::conj(bp[K21]);
    b[K00] -= sqrt(2.)*ri * bi[K10] * std::conj(bp[K21]);
    b[K32] += sqrt(3.) * bi[K10] * bp[K22];
    b[K21] -= sqrt(2.)*(1.-2.*ri) * bi[K10] * bp[K22];
    b[K10] -= (2.-ri)*ri * bi[K10] * bp[K22];

    b[K11] += std::conj(bi[K10]) * bp[K10];
    b[K00] -= std::conj(bi[K10]) * bp[K10];
    b[K21] += std::conj(bi[K10]) * bp[K20];
    b[K10] -= sqrt(2.) * std::conj(bi[K10]) * bp[K20];
    b[K31] += std::conj(bi[K10]) * bp[K30];
    b[K20] -= sqrt(3.) * std::conj(bi[K10]) * bp[K30];
    b[K41] += std::conj(bi[K10]) * bp[K40];
    b[K30] -= 2. * std::conj(bi[K10]) * bp[K40];
    b[K22] += sqrt(2.) * std::conj(bi[K10]) * bp[K21];
    b[K11] -= sqrt(2.)*(1.-ri) * std::conj(bi[K10]) * bp[K21];
    b[K00] -= sqrt(2.)*ri * std::conj(bi[K10]) * bp[K21];
    b[K32] += sqrt(2.) * std::conj(bi[K10]) * bp[K31];
    b[K21] -= sqrt(3.)*(1.-ri) * std::conj(bi[K10]) * bp[K31];
    b[K10] -= sqrt(6.)*ri * std::conj(bi[K10]) * bp[K31];

    b[K20] += bi[K20] * bp[K00];
    b[K30] += sqrt(3.) * bi[K20] * bp[K10];
    b[K40] += sqrt(6.) * bi[K20] * bp[K20];
    b[K50] += sqrt(10.) * bi[K20] * bp[K30];
    b[K60] += sqrt(15.) * bi[K20] * bp[K40];
    b[K21] += bi[K20] * std::conj(bp[K10]);
    b[K10] -= sqrt(2.) * bi[K20] * std::conj(bp[K10]);
    b[K31] += sqrt(3.) * bi[K20] * bp[K11];
    b[K20] -= (2.-ri) * bi[K20] * bp[K11];
    b[K41] += sqrt(6.) * bi[K20] * bp[K21];
    b[K30] -= sqrt(6.)*(1.-ri) * bi[K20] * bp[K21];
    b[K51] += sqrt(10.) * bi[K20] * bp[K31];
    b[K40] -= sqrt(2.)*(2.-3.*ri) * bi[K20] * bp[K31];
    b[K22] += bi[K20] * std::conj(bp[K20]);
    b[K11] -= 2. * bi[K20] * std::conj(bp[K20]);
    b[K00] += bi[K20] * std::conj(bp[K20]);
    b[K32] += sqrt(3.) * bi[K20] * std::conj(bp[K21]);
    b[K21] -= sqrt(2.)*(2.-ri) * bi[K20] * std::conj(bp[K21]);
    b[K10] += (1.-2.*ri) * bi[K20] * std::conj(bp[K21]);
    b[K42] += sqrt(6.) * bi[K20] * bp[K22];
    b[K31] -= 2.*sqrt(3.)*(1.-ri) * bi[K20] * bp[K22];
    b[K20] += (1.-4.*ri+ri2) * bi[K20] * bp[K22];
    b[K33] += sqrt(3.) * bi[K20] * std::conj(bp[K31]);
    b[K22] -= sqrt(3.)*(2.-ri) * bi[K20] * std::conj(bp[K31]);
    b[K11] += sqrt(3.)*(1.-2.*ri) * bi[K20] * std::conj(bp[K31]);
    b[K00] += sqrt(3.)*ri * bi[K20] * std::conj(bp[K31]);

    b[K22] += std::conj(bi[K20]) * bp[K20];
    b[K11] -= 2. * std::conj(bi[K20]) * bp[K20];
    b[K00] += std::conj(bi[K20]) * bp[K20];
    b[K32] += std::conj(bi[K20]) * bp[K30];
    b[K21] -= sqrt(6.) * std::conj(bi[K20]) * bp[K30];
    b[K10] += sqrt(3.) * std::conj(bi[K20]) * bp[K30];
    b[K42] += std::conj(bi[K20]) * bp[K40];
    b[K31] -= 2.*sqrt(2.) * std::conj(bi[K20]) * bp[K40];
    b[K20] += sqrt(6.) * std::conj(bi[K20]) * bp[K40];
    b[K33] += sqrt(3.) * std::conj(bi[K20]) * bp[K31];
    b[K22] -= sqrt(3.)*(2.-ri) * std::conj(bi[K20]) * bp[K31];
    b[K11] += sqrt(3.)*(1.-2.*ri) * std::conj(bi[K20]) * bp[K31];
    b[K00] += sqrt(3.)*ri * std::conj(bi[K20]) * bp[K31];

    b[K30] += bi[K30] * bp[K00];
    b[K40] += 2. * bi[K30] * bp[K10];
    b[K50] += sqrt(10.) * bi[K30] * bp[K20];
    b[K60] += 2.*sqrt(5.) * bi[K30] * bp[K30];
    b[K70] += sqrt(35.) * bi[K30] * bp[K40];
    b[K31] += bi[K30] * std::conj(bp[K10]);
    b[K20] -= sqrt(3.) * bi[K30] * std::conj(bp[K10]);
    b[K41] += 2. * bi[K30] * bp[K11];
    b[K30] -= (3.-ri) * bi[K30] * bp[K11];
    b[K51] += sqrt(10.) * bi[K30] * bp[K21];
    b[K40] -= sqrt(2.)*(3.-2.*ri) * bi[K30] * bp[K21];
    b[K61] += 2.*sqrt(5.) * bi[K30] * bp[K31];
    b[K50] -= sqrt(30.)*(1.-ri) * bi[K30] * bp[K31];
    b[K32] += bi[K30] * std::conj(bp[K20]);
    b[K21] -= sqrt(6.) * bi[K30] * std::conj(bp[K20]);
    b[K10] += sqrt(3.) * bi[K30] * std::conj(bp[K20]);
    b[K42] += 2. * bi[K30] * std::conj(bp[K21]);
    b[K31] -= sqrt(2.)*(3.-ri) * bi[K30] * std::conj(bp[K21]);
    b[K20] += sqrt(6.)*(1.-ri) * bi[K30] * std::conj(bp[K21]);
    b[K52] += sqrt(10.) * bi[K30] * bp[K22];
    b[K41] -= 2.*(3.-2.*ri) * bi[K30] * bp[K22];
    b[K30] += (3.-6.*ri+ri2) * bi[K30] * bp[K22];
    b[K33] += bi[K30] * std::conj(bp[K30]);
    b[K22] -= 3.* bi[K30] * std::conj(bp[K30]);
    b[K11] += 3.* bi[K30] * std::conj(bp[K30]);
    b[K00] -= bi[K30] * std::conj(bp[K30]);
    b[K43] += 2. * bi[K30] * std::conj(bp[K31]);
    b[K32] -= sqrt(3.)*(3.-ri) * bi[K30] * std::conj(bp[K31]);
    b[K21] += 3.*sqrt(2.)*(1.-ri) * bi[K30] * std::conj(bp[K31]);
    b[K10] -= (1.-3.*ri) * bi[K30] * std::conj(bp[K31]);

    b[K33] += std::conj(bi[K30]) * bp[K30];
    b[K22] -= 3. * std::conj(bi[K30]) * bp[K30];
    b[K11] += 3. * std::conj(bi[K30]) * bp[K30];
    b[K00] -= std::conj(bi[K30]) * bp[K30];
    b[K43] += std::conj(bi[K30]) * bp[K40];
    b[K32] -= 2.*sqrt(3.) * std::conj(bi[K30]) * bp[K40];
    b[K21] += 3.*sqrt(2.) * std::conj(bi[K30]) * bp[K40];
    b[K10] -= 2. * std::conj(bi[K30]) * bp[K40];

    b[K40] += bi[K40] * bp[K00];
    b[K50] += sqrt(5.) * bi[K40] * bp[K10];
    b[K60] += sqrt(15.) * bi[K40] * bp[K20];
    b[K70] += sqrt(35.) * bi[K40] * bp[K30];
    b[K80] += sqrt(70.) * bi[K40] * bp[K40];
    b[K41] += bi[K40] * std::conj(bp[K10]);
    b[K30] -= 2. * bi[K40] * std::conj(bp[K10]);
    b[K51] += sqrt(5.) * bi[K40] * bp[K11];
    b[K40] -= (4.-ri) * bi[K40] * bp[K11];
    b[K61] += sqrt(15.) * bi[K40] * bp[K21];
    b[K50] -= sqrt(10.)*(2.-ri) * bi[K40] * bp[K21];
    b[K71] += sqrt(35.) * bi[K40] * bp[K31];
    b[K60] -= sqrt(5.)*(4.-3.*ri) * bi[K40] * bp[K31];
    b[K42] += bi[K40] * std::conj(bp[K20]);
    b[K31] -= 2.*sqrt(2.) * bi[K40] * std::conj(bp[K20]);
    b[K20] += sqrt(6.) * bi[K40] * std::conj(bp[K20]);
    b[K52] += sqrt(5.) * bi[K40] * std::conj(bp[K21]);
    b[K41] -= sqrt(2.)*(4.-ri) * bi[K40] * std::conj(bp[K21]);
    b[K30] += sqrt(2.)*(3.-2.*ri) * bi[K40] * std::conj(bp[K21]);
    b[K62] += sqrt(15.) * bi[K40] * bp[K22];
    b[K51] -= 2.*sqrt(5.)*(2.-ri) * bi[K40] * bp[K22];
    b[K40] += (6.-8.*ri+ri2) * bi[K40] * bp[K22];
    b[K43] += bi[K40] * std::conj(bp[K30]);
    b[K32] -= 2.*sqrt(3.) * bi[K40] * std::conj(bp[K30]);
    b[K21] += 3.*sqrt(2.) * bi[K40] * std::conj(bp[K30]);
    b[K10] -= 2.*bi[K40] * std::conj(bp[K30]);
    b[K53] += sqrt(5.) * bi[K40] * std::conj(bp[K31]);
    b[K42] -= sqrt(3.)*(4.-ri) * bi[K40] * std::conj(bp[K31]);
    b[K31] += sqrt(6.)*(3.-2.*ri) * bi[K40] * std::conj(bp[K31]);
    b[K20] -= sqrt(2.)*(2.-3.*ri) * bi[K40] * std::conj(bp[K31]);
    b[K44] += bi[K40] * std::conj(bp[K40]);
    b[K33] -= 4. * bi[K40] * std::conj(bp[K40]);
    b[K22] += 6. * bi[K40] * std::conj(bp[K40]);
    b[K11] -= 4. * bi[K40] * std::conj(bp[K40]);
    b[K00] += bi[K40] * std::conj(bp[K40]);

    b[K44] += std::conj(bi[K40]) * bp[K40];
    b[K33] -= 4. * std::conj(bi[K40]) * bp[K40];
    b[K22] += 6. * std::conj(bi[K40]) * bp[K40];
    b[K11] -= 4. * std::conj(bi[K40]) * bp[K40];
    b[K00] += std::conj(bi[K40]) * bp[K40];

    b[K11] += bi[K11] * bp[K00];
    b[K00] += rp * bi[K11] * bp[K00];
    b[K21] += sqrt(2.) * bi[K11] * bp[K10];
    b[K10] += (rp-1.) * bi[K11] * bp[K10];
    b[K31] += sqrt(3.) * bi[K11] * bp[K20];
    b[K20] += (rp-2.) * bi[K11] * bp[K20];
    b[K41] += 2. * bi[K11] * bp[K30];
    b[K30] += (rp-3.) * bi[K11] * bp[K30];
    b[K51] += sqrt(5.) * bi[K11] * bp[K40];
    b[K40] += (rp-4.) * bi[K11] * bp[K40];
    b[K22] += 2. * bi[K11] * bp[K11];
    b[K11] += (rp-2.+ri) * bi[K11] * bp[K11];
    b[K00] += 2. * bi[K11] * bp[K11];
    b[K32] += sqrt(6.) * bi[K11] * bp[K21];
    b[K21] += (rp-3.+2.*ri) * bi[K11] * bp[K21];
    b[K10] += sqrt(2.)*(2.-ri) * bi[K11] * bp[K21];
    b[K42] += 2.*sqrt(2.) * bi[K11] * bp[K31];
    b[K31] += (rp-4.+3.*ri) * bi[K11] * bp[K31];
    b[K20] += 2.*sqrt(3.)*(1.-ri) * bi[K11] * bp[K31];
    b[K33] += 3. * bi[K11] * bp[K22];
    b[K22] += (rp-4.+4.*ri) * bi[K11] * bp[K22];
    b[K11] += (4.-4.*ri+ri2) * bi[K11] * bp[K22];
    b[K00] += 3.*ri * bi[K11] * bp[K22];

    b[K21] += bi[K21] * bp[K00];
    b[K10] += sqrt(2.)*rp * bi[K21] * bp[K00];
    b[K31] += sqrt(3.) * bi[K21] * bp[K10];
    b[K20] += (2.*rp-1.) * bi[K21] * bp[K10];
    b[K41] += sqrt(6.) * bi[K21] * bp[K20];
    b[K30] += sqrt(6.)*(rp-1.) * bi[K21] * bp[K20];
    b[K51] += sqrt(10.) * bi[K21] * bp[K30];
    b[K40] += sqrt(2.)*(2.*rp-3.) * bi[K21] * bp[K30];
    b[K61] += sqrt(15.) * bi[K21] * bp[K40];
    b[K50] += sqrt(10.)*(rp-2.) * bi[K21] * bp[K40];
    b[K22] += sqrt(2.) * bi[K21] * std::conj(bp[K10]);
    b[K11] += sqrt(2.)*(rp-1.) * bi[K21] * std::conj(bp[K10]);
    b[K00] -= sqrt(2.)*rp * bi[K21] * std::conj(bp[K10]);
    b[K32] += sqrt(6.) * bi[K21] * bp[K11];
    b[K21] += (2.*rp-3.+ri) * bi[K21] * bp[K11];
    b[K10] -= sqrt(2.)*(rp-2.) * bi[K21] * bp[K11];
    b[K42] += 2.*sqrt(3.) * bi[K21] * bp[K21];
    b[K31] += sqrt(6.)*(rp-2.+ri) * bi[K21] * bp[K21];
    b[K20] -= sqrt(2.)*(rp-4.+ri) * bi[K21] * bp[K21];
    b[K52] += 2.*sqrt(5.) * bi[K21] * bp[K31];
    b[K41] += sqrt(2.)*(2.*rp-5.+3.*ri) * bi[K21] * bp[K31];
    b[K30] -= sqrt(2.)*(rp-6.+3.*ri) * bi[K21] * bp[K31];
    b[K33] += 3. * bi[K21] * std::conj(bp[K21]);
    b[K22] += (2.*rp-5.+2.*ri) * bi[K21] * std::conj(bp[K21]);
    b[K11] -= (2.*rp-5.+2.*ri) * bi[K21] * std::conj(bp[K21]);
    b[K00] -= 3. * bi[K21] * std::conj(bp[K21]);
    b[K43] += 3.*sqrt(2.) * bi[K21] * bp[K22];
    b[K32] += sqrt(6.)*(rp-3.+2.*ri) * bi[K21] * bp[K22];
    b[K21] -= (2.*rp-9.+6.*ri-ri2) * bi[K21] * bp[K22];
    b[K10] -= 3.*sqrt(2.)*(1.-ri) * bi[K21] * bp[K22];

    b[K22] += sqrt(2.) * std::conj(bi[K21]) * bp[K10];
    b[K11] += sqrt(2.)*(rp-1.) * std::conj(bi[K21]) * bp[K10];
    b[K00] -= sqrt(2.)*rp * std::conj(bi[K21]) * bp[K10];
    b[K32] += sqrt(3.) * std::conj(bi[K21]) * bp[K20];
    b[K21] += sqrt(2.)*(rp-2.) * std::conj(bi[K21]) * bp[K20];
    b[K10] -= (2.*rp-1.) * std::conj(bi[K21]) * bp[K20];
    b[K42] += 2. * std::conj(bi[K21]) * bp[K30];
    b[K31] += sqrt(2.)*(rp-3.) * std::conj(bi[K21]) * bp[K30];
    b[K20] -= sqrt(6.)*(rp-1.) * std::conj(bi[K21]) * bp[K30];
    b[K52] += sqrt(5.) * std::conj(bi[K21]) * bp[K40];
    b[K41] += sqrt(2.)*(rp-4.) * std::conj(bi[K21]) * bp[K40];
    b[K30] -= sqrt(2.)*(2.*rp-3.) * std::conj(bi[K21]) * bp[K40];
    b[K33] += 3. * std::conj(bi[K21]) * bp[K21];
    b[K22] += (2.*rp-5.+2.*ri) * std::conj(bi[K21]) * bp[K21];
    b[K11] -= (2.*rp-5.+2.*ri) * std::conj(bi[K21]) * bp[K21];
    b[K00] -= 3. * std::conj(bi[K21]) * bp[K21];
    b[K43] += 2.*sqrt(3.) * std::conj(bi[K21]) * bp[K31];
    b[K32] += (2.*rp-7.+3.*ri) * std::conj(bi[K21]) * bp[K31];
    b[K21] -= sqrt(6.)*(rp-3.+2.*ri) * std::conj(bi[K21]) * bp[K31];
    b[K10] -= sqrt(3.)*(3.-ri) * std::conj(bi[K21]) * bp[K31];

    b[K31] += bi[K31] * bp[K00];
    b[K20] += sqrt(3.)*rp * bi[K31] * bp[K00];
    b[K41] += 2. * bi[K31] * bp[K10];
    b[K30] += (3.*rp-1.) * bi[K31] * bp[K10];
    b[K51] += sqrt(10.) * bi[K31] * bp[K20];
    b[K40] += sqrt(2.)*(3.*rp-2.) * bi[K31] * bp[K20];
    b[K61] += 2.*sqrt(5.) * bi[K31] * bp[K30];
    b[K50] += sqrt(30.)*(rp-1.) * bi[K31] * bp[K30];
    b[K71] += sqrt(35.) * bi[K31] * bp[K40];
    b[K60] += sqrt(5.)*(3.*rp-4.) * bi[K31] * bp[K40];
    b[K32] += sqrt(2.) * bi[K31] * std::conj(bp[K10]);
    b[K21] += sqrt(3.)*(rp-1.) * bi[K31] * std::conj(bp[K10]);
    b[K10] -= sqrt(6.)*rp * bi[K31] * std::conj(bp[K10]);
    b[K42] += 2.*sqrt(2.) * bi[K31] * bp[K11];
    b[K31] += (3.*rp-4.+ri) * bi[K31] * bp[K11];
    b[K20] -= 2.*sqrt(3.)*(rp-1.) * bi[K31] * bp[K11];
    b[K52] += 2.*sqrt(5.) * bi[K31] * bp[K21];
    b[K41] += sqrt(2.)*(3.*rp-5.+2.*ri) * bi[K31] * bp[K21];
    b[K30] -= sqrt(2.)*(3.*rp-6.+ri) * bi[K31] * bp[K21];
    b[K62] += 2.*sqrt(10.) * bi[K31] * bp[K31];
    b[K51] += sqrt(30.)*(rp-2.+ri) * bi[K31] * bp[K31];
    b[K40] -= 2.*sqrt(6.)*(rp-3.+ri) * bi[K31] * bp[K31];
    b[K33] += sqrt(3.) * bi[K31] * std::conj(bp[K20]);
    b[K22] += sqrt(3.)*(rp-2.) * bi[K31] * std::conj(bp[K20]);
    b[K11] -= sqrt(3.)*(2.*rp-1.) * bi[K31] * std::conj(bp[K20]);
    b[K00] += sqrt(3.)*rp * bi[K31] * std::conj(bp[K20]);
    b[K43] += 2.*sqrt(3.) * bi[K31] * std::conj(bp[K21]);
    b[K32] += (3.*rp-7.+2.*ri) * bi[K31] * std::conj(bp[K21]);
    b[K21] -= sqrt(6.)*(2.*rp-3.+ri) * bi[K31] * std::conj(bp[K21]);
    b[K10] += sqrt(3.)*(rp-3.) * bi[K31] * std::conj(bp[K21]);
    b[K53] += sqrt(30.) * bi[K31] * bp[K22];
    b[K42] += sqrt(2.)*(3.*rp-8.+4.*ri) * bi[K31] * bp[K22];
    b[K31] -= (6.*rp-15.+8.*ri-ri2) * bi[K31] * bp[K22];
    b[K20] += sqrt(3.)*(rp-6.+3.*ri) * bi[K31] * bp[K22];
    b[K44] += 4. * bi[K31] * std::conj(bp[K31]);
    b[K33] += (3.*rp-10.+3.*ri) * bi[K31] * std::conj(bp[K31]);
    b[K22] -= 6.*(rp-2.+ri) * bi[K31] * std::conj(bp[K31]);
    b[K11] += (3.*rp-10.+3.*ri) * bi[K31] * std::conj(bp[K31]);
    b[K00] += 4. * bi[K31] * std::conj(bp[K31]);

    b[K33] += sqrt(3.) * std::conj(bi[K31]) * bp[K20];
    b[K22] += sqrt(3.)*(rp-2.) * std::conj(bi[K31]) * bp[K20];
    b[K11] -= sqrt(3.)*(2.*rp-1.) * std::conj(bi[K31]) * bp[K20];
    b[K00] += sqrt(3.)*rp * std::conj(bi[K31]) * bp[K20];
    b[K43] += 2. * std::conj(bi[K31]) * bp[K30];
    b[K32] += sqrt(3.)*(rp-3.) * std::conj(bi[K31]) * bp[K30];
    b[K21] -= 3.*sqrt(2.)*(rp-1.) * std::conj(bi[K31]) * bp[K30];
    b[K10] += (3.*rp-1.) * std::conj(bi[K31]) * bp[K30];
    b[K53] += sqrt(5.) * std::conj(bi[K31]) * bp[K40];
    b[K42] += sqrt(3.)*(rp-4.) * std::conj(bi[K31]) * bp[K40];
    b[K31] -= sqrt(6.)*(2.*rp-3.) * std::conj(bi[K31]) * bp[K40];
    b[K20] += sqrt(2.)*(3.*rp-2.) * std::conj(bi[K31]) * bp[K40];
    b[K44] += 4. * std::conj(bi[K31]) * bp[K31];
    b[K33] += (3.*rp-10.+3.*ri) * std::conj(bi[K31]) * bp[K31];
    b[K22] -= 6.*(rp-2.+ri) * std::conj(bi[K31]) * bp[K31];
    b[K11] += (3.*rp-10.+3.*ri) * std::conj(bi[K31]) * bp[K31];
    b[K00] += 4. * std::conj(bi[K31]) * bp[K31];

    b[K22] += bi[K22] * bp[K00];
    b[K11] += 2.*rp * bi[K22] * bp[K00];
    b[K00] += rp2 * bi[K22] * bp[K00];
    b[K32] += sqrt(3.) * bi[K22] * bp[K10];
    b[K21] += sqrt(2.)*(2.*rp-1.) * bi[K22] * bp[K10];
    b[K10] += (rp-2.)*rp * bi[K22] * bp[K10];
    b[K42] += sqrt(6.) * bi[K22] * bp[K20];
    b[K31] += 2.*sqrt(3.)*(rp-1.) * bi[K22] * bp[K20];
    b[K20] += (rp2-4.*rp+1.) * bi[K22] * bp[K20];
    b[K52] += sqrt(10.) * bi[K22] * bp[K30];
    b[K41] += 2.*(2.*rp-3.) * bi[K22] * bp[K30];
    b[K30] += (rp2-6.*rp+3.) * bi[K22] * bp[K30];
    b[K62] += sqrt(15.) * bi[K22] * bp[K40];
    b[K51] += 2.*sqrt(5.)*(rp-2.) * bi[K22] * bp[K40];
    b[K40] += (rp2-8.*rp+6) * bi[K22] * bp[K40];
    b[K33] += 3. * bi[K22] * bp[K11];
    b[K22] += (4.*rp-4.+ri) * bi[K22] * bp[K11];
    b[K11] += (rp2-4.*rp+4.) * bi[K22] * bp[K11];
    b[K00] += 3.*rp * bi[K22] * bp[K11];
    b[K43] += 3.*sqrt(2.) * bi[K22] * bp[K21];
    b[K32] += sqrt(6.)*(2.*rp-3.+ri) * bi[K22] * bp[K21];
    b[K21] += (rp2-6.*rp+9.-2.*ri) * bi[K22] * bp[K21];
    b[K10] += 3.*sqrt(2.)*(rp-1.) * bi[K22] * bp[K21];
    b[K53] += sqrt(30.) * bi[K22] * bp[K31];
    b[K42] += sqrt(2.)*(4.*rp-8.+3.*ri) * bi[K22] * bp[K31];
    b[K31] += (rp2-8.*rp+15.-6.*ri) * bi[K22] * bp[K31];
    b[K20] += sqrt(3.)*(3.*rp-6.+ri) * bi[K22] * bp[K31];
    b[K44] += 6. * bi[K22] * bp[K22];
    b[K33] += 6.*(rp-2.+ri) * bi[K22] * bp[K22];
    b[K22] += (rp2-8.*rp+18.-8.*ri+ri2) * bi[K22] * bp[K22];
    b[K11] += 6.*(rp-2.+ri) * bi[K22] * bp[K22];
    b[K00] += 6. * bi[K22] * bp[K22];

    b *= twosqrtpi;
    CtoR(b,rb);
}

#define EPS (stype == 0 ? 1.e-8 : 3.e-6)

#define xdbgout (XDEBUG ? dbgout : 0)

int main(int argc, char **argv) try 
{
#if 1
    std::string dbgfile = "testwl.debug";
    dbgout = new std::ofstream(dbgfile.c_str());
#ifdef _OPENMP
    omp_set_dynamic(0); 
#pragma omp parallel copyin(dbgout, XDEBUG)
    {
        int threadnum = omp_get_thread_num();
        std::stringstream ss;
        ss << threadnum;
        std::string dbgfile2 = dbgfile + "_" + ss.str();
        if (threadnum > 0) dbgout = new std::ofstream(dbgfile2.c_str());
    }
#endif
#else
    dbgout = &std::cout;
#endif

    tmv::WriteWarningsTo(dbgout);

    // First check that the functions work very well in the limit of
    // well sampled (pixScale = 0.1), large aperture (aperture = 15.)
    // data with exact shapelet intensity patterns.
    // Read data will of course fail to work this well, but the answers
    // will nonetheless be as accurate as possible give the data available.

    tmv::SmallMatrix<double,2,2,tmv::RowMajor> D;
#if 1
    D = tmv::ListInit,
      0.08, -0.02,
      -0.015, 0.11;
    double pixScale = sqrt(D.Det());
#else
    double pixScale = 0.1;
    D.SetToIdentity(pixScale);
#endif
    dbg<<"pixscale = "<<pixScale<<std::endl;

#ifdef TEST12
    const int NB = 2;
    double b_vecs[NB][15] = {
        // 00  10  10  20  20  11  30  30  21  21  40  40  31  31  22
        //
        // {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        // {1.0,0.5,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
        // {1.0,0.0,0.0,0.5,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
        // {1.0,0.0,0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
        // {1.0,0.0,0.0,0.0,0.0,0.0,0.5,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
        // {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.3,0.0,0.0,0.0,0.0,0.0},
        // {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.3,0.0,0.0,0.0},
        // {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.3,0.0},
        // {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5},
        
        {10.0,0.5,0.1,0.2,0.3,0.4,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1},
        {2860,-143,286,-151,-202,358,160,-29,-126,253,-25,87,-109,-146,223} 
    };
#endif

#ifdef TEST1
    const int NLONG = 2;
    double blong_vec[NLONG][45] = {
        {10.0,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,
            0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,
            0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4},
        {8120,-122,331,224,-109,412,-142,-100,-399,200,103,-193,89,-407,142,
            -382,451,110,-538,287,-294,377,295,-275,165,-198,163,528,572,-312,
            -361,-189,140,226,155,-175,-372,238,210,-528,309,-241,314,-252,-111}
    };
#endif

#ifdef TEST23
    const int NPSF = 2;
    double bpsf_vecs[NPSF][15] = {
        //{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        //{1,0.1,0.2,0,0,0,0,0,0,0,0,0,0,0,0},
        //{1,0,0,0.1,0.2,0,0,0,0,0,0,0,0,0,0},
        //{1,0,0,0,0,0.2,0,0,0,0,0,0,0,0,0},
        //{1,0,0,0,0,0,0.1,0.2,0,0,0,0,0,0,0},
        //{1,0,0,0,0,0,0,0,0.1,0.2,0,0,0,0,0},
        //{1,0,0,0,0,0,0,0,0,0,0.1,0.2,0,0,0},
        //{1,0,0,0,0,0,0,0,0,0,0,0,0.1,0.2,0},
        //{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2},
        {1.0,0.0,0.0,-0.2,0.1,0.1,-0.2,0.3,0.1,-0.2,0.2,0.1,-0.1,-0.2,-0.1},
        {1.0,0.22,0.31,-0.28,0.12,0.11,-0.18,0.30,0.09,-0.28,0.23,0.14,-0.12,-0.20,-0.09}
    };
    double sigmaPsf = 1.1;
#endif

#ifdef TEST234
    double f_psf = 1.0;
#endif

#ifdef TEST45
    const int NPSFx = 7;
    double bpsfx_vecs[NPSFx][15] = {
        {1.0,0.00,0.00,-0.28,0.12,0.11,-0.18,0.30,0.09,-0.28,0.23,0.14,-0.12,-0.20,-0.09},
        {1.0,0.00,0.00,0.18,-0.29,0.41,0.33,-0.41,0.19,0.47,0.11,-0.18,-0.31,0.15,-0.38},
        {1.0,0.00,0.00,-0.06,0.01,0.00,-0.01,-0.03,0.06,-0.07,-0.07,0.06,0.09,-0.07,0.05},
        {1.0,0.00,0.00,0.19,-0.12,-0.17,0.17,0.12,-0.16,0.18,-0.18,0.14,-0.18,0.14,-0.10},
        {1.0,0.00,0.00,0.31,-0.34,-0.32,-0.35,-0.33,0.36,-0.30,0.34,0.33,-0.34,-0.38,0.38},
        {1.0,0.00,0.00,0.30,-0.43,0.23,0.04,-0.28,-0.37,-0.02,-0.21,0.32,0.21,0.44,-0.26},
        {1.0,0.00,0.00,-0.24,0.28,0.29,-0.29,0.27,0.24,-0.24,0.23,-0.28,0.20,0.22,0.22}
    };
    double xcenx[NPSFx] = { 0.7, -0.4, 0.1, -0.2, -0.3, 0.9, -0.8  };
    double ycenx[NPSFx] = { -0.2, -0.6, 0.4, 0.0, 0.5, -0.1, 0.3 };
    double sigmaPsfx[NPSFx] = { 1.3, 1.0, 1.1, 0.9, 1.1, 0.7, 1.2 };
#endif

#ifdef TEST12345
    double aperture = 15.;
    const int NELL = 2;
    double ell_vecs[NELL][5] = {
        //{0.0,0.0,0.0,0.0,0.0},
        //{0.3,0.0,0.0,0.0,0.0},
        //{0.0,0.3,0.0,0.0,0.0},
        //{0.0,0.0,0.3,0.0,0.0},
        //{0.0,0.0,0.0,0.3,0.0},
        //{0.0,0.0,0.0,0.0,0.3},
        {0.1,-0.2,0.1,0.2,-0.2},
        {-0.11,0.24,-0.18,0.15,-0.14}
    };

    double sigma_i = 1.7;
    Ellipse e0; // all 0's.

    std::auto_ptr<BaseEllipseSolver> s;
    std::vector<PixelList> allpix(1);
#endif

#ifdef TEST1
    // First test the applyZ, G, Mu routines along with basic b vector 
    // measurements.
    for(int stype = 0; stype<=1; ++stype) {
        dbg<<"stype = "<<stype<<std::endl;
        for(int ie = 0; ie < NELL; ++ie) {
            Ellipse ell(ell_vecs[ie]);
            dbg<<"ell = "<<ell<<std::endl;
            //std::complex<double> z0 = ell.getCen();
            std::complex<double> g0 = ell.getGamma();
            double absg0 = std::abs(g0);
            double normg0 = std::norm(g0);
            std::complex<double> argg0 = absg0 == 0. ? 
                std::complex<double>(0.6,0.8) : g0 / absg0;
            double eta0 = absg0 == 0. ? 0. : 2.*atanh(absg0);
            double m0 = ell.getMu();
            for(int ib = 0; ib < NB; ++ib) {
                BVec bin(4,sigma_i,b_vecs[ib]);
                dbg<<"bin = "<<bin.vec()<<std::endl;
                int order = 8;
                double xcen = 0.709;
                double ycen = 0.243;
                tmv::Vector<double> x0(5,ell_vecs[ie]);
                tmv::Vector<double> xe(5,ell_vecs[ie]);
                tmv::Vector<double> f(5);

                allpix.resize(1);
                allpix[0].clear();
                getFakePixList(allpix[0],xcen,ycen,D,aperture,bin,ell);

                if (stype == 0) s.reset(
                    new EllipseSolver(allpix,order,sigma_i));
                else s.reset(
                    new EllipseSolver2(allpix,order,sigma_i,pixScale));
                s->calculateF(x0,f);
                BVec b0 = s->getB();
                dbg<<"measured b = "<<b0.vec()<<std::endl;
                BVec b0save = b0;
                test(Norm(b0.vec().SubVector(0,15)-bin.vec()) < 
                     EPS*Norm(bin.vec()),
                     "Basic B Measurement");
                test(Norm(b0.vec().SubVector(15,b0.size())) < 
                     EPS*Norm(bin.vec()),
                     "Basic B Measurement 2");
#ifdef TESTJ
                test(s->testJ(x0,f,xdbgout,1.e-6),"TestJ - Basic");
#endif

                std::complex<double> z1(0.2,0.3);
                tmv::Vector<double> x1 = x0;
                //x1[0] += std::real(z1);
                //x1[1] += std::imag(z1);
                // The above is right if ell = (0,0,0,0,0)
                // But if not, we need the following adjustment to account 
                // for the fact that the shift happens after having already 
                // been transformed:
                x1[0] += std::real(
                    z1 + g0*std::conj(z1))*exp(m0)/sqrt(1.-normg0);
                x1[1] += std::imag(
                    z1 + g0*std::conj(z1))*exp(m0)/sqrt(1.-normg0);
                s->calculateF(x1,f);
                BVec b1 = s->getB();
                xdbg<<"b(z1 = "<<z1<<") = "<<b1.vec()<<std::endl;
                applyZ(z1,b0);
                xdbg<<"b(predicted) = "<<b0.vec()<<std::endl;
                dbg<<"diff for z = "<<z1<<" = "<<
                    Norm(b0.vec()-b1.vec())<<std::endl;
                test(Norm(b1.vec()-b0.vec()) < EPS*Norm(b0.vec()),"ApplyZ 1");
                x1 = x0;
                x1[0] -= std::real(
                    z1 + g0*std::conj(z1))*exp(m0)/sqrt(1.-normg0);
                x1[1] -= std::imag(
                    z1 + g0*std::conj(z1))*exp(m0)/sqrt(1.-normg0);
                s->calculateF(x1,f);
                b1 = s->getB();
                applyZ(-z1,b0=b0save);
                dbg<<"diff for z = "<<-z1<<" = "<<
                    Norm(b0.vec()-b1.vec())<<std::endl;
                test(Norm(b1.vec()-b0.vec()) < EPS*Norm(b0.vec()),"ApplyZ 2");

                double m1(0.2);
                x1 = x0;
                x1[4] += m1;
                s->calculateF(x1,f);
                b1 = s->getB();
                xdbg<<"b(mu = "<<m1<<") = "<<b1.vec()<<std::endl;
                applyMu(m1,b0=b0save);
                xdbg<<"b(predicted) = "<<b0.vec()<<std::endl;
                xdbg<<"diff = "<<b1.vec()-b0.vec()<<std::endl;
                dbg<<"diff for mu = "<<m1<<" = "<<
                    Norm(b0.vec()-b1.vec())<<std::endl;
                dbg<<"Norm(b0) = "<<Norm(b0.vec())<<std::endl;
#if 1
                test(Norm(b1.vec()-b0.vec()) < 100.*EPS*Norm(b0.vec()),
                     "ApplyMu 1");
#else
                test(Norm(b1.vec().SubVector(0,15)-b0.vec().SubVector(0,15)) <
                     EPS*Norm(b0.vec()),
                     "ApplyMu 1");
                test(Norm(b1.vec().SubVector(15,b1.size()) -
                          b0.vec().SubVector(15,b1.size())) < 
                     10.*b1.size()*EPS*Norm(b0.vec()),
                     "ApplyMu 2");
#endif
                x1 = x0;
                x1[4] += -m1;
                s->calculateF(x1,f);
                b1 = s->getB();
                applyMu(-m1,b0=b0save);
                dbg<<"diff for mu = "<<-m1<<" = "<<
                    Norm(b0.vec()-b1.vec())<<std::endl;
#if 1
                test(Norm(b1.vec()-b0.vec()) < 100.*EPS*Norm(b0.vec()),
                     "ApplyMu 2");
#else
                test(Norm(b1.vec().SubVector(0,15)-b0.vec().SubVector(0,15)) < 
                     EPS*Norm(b0.vec()),
                     "ApplyMu 3");
                test(Norm(b1.vec().SubVector(15,b1.size()) -
                          b0.vec().SubVector(15,b1.size())) < 
                     10.*b1.size()*EPS*Norm(b0.vec()),
                     "ApplyMu 4");
#endif

                double e1(0.3);
                double neweta = e1 + eta0;
                double g1 = tanh(e1/2);
                x1 = x0;
                // Need to apply eta in the same direction as old eta.
                // Otherwise the equations get really hairy.
                x1[2] = tanh(neweta/2.)*std::real(argg0);
                x1[3] = tanh(neweta/2.)*std::imag(argg0);
                s->calculateF(x1,f);
                b1 = s->getB();
                xdbg<<"b(gamma = "<<g1<<") = "<<b1.vec()<<std::endl;
                applyG(g1*argg0,b0=b0save);
                xdbg<<"b(predicted) = "<<b0.vec()<<std::endl;
                dbg<<"diff for gamma = "<<g1<<" = "<<
                    Norm(b0.vec()-b1.vec())<<std::endl;
#if 1
                test(Norm(b1.vec()-b0.vec()) < 10.*EPS*Norm(b0.vec()),
                     "ApplyG 1");
#else
                test(Norm(b1.vec().SubVector(0,15)-b0.vec().SubVector(0,15)) < 
                     EPS*Norm(b0.vec()),
                     "ApplyG 1");
                test(Norm(b1.vec().SubVector(15,b1.size()) -
                          b0.vec().SubVector(15,b1.size())) < 
                     b1.size()*EPS*Norm(b0.vec()),
                     "ApplyG 2");
#endif
                x1 = x0;
                neweta = -e1 + eta0;
                x1[2] = tanh(neweta/2.)*std::real(argg0);
                x1[3] = tanh(neweta/2.)*std::imag(argg0);
                s->calculateF(x1,f);
                b1 = s->getB();
                applyG(-g1*argg0,b0=b0save);
                dbg<<"diff for gamma = "<<-g1<<" = "<<
                    Norm(b0.vec()-b1.vec())<<std::endl;
#if 1
                test(Norm(b1.vec()-b0.vec()) < 10.*EPS*Norm(b0.vec()),
                     "ApplyG 2");
#else
                test(Norm(b1.vec().SubVector(0,15)-b0.vec().SubVector(0,15)) < 
                     EPS*Norm(b0.vec()),
                     "ApplyG 3");
                test(Norm(b1.vec().SubVector(15,b1.size()) -
                          b0.vec().SubVector(15,b1.size())) < 
                     b1.size()*EPS*Norm(b0.vec()),
                     "ApplyG 4");
#endif

                if (ib < NLONG) {
                    BVec blong(8,sigma_i,blong_vec[ib]);

                    xcen = 0.709;
                    ycen = 0.243;
                    allpix.resize(1);
                    allpix[0].clear();
                    getFakePixList(allpix[0],xcen,ycen,D,aperture,blong,ell);

                    order = 10;
                    if (stype == 0) s.reset(
                        new EllipseSolver(allpix,order,sigma_i));
                    else s.reset(
                        new EllipseSolver2(allpix,order,sigma_i,pixScale));
                    tmv::Vector<double> x2(5,ell_vecs[ie]);
                    s->calculateF(x2,f);
                    BVec b = s->getB();
                    dbg<<"blong = "<<blong.vec()<<std::endl;
                    dbg<<"b = "<<b.vec()<<std::endl;
                    test(Norm(b.vec().SubVector(0,45)-blong.vec()) < 
                         EPS*Norm(blong.vec()),
                         "Basic Long B Measurement");
                    test(Norm(b.vec().SubVector(45,b.size())) < 
                         EPS*Norm(blong.vec()),
                         "Basic Long B Measurement (2)");
#ifdef TESTJ
                    test(s->testJ(x2,f,xdbgout,1.e-6),"TestJ - Long 1");
#endif

                    double dmu = 0.01;
                    double sigma_ix = exp(dmu) * sigma_i;
                    if (stype == 0) s.reset(
                        new EllipseSolver(allpix,order,sigma_ix));
                    else s.reset(
                        new EllipseSolver2(allpix,order,sigma_ix,pixScale));
                    x2[4] -= dmu;
                    s->calculateF(x2,f);
                    BVec bx = s->getB();
                    dbg<<"blong = "<<blong.vec()<<std::endl;
                    dbg<<"bx = "<<bx.vec()<<std::endl;
                    test(Norm(bx.vec().SubVector(0,45)-blong.vec()) < 
                         EPS*Norm(blong.vec()),
                         "Basic Long B Measurement");
                    test(Norm(bx.vec().SubVector(45,bx.size())) < 
                         EPS*Norm(blong.vec()),
                         "Basic Long B Measurement (2)");
#ifdef TESTJ
                    test(s->testJ(x2,f,xdbgout,1.e-6),"TestJ - Long 2");
#endif
                }
            }
        }
    }
    std::cout<<"Passed tests of ApplyZ, G, Mu and basic measurements\n";
#endif

#ifdef TEST2
    // Next test calculatePSFConvolve
    for(int stype = 0; stype<=1; ++stype) {
        dbg<<"Start stype = "<<stype<<std::endl;
        for(int ie = 0; ie < NELL; ++ie) {
            dbg<<"Start ie = "<<ie<<std::endl;
            Ellipse ell(ell_vecs[ie]);
            dbg<<"ell = "<<ell<<std::endl;
            std::complex<double> z0 = ell.getCen();
            std::complex<double> g0 = ell.getGamma();
            double m0 = ell.getMu();
            for(int ib = 0; ib < NB; ++ib) {
                dbg<<"Start ib = "<<ib<<std::endl;
                BVec bin(4,sigma_i,b_vecs[ib]);
                dbg<<"bin = "<<bin.vec()<<std::endl;
                int order = 8;
                double xcen = 0.709;
                double ycen = 0.243;
                tmv::Vector<double> xe(5,ell_vecs[ie]);
                tmv::Vector<double> f(5);

                for(int ip = 0; ip < NPSF; ++ip) {
                    dbg<<"Start ip = "<<ip<<std::endl;
                    // b0 is the galaxy shape in the e0 frame.
                    BVec b0 = bin;

                    // bpsf is the psf in the e0 frame.
                    BVec bpsf(4,sigmaPsf,bpsf_vecs[ip]);
                    dbg<<"bpsf = "<<bpsf.vec()<<std::endl;
                    xdbg<<"ell = "<<ell<<std::endl;
                    xdbg<<"b0 = "<<b0.vec()<<std::endl;

                    // c0 = C(bpsf) * b0
                    // This is one estimate of the observed galaxy in the e0
                    // frame.
                    double sigma_o = sqrt(pow(sigma_i,2)+pow(sigmaPsf,2));
                    tmv::Matrix<double> C(45,45);
                    calculatePsfConvolve(bpsf,8,sigma_i,C);
                    BVec c0(8,sigma_o);
                    c0.vec() = C.Cols(0,15) * b0.vec();
                    xdbg<<"c0 = "<<c0.vec()<<std::endl;

                    // c0x is another estimate of the same thing based on a 
                    // direct convolution from Maple calculations rather than
                    // the BJ02 recursion formulae.
                    BVec c0x(8,sigma_o);
                    DirectConvolveB(b0,bpsf,c0x);
                    test(Norm(c0.vec()-c0x.vec()) < EPS*Norm(c0.vec()),
                         "PSFConvolve");

                    // be is the galaxy in ell fram:
                    allpix.resize(1);
                    allpix[0].clear();
                    getFakePixList(allpix[0],xcen,ycen,D,aperture,b0,e0);
                    order = 12;
                    if (stype == 0) s.reset(
                        new EllipseSolver(allpix,order,sigma_i));
                    else s.reset(
                        new EllipseSolver2(allpix,order,sigma_i,pixScale));
                    s->calculateF(xe,f);
                    BVec be = s->getB();
                    xdbg<<"be = "<<be.vec()<<std::endl;
#ifdef TESTJ
                    test(s->testJ(xe,f,xdbgout,1.e-6),
                         "TestJ - be in ell frame");
#endif

                    // bex is another estimate of the galaxy in the ell frame:
                    BVec bex(order,sigma_i);
                    xdbg<<"b0 = "<<b0.vec()<<std::endl;
                    xdbg<<"bex = "<<bex.vec()<<std::endl;
                    bex = b0;
                    xdbg<<"after bex = b0: "<<bex.vec()<<std::endl;
                    applyZ(z0,bex);
                    xdbg<<"after z: bex = "<<bex.vec()<<std::endl;
                    applyG(g0,bex);
                    xdbg<<"after g: bex = "<<bex.vec()<<std::endl;
                    applyMu(m0,bex);
                    xdbg<<"after mu: bex = "<<bex.vec()<<std::endl;
                    // The higher order terms don't mesh perfectly,
                    // because of the finite order of the intermediate 
                    // transformations.
                    // But through 4th order, the match should be quite good.
                    test(Norm(bex.vec().SubVector(0,15) - 
                              be.vec().SubVector(0,15)) <
                         3.e-5*(Norm(be.vec())+Norm(b0.vec())),
                         "View B in ell frame");
                    test(Norm(bex.vec().SubVector(0,45) - 
                              be.vec().SubVector(0,45)) <
                         1.e-3*(Norm(be.vec())+Norm(b0.vec())),
                         "View B in ell frame");
                    test(Norm(bex.vec()-be.vec()) <
                         3.e-2*(Norm(be.vec())+Norm(b0.vec())),
                         "View B in ell frame");

                    // ce is the convolved galaxy in the ell frame:
                    allpix.resize(1);
                    allpix[0].clear();
                    getFakePixList(allpix[0],xcen,ycen,D,aperture,c0,e0);
                    if (stype == 0) s.reset(
                        new EllipseSolver(allpix,order,sigma_o));
                    else s.reset(
                        new EllipseSolver2(allpix,order,sigma_o,pixScale));
                    s->calculateF(xe,f);
                    BVec ce = s->getB();
#ifdef TESTJ
                    test(s->testJ(xe,f,xdbgout,1.e-5),
                         "TestJ - ce in ell frame");
#endif
                    // cex is another estimate of the same: transform c0
                    BVec cex(order,sigma_o);
                    cex = c0;
                    applyZ(z0,cex);
                    applyG(g0,cex);
                    applyMu(m0,cex);
                    test(Norm(cex.vec().SubVector(0,15) - 
                              ce.vec().SubVector(0,15)) <
                         3.e-5*(Norm(ce.vec())+Norm(b0.vec())),
                         "Convolved B in ell frame");
                    test(Norm(cex.vec().SubVector(0,45) - 
                              ce.vec().SubVector(0,45)) <
                         1.e-3*(Norm(ce.vec())+Norm(b0.vec())),
                         "Convolved B in ell frame");
                    test(Norm(cex.vec()-ce.vec()) <
                         3.e-2*(Norm(ce.vec())+Norm(c0.vec())),
                         "Convolved B in ell frame");

                    // cey is yet another estimate of the same: convolve be
                    BVec cey(order,sigma_o);
                    tmv::Matrix<double> Ce(cey.size(),cey.size());
                    BVec bpsfe(order,sigmaPsf);
                    bpsfe = bpsf;
                    //applyZ(z0,bpsfe);
                    applyG(g0,bpsfe);
                    applyMu(m0,bpsfe);
                    calculatePsfConvolve(bpsfe,order,sigma_i,Ce);
                    cey.vec() = Ce * be.vec();
                    cey.vec() *= exp(2.*m0);
                    test(Norm(cey.vec()-ce.vec()) <
                         3.e-2*(Norm(ce.vec())+Norm(c0.vec())),
                         "Convolved B in ell frame");

                    // Next we test the deconovolving shape estimator.
                    // Since the cey calculation doesn't really match ce 
                    // very well due to numerical issues from the finite 
                    // order of the calculations, we recast the problem so 
                    // b0 is the shape in the ell frame.
                    BVec bpsf2(order,sigmaPsf);
                    bpsf2 = bpsf;
#if 1
                    tmv::Matrix<double> S1(bpsf2.size(),bpsf2.size());
                    calculateGTransform(g0,bpsf2.getOrder(),S1);
                    tmv::Matrix<double> D1(bpsf2.size(),bpsf2.size());
                    calculateMuTransform(m0,bpsf2.getOrder(),D1);
                    bpsf2.vec() /= D1;
                    bpsf2.vec() /= S1;
#else
                    // With this version, the numerical errors mean that
                    // the EPS in the two Test statements below need to 
                    // be changed to around 1.e-2.
                    tmv::Matrix<double> S1(bpsf2.size(),bpsf2.size());
                    calculateGTransform(-g0,bpsf2.getOrder(),S1);
                    tmv::Matrix<double> D1(bpsf2.size(),bpsf2.size());
                    calculateMuTransform(-m0,bpsf2.getOrder(),D1);
                    bpsf2 = D1*bpsf2;
                    bpsf2 = S1*bpsf2;
#endif
                    allpix.resize(1);
                    allpix[0].clear();
                    std::vector<BVec> allpsf(1,bpsf2);
                    getFakePixList(allpix[0],xcen,ycen,D,aperture,c0,ell);

                    if (stype == 0) s.reset(
                        new EllipseSolver(allpix,allpsf,f_psf,order,sigma_i));
                    else s.reset(
                        new EllipseSolver2(allpix,allpsf,f_psf,order,sigma_i,
                                           pixScale));
                    s->calculateF(xe,f);
                    BVec bes = s->getB();
                    xdbg<<"bes = "<<bes.vec()<<std::endl;
                    xdbg<<"b0 = "<<b0.vec()<<std::endl;
                    xdbg<<"Norm(bes-b0) = "<<
                        Norm(bes.vec().SubVector(0,15)-b0.vec())<<std::endl;
                    xdbg<<"Norm(bes(15:)) = "<<
                        Norm(bes.vec().SubVector(15,bes.size()))<<std::endl;
                    xdbg<<"cf "<<EPS<<" * "<<Norm(b0.vec())+Norm(c0.vec())<<std::endl;
                    test(Norm(bes.vec().SubVector(0,15)-b0.vec()) < 
                         EPS*(Norm(b0.vec())+Norm(c0.vec())),
                         "Deconvolving solver BE Measurement");
                    test(Norm(bes.vec().SubVector(15,bes.size())) <
                         3.*EPS*(Norm(b0.vec())+Norm(c0.vec())),
                         "Deconvolving solver BE Measurement");
#ifdef TESTJ
                    test(s->testJ(xe,f,xdbgout,1.e-6),
                         "TestJ - Deconvolved b0 in ell frame");
#endif

#if 0
                    // Next do it with a different sigma for the galaxy
                    double dmu = 0.002;
                    double sigma_ix = exp(dmu) * sigma_i;
                    double sigma_ox = 
                        sqrt(sigma_ix*sigma_ix + sigmaPsf*sigmaPsf);
                    double dmu2 = dmu + 0.0*dmu*dmu - 0.0*dmu*dmu*dmu;
                    xdbg<<"sigma_i = "<<sigma_i<<std::endl;
                    xdbg<<"sigma_o = "<<sigma_o<<std::endl;
                    xdbg<<"dmu = "<<dmu<<std::endl;
                    xdbg<<"sigma_ix = "<<sigma_ix<<std::endl;
                    xdbg<<"sigma_ox = "<<sigma_ox<<std::endl;
                    xdbg<<"dmu2 = "<<dmu2<<std::endl;
                    if (stype == 0) s.reset(
                        new EllipseSolver(allpix,allpsf,f_psf,order,sigma_ix));
                    else s.reset(
                        new EllipseSolver2(allpix,allpsf,f_psf,order,sigma_ix,
                                           pixScale));
                    xe[4] -= dmu2;
                    s->calculateF(xe,f);
                    BVec besx = s->getB();
                    xdbg<<"besx = "<<besx.vec()<<std::endl;
                    xdbg<<"b0 = "<<b0.vec()<<std::endl;
                    besx.vec() *= b0(0)/besx(0);
                    xdbg<<"besx => "<<besx.vec()<<std::endl;
                    xdbg<<"Norm(besx-b0) = "<<
                        Norm(besx.vec().SubVector(0,15)-b0.vec())<<std::endl;
                    xdbg<<"Norm(bes(15:)) = "<<
                        Norm(besx.vec().SubVector(15,besx.size()))<<std::endl;
                    test(Norm(besx.vec().SubVector(0,15)-b0.vec()) < 
                         EPS*(Norm(b0.vec())+Norm(c0.vec())),
                         "Deconvolving solver BE Measurement #2");
                    test(Norm(besx.vec().SubVector(15,besx.size())) <
                         3.*EPS*(Norm(b0.vec())+Norm(c0.vec())),
                         "Deconvolving solver BE Measurement #2");
#ifdef TESTJ
                    test(s->testJ(xe,f,xdbgout,1.e-6),
                         "TestJ - Deconvolved b0 in ell frame #2");
#endif
                    xe[4] += dmu2;
#endif
                }
            }
        }
    }
    std::cout<<"Passed tests of PSFConvolve and deconvolving measurements.\n";
#endif

#ifdef TEST3
    // Next test Solve for the coordinates in which a Gaussian galaxy is round
    for(int stype = 0; stype<=1; ++stype) {
        dbg<<"Start stype = "<<stype<<std::endl;
        for(int ie = 0; ie < NELL; ++ie) {
            dbg<<"Start ie = "<<ie<<std::endl;
            Ellipse ell(ell_vecs[ie]);
            dbg<<"ell = "<<ell<<std::endl;
            std::complex<double> g0 = ell.getGamma();
            double m0 = ell.getMu();

            BVec bin(4,sigma_i);
            bin.vec().Zero();
            bin(0) = 1.;
            dbg<<"bin = "<<bin.vec()<<std::endl;
            int order = 8;
            double xcen = 0.709;
            double ycen = 0.243;
            tmv::Vector<double> xe(5,ell_vecs[ie]);

            // b0 is the galaxy shape in the ell frame.
            BVec b0 = bin;
            allpix.resize(1);
            allpix[0].clear();
            getFakePixList(allpix[0],xcen,ycen,D,aperture,b0,ell);
            if (stype == 0) s.reset(
                new EllipseSolver(allpix,order,sigma_i));
            else s.reset(
                new EllipseSolver2(allpix,order,sigma_i,pixScale));
            tmv::Vector<double> x(5,0.);
            tmv::Vector<double> f(5);
            if (XDEBUG) s->setOutput(*dbgout);
            s->useDogleg();
            s->setDelta0(0.1);
            s->setFTol(1.e-6);
#ifdef TESTJ
            test(s->testJ(x,f,xdbgout,1.e-6),
                 "TestJ - Single Image no PSF");
#endif
            bool success = s->solve(x,f);
            test(success,"No PSF Solve - success");
            test(Norm(x-xe) < 1.e-5,"No PSF Solve - x == xe");

            for(int ip = 0; ip < NPSF; ++ip) {
                dbg<<"Start ip = "<<ip<<std::endl;

                // bpsf is the psf in the ell frame.
                BVec bpsf(4,sigmaPsf,bpsf_vecs[ip]);
                dbg<<"bpsf = "<<bpsf.vec()<<std::endl;
                xdbg<<"ell = "<<ell<<std::endl;
                xdbg<<"b0 = "<<b0.vec()<<std::endl;

                // bpsf2 is the psf in the e0 frame.
                BVec bpsf2(order+4,sigmaPsf);
                bpsf2 = bpsf;
                tmv::Matrix<double> S1(bpsf2.size(),bpsf2.size());
                calculateGTransform(g0,bpsf2.getOrder(),S1);
                tmv::Matrix<double> D1(bpsf2.size(),bpsf2.size());
                calculateMuTransform(m0,bpsf2.getOrder(),D1);
                bpsf2.vec() /= D1;
                bpsf2.vec() /= S1;

                // c0 is the convolved galaxy shape in the ell frame.
                double sigma_o = sqrt(pow(sigma_i,2)+pow(sigmaPsf,2));
                tmv::Matrix<double> C(45,45);
                calculatePsfConvolve(bpsf,8,sigma_i,C);
                BVec c0(8,sigma_o);
                c0.vec() = C.Cols(0,15) * b0.vec();
                xdbg<<"c0 = "<<c0.vec()<<std::endl;
                allpix.resize(1);
                allpix[0].clear();
                std::vector<BVec> allpsf(1,bpsf2);
                getFakePixList(allpix[0],xcen,ycen,D,aperture,c0,ell);

                if (stype == 0) s.reset(
                    new EllipseSolver(allpix,allpsf,f_psf,order,sigma_i));
                else s.reset(
                    new EllipseSolver2(allpix,allpsf,f_psf,order,sigma_i,
                                       pixScale));
                s->calculateF(xe,f);
                BVec be = s->getB();
                xdbg<<"be = "<<be.vec()<<std::endl;
                test(Norm(be.vec().SubVector(0,15)-bin.vec()) < 
                     EPS*Norm(bin.vec()),
                     "B Measurement for Single Image w/ PSF");
                test(Norm(be.vec().SubVector(15,be.size())) < 
                     3.*EPS*Norm(bin.vec()),
                     "B Measurement 2 for Single Image w/ PSF");
                x.Zero();
                if (XDEBUG) s->setOutput(*dbgout);
                s->useDogleg();
                s->setDelta0(0.1);
                s->setFTol(1.e-6);
#ifdef TESTJ
                test(s->testJ(x,f,xdbgout,1.e-6),
                     "TestJ - Single Image w/ PSF");
#endif
                success = s->solve(x,f);
#ifdef TESTJ
                test(s->testJ(x,f,xdbgout,1.e-6),
                     "TestJ - Single Image w/ PSF");
#endif
                test(success,"PSF Solve - success");
                test(Norm(x-xe) < 1.e-5,"PSF Solve - x == xe");

#if 0
                double dmu = 0.1;
                double sigma_ix = exp(dmu) * sigma_i;
                if (stype == 0) s.reset(
                    new EllipseSolver(allpix,allpsf,f_psf,order,sigma_ix));
                else s.reset(
                    new EllipseSolver2(allpix,allpsf,f_psf,order,sigma_ix,
                                       pixScale));
                xe[4] -= dmu;
                s->calculateF(xe,f);
                BVec bex = s->getB();
                xdbg<<"bex = "<<bex.vec()<<std::endl;
                test(Norm(bex.vec().SubVector(0,15)-bin.vec()) < 
                     EPS*Norm(bin.vec()),
                     "B Measurement for Single Image w/ PSF");
                test(Norm(bex.vec().SubVector(15,bex.size())) < 
                     3.*EPS*Norm(bin.vec()),
                     "B Measurement 2 for Single Image w/ PSF");
                x.Zero();
                if (XDEBUG) s->setOutput(*dbgout);
                s->useDogleg();
                s->setDelta0(0.1);
                s->setFTol(1.e-6);
                success = s->solve(x,f);
                test(success,"PSF Solve - success");
                test(Norm(x-xe) < 1.e-5,"PSF Solve - x == xe");
                xe[4] += dmu;
#endif
            }
        }
    }
    std::cout<<"Passed tests of solve for Ellipse.\n";
#endif

#ifdef TEST4
    // Next test Solve with multiple PSFs
    for(int stype = 0; stype<=1; ++stype) {
        dbg<<"Start stype = "<<stype<<std::endl;
        for(int ie = 0; ie < NELL; ++ie) {
            dbg<<"Start ie = "<<ie<<std::endl;
            Ellipse ell(ell_vecs[ie]);
            dbg<<"ell = "<<ell<<std::endl;
            std::complex<double> g0 = ell.getGamma();
            double m0 = ell.getMu();

            BVec bin(4,sigma_i);
            bin.vec().Zero();
            bin(0) = 1.;
            dbg<<"bin = "<<bin.vec()<<std::endl;
            int order = 8;
            tmv::Vector<double> xe(5,ell_vecs[ie]);

            // b0 is the galaxy shape in the ell frame.
            BVec b0 = bin;
            allpix.resize(NPSFx);

#if 1
            {
                // First do the test without the PSFs
                for(int k=0;k<NPSFx;++k) {
                    allpix[k].clear();
                    getFakePixList(
                        allpix[k],xcenx[k],ycenx[k],D,aperture,b0,ell);
                }
                if (stype == 0) s.reset(
                    new EllipseSolver(allpix,order,sigma_i));
                else s.reset(
                    new EllipseSolver2(allpix,order,sigma_i,pixScale));
                tmv::Vector<double> x(5,0.);
                tmv::Vector<double> f(5);
                s->calculateF(xe,f);
                BVec be = s->getB();
                xdbg<<"be = "<<be.vec()<<std::endl;
                xdbg<<"Norm(be(0:15)-bin) = "<<
                    Norm(be.vec().SubVector(0,15)-bin.vec())<<std::endl;
                xdbg<<"Norm(be(15:)) = "<<
                    Norm(be.vec().SubVector(15,be.size()))<<std::endl;
                test(Norm(be.vec().SubVector(0,15)-bin.vec()) < 
                     EPS*Norm(bin.vec()),
                     "B Measurement for Multi Image");
                test(Norm(be.vec().SubVector(15,be.size())) < 
                     3.*EPS*Norm(bin.vec()),
                     "B Measurement 2 for Multi Image");

                if (XDEBUG) s->setOutput(*dbgout);
                s->useDogleg();
                s->setDelta0(0.1);
                s->setFTol(1.e-6);
#ifdef TESTJ
                test(s->testJ(x,f,xdbgout,1.e-6),
                     "TestJ - Multi Image no PSF 1");
#endif
                bool success = s->solve(x,f);
                xdbg<<"x = "<<x<<std::endl;
                xdbg<<"xe = "<<xe<<std::endl;
                xdbg<<"f = "<<f<<std::endl;
                xdbg<<"success = "<<success<<std::endl;
                test(success,"No PSF Solve - success");
                test(Norm(x-xe) < 1.e-5,"No PSF Solve - x == xe");
#ifdef TESTJ
                test(s->testJ(x,f,xdbgout,1.e-6),
                     "TestJ - Multi Image no PSF 2");
#endif
            }
#endif


            {
                // Next do the test with the PSFs
                std::vector<BVec> allpsf;
                std::vector<BVec> allpsf2;

                for(int k=0;k<NPSFx;++k) {
                    // allpsf are the psfs in the ell frame.
                    allpsf.push_back(BVec(4,sigmaPsfx[k],bpsfx_vecs[k]));
                    dbg<<"allpsf["<<k<<"] = "<<allpsf[k].vec()<<std::endl;

                    // allpsf2 is the psf in the e0 frame.
                    allpsf2.push_back(BVec(order+4,sigmaPsfx[k]));
                    allpsf2[k] = allpsf[k];
                    tmv::Matrix<double> S1(allpsf2[k].size(),allpsf2[k].size());
                    calculateGTransform(g0,allpsf2[k].getOrder(),S1);
                    tmv::Matrix<double> D1(allpsf2[k].size(),allpsf2[k].size());
                    calculateMuTransform(m0,allpsf2[k].getOrder(),D1);
                    allpsf2[k].vec() /= D1;
                    allpsf2[k].vec() /= S1;
                    xdbg<<"allpsf2[k] = "<<allpsf2[k].vec()<<std::endl;

                    // c0 is the convolved galaxy shape in the ell frame.
                    double sigma_o = sqrt(pow(sigma_i,2)+pow(sigmaPsfx[k],2));
                    tmv::Matrix<double> C(45,45);
                    calculatePsfConvolve(allpsf[k],8,sigma_i,C);
                    BVec c0(8,sigma_o);
                    c0.vec() = C.Cols(0,15) * b0.vec();
                    xdbg<<"c0 = "<<c0.vec()<<std::endl;

                    allpix[k].clear();
                    getFakePixList(
                        allpix[k],xcenx[k],ycenx[k],D,aperture,c0,ell);
                }

                if (stype == 0) s.reset(
                    new EllipseSolver(allpix,allpsf2,f_psf,order,sigma_i));
                else s.reset(
                    new EllipseSolver2(allpix,allpsf2,f_psf,order,sigma_i,
                                       pixScale));

                tmv::Vector<double> x(5,0.);
                tmv::Vector<double> f(5);
                s->calculateF(xe,f);
                BVec be = s->getB();
                xdbg<<"be = "<<be.vec()<<std::endl;
                xdbg<<"Norm(be(0:15)-bin) = "<<
                    Norm(be.vec().SubVector(0,15)-bin.vec())<<std::endl;
                xdbg<<"Norm(be(15:)) = "<<
                    Norm(be.vec().SubVector(15,be.size()))<<std::endl;
                test(Norm(be.vec().SubVector(0,15)-bin.vec()) < 
                     EPS*Norm(bin.vec()),
                     "B Measurement for Multi Image w/ PSF");
                test(Norm(be.vec().SubVector(15,be.size())) < 
                     3.*EPS*Norm(bin.vec()),
                     "B Measurement 2 for Multi Image w/ PSF");

                x.Zero();
                if (XDEBUG) s->setOutput(*dbgout);
                s->useDogleg();
                s->setDelta0(0.1);
                s->setFTol(1.e-6);
#ifdef TESTJ
                test(s->testJ(x,f,xdbgout,1.e-6),
                     "TestJ - Multi Image w/ PSF 1");
#endif
                bool success = s->solve(x,f);
                xdbg<<"x = "<<x<<std::endl;
                xdbg<<"xe = "<<xe<<std::endl;
                xdbg<<"f = "<<f<<std::endl;
                xdbg<<"success = "<<success<<std::endl;
#ifdef TESTJ
                test(s->testJ(x,f,xdbgout,1.e-6),
                     "TestJ - Multi Image w/ PSF 2");
#endif
                test(success,"PSF Solve - success");
                test(Norm(x-xe) < 1.e-5,"PSF Solve - x == xe");

#if 0
                double dmu = 0.1;
                double sigma_ix = exp(dmu) * sigma_i;
                if (stype == 0) s.reset(
                    new EllipseSolver(allpix,allpsf2,f_psf,order,sigma_ix));
                else s.reset(
                    new EllipseSolver2(allpix,allpsf2,f_psf,order,sigma_ix,
                                       pixScale));

                xe[4] -= dmu;
                s->calculateF(xe,f);
                BVec bex = s->getB();
                xdbg<<"bex = "<<bex.vec()<<std::endl;
                xdbg<<"Norm(bex(0:15)-bin.vec()) = "<<
                    Norm(bex.vec().SubVector(0,15)-bin.vec())<<std::endl;
                xdbg<<"Norm(bex(15:)) = "<<
                    Norm(bex.vec().SubVector(15,bex.size()))<<std::endl;
                test(Norm(bex.vec().SubVector(0,15)-bin.vec()) < 
                     EPS*Norm(bin.vec()),
                     "B Measurement for Multi Image w/ PSF");
                test(Norm(bex.vec().SubVector(15,bex.size())) < 
                     3.*EPS*Norm(bin.vec()),
                     "B Measurement 2 for Multi Image w/ PSF");

                x.Zero();
                if (XDEBUG) s->setOutput(*dbgout);
                s->useDogleg();
                s->setDelta0(0.1);
                s->setFTol(1.e-6);
#ifdef TESTJ
                test(s->testJ(x,f,xdbgout,1.e-6),
                     "TestJ - Multi Image w/ PSF 1");
#endif
                success = s->solve(x,f);
                xdbg<<"x = "<<x<<std::endl;
                xdbg<<"xe = "<<xe<<std::endl;
                xdbg<<"f = "<<f<<std::endl;
                xdbg<<"success = "<<success<<std::endl;
                test(success,"PSF Solve - success");
                test(Norm(x-xe) < 1.e-5,"PSF Solve - x == xe");
#ifdef TESTJ
                test(s->testJ(x,f,xdbgout,1.e-6),
                     "TestJ - Multi Image w/ PSF 2");
#endif
                xe[4] += dmu;
#endif
            }
        }
    }
    std::cout<<"Passed tests of solve for Ellipse with multiple images.\n";
#endif

#ifdef TEST5
    // Next test direct Ellipse class Measure method
    for(int ie = 0; ie < NELL; ++ie) {
        dbg<<"Start ie = "<<ie<<std::endl;
        Ellipse ell(ell_vecs[ie]);
        dbg<<"ell = "<<ell<<std::endl;
        std::complex<double> g0 = ell.getGamma();
        double m0 = ell.getMu();

        // First do the test without the PSFs
        BVec bin(4,sigma_i);
        bin.vec().Zero();
        bin(0) = 1.;
        xdbg<<"bin = "<<bin.vec()<<std::endl;
        int order = 6;
        tmv::Vector<double> xe(5,ell_vecs[ie]);

        // b0 is the galaxy shape in the ell frame.
        BVec b0 = bin;
        allpix.resize(NPSFx,0);

        // Test CrudeMeasure:
        // CrudeMeasure only works precisely for gamma = 0 shapes.
        Ellipse ell_nogamma = ell;
        ell_nogamma.setGamma(0.);
        for(int k=0;k<NPSFx;++k) {
            allpix[k].clear();
            getFakePixList(allpix[k],xcenx[k],ycenx[k],D,aperture,b0,
                           ell_nogamma);
        }
        Ellipse e_0;
        e_0.crudeMeasure(allpix[0],sigma_i);
        xdbg<<"Done e_0 Measure\n";
        xdbg<<"e_0 = "<<e_0<<std::endl;
        test(std::abs(e_0.getCen()-ell.getCen()) < 1.e-3, "CrudeMeasure - cen");
        test(std::abs(e_0.getGamma()) < 1.e-3, "CrudeMeasure - gam");
        test(std::abs(e_0.getMu()-ell.getMu()) < 1.e-3, "CrudeMeasure - mu");

        if (ell.getGamma() != 0.) {
            for(int k=0;k<NPSFx;++k) {
                allpix[k].clear();
                getFakePixList(allpix[k],xcenx[k],ycenx[k],D,aperture,b0,ell);
            }
        }

        // Now the regular Measure
        Ellipse e_1;
        BVec b_1(order,sigma_i);
        long flag;
        bool success = e_1.measure(allpix,order,sigma_i,true,flag,0,&b_1);
        xdbg<<"e_1 Measure: success = "<<success<<std::endl;
        xdbg<<"e_1 = "<<e_1<<std::endl;
        xdbg<<"b_1 = "<<b_1.vec()<<std::endl;
        test(success,"Measure for Multi Image w/ PSF - success");
        test(std::abs(e_1.getCen()-ell.getCen()) < 1.e-5,
             "Measure for Multi Image w/ PSF - cen");
        test(std::abs(e_1.getGamma()-ell.getGamma()) < 1.e-5,
             "Measure for Multi Image w/ PSF - gam");
        test(std::abs(e_1.getMu()-ell.getMu()) < 1.e-5,
             "Measure for Multi Image w/ PSF - mu");
        test(Norm(b_1.vec().SubVector(0,15)-bin.vec()) < 1.e-5*Norm(bin.vec()),
             "Measure for Multi Image w/ PSF - b");
        test(Norm(b_1.vec().SubVector(15,b_1.size())) < 1.e-5*Norm(bin.vec()),
             "Measure 2 for Multi Image w/ PSF - b");


        // Next do the test with the PSFs

        std::vector<BVec> allpsf;
        std::vector<BVec> allpsf2;
        xdbg<<"ell = "<<ell<<std::endl;
        xdbg<<"b0 = "<<b0.vec()<<std::endl;

        for(int k=0;k<NPSFx;++k) {
            // allpsf are the psfs in the ell frame.
            allpsf.push_back(BVec(4,sigmaPsfx[k],bpsfx_vecs[k]));
            dbg<<"allpsf["<<k<<"] = "<<allpsf[k].vec()<<std::endl;

            // allpsf2 is the psf in the e0 frame.
            allpsf2.push_back(BVec(order+4,sigmaPsfx[k]));
            allpsf2[k] = allpsf[k];
            tmv::Matrix<double> S1(allpsf2[k].size(),allpsf2[k].size());
            calculateGTransform(g0,allpsf2[k].getOrder(),S1);
            tmv::Matrix<double> D1(allpsf2[k].size(),allpsf2[k].size());
            calculateMuTransform(m0,allpsf2[k].getOrder(),D1);
            allpsf2[k].vec() /= D1;
            allpsf2[k].vec() /= S1;
            xdbg<<"allpsf2[k] = "<<allpsf2[k].vec()<<std::endl;

            // c0 is the convolved galaxy shape in the ell frame.
            double sigma_o = sqrt(pow(sigma_i,2)+pow(sigmaPsfx[k],2));
            tmv::Matrix<double> C(45,45);
            calculatePsfConvolve(allpsf[k],8,sigma_i,C);
            BVec c0(8,sigma_o);
            c0.vec() = C.Cols(0,15) * b0.vec();
            xdbg<<"c0 = "<<c0.vec()<<std::endl;

            allpix[k].clear();
            getFakePixList(allpix[k],xcenx[k],ycenx[k],D,aperture,c0,ell);
        }

        Ellipse e_2;
        BVec b_2(order,sigma_i);
        success = e_2.measure(allpix,allpsf2,order,sigma_i,true,flag,0,&b_2);
        xdbg<<"sigma_i = "<<sigma_i<<std::endl;
        xdbg<<"e_2 = "<<e_2<<std::endl;
        xdbg<<"b_2 = "<<b_2.vec()<<std::endl;
        xdbg<<"flag = "<<flag<<std::endl;
        test(success,"Measure for Multi Image w/ PSF - success");
        test(std::abs(e_2.getCen()-ell.getCen()) < 1.e-5,
             "Measure for Multi Image w/ PSF - cen");
        test(std::abs(e_2.getGamma()-ell.getGamma()) < 1.e-5,
             "Measure for Multi Image w/ PSF - gam");
        test(std::abs(e_2.getMu()-ell.getMu()) < 1.e-5,
             "Measure for Multi Image w/ PSF - mu");
        test(Norm(b_2.vec().SubVector(0,15)-bin.vec()) < 1.e-5*Norm(bin.vec()),
             "Measure for Multi Image w/ PSF - b");
        test(Norm(b_2.vec().SubVector(15,b_2.size())) < 1.e-5*Norm(bin.vec()),
             "Measure 2 for Multi Image w/ PSF - b");

#if 0
        double dmu = 0.1;
        double sigma_ix = exp(dmu)*sigma_i;
        Ellipse e_3;
        BVec b_3(order,sigma_ix);
        success = e_3.measure(allpix,allpsf2,order,sigma_ix,true,flag,0,&b_3);
        xdbg<<"sigma_ix = "<<sigma_ix<<std::endl;
        xdbg<<"e_3 = "<<e_3<<std::endl;
        xdbg<<"b_3 = "<<b_3.vec()<<std::endl;
        test(success,"Measure for Multi Image w/ PSF #2 - success");
        test(std::abs(e_3.getCen()-ell.getCen()) < 1.e-5,
             "Measure for Multi Image #2 w/ PSF - cen");
        test(std::abs(e_3.getGamma()-ell.getGamma()) < 1.e-5,
             "Measure for Multi Image #2 w/ PSF - gam");
        test(std::abs(e_3.getMu()-ell.getMu()) < 1.e-5,
             "Measure for Multi Image #2 w/ PSF - mu");
        test(Norm(b_3.vec().SubVector(0,15)-bin.vec()) < 1.e-5*Norm(bin.vec()),
             "Measure for Multi Image #2 w/ PSF - b");
        test(Norm(b_3.vec().SubVector(15,b_3.size())) < 1.e-5*Norm(bin.vec()),
             "Measure 2 for Multi Image #2 w/ PSF - b");
#endif
    }
    std::cout<<"Passed tests of Ellipse Measure() method\n";
#endif

#ifdef TEST6
    ConfigFile params("testwl.config");
    params.load("findstars.params.btc","\t");
    params.load("fitsparams.config");

    // Read input files
    std::auto_ptr<Image<double> > weight_im;
    Image<double> im(params,weight_im);
    Transformation trans(params);
    InputCatalog incat(params);

    // Find stars
    StarCatalog starCat(incat,params,"");
    //starCat.calcSizes(im,weight_im.get(),trans);
    FindStarsLog fslog(params,"testfs.log");
    int nstars = starCat.findStars(fslog);
    dbg<<fslog<<std::endl;
    test(nstars >= 240,"Find Stars");

    // Test I/O
    starCat.write();
    StarCatalog starCat2(params,"");
    std::vector<std::string> all_ext = params["stars_ext"];
    for(size_t k=0;k<all_ext.size();++k) {
        dbg<<"Test I/O with extension "<<all_ext[k]<<std::endl;
        params["stars_ext"] = all_ext[k];
        if (k > 0) starCat2.read();
        test(starCat.size() == starCat2.size(),"starcat size I/O");
        for(size_t i=0;i<starCat.size();++i) {
            test(starCat.getId(i) == starCat2.getId(i), "starcat id I/O");
            test(std::abs(starCat.getPos(i) - starCat2.getPos(i)) <= 0.01,
                 "starcat pos I/O");
            test(std::abs(starCat.getSky(i) - starCat2.getSky(i)) <= 0.01,
                 "starcat sky I/O");
            test(std::abs(starCat.getNoise(i) - starCat2.getNoise(i)) <= 0.01,
                 "starcat noise I/O");
            test(starCat.getFlags(i) == starCat2.getFlags(i),
                 "starcat flags I/O");
            dbg<<"starcat.mag[i] = "<<starCat.getMag(i)<<
                ", starcat2.mag[i] = "<<starCat2.getMag(i)<<std::endl;
            dbg<<"abs(diff) = "<<
                std::abs(starCat.getMag(i) - starCat2.getMag(i))<<std::endl;
            test(std::abs(starCat.getMag(i) - starCat2.getMag(i)) <= 0.01,
                 "starcat mag I/O");
            test(std::abs(starCat.getObjSize(i)-starCat2.getObjSize(i)) <= 0.01,
                 "starcat objsize I/O");
            test(starCat.isAStar(i) == starCat2.isAStar(i),
                 "starcat isastar I/O");
        }
    }
    params["stars_ext"] = all_ext;

    // Measure PSF
    PsfCatalog psfcat(starCat,params);
    double sigmaP = psfcat.estimateSigma(im,weight_im.get(),trans);
    PsfLog psflog(params,"testpsf.log");
    int npsf = psfcat.measurePsf(im,weight_im.get(),trans,sigmaP,psflog);
    dbg<<psflog<<std::endl;
    // There are 244 stars in the file, but depending on the exact parameters,
    // one or two cross the edge, which gives an error flag.
    // The rest should all be measured successfully.
    test(npsf >= 240,"Measure PSF");

    // Test I/O
    psfcat.write();
    PsfCatalog psfcat2(params);
    all_ext = params["psf_ext"];
    for(size_t k=0;k<all_ext.size();++k) {
        dbg<<"Test I/O with extension "<<all_ext[k]<<std::endl;
        params["psf_ext"] = all_ext[k];
        if (k > 0) psfcat2.read();
        test(psfcat.size() == psfcat2.size(),"psfcat size I/O");
        for(size_t i=0;i<psfcat.size();++i) {
            test(psfcat.getId(i) == psfcat2.getId(i),"psfcat id I/O");
            test(std::abs(psfcat.getPos(i) - psfcat2.getPos(i)) <= 0.01,
                 "psfcat pos I/O");
            test(std::abs(psfcat.getSky(i) - psfcat2.getSky(i)) <= 0.01,
                 "psfcat sky I/O");
            test(std::abs(psfcat.getNoise(i) - psfcat2.getNoise(i)) <= 0.01,
                 "psfcat noise I/O");
            test(psfcat.getFlags(i) == psfcat2.getFlags(i),"psfcat flags I/O");
            test(std::abs(psfcat.getNu(i) - psfcat2.getNu(i)) <= 0.01,
                 "psfcat nu I/O");
            test(psfcat.getPsf(i).size() == psfcat2.getPsf(i).size(),
                 "psfcat psf.size I/O");
            test(std::abs(psfcat.getPsf(i).getSigma() -
                          psfcat2.getPsf(i).getSigma()) 
                 <= 0.01, "psfcat psf.getSigma I/O");
            test(Norm(psfcat.getPsf(i).vec() - psfcat2.getPsf(i).vec()) <= 
                 1.e-4*Norm(psfcat.getPsf(i).vec()), "psfcat psf vector I/O");
        }
    }
    params["psf_ext"] = all_ext;

    // Fit PSF
    FittedPsf fitpsf(psfcat,params);
    double rms = 0.; 
    int count = 0;
    for(int i=0;i<nstars;++i) if (!psfcat.getFlags(i)) {
        BVec checkpsf(fitpsf.getPsfOrder(),fitpsf.getSigma());
        checkpsf = fitpsf(starCat.getPos(i));
        double normsqdiff = NormSq(psfcat.getPsf(i).vec()-checkpsf.vec());
        rms += normsqdiff;
        ++count;
    }
    rms /= count;
    rms = sqrt(rms);
    dbg<<"fitpsf rms = "<<rms<<std::endl;
    test(rms < double(fitpsf.getPsfSize())/double(fitpsf.getFitSize()),
         "Fit PSF rms");

    // Test I/O
    fitpsf.write();
    FittedPsf fitpsf2(params);
    all_ext = params["fitpsf_ext"];
    for(size_t k=0;k<all_ext.size();++k) {
        dbg<<"Test I/O with extension "<<all_ext[k]<<std::endl;
        params["fitpsf_ext"] = all_ext[k];
        if (k > 0) fitpsf2.read();
        rms = 0.; 
        count = 0;
        test(fitpsf2.getPsfOrder() == fitpsf.getPsfOrder(),
             "FittedPSF I/O: order");
        test(std::abs(fitpsf2.getSigma() - fitpsf.getSigma()) < 0.01, 
             "FittedPSF I/O: sigma");
        for(int i=0;i<nstars;++i) if (!psfcat.getFlags(i)) {
            BVec checkpsf(fitpsf.getPsfOrder(),fitpsf.getSigma());
            checkpsf = fitpsf(starCat.getPos(i));
            BVec checkpsf2(fitpsf2.getPsfOrder(),fitpsf2.getSigma());
            checkpsf2 = fitpsf2(starCat.getPos(i));
            double normsqdiff = NormSq(checkpsf2.vec()-checkpsf.vec());
            rms += normsqdiff;
            ++count;
        }
        rms /= count;
        rms = sqrt(rms);
        dbg<<"fitpsf I/O rms = "<<rms<<std::endl;
        test(rms < 1.e-4,"Fit PSF I/O");
    }
    params["fitpsf_ext"] = all_ext;

    // Measure shears
    ShearCatalog shearcat(incat,trans,params);
    ShearLog shearlog(params,"testshear.log");
    int nshear = shearcat.measureShears(
        im,weight_im.get(),trans,fitpsf,shearlog);
    dbg<<shearlog<<std::endl;
    // There are 4557 galaxies in the file without error codes.
    // The code currently converges on more than 3000 of them,
    // although only about 2300 or so have no error flag.
    test(nshear >= 3000,"Measure Shear");

    // Test I/O
    shearcat.write();
    ShearCatalog shearcat2(params);
    all_ext = params["shear_ext"];
    for(size_t k=0;k<all_ext.size();++k) {
        dbg<<"Test I/O with extension "<<all_ext[k]<<std::endl;
        params["shear_ext"] = all_ext[k];
        if (k > 0) shearcat2.read();
        test(shearcat.size() == shearcat2.size(),"shearcat size I/O");
        for(size_t i=0;i<shearcat.size();++i) {
            test(shearcat.getId(i) == shearcat2.getId(i),"shearcat id I/O");
            test(std::abs(shearcat.getPos(i) - shearcat2.getPos(i)) <= 0.01,
                 "shearcat pos I/O");
            test(std::abs(shearcat.getSky(i) - shearcat2.getSky(i)) <= 0.01,
                 "shearcat sky I/O");
            test(std::abs(shearcat.getNoise(i) - shearcat2.getNoise(i)) <= 0.01,
                 "shearcat noise I/O");
            test(shearcat.getFlags(i) == shearcat2.getFlags(i),
                 "shearcat flags I/O");
            dbg<<"skypos = "<<shearcat.getSkyPos(i)<<"  "<<shearcat2.getSkyPos(i)<<std::endl;
            test(std::abs(shearcat.getSkyPos(i) - shearcat2.getSkyPos(i)) 
                 <= 0.01,
                 "shearcat skypos I/O");
            test(std::abs(shearcat.getShear(i) - shearcat2.getShear(i)) <= 0.01,
                 "shearcat shear I/O");
            test(std::abs(shearcat.getNu(i) - shearcat2.getNu(i)) <= 0.01,
                 "shearcat nu I/O");
            test(Norm(tmv::Matrix<double>(shearcat.getCov(i) - 
                                          shearcat2.getCov(i)))
                 <= 1.e-4*Norm(shearcat.getCov(i)),
                 "shearcat cov I/O");
            test(shearcat.getShape(i).size() == shearcat2.getShape(i).size(),
                 "shearcat shape.size I/O");
            test(std::abs(shearcat.getShape(i).getSigma() -
                          shearcat2.getShape(i).getSigma()) <= 0.01,
                 "shearcat shape.getSigma I/O");
            test(Norm(shearcat.getShape(i).vec() - shearcat2.getShape(i).vec()) 
                 <= 1.e-4*Norm(shearcat.getShape(i).vec()),
                 "shearcat shape vector I/O");
        }
    }
    params["shear_ext"] = all_ext;

    std::cout<<"Passed full pipeline\n";
#endif

    if (dbgout && dbgout != &std::cout) {delete dbgout; dbgout=0;}
    return 0;
}
#if 0
// Change to 1 to let gdb see where the program bombed out.
catch(int) {}
#else
catch (tmv::Error& e)
{
    std::cerr<<"Caught "<<e<<std::endl;
    return 1;
}
catch (std::exception& e)
{
    std::cerr<<"Caught std::exception: "<<e.what()<<std::endl;
}
catch (...)
{
    std::cerr<<"Caught Unknown error\n";
}
#endif
