
#include <valarray>
#include <fstream>
#include <iostream>
#include <cmath>

#ifdef _OPENMP
#include "omp.h"
#endif

#include "dbg.h"
#include "BVec.h"
#include "Ellipse.h"
#include "EllipseSolver.h"
#include "TestHelper.h"
#include "Log.h"
#include "InputCatalog.h"
#include "StarCatalog.h"
#include "PsfCatalog.h"
#include "FittedPsf.h"
#include "ShearCatalog.h"
#include "PsiHelper.h"
#include "BinomFact.h"

std::ostream* dbgout = 0;
bool XDEBUG = true;
//bool XDEBUG = false;

bool show_tests = false;
bool should_throw = true;
std::string last_success = "";
std::ostream* testout = &std::cout;

//#define MIN_TESTS  // Only do one of each testing loop.

//#define TEST1  // Basic measurements, ApplyZ/G/Mu
#define TEST2  // Deconvolving measurements, PsfConvolve
//#define TEST3  // Solve for Ellipse -- single image
//#define TEST4  // Solve for Ellipse -- multiple images
//#define TEST5  // Use Ellipse's measure function
//#define TEST6  // Test on real data images
//#define TEST7  // Compare with Gary's shapelet code

#ifdef TEST1
#define TEST12
#define TEST123
#define TEST12345
#endif
#ifdef TEST2
#define TEST12
#define TEST123
#define TEST237
#define TEST12345
#endif
#ifdef TEST3
#define TEST123
#define TEST34
#define TEST237
#define TEST12345
#endif
#ifdef TEST4
#define TEST34
#define TEST45
#define TEST12345
#endif
#ifdef TEST5
#define TEST45
#define TEST12345
#endif
#ifdef TEST7
#define TEST237
#endif

// The tests of the jacobians are pretty time consuming.
// They are worth testing, but once the code is working, I usually turn
// turn them off (by commenting out the next line) to save time.
//#define TESTJ

// Usually, we use the method measureShapelet to measure the shapelet
// decompositions.  Uncomment the next line to use altMeasureShapelet instead.
//#define USE_ALT

#ifdef TEST12345
inline void GetFakePixList(
    PixelList& pix,
    double xcen, double ycen,
    const DSmallMatrix22& D,
    double aperture, const BVec& b, const Ellipse& ell)
{
    dbg<<"Start GetFakePixList:\n";
    dbg<<"cen = "<<xcen<<','<<ycen<<std::endl;
    dbg<<"D = "<<D<<std::endl;
    dbg<<"ap = "<<aperture<<std::endl;
    dbg<<"b = "<<b<<std::endl;
    dbg<<"ell = "<<ell<<std::endl;
    const int border = b.getOrder();
    //Assert(border <= 8); // This is as high as I have implemented.

    double det = std::abs(D.TMV_det());
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
    // things like ApplyZ and such.  Not sure why yet.
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
            std::vector<std::complex<double> > ztothem(border+1);
            ztothem[0] = 1.;
            if (border > 0) ztothem[1] = z;
            for(int m=2;m<=border;++m) ztothem[m] = ztothem[m-1] * z;

            // The factors of 2 and -2 for m!=0 terms below are due to the fact
            // that we want real I values, so we assume b(q,p) = b*(p,q).
            // So I_pq = b(p,q) psi(p,q) + b(q,p) psi(q,p) 
            //         = 2 * Re(b(p,q) psi(p,q))
            //         = 2 * Re(b(p,q)) Re(psi(p,q)) - 2 * Im(b(p,q)) Im(psi(p,q))

#ifdef USE_TMV
            DVector::const_iterator bi = b.vec().begin();
#else
            const double* bi = b.vec().data();
#endif

#if 1
            double I = 0.;
            for(int n=0;n<=border;++n) {
                for(int m=n;m>=0;m-=2) {
                    int p = (n+m)/2;
                    int q = (n-m)/2;
                    //xdbg<<"n,m,p,q = "<<n<<','<<m<<','<<p<<','<<q<<std::endl;
                    //xdbg<<"b_pq = "<<*bi<<std::endl;

                    double laguerre = 0.;
                    double temp = 1.;
                    for(int k=m;k<=p;++k, temp*=-rsq) {
                        laguerre += binom(p,k) * temp / fact(k-m);
                    }
                    //xdbg<<"laguerre = "<<laguerre<<std::endl;
                    // So far, this is L_q^(m)(rsq)
                    // Still need (-)^q * sqrt(q!/p!) z^m
                    if (q % 2 == 1) laguerre = -laguerre;
                    laguerre *= sqrtfact(q)/sqrtfact(p);
                    //xdbg<<"laguerre => "<<laguerre<<std::endl;

                    if (m == 0) {
                        I += (*bi++) * laguerre;
                    } else {
                        I += 2. * (*bi++) * laguerre * real(ztothem[m]);
                        I -= 2. * (*bi++) * laguerre * imag(ztothem[m]);
                    }
                }
            }
#else
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
#endif

            I *= expfactor;

            pix.push_back(Pixel(std::real(globz),std::imag(globz),I,var));
        }
    }
}
#endif

#ifdef TEST2
inline void RtoC(const BVec& b, CDVector& bc)
{
    const int border = b.getOrder();
    Assert(int(bc.size()) == (border+2)*(border+2)/4);

    int kr=0, kc=0;
    for(int ppq=0;ppq<=border;++ppq) {
        for(int p=ppq,q=0;p>=q;--p,++q,++kc,++kr) {
            if (p>q) {
                bc(kc) = std::complex<double>(b(kr),b(kr+1)); ++kr;
            } else {
                bc(kc) = b(kr);
            }
        }
    }
    Assert(kr == int(b.size()));
    Assert(kc == int(bc.size()));
}

inline void CtoR(const CDVector& bc, BVec& b)
{
    const int border = b.getOrder();
    Assert(int(bc.size()) == (border+2)*(border+2)/4);

    int kr=0, kc=0;
    for(int ppq=0;ppq<=border;++ppq) {
        for(int p=ppq,q=0;p>=q;--p,++q,++kc,++kr) {
            b(kr) = std::real(bc[kc]);
            if (p>q) b(++kr) = std::imag(bc[kc]);
            else Assert(std::abs(std::imag(bc[kc])) < 1.e-5);
        }
    }
    Assert(kr == int(b.size()));
    Assert(kc == int(bc.size()));
}

inline void DirectConvolveB(const BVec& rbi, const BVec& rbp, BVec& rb)
    // This is based on direct convolutions of the various psi_pq o psi_st
    // combinations in Maple.  So it is independent of the formulae in 
    // BJ02.
{
    //const double twosqrtpi = 2.*sqrtpi;
    Assert(rbi.getOrder() == 4);
    Assert(rbp.getOrder() == 4);
    Assert(rb.getOrder() == 8);

    CDVector bi(9);
    CDVector bp(9);
    CDVector b(25);
    RtoC(rbi,bi);
    RtoC(rbp,bp);

    const int K00 = 0;
    const int K10 = 1;
    const int K20 = 2, K11 = 3;
    const int K30 = 4,  K21 = 5;
    const int K40 = 6, K31 = 7, K22 = 8;
    const int K50 = 9, K41 = 10, K32 = 11;
    const int K60 = 12, K51 = 13, K42 = 14, K33 = 15;
    const int K70 = 16, K61 = 17, K52 = 18, K43 = 19;
    const int K80 = 20, K71 = 21, K62 = 22, K53 = 23, K44 = 24;

    double sigma = rb.getSigma();
    double sigma_i = rbi.getSigma();
    double sigma_psf = rbp.getSigma();
    Assert(std::abs(sigma-sqrt(pow(sigma_i,2) + pow(sigma_psf,2)))<1.e-5);

    double fact = 1.;
    for(int ppq=0,k=0;ppq<=4;++ppq) {
        if (ppq>0) fact *= sigma_i/sigma;
        for(int p=ppq,q=0;p>=q;--p,++q,++k)
            bi[k] *= fact;
    }
    fact = 1.;
    for(int spt=0,m=0;spt<=4;++spt) {
        if (spt>0) fact *= sigma_psf/sigma;
        for(int s=spt,t=0;s>=t;--s,++t,++m)
            bp[m] *= fact;
    }

    double ri = sigma_i/sigma_psf; ri *= ri;
    double ri2 = ri*ri;
    double rp = sigma_psf/sigma_i; rp *= rp;
    double rp2 = rp*rp;

    b.setZero();

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

    b /= sqrtpi;
    CtoR(b,rb);
}
#endif

#ifdef TEST7
#include "gary/Laguerre.h"

template <class T>
TVector(T) convertVector(const mv::Vector<T>& v1)
{
    TVector(T) v2(int(v1.size()));
    for(int i=0;i<int(v1.size());++i)
        v2(i) = v1[i];
    return v2;
}
template <class T>
TMatrix(T) convertMatrix(const mv::Matrix<T>& m1)
{
    TMatrix(T) m2(int(m1.getM()),int(m1.getN()));
    for(int i=0;i<int(m1.getM());++i)
        for(int j=0;j<int(m1.getN());++j) 
            m2(i,j) = m1(i,j);
    return m2;
}
#endif

#define xdbgout (XDEBUG ? dbgout : 0)

int main(int argc, char **argv) try 
{
#if 1
    std::string dbgfile = "testwl.debug";
    dbgout = new std::ofstream(dbgfile.c_str());
#ifdef _OPENMP
    omp_set_dynamic(0); 
#ifndef __PGI
#pragma omp parallel copyin(dbgout, XDEBUG)
    {
        int threadnum = omp_get_thread_num();
        std::stringstream ss;
        ss << threadnum;
        std::string dbgfile2 = dbgfile + "_" + ss.str();
        if (threadnum > 0) dbgout = new std::ofstream(dbgfile2.c_str());
    }
#endif
#endif
#else
    dbgout = &std::cout;
#endif

#ifdef USE_TMV
    tmv::WriteWarningsTo(dbgout);
#endif

    // First check that the functions work very well in the limit of
    // well sampled (pixScale = 0.1), large aperture (aperture = 15.)
    // data with exact shapelet intensity patterns.
    // Read data will of course fail to work this well, but the answers
    // will nonetheless be as accurate as possible give the data available.

    DSmallMatrix22 D;
#ifdef USE_ALT
    double pixScale = 0.1;
    D.setToIdentity(pixScale);
#else
    D << 0.08, -0.02, -0.015, 0.11;
    double pixScale = sqrt(D.TMV_det());
#endif
    dbg<<"pixscale = "<<pixScale<<std::endl;

#ifdef TEST12
#ifdef MIN_TESTS
    const int NB = 1;
#else
    const int NB = 11;
#endif
    const int FIRSTB = 0;
    double b_vecs[NB][15] = {
        // 00  10  10  20  20  11  30  30  21  21  40  40  31  31  22
        //
#ifndef MIN_TESTS
        {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {1.0,0.5,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
        {1.0,0.0,0.0,0.5,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
        {1.0,0.0,0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
        {1.0,0.0,0.0,0.0,0.0,0.0,0.5,0.3,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
        {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.3,0.0,0.0,0.0,0.0,0.0},
        {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.3,0.0,0.0,0.0},
        {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.3,0.0},
        {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5},
        {10.0,0.5,0.1,0.2,0.3,0.4,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1},
#endif
        {2860,-143,286,-151,-202,358,160,-29,-126,253,-25,87,-109,-146,223} 
    };
#endif

#ifdef TEST1
#ifdef MIN_TESTS
    const int NLONG = 1;
#else
    const int NLONG = 2;
#endif
    double blong_vec[NLONG][45] = {
#ifndef MIN_TESTS
        {10.0,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,
            0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,
            0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.4},
#endif
        {8120,-122,331,224,-109,412,-142,-100,-399,200,103,-193,89,-407,142,
            -382,451,110,-538,287,-294,377,295,-275,165,-198,163,528,572,-312,
            -361,-189,140,226,155,-175,-372,238,210,-528,309,-241,314,-252,-111}
    };
#endif

#ifdef TEST237
#ifdef MIN_TESTS
    const int NPSF = 1;
#else
    const int NPSF = 11;
#endif
    const int FIRSTPSF = 0;
    double bpsf_vecs[NPSF][15] = {
#ifndef MIN_TESTS
        {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {1,0.1,0.2,0,0,0,0,0,0,0,0,0,0,0,0},
        {1,0,0,0.1,0.2,0,0,0,0,0,0,0,0,0,0},
        {1,0,0,0,0,0.2,0,0,0,0,0,0,0,0,0},
        {1,0,0,0,0,0,0.1,0.2,0,0,0,0,0,0,0},
        {1,0,0,0,0,0,0,0,0.1,0.2,0,0,0,0,0},
        {1,0,0,0,0,0,0,0,0,0,0.1,0.2,0,0,0},
        {1,0,0,0,0,0,0,0,0,0,0,0,0.1,0.2,0},
        {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2},
        {1.0,0.0,0.0,-0.2,0.1,0.1,-0.2,0.3,0.1,-0.2,0.2,0.1,-0.1,-0.2,-0.1},
#endif
        {1.0,0.22,0.31,-0.28,0.12,0.11,-0.18,0.30,0.09,-0.28,0.23,0.14,-0.12,-0.20,-0.09}
    };
    double sigma_psf = 1.5;
#endif

#ifdef TEST34
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
    double sigma_psfx[NPSFx] = { 1.3, 1.0, 1.1, 0.9, 1.1, 0.7, 1.2 };
#endif

#ifdef TEST123
    const double xcen = 0.709;
    const double ycen = 0.243;
#endif

#ifdef TEST12345
    double aperture = 20.;
#ifdef MIN_TESTS
    const int NELL = 1;
#else
    const int NELL = 8;
#endif
    const int FIRSTELL = 0;
    double ell_vecs[NELL][5] = {
#ifndef MIN_TESTS
        {0.0,0.0,0.0,0.0,0.0},
        {0.3,0.0,0.0,0.0,0.0},
        {0.0,0.3,0.0,0.0,0.0},
        {0.0,0.0,0.3,0.0,0.0},
        {0.0,0.0,0.0,0.3,0.0},
        {0.0,0.0,0.0,0.0,0.3},
        {0.1,-0.2,0.1,0.2,-0.2},
#endif
        {-0.11,0.24,-0.18,0.15,-0.14}
    };

    double sigma_i = 1.7;
    Ellipse e0; // all 0's.

    std::vector<PixelList> allpix(1);
#endif

#ifdef TEST1
    // First test the Ellipse preShiftBy function:
    for(int ie2 = 0; ie2 < NELL; ++ie2) {
        dbg<<"ie2 = "<<ie2<<std::endl;
        Ellipse e2(ell_vecs[ie2]);
        dbg<<"e2 = "<<e2<<std::endl;
        std::complex<double> c2 = e2.getCen();
        std::complex<double> g2 = e2.getGamma();
        double m2 = e2.getMu();
        for(int ie1 = 0; ie1 < NELL; ++ie1) {
            dbg<<"ie1 = "<<ie1<<std::endl;
            Ellipse e1(ell_vecs[ie1]);
            dbg<<"e1 = "<<e1<<std::endl;
            std::complex<double> c1 = e1.getCen();
            std::complex<double> g1 = e1.getGamma();
            double m1 = e1.getMu();
            Ellipse e3 = e2;
            e3.preShiftBy(c1,g1,m1);
            dbg<<"e3 -> "<<e3<<std::endl;
            std::complex<double> c3 = e3.getCen();
            std::complex<double> g3 = e3.getGamma();
            std::complex<double> m3(e3.getMu(),e3.getTheta());
            // The moments produced from the final Ellipse
            // should match the moments from doing one Ellipse at a time.
            double Ix=0., Iy=0., Ixx=0., Ixy=0., Iyy=0.;
            double Ix2=0., Iy2=0., Ixx2=0., Ixy2=0., Iyy2=0.;
            double Ix3=0., Iy3=0., Ixx3=0., Ixy3=0., Iyy3=0.;
            dbg<<"z      z1,z2       z3\n";
            for(double t = 0.; t < 360.; t+= 1.) {
                double t_rad = t * 3.141592653589793 / 180.;
                double x = cos(t_rad);
                double y = sin(t_rad);
                std::complex<double> z(x,y);
                std::complex<double> z1 = exp(-m1)/sqrt(1.-norm(g1)) *
                    ( z-c1 - g1 * conj(z-c1) );
                std::complex<double> z2 = exp(-m2)/sqrt(1.-norm(g2)) *
                    ( z1-c2 - g2 * conj(z1-c2) );
                std::complex<double> z3 = exp(-m3)/sqrt(1.-norm(g3)) *
                    ( z-c3 - g3 * conj(z-c3) );
                dbg<<z<<"  "<<z2<<"  "<<z3<<std::endl;
                Ix += x;
                Iy += y;
                Ixx += x*x;
                Ixy += x*y;
                Iyy += y*y;
                x = real(z2); y = imag(z2);
                Ix2 += x;
                Iy2 += y;
                Ixx2 += x*x;
                Ixy2 += x*y;
                Iyy2 += y*y;
                x = real(z3); y = imag(z3);
                Ix3 += x;
                Iy3 += y;
                Ixx3 += x*x;
                Ixy3 += x*y;
                Iyy3 += y*y;
            }
            dbg<<"Ix = "<<Ix<<"  "<<Ix2<<"  "<<Ix3<<std::endl;
            dbg<<"Iy = "<<Iy<<"  "<<Iy2<<"  "<<Iy3<<std::endl;
            dbg<<"Ixx = "<<Ixx<<"  "<<Ixx2<<"  "<<Ixx3<<std::endl;
            dbg<<"Ixy = "<<Ixy<<"  "<<Ixy2<<"  "<<Ixy3<<std::endl;
            dbg<<"Iyy = "<<Iyy<<"  "<<Iyy2<<"  "<<Iyy3<<std::endl;
            Test(std::abs(Ix2-Ix3) < 1.e-3, "shiftBy: Ix");
            Test(std::abs(Iy2-Iy3) < 1.e-3, "shiftBy: Iy");
            Test(std::abs(Ixx2-Ixx3) < 1.e-3, "shiftBy: Ixx");
            Test(std::abs(Ixy2-Ixy3) < 1.e-3, "shiftBy: Ixy");
            Test(std::abs(Iyy2-Iyy3) < 1.e-3, "shiftBy: Iyy");
        }
    }
    // Now test the ApplyZ, G, Mu routines along with basic b vector 
    // measurements.
    for(int ie = FIRSTELL; ie < NELL; ++ie) {
        Ellipse ell(ell_vecs[ie]);
        dbg<<"ell = "<<ell<<std::endl;
        for(int ib = FIRSTB; ib < NB; ++ib) {
            BVec bin(4,sigma_i,b_vecs[ib]);
            dbg<<"bin = "<<bin.vec()<<std::endl;
            int order = 12;
            int bsize = (order+1)*(order+2)/2;
            int order2 = 40;
            int maxm = order;

            allpix.resize(1);
            allpix[0].clear();
            GetFakePixList(allpix[0],xcen,ycen,D,aperture,bin,ell);

            BVec b0(order,sigma_i);
#ifdef USE_ALT
            ell.altMeasureShapelet(allpix,b0,order,order2,pixScale);
#else
            ell.measureShapelet(allpix,b0,order,order2,maxm);
#endif
            dbg<<"measured b = "<<b0.vec()<<std::endl;
            dbg<<ie<<"  "<<ib<<"  basic  ";
            dbg<<"NormInf(diff) = "<<
                (b0.vec().TMV_subVector(15,bsize)).TMV_normInf()<<std::endl;
            Test((b0.vec().TMV_subVector(0,15)-bin.vec()).TMV_normInf() < 
                 1.e-6*bin.vec().norm(),
                 "Basic B Measurement");
            Test((b0.vec().TMV_subVector(15,bsize)).TMV_normInf() < 
                 1.e-6*bin.vec().norm(),
                 "Basic B Measurement 2");

            std::complex<double> z1(0.2,0.3);
            Ellipse e1 = ell;
            e1.postShiftBy(z1,0.,0.);
            BVec b1(order,sigma_i);
#ifdef USE_ALT
            e1.altMeasureShapelet(allpix,b1,order,order2,pixScale);
#else
            e1.measureShapelet(allpix,b1,order,order2,maxm);
#endif
            xdbg<<"b(z1 = "<<z1<<") = "<<b1.vec()<<std::endl;
            BVec b2(order,sigma_i);
            ApplyZ(z1,b2=b0);
            xdbg<<"b(predicted) = "<<b2.vec()<<std::endl;
            dbg<<"diff for z = "<<z1<<" = "<<b2.vec()-b1.vec()<<std::endl;
            dbg<<ie<<"  "<<ib<<"  Z1  ";
            dbg<<"NormInf(diff) = "<<(b1.vec()-b2.vec()).TMV_normInf()<<std::endl;
            Test((b1.vec()-b2.vec()).TMV_subVector(0,15).TMV_normInf() < 
                 1.e-6*b2.vec().norm(),"ApplyZ 1");
            Test((b1.vec()-b2.vec()).TMV_normInf() < 1.e-4*b2.vec().norm(),
                 "ApplyZ 1b");
            e1 = ell;
            e1.postShiftBy(-z1,0.,0.);
#ifdef USE_ALT
            e1.altMeasureShapelet(allpix,b1,order,order2,pixScale);
#else
            e1.measureShapelet(allpix,b1,order,order2,maxm);
#endif
            ApplyZ(-z1,b2=b0);
            dbg<<"diff for z = "<<-z1<<" = "<<b2.vec()-b1.vec()<<std::endl;
            dbg<<ie<<"  "<<ib<<"  Z2  ";
            dbg<<"NormInf(diff) = "<<(b1.vec()-b2.vec()).TMV_normInf()<<std::endl;
            Test((b1.vec()-b2.vec()).TMV_subVector(0,15).TMV_normInf() < 
                 1.e-6*b2.vec().norm(),"ApplyZ 2");
            Test((b1.vec()-b2.vec()).TMV_normInf() < 1.e-4*b2.vec().norm(),
                 "ApplyZ 2b");

            double m1(0.2);
            e1 = ell;
            e1.postShiftBy(0.,0.,m1);
#ifdef USE_ALT
            e1.altMeasureShapelet(allpix,b1,order,order2,pixScale);
#else
            e1.measureShapelet(allpix,b1,order,order2,maxm);
#endif
            xdbg<<"b(mu = "<<m1<<") = "<<b1.vec()<<std::endl;
            ApplyMu(m1,b2=b0);
            xdbg<<"b(predicted) = "<<b2.vec()<<std::endl;
            dbg<<"diff for mu = "<<m1<<" = "<<b2.vec()-b1.vec()<<std::endl;
            dbg<<ie<<"  "<<ib<<"  Mu1  ";
            dbg<<"NormInf(diff) = "<<(b1.vec()-b2.vec()).TMV_normInf()<<std::endl;
            Test((b1.vec()-b2.vec()).TMV_subVector(0,15).TMV_normInf() < 
                 1.e-6*b2.vec().norm(),"ApplyMu 1");
            Test((b1.vec()-b2.vec()).TMV_normInf() < 1.e-4*b2.vec().norm(),
                 "ApplyMu 1b");

            e1 = ell;
            e1.postShiftBy(0.,0.,-m1);
#ifdef USE_ALT
            e1.altMeasureShapelet(allpix,b1,order,order2,pixScale);
#else
            e1.measureShapelet(allpix,b1,order,order2,maxm);
#endif
            ApplyMu(-m1,b2=b0);
            dbg<<"diff for mu = "<<-m1<<" = "<<b2.vec()-b1.vec()<<std::endl;
            dbg<<ie<<"  "<<ib<<"  Mu2  ";
            dbg<<"NormInf(diff) = "<<(b1.vec()-b2.vec()).TMV_normInf()<<std::endl;
            Test((b1.vec()-b2.vec()).TMV_subVector(0,15).TMV_normInf() < 
                 1.e-6*b2.vec().norm(),"ApplyMu 2");
            Test((b1.vec()-b2.vec()).TMV_normInf() < 1.e-4*b2.vec().norm(),
                 "ApplyMu 2b");

            std::complex<double> g1(0.2,0.1);
            e1 = ell;
            e1.postShiftBy(0.,g1,0.);
#ifdef USE_ALT
            e1.altMeasureShapelet(allpix,b1,order,order2,pixScale);
#else
            e1.measureShapelet(allpix,b1,order,order2,maxm);
#endif
            ApplyG(g1,b2=b0);
            xdbg<<"b(predicted) = "<<b2.vec()<<std::endl;
            dbg<<"diff for gamma = "<<g1<<" = "<<b2.vec()-b1.vec()<<std::endl;
            dbg<<ie<<"  "<<ib<<"  G1  ";
            dbg<<"NormInf(diff) = "<<(b1.vec()-b2.vec()).TMV_normInf()<<std::endl;
            Test((b1.vec()-b2.vec()).TMV_subVector(0,15).TMV_normInf() < 
                 1.e-6*b2.vec().norm(),"ApplyG 1");
            Test((b1.vec()-b2.vec()).TMV_normInf() < 1.e-4*b2.vec().norm(),
                 "ApplyG 1b");

            e1 = ell;
            e1.postShiftBy(0.,-g1,0.);
#ifdef USE_ALT
            e1.altMeasureShapelet(allpix,b1,order,order2,pixScale);
#else
            e1.measureShapelet(allpix,b1,order,order2,maxm);
#endif
            ApplyG(-g1,b2=b0);
            dbg<<"diff for gamma = "<<-g1<<" = "<<b2.vec()-b1.vec()<<std::endl;
            dbg<<ie<<"  "<<ib<<"  G2  ";
            dbg<<"NormInf(diff) = "<<(b1.vec()-b2.vec()).TMV_normInf()<<std::endl;
            Test((b1.vec()-b2.vec()).TMV_subVector(0,15).TMV_normInf() < 
                 1.e-6*b2.vec().norm(),"ApplyG 2");
            Test((b1.vec()-b2.vec()).TMV_normInf() < 1.e-4*b2.vec().norm(),
                 "ApplyG 2b");

            if (ib < NLONG) {
                BVec blong(8,sigma_i,blong_vec[ib]);

                allpix.resize(1);
                allpix[0].clear();
                GetFakePixList(allpix[0],xcen,ycen,D,aperture,blong,ell);

                BVec b(order,sigma_i);
#ifdef USE_ALT
                ell.altMeasureShapelet(allpix,b,order,order2,pixScale);
#else
                ell.measureShapelet(allpix,b,order,order2,maxm);
#endif
                dbg<<"blong = "<<blong.vec()<<std::endl;
                dbg<<"b = "<<b.vec()<<std::endl;
                Test((b.vec().TMV_subVector(0,45)-blong.vec()).TMV_normInf() < 
                     1.e-6*blong.vec().norm(),
                     "Basic Long B Measurement");
                Test((b.vec().TMV_subVector(45,b.size())).TMV_normInf() < 
                     1.e-6*blong.vec().norm(),
                     "Basic Long B Measurement (2)");

                double dmu = 0.01;
                double sigma_ix = exp(dmu) * sigma_i;
                BVec bx(order,sigma_ix);
                Ellipse e1 = ell;
                e1.postShiftBy(0.,0.,-dmu);
#ifdef USE_ALT
                e1.altMeasureShapelet(allpix,bx,order,order2,pixScale);
#else
                e1.measureShapelet(allpix,bx,order,order2,maxm);
#endif
                dbg<<"blong = "<<blong.vec()<<std::endl;
                dbg<<"bx = "<<bx.vec()<<std::endl;
                Test((bx.vec().TMV_subVector(0,45)-blong.vec()).TMV_normInf() < 
                     1.e-6*blong.vec().norm(),
                     "Basic Long B Measurement");
                Test((bx.vec().TMV_subVector(45,bx.size())).TMV_normInf() < 
                     1.e-6*blong.vec().norm(),
                     "Basic Long B Measurement (2)");
            }
        }
    }
    std::cout<<"Passed tests of ApplyZ, G, Mu and basic measurements\n";
#endif

#ifdef TEST2
    // Next test CalculatePsfConvolve
    for(int ie = FIRSTELL; ie < NELL; ++ie) {
        dbg<<"Start ie = "<<ie<<std::endl;
        Ellipse ell(ell_vecs[ie]);
        dbg<<"ell = "<<ell<<std::endl;
        std::complex<double> z0 = ell.getCen();
        std::complex<double> g0 = ell.getGamma();
        double m0 = ell.getMu();
        for(int ib = FIRSTB; ib < NB; ++ib) {
            dbg<<"Start ib = "<<ib<<std::endl;
            BVec bin(4,sigma_i,b_vecs[ib]);
            dbg<<"bin = "<<bin.vec()<<std::endl;
            for(int ip = FIRSTPSF; ip < NPSF; ++ip) {
                dbg<<"Start ip = "<<ip<<std::endl;
                // b0 is the galaxy shape in the e0 frame.
                BVec b0 = bin;

                int order = 12;
                int bsize = (order+1)*(order+2)/2;
                int order2 = 40;
                int maxm = order;
                int bsize2 = (order2+1)*(order2+2)/2;

                // bpsf is the psf in the e0 frame.
                BVec bpsf(4,sigma_psf,bpsf_vecs[ip]);
                dbg<<"bpsf = "<<bpsf.vec()<<std::endl;
                xdbg<<"ell = "<<ell<<std::endl;
                xdbg<<"b0 = "<<b0.vec()<<std::endl;

                // c0 = C(bpsf) * b0
                // This is one estimate of the observed galaxy in the e0
                // frame.
                double sigma_o = sqrt(pow(sigma_i,2)+pow(sigma_psf,2));
                DMatrix C(bsize2,bsize2);
                CalculatePsfConvolve(bpsf,order2,sigma_i,C);
                BVec c0(order2,sigma_o);
                //xdbg<<"C = "<<TMV_colRange(C,0,15)<<std::endl;
                c0.vec() = TMV_colRange(C,0,15) * b0.vec();
                xdbg<<"c0 = "<<c0.vec()<<std::endl;

                // c0x is another estimate of the same thing based on a 
                // direct convolution from Maple calculations rather than
                // the BJ02 recursion formulae.
                BVec c0x(8,sigma_o); 
                // This order must be 8.  Indeed all higher terms are == 0.
                DirectConvolveB(b0,bpsf,c0x);
                xdbg<<"c0x = "<<c0x<<std::endl;
                Test((c0.vec().TMV_subVector(0,c0x.size())-
                      c0x.vec()).TMV_normInf() < 1.e-6*c0.vec().norm(),
                     "PsfConvolve");
                Test(c0.vec().TMV_subVector(c0x.size(),c0.size()).TMV_normInf() <
                     1.e-6*c0.vec().norm(),"PsfConvolve");

                // be is the galaxy in ell fram:
                allpix.resize(1);
                allpix[0].clear();
                GetFakePixList(allpix[0],xcen,ycen,D,aperture,b0,e0);
                BVec be(order,sigma_i);
#ifdef USE_ALT
                ell.altMeasureShapelet(allpix,be,order,order2,pixScale);
#else
                ell.measureShapelet(allpix,be,order,order2,maxm);
#endif
                xdbg<<"be = "<<be.vec()<<std::endl;

                // bex is another estimate of the galaxy in the ell frame:
                BVec bex(order2,sigma_i);
                bex = b0;
                ApplyZ(z0,bex);
                ApplyG(g0,bex);
                ApplyMu(m0,bex);
                xdbg<<"bex = "<<bex.vec()<<std::endl;
                xdbg<<"NormInf(bex-be) = "<<
                    (bex.vec().TMV_subVector(0,bsize) - be.vec()).TMV_normInf()<<std::endl;
                Test((bex.vec().TMV_subVector(0,bsize) - be.vec()).TMV_normInf() <
                    1.e-6*(be.vec().norm()+b0.vec().norm()),
                    "View B in ell frame");

                // ce is the convolved galaxy in the ell frame:
                allpix.resize(1);
                allpix[0].clear();
                GetFakePixList(allpix[0],xcen,ycen,D,aperture,c0,e0);
                BVec ce(order,sigma_o);
#ifdef USE_ALT
                ell.altMeasureShapelet(allpix,ce,order,order2,pixScale);
#else
                ell.measureShapelet(allpix,ce,order,order2,maxm);
#endif
                xdbg<<"ce = "<<ce<<std::endl;
                xdbg<<"Norm(ce)+Norm(c0) = "<<
                    ce.vec().norm()+c0.vec().norm()<<std::endl;

                // cex is another estimate of the same: transform c0
                BVec cex(order2,sigma_o);
                cex = c0;
                ApplyZ(z0,cex);
                ApplyG(g0,cex);
                ApplyMu(m0,cex);
                xdbg<<"cex = "<<cex<<std::endl;
                xdbg<<"diff = "<<
                    cex.vec().TMV_subVector(0,bsize) - ce.vec()<<std::endl;
                xdbg<<"NormInf(cex-ce)(0:15) = "<<
                    (cex.vec().TMV_subVector(0,15) -
                     ce.vec().TMV_subVector(0,15)).TMV_normInf()<<std::endl;
                xdbg<<"NormInf(cex-ce)(0:28) = "<<
                    (cex.vec().TMV_subVector(0,28) -
                     ce.vec().TMV_subVector(0,28)).TMV_normInf()<<std::endl;
                xdbg<<"NormInf(cex-ce) = "<<
                    (cex.vec().TMV_subVector(0,bsize) - ce.vec()).TMV_normInf()<<std::endl;
                Test((cex.vec().TMV_subVector(0,15) -
                      ce.vec().TMV_subVector(0,15)).TMV_normInf() <
                     1.e-6*(ce.vec().norm()+c0.vec().norm()),
                     "Convolved B in ell frame");
                Test((cex.vec().TMV_subVector(0,28) -
                      ce.vec().TMV_subVector(0,28)).TMV_normInf() <
                     1.e-6*(ce.vec().norm()+c0.vec().norm()),
                     "Convolved B in ell frame");
                Test((cex.vec().TMV_subVector(0,bsize) - ce.vec()).TMV_normInf() <
                     1.e-6*(ce.vec().norm()+c0.vec().norm()),
                     "Convolved B in ell frame");

                // cey is yet another estimate of the same: convolve be
                BVec cey(order2,sigma_o);
                DMatrix Ce(bsize2,bsize2);
                BVec bpsfe(order2,sigma_psf);
                bpsfe = bpsf;
                ApplyG(g0,bpsfe);
                ApplyMu(m0,bpsfe);
                xdbg<<"psfe = "<<bpsfe<<std::endl;
                bpsfe.vec() *= exp(2.*m0);
                xdbg<<"psfe => "<<bpsfe<<std::endl;
                CalculatePsfConvolve(bpsfe,order2,sigma_i,Ce);
                cey.vec() = Ce * bex.vec();
                xdbg<<"cey = "<<cey<<std::endl;
                xdbg<<"ce = "<<ce<<std::endl;
                xdbg<<"diff = "<<(
                    cey.vec().TMV_subVector(0,bsize)-ce.vec())<<std::endl;
                xdbg<<"NormInf(cey-ce)(0:15) = "<<
                    (cey.vec().TMV_subVector(0,15) -
                     ce.vec().TMV_subVector(0,15)).TMV_normInf()<<std::endl;
                xdbg<<"NormInf(cey-ce)(0:28) = "<<
                    (cey.vec().TMV_subVector(0,28) -
                     ce.vec().TMV_subVector(0,28)).TMV_normInf()<<std::endl;
                xdbg<<"NormInf(cey-ce) = "<<
                    (cey.vec().TMV_subVector(0,bsize) - ce.vec()).TMV_normInf()<<std::endl;
                Test((cey.vec().TMV_subVector(0,15) -
                      ce.vec().TMV_subVector(0,15)).TMV_normInf() <
                     1.e-6*(ce.vec().norm()+c0.vec().norm()),
                     "Convolved B in ell frame");
                Test((cey.vec().TMV_subVector(0,28) -
                      ce.vec().TMV_subVector(0,28)).TMV_normInf() <
                     1.e-6*(ce.vec().norm()+c0.vec().norm()),
                     "Convolved B in ell frame");
                Test((cey.vec().TMV_subVector(0,bsize) - ce.vec()).TMV_normInf() <
                     1.e-6*(ce.vec().norm()+c0.vec().norm()),
                     "Convolved B in ell frame");

                // Since the cey calculation doesn't really match ce 
                // very well due to numerical issues from the finite 
                // order of the calculations, we recast the problem so 
                // b0 is the shape in the ell frame.
                allpix.resize(1);
                allpix[0].clear();
                std::vector<BVec> allpsf(1,bpsfe);
                GetFakePixList(allpix[0],xcen,ycen,D,aperture,c0,ell);

                BVec bes(order,sigma_i);
#ifdef USE_ALT
                ell.altMeasureShapelet(allpix,allpsf,bes,order,order2,pixScale);
#else
                ell.measureShapelet(allpix,allpsf,bes,order,order2,maxm);
#endif
                xdbg<<"bes = "<<bes.vec()<<std::endl;
                xdbg<<"b0 = "<<b0.vec()<<std::endl;
                xdbg<<"NormInf(bes-b0) = "<<
                    (bes.vec().TMV_subVector(0,15)-b0.vec()).TMV_normInf()<<std::endl;
                xdbg<<"NormInf(bes(15:)) = "<<
                    bes.vec().TMV_subVector(15,bes.size()).TMV_normInf()<<std::endl;
                Test((bes.vec().TMV_subVector(0,15)-b0.vec()).TMV_normInf() < 
                     1.e-6*(b0.vec().norm()+c0.vec().norm()),
                     "Deconvolving solver BE Measurement");
                Test(bes.vec().TMV_subVector(15,bes.size()).TMV_normInf() <
                     1.e-6*(b0.vec().norm()+c0.vec().norm()),
                     "Deconvolving solver BE Measurement");

                // Now test the deconvolving measurement after a shift
                double eps = 1.e-3 * (b0.vec().norm() + bpsfe.vec().norm());

                std::complex<double> z1(0.2,0.3);
                Ellipse e1 = ell;
                e1.postShiftBy(z1,0.,0.);
                BVec b1(order,sigma_i);
#ifdef USE_ALT
                e1.altMeasureShapelet(allpix,allpsf,b1,order,order2,pixScale);
#else
                e1.measureShapelet(allpix,allpsf,b1,order,order2,maxm);
#endif
                xdbg<<"b(z1 = "<<z1<<") = "<<b1.vec()<<std::endl;
                BVec b2(order,sigma_i);
                ApplyZ(z1,b2=b0);
                xdbg<<"b(predicted) = "<<b2.vec()<<std::endl;
                dbg<<"diff for z = "<<z1<<" = "<<b2.vec()-b1.vec()<<std::endl;
                dbg<<ie<<"  "<<ib<<"  "<<ip<<"  Z1  ";
                dbg<<"NormInf(diff) = "<<(b1.vec()-b2.vec()).TMV_normInf()<<std::endl;
                dbg<<"eps = "<<eps<<std::endl;
                Test((b1.vec()-b2.vec()).TMV_subVector(0,15).TMV_normInf() < eps,
                     "Deconvolving ApplyZ 1");
                Test((b1.vec()-b2.vec()).TMV_normInf() < 10.*eps,
                     "Deconvolving ApplyZ 1b");
                e1 = ell;
                e1.postShiftBy(-z1,0.,0.);
#ifdef USE_ALT
                e1.altMeasureShapelet(allpix,allpsf,b1,order,order2,pixScale);
#else
                e1.measureShapelet(allpix,allpsf,b1,order,order2,maxm);
#endif
                xdbg<<"b(z1 = "<<-z1<<") = "<<b1.vec()<<std::endl;
                ApplyZ(-z1,b2=b0);
                xdbg<<"b(predicted) = "<<b2.vec()<<std::endl;
                dbg<<"diff for z = "<<-z1<<" = "<<b2.vec()-b1.vec()<<std::endl;
                dbg<<ie<<"  "<<ib<<"  "<<ip<<"  Z2  ";
                dbg<<"NormInf(diff) = "<<(b1.vec()-b2.vec()).TMV_normInf()<<std::endl;
                dbg<<"eps = "<<eps<<std::endl;
                Test((b1.vec()-b2.vec()).TMV_subVector(0,15).TMV_normInf() < eps,
                     "Deconvolving ApplyZ 2");
                Test((b1.vec()-b2.vec()).TMV_normInf() < 10.*eps,
                     "Deconvolving ApplyZ 2b");

                double m1(0.2);
                e1 = ell;
                e1.postShiftBy(0.,0.,m1);
#ifdef USE_ALT
                e1.altMeasureShapelet(allpix,allpsf,b1,order,order2,pixScale);
#else
                e1.measureShapelet(allpix,allpsf,b1,order,order2,maxm);
#endif
                xdbg<<"b(mu = "<<m1<<") = "<<b1.vec()<<std::endl;
                ApplyMu(m1,b2=b0);
                xdbg<<"b(predicted) = "<<b2.vec()<<std::endl;
                dbg<<"diff for mu = "<<m1<<" = "<<b2.vec()-b1.vec()<<std::endl;
                dbg<<ie<<"  "<<ib<<"  "<<ip<<"  Mu1  ";
                dbg<<"NormInf(diff) = "<<(b1.vec()-b2.vec()).TMV_normInf()<<std::endl;
                dbg<<"eps = "<<eps<<std::endl;
                Test((b1.vec()-b2.vec()).TMV_subVector(0,15).TMV_normInf() < eps,
                     "Deconvolving ApplyMu 1");
                Test((b1.vec()-b2.vec()).TMV_normInf() < 10.*eps,
                     "Deconvolving ApplyMu 1b");

                e1 = ell;
                e1.postShiftBy(0.,0.,-m1);
#ifdef USE_ALT
                e1.altMeasureShapelet(allpix,allpsf,b1,order,order2,pixScale);
#else
                e1.measureShapelet(allpix,allpsf,b1,order,order2,maxm);
#endif
                xdbg<<"b(mu = "<<-m1<<") = "<<b1.vec()<<std::endl;
                ApplyMu(-m1,b2=b0);
                xdbg<<"b(predicted) = "<<b2.vec()<<std::endl;
                dbg<<"diff for mu = "<<-m1<<" = "<<b2.vec()-b1.vec()<<std::endl;
                dbg<<ie<<"  "<<ib<<"  "<<ip<<"  Mu2  ";
                dbg<<"NormInf(diff) = "<<(b1.vec()-b2.vec()).TMV_normInf()<<std::endl;
                dbg<<"eps = "<<eps<<std::endl;
                Test((b1.vec()-b2.vec()).TMV_subVector(0,15).TMV_normInf() < eps,
                     "Deconvolving ApplyMu 2");
                Test((b1.vec()-b2.vec()).TMV_normInf() < 10.*eps,
                     "Deconvolving ApplyMu 2b");

                std::complex<double> g1(0.2,0.1);
                e1 = ell;
                e1.postShiftBy(0.,g1,0.);
#ifdef USE_ALT
                e1.altMeasureShapelet(allpix,allpsf,b1,order,order2,pixScale);
#else
                e1.measureShapelet(allpix,allpsf,b1,order,order2,maxm);
#endif
                xdbg<<"b(g = "<<g1<<") = "<<b1.vec()<<std::endl;
                ApplyG(g1,b2=b0);
                xdbg<<"b(predicted) = "<<b2.vec()<<std::endl;
                dbg<<"diff for gamma = "<<g1<<" = "<<b2.vec()-b1.vec()<<std::endl;
                dbg<<ie<<"  "<<ib<<"  "<<ip<<"  G1  ";
                dbg<<"NormInf(diff) = "<<(b1.vec()-b2.vec()).TMV_normInf()<<std::endl;
                dbg<<"eps = "<<eps<<std::endl;
                Test((b1.vec()-b2.vec()).TMV_subVector(0,15).TMV_normInf() < eps,
                     "Deconvolving ApplyG 1");
                Test((b1.vec()-b2.vec()).TMV_normInf() < 10.*eps,
                     "Deconvolving ApplyG 1b");

                e1 = ell;
                e1.postShiftBy(0.,-g1,0.);
#ifdef USE_ALT
                e1.altMeasureShapelet(allpix,allpsf,b1,order,order2,pixScale);
#else
                e1.measureShapelet(allpix,allpsf,b1,order,order2,maxm);
#endif
                xdbg<<"b(g = "<<-g1<<") = "<<b1.vec()<<std::endl;
                ApplyG(-g1,b2=b0);
                xdbg<<"b(predicted) = "<<b2.vec()<<std::endl;
                dbg<<"diff for gamma = "<<-g1<<" = "<<b2.vec()-b1.vec()<<std::endl;
                dbg<<ie<<"  "<<ib<<"  "<<ip<<"  G2  ";
                dbg<<"NormInf(diff) = "<<(b1.vec()-b2.vec()).TMV_normInf()<<std::endl;
                dbg<<"eps = "<<eps<<std::endl;
                Test((b1.vec()-b2.vec()).TMV_subVector(0,15).TMV_normInf() < eps,
                     "Deconvolving ApplyG 2");
                Test((b1.vec()-b2.vec()).TMV_normInf() < 10.*eps,
                     "Deconvolving ApplyG 2b");
                for(int ie2 = 0; ie2 < NELL; ++ie2) {
                    dbg<<"Start ie2 = "<<ie2<<std::endl;
                    Ellipse e2(ell_vecs[ie2]);
                    std::complex<double> z2 = e2.getCen();
                    std::complex<double> g2 = e2.getGamma();
                    double m2 = e2.getMu();
                    dbg<<"z2 = "<<z2<<std::endl;
                    dbg<<"g2 = "<<g2<<std::endl;
                    dbg<<"m2 = "<<m2<<std::endl;
                    e1 = ell;
                    e1.postShiftBy(z2,g2,m2);
#ifdef USE_ALT
                    e1.altMeasureShapelet(allpix,allpsf,b1,order,order2,pixScale);
#else
                    e1.measureShapelet(allpix,allpsf,b1,order,order2,maxm);
#endif
                    dbg<<"b1 = "<<b1<<std::endl;
                    BVec b3(order2,sigma_i);
                    ApplyZ(z2,b3=b0);
                    ApplyG(g2,b3);
                    ApplyMu(m2,b3);
                    dbg<<"b3 = "<<b3<<std::endl;
                    dbg<<"diff = "<<
                        b3.vec().TMV_subVector(0,bsize)-b1.vec()<<std::endl;
                    dbg<<ie<<"  "<<ib<<"  "<<ip<<"  "<<ie2<<"  ";
                    dbg<<"NormInf(diff) = "<<
                        (b3.vec().TMV_subVector(0,bsize)-b1.vec()).TMV_normInf()<<std::endl;
                    dbg<<"eps = "<<eps<<std::endl;
                    Test((b3.vec().TMV_subVector(0,15) -
                          b1.vec().TMV_subVector(0,15)).TMV_normInf() < 100.*eps,
                         "Deconvolving shift a");
                    Test((b3.vec().TMV_subVector(0,bsize) -
                          b1.vec()).TMV_normInf() < 1000.*eps,
                         "Deconvolving shift b");
                }
            }
        }
    }
    std::cout<<"Passed tests of PsfConvolve and deconvolving measurements.\n";
#endif

#ifdef TEST3
    // Next test Solve for the coordinates in which a Gaussian galaxy is round
    for(int ie = 0; ie < NELL; ++ie) {
        dbg<<"Start ie = "<<ie<<std::endl;
        Ellipse ell(ell_vecs[ie]);
        dbg<<"ell = "<<ell<<std::endl;
        std::complex<double> g0 = ell.getGamma();
        double m0 = ell.getMu();

        BVec bin(4,sigma_i);
        bin.vec().setZero();
        bin(0) = 1.;
        dbg<<"bin = "<<bin.vec()<<std::endl;
        int order = 8;
#ifdef USE_TMV
        DVector xe(5,ell_vecs[ie]);
#else
        DVector xe(5);
        for(int i=0;i<5;++i) xe(i) = ell_vecs[ie][i];
#endif

        // b0 is the galaxy shape in the ell frame.
        BVec b0 = bin;
        allpix.resize(1);
        allpix[0].clear();
        GetFakePixList(allpix[0],xcen,ycen,D,aperture,b0,ell);
        if (stype == 0) s.reset(
            new EllipseSolver(allpix,order,sigma_i));
        else s.reset(
            new EllipseSolver2(allpix,order,sigma_i,pixScale));
        DVector x(5);
        x.setZero();
        DVector f(5);
        if (XDEBUG) s->setOutput(*dbgout);
        s->useDogleg();
        s->setDelta0(0.1);
        s->setFTol(1.e-6);
#ifdef TESTJ
        Test(s->testJ(x,f,xdbgout,1.e-6),
             "TestJ - Single Image no PSF");
#endif
        bool success = s->solve(x,f);
        Test(success,"No PSF Solve - success");
        Test((x-xe).TMV_normInf() < 1.e-5,"No PSF Solve - x == xe");

        for(int ip = 0; ip < NPSF; ++ip) {
            dbg<<"Start ip = "<<ip<<std::endl;

            // bpsf is the psf in the ell frame.
            BVec bpsf(4,sigma_psf,bpsf_vecs[ip]);
            dbg<<"bpsf = "<<bpsf.vec()<<std::endl;
            xdbg<<"ell = "<<ell<<std::endl;
            xdbg<<"b0 = "<<b0.vec()<<std::endl;

            // bpsf2 is the psf in the e0 frame.
            BVec bpsf2(order+4,sigma_psf);
            bpsf2 = bpsf;
            DMatrix S1(bpsf2.size(),bpsf2.size());
            CalculateGTransform(g0,bpsf2.getOrder(),S1);
            DMatrix D1(bpsf2.size(),bpsf2.size());
            CalculateMuTransform(m0,bpsf2.getOrder(),D1);
#ifdef USE_TMV
            bpsf2.vec() /= D1;
            bpsf2.vec() /= S1;
#else
            D1.lu().solve(bpsf2.vec(),&bpsf2.vec());
            S1.lu().solve(bpsf2.vec(),&bpsf2.vec());
#endif

            // c0 is the convolved galaxy shape in the ell frame.
            double sigma_o = sqrt(pow(sigma_i,2)+pow(sigma_psf,2));
            DMatrix C(45,45);
            CalculatePsfConvolve(bpsf,8,sigma_i,C);
            BVec c0(8,sigma_o);
            c0.vec() = TMV_colRange(C,0,15) * b0.vec();
            xdbg<<"c0 = "<<c0.vec()<<std::endl;
            allpix.resize(1);
            allpix[0].clear();
            std::vector<BVec> allpsf(1,bpsf2);
            GetFakePixList(allpix[0],xcen,ycen,D,aperture,c0,ell);

            if (stype == 0) s.reset(
                new EllipseSolver(allpix,allpsf,f_psf,order,sigma_i));
            else s.reset(
                new EllipseSolver2(allpix,allpsf,f_psf,order,sigma_i,
                                   pixScale));
            s->calculateF(xe,f);
            BVec be = s->getB();
            xdbg<<"be = "<<be.vec()<<std::endl;
            Test((be.vec().TMV_subVector(0,15)-bin.vec()).TMV_normInf() < 
                 1.e-6*bin.vec().norm(),
                 "B Measurement for Single Image w/ PSF");
            Test((be.vec().TMV_subVector(15,be.size())).TMV_normInf() < 
                 3.e-6*bin.vec().norm(),
                 "B Measurement 2 for Single Image w/ PSF");
            x.setZero();
            if (XDEBUG) s->setOutput(*dbgout);
            s->useDogleg();
            s->setDelta0(0.1);
            s->setFTol(1.e-6);
#ifdef TESTJ
            Test(s->testJ(x,f,xdbgout,1.e-6),
                 "TestJ - Single Image w/ PSF");
#endif
            success = s->solve(x,f);
#ifdef TESTJ
            Test(s->testJ(x,f,xdbgout,1.e-6),
                 "TestJ - Single Image w/ PSF");
#endif
            Test(success,"PSF Solve - success");
            Test((x-xe).TMV_normInf() < 1.e-5,"PSF Solve - x == xe");
        }
    }
    std::cout<<"Passed tests of solve for Ellipse.\n";
#endif

#ifdef TEST4
    // Next test Solve with multiple PSFs
    for(int ie = 0; ie < NELL; ++ie) {
        dbg<<"Start ie = "<<ie<<std::endl;
        Ellipse ell(ell_vecs[ie]);
        dbg<<"ell = "<<ell<<std::endl;
        std::complex<double> g0 = ell.getGamma();
        double m0 = ell.getMu();

        BVec bin(4,sigma_i);
        bin.vec().setZero();
        bin(0) = 1.;
        dbg<<"bin = "<<bin.vec()<<std::endl;
        int order = 8;
#ifdef USE_TMV
        DVector xe(5,ell_vecs[ie]);
#else
        DVector xe(5);
        for(int i=0;i<5;++i) xe(i) = ell_vecs[ie][i];
#endif

        // b0 is the galaxy shape in the ell frame.
        BVec b0 = bin;
        allpix.resize(NPSFx);

#if 1
        {
            // First do the test without the PSFs
            for(int k=0;k<NPSFx;++k) {
                allpix[k].clear();
                GetFakePixList(
                    allpix[k],xcenx[k],ycenx[k],D,aperture,b0,ell);
            }
            if (stype == 0) s.reset(
                new EllipseSolver(allpix,order,sigma_i));
            else s.reset(
                new EllipseSolver2(allpix,order,sigma_i,pixScale));
            DVector x(5);
            x.setZero();
            DVector f(5);
            s->calculateF(xe,f);
            BVec be = s->getB();
            xdbg<<"be = "<<be.vec()<<std::endl;
            xdbg<<"NormInf(be(0:15)-bin) = "<<
                (be.vec().TMV_subVector(0,15)-bin.vec()).TMV_normInf()<<std::endl;
            xdbg<<"NormInf(be(15:)) = "<<
                be.vec().TMV_subVector(15,be.size()).TMV_normInf()<<std::endl;
            Test((be.vec().TMV_subVector(0,15)-bin.vec()).TMV_normInf() < 
                 1.e-6*bin.vec().norm(),
                 "B Measurement for Multi Image");
            Test(be.vec().TMV_subVector(15,be.size()).TMV_normInf() < 
                 3.e-6*bin.vec().norm(),
                 "B Measurement 2 for Multi Image");

            if (XDEBUG) s->setOutput(*dbgout);
            s->useDogleg();
            s->setDelta0(0.1);
            s->setFTol(1.e-6);
#ifdef TESTJ
            Test(s->testJ(x,f,xdbgout,1.e-6),
                 "TestJ - Multi Image no PSF 1");
#endif
            bool success = s->solve(x,f);
            xdbg<<"x = "<<x<<std::endl;
            xdbg<<"xe = "<<xe<<std::endl;
            xdbg<<"f = "<<f<<std::endl;
            xdbg<<"success = "<<success<<std::endl;
            Test(success,"No PSF Solve - success");
            Test((x-xe).TMV_normInf() < 1.e-5,"No PSF Solve - x == xe");
#ifdef TESTJ
            Test(s->testJ(x,f,xdbgout,1.e-6),
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
                allpsf.push_back(BVec(4,sigma_psfx[k],bpsfx_vecs[k]));
                dbg<<"allpsf["<<k<<"] = "<<allpsf[k].vec()<<std::endl;

                // allpsf2 is the psf in the e0 frame.
                allpsf2.push_back(BVec(order+4,sigma_psfx[k]));
                allpsf2[k] = allpsf[k];
                DMatrix S1(allpsf2[k].size(),allpsf2[k].size());
                CalculateGTransform(g0,allpsf2[k].getOrder(),S1);
                DMatrix D1(allpsf2[k].size(),allpsf2[k].size());
                CalculateMuTransform(m0,allpsf2[k].getOrder(),D1);
#ifdef USE_TMV
                allpsf2[k].vec() /= D1;
                allpsf2[k].vec() /= S1;
#else
                D1.lu().solve(allpsf2[k].vec(),&allpsf2[k].vec());
                S1.lu().solve(allpsf2[k].vec(),&allpsf2[k].vec());
#endif
                xdbg<<"allpsf2[k] = "<<allpsf2[k].vec()<<std::endl;

                // c0 is the convolved galaxy shape in the ell frame.
                double sigma_o = sqrt(pow(sigma_i,2)+pow(sigma_psfx[k],2));
                DMatrix C(45,45);
                CalculatePsfConvolve(allpsf[k],8,sigma_i,C);
                BVec c0(8,sigma_o);
                c0.vec() = TMV_colRange(C,0,15) * b0.vec();
                xdbg<<"c0 = "<<c0.vec()<<std::endl;

                allpix[k].clear();
                GetFakePixList(
                    allpix[k],xcenx[k],ycenx[k],D,aperture,c0,ell);
            }

            if (stype == 0) s.reset(
                new EllipseSolver(allpix,allpsf2,f_psf,order,sigma_i));
            else s.reset(
                new EllipseSolver2(allpix,allpsf2,f_psf,order,sigma_i,
                                   pixScale));

            DVector x(5);
            x.setZero();
            DVector f(5);
            s->calculateF(xe,f);
            BVec be = s->getB();
            xdbg<<"be = "<<be.vec()<<std::endl;
            xdbg<<"NormInf(be(0:15)-bin) = "<<
                (be.vec().TMV_subVector(0,15)-bin.vec()).TMV_normInf()<<std::endl;
            xdbg<<"NormInf(be(15:)) = "<<
                be.vec().TMV_subVector(15,be.size()).TMV_normInf()<<std::endl;
            Test((be.vec().TMV_subVector(0,15)-bin.vec()).TMV_normInf() < 
                 1.e-6*bin.vec().norm(),
                 "B Measurement for Multi Image w/ PSF");
            Test(be.vec().TMV_subVector(15,be.size()).TMV_normInf() < 
                 3.e-6*bin.vec().norm(),
                 "B Measurement 2 for Multi Image w/ PSF");

            x.setZero();
            if (XDEBUG) s->setOutput(*dbgout);
            s->useDogleg();
            s->setDelta0(0.1);
            s->setFTol(1.e-6);
#ifdef TESTJ
            Test(s->testJ(x,f,xdbgout,1.e-6),
                 "TestJ - Multi Image w/ PSF 1");
#endif
            bool success = s->solve(x,f);
            xdbg<<"x = "<<x<<std::endl;
            xdbg<<"xe = "<<xe<<std::endl;
            xdbg<<"f = "<<f<<std::endl;
            xdbg<<"success = "<<success<<std::endl;
#ifdef TESTJ
            Test(s->testJ(x,f,xdbgout,1.e-6),
                 "TestJ - Multi Image w/ PSF 2");
#endif
            Test(success,"PSF Solve - success");
            Test((x-xe).TMV_normInf() < 1.e-5,"PSF Solve - x == xe");
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
        bin.vec().setZero();
        bin(0) = 1.;
        xdbg<<"bin = "<<bin.vec()<<std::endl;
        int order = 6;
#ifdef USE_TMV
        DVector xe(5,ell_vecs[ie]);
#else
        DVector xe(5);
        for(int i=0;i<5;++i) xe(i) = ell_vecs[ie][i];
#endif

        // b0 is the galaxy shape in the ell frame.
        BVec b0 = bin;
        allpix.resize(NPSFx,0);

        // Test CrudeMeasure:
        // CrudeMeasure only works precisely for gamma = 0 shapes.
        Ellipse ell_nogamma = ell;
        ell_nogamma.setGamma(0.);
        for(int k=0;k<NPSFx;++k) {
            allpix[k].clear();
            GetFakePixList(allpix[k],xcenx[k],ycenx[k],D,aperture,b0,
                           ell_nogamma);
        }
        Ellipse e_0;
        e_0.crudeMeasure(allpix[0],sigma_i);
        xdbg<<"Done e_0 Measure\n";
        xdbg<<"e_0 = "<<e_0<<std::endl;
        Test(std::abs(e_0.getCen()-ell.getCen()) < 1.e-3, "CrudeMeasure - cen");
        Test(std::abs(e_0.getGamma()) < 1.e-3, "CrudeMeasure - gam");
        Test(std::abs(e_0.getMu()-ell.getMu()) < 1.e-3, "CrudeMeasure - mu");

        if (ell.getGamma() != 0.) {
            for(int k=0;k<NPSFx;++k) {
                allpix[k].clear();
                GetFakePixList(allpix[k],xcenx[k],ycenx[k],D,aperture,b0,ell);
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
        Test(success,"Measure for Multi Image w/ PSF - success");
        Test(std::abs(e_1.getCen()-ell.getCen()) < 1.e-5,
             "Measure for Multi Image w/ PSF - cen");
        Test(std::abs(e_1.getGamma()-ell.getGamma()) < 1.e-5,
             "Measure for Multi Image w/ PSF - gam");
        Test(std::abs(e_1.getMu()-ell.getMu()) < 1.e-5,
             "Measure for Multi Image w/ PSF - mu");
        Test((b_1.vec().TMV_subVector(0,15)-bin.vec()).TMV_normInf() < 1.e-5*bin.vec().norm(),
             "Measure for Multi Image w/ PSF - b");
        Test(b_1.vec().TMV_subVector(15,b_1.size()).TMV_normInf() < 1.e-5*bin.vec().norm(),
             "Measure 2 for Multi Image w/ PSF - b");


        // Next do the test with the PSFs

        std::vector<BVec> allpsf;
        std::vector<BVec> allpsf2;
        xdbg<<"ell = "<<ell<<std::endl;
        xdbg<<"b0 = "<<b0.vec()<<std::endl;

        for(int k=0;k<NPSFx;++k) {
            // allpsf are the psfs in the ell frame.
            allpsf.push_back(BVec(4,sigma_psfx[k],bpsfx_vecs[k]));
            dbg<<"allpsf["<<k<<"] = "<<allpsf[k].vec()<<std::endl;

            // allpsf2 is the psf in the e0 frame.
            allpsf2.push_back(BVec(order+4,sigma_psfx[k]));
            allpsf2[k] = allpsf[k];
            DMatrix S1(allpsf2[k].size(),allpsf2[k].size());
            CalculateGTransform(g0,allpsf2[k].getOrder(),S1);
            DMatrix D1(allpsf2[k].size(),allpsf2[k].size());
            CalculateMuTransform(m0,allpsf2[k].getOrder(),D1);
#ifdef USE_TMV
            allpsf2[k].vec() /= D1;
            allpsf2[k].vec() /= S1;
#else
            D1.lu().solve(allpsf2[k].vec(),&allpsf2[k].vec());
            S1.lu().solve(allpsf2[k].vec(),&allpsf2[k].vec());
#endif
            xdbg<<"allpsf2[k] = "<<allpsf2[k].vec()<<std::endl;

            // c0 is the convolved galaxy shape in the ell frame.
            double sigma_o = sqrt(pow(sigma_i,2)+pow(sigma_psfx[k],2));
            DMatrix C(45,45);
            CalculatePsfConvolve(allpsf[k],8,sigma_i,C);
            BVec c0(8,sigma_o);
            c0.vec() = TMV_colRange(C,0,15) * b0.vec();
            xdbg<<"c0 = "<<c0.vec()<<std::endl;

            allpix[k].clear();
            GetFakePixList(allpix[k],xcenx[k],ycenx[k],D,aperture,c0,ell);
        }

        Ellipse e_2;
        BVec b_2(order,sigma_i);
        success = e_2.measure(allpix,allpsf2,order,sigma_i,true,flag,0,&b_2);
        xdbg<<"sigma_i = "<<sigma_i<<std::endl;
        xdbg<<"e_2 = "<<e_2<<std::endl;
        xdbg<<"b_2 = "<<b_2.vec()<<std::endl;
        xdbg<<"flag = "<<flag<<std::endl;
        Test(success,"Measure for Multi Image w/ PSF - success");
        Test(std::abs(e_2.getCen()-ell.getCen()) < 1.e-5,
             "Measure for Multi Image w/ PSF - cen");
        Test(std::abs(e_2.getGamma()-ell.getGamma()) < 1.e-5,
             "Measure for Multi Image w/ PSF - gam");
        Test(std::abs(e_2.getMu()-ell.getMu()) < 1.e-5,
             "Measure for Multi Image w/ PSF - mu");
        Test((b_2.vec().TMV_subVector(0,15)-bin.vec()).TMV_normInf() < 1.e-5*bin.vec().norm(),
             "Measure for Multi Image w/ PSF - b");
        Test(b_2.vec().TMV_subVector(15,b_2.size()).TMV_normInf() < 1.e-5*bin.vec().norm(),
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
        Test(success,"Measure for Multi Image w/ PSF #2 - success");
        Test(std::abs(e_3.getCen()-ell.getCen()) < 1.e-5,
             "Measure for Multi Image #2 w/ PSF - cen");
        Test(std::abs(e_3.getGamma()-ell.getGamma()) < 1.e-5,
             "Measure for Multi Image #2 w/ PSF - gam");
        Test(std::abs(e_3.getMu()-ell.getMu()) < 1.e-5,
             "Measure for Multi Image #2 w/ PSF - mu");
        Test((b_3.vec().TMV_subVector(0,15)-bin.vec()).TMV_normInf() < 1.e-5*bin.vec().norm(),
             "Measure for Multi Image #2 w/ PSF - b");
        Test(b_3.vec().TMV_subVector(15,b_3.size()).TMV_normInf() < 1.e-5*bin.vec().norm(),
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
    incat.read();

    // Find stars
    StarCatalog starcat(incat,params,"");
    //starcat.calcSizes(im,weight_im.get(),trans);
    FindStarsLog fslog(params,"testfs.log");
    int nstars = starcat.findStars(fslog);
    dbg<<fslog<<std::endl;
    Test(nstars >= 240,"Find Stars");

    // Test I/O
    starcat.write();
    StarCatalog starcat2(params,"");
    starcat2.read();
    std::vector<std::string> all_ext = params["stars_ext"];
    for(int k=0;k<int(all_ext.size());++k) {
        dbg<<"Test I/O with extension "<<all_ext[k]<<std::endl;
        params["stars_ext"] = all_ext[k];
        if (k > 0) starcat2.read();
        Test(starcat.size() == starcat2.size(),"starcat size I/O");
        for(int i=0;i<starcat.size();++i) {
            Test(starcat.getId(i) == starcat2.getId(i), "starcat id I/O");
            Test(std::abs(starcat.getPos(i) - starcat2.getPos(i)) <= 0.01,
                 "starcat pos I/O");
            Test(std::abs(starcat.getSky(i) - starcat2.getSky(i)) <= 0.01,
                 "starcat sky I/O");
            Test(std::abs(starcat.getNoise(i) - starcat2.getNoise(i)) <= 0.01,
                 "starcat noise I/O");
            Test(starcat.getFlags(i) == starcat2.getFlags(i),
                 "starcat flags I/O");
            dbg<<"starcat.mag[i] = "<<starcat.getMag(i)<<
                ", starcat2.mag[i] = "<<starcat2.getMag(i)<<std::endl;
            dbg<<"abs(diff) = "<<
                std::abs(starcat.getMag(i) - starcat2.getMag(i))<<std::endl;
            Test(std::abs(starcat.getMag(i) - starcat2.getMag(i)) <= 0.01,
                 "starcat mag I/O");
            Test(std::abs(starcat.getObjSize(i)-starcat2.getObjSize(i)) <= 0.01,
                 "starcat objsize I/O");
            Test(starcat.getIsAStar(i) == starcat2.getIsAStar(i),
                 "starcat isastar I/O");
        }
    }
    params["stars_ext"] = all_ext;

    // Measure PSF
    PsfCatalog psfcat(starcat,params);
    double sigma_p = psfcat.estimateSigma(im,weight_im.get(),trans);
    PsfLog psflog(params,"testpsf.log");
    int npsf = psfcat.measurePsf(im,weight_im.get(),trans,sigma_p,psflog);
    dbg<<psflog<<std::endl;
    // There are 244 stars in the file, but depending on the exact parameters,
    // one or two cross the edge, which gives an error flag.
    // The rest should all be measured successfully.
    Test(npsf >= 240,"Measure PSF");

    // Test I/O
    psfcat.write();
    PsfCatalog psfcat2(params);
    psfcat2.read();
    all_ext = params["psf_ext"];
    for(int k=0;k<int(all_ext.size());++k) {
        dbg<<"Test I/O with extension "<<all_ext[k]<<std::endl;
        params["psf_ext"] = all_ext[k];
        if (k > 0) psfcat2.read();
        Test(psfcat.size() == psfcat2.size(),"psfcat size I/O");
        for(int i=0;i<psfcat.size();++i) {
            Test(psfcat.getId(i) == psfcat2.getId(i),"psfcat id I/O");
            Test(std::abs(psfcat.getPos(i) - psfcat2.getPos(i)) <= 0.01,
                 "psfcat pos I/O");
            Test(std::abs(psfcat.getSky(i) - psfcat2.getSky(i)) <= 0.01,
                 "psfcat sky I/O");
            Test(std::abs(psfcat.getNoise(i) - psfcat2.getNoise(i)) <= 0.01,
                 "psfcat noise I/O");
            Test(psfcat.getFlags(i) == psfcat2.getFlags(i),"psfcat flags I/O");
            Test(std::abs(psfcat.getNu(i) - psfcat2.getNu(i)) <= 0.01,
                 "psfcat nu I/O");
            Test(psfcat.getPsf(i).size() == psfcat2.getPsf(i).size(),
                 "psfcat psf.size I/O");
            Test(std::abs(psfcat.getPsf(i).getSigma() -
                          psfcat2.getPsf(i).getSigma()) 
                 <= 0.01, "psfcat psf.getSigma I/O");
            Test((psfcat.getPsf(i).vec() - psfcat2.getPsf(i).vec()).TMV_normInf() <= 
                 1.e-4*psfcat.getPsf(i).vec().norm(), "psfcat psf vector I/O");
        }
    }
    params["psf_ext"] = all_ext;

    // Fit PSF
    FittedPsf fitpsf(psfcat,params,psflog);
    double rms = 0.; 
    int count = 0;
    for(int i=0;i<nstars;++i) if (!psfcat.getFlags(i)) {
        BVec checkpsf(fitpsf.getPsfOrder(),fitpsf.getSigma());
        checkpsf = fitpsf(starcat.getPos(i));
        double normsqdiff = (psfcat.getPsf(i).vec()-checkpsf.vec()).TMV_normSq();
        rms += normsqdiff;
        ++count;
    }
    rms /= count;
    rms = sqrt(rms);
    dbg<<"fitpsf rms = "<<rms<<std::endl;
    Test(rms < double(fitpsf.getPsfSize())/double(fitpsf.getFitSize()),
         "Fit PSF rms");

    // Test I/O
    fitpsf.write();
    FittedPsf fitpsf2(params);
    fitpsf2.read();
    all_ext = params["fitpsf_ext"];
    for(int k=0;k<int(all_ext.size());++k) {
        dbg<<"Test I/O with extension "<<all_ext[k]<<std::endl;
        params["fitpsf_ext"] = all_ext[k];
        if (k > 0) fitpsf2.read();
        rms = 0.; 
        count = 0;
        Test(fitpsf2.getPsfOrder() == fitpsf.getPsfOrder(),
             "FittedPSF I/O: order");
        Test(std::abs(fitpsf2.getSigma() - fitpsf.getSigma()) < 0.01, 
             "FittedPSF I/O: sigma");
        for(int i=0;i<nstars;++i) if (!psfcat.getFlags(i)) {
            BVec checkpsf(fitpsf.getPsfOrder(),fitpsf.getSigma());
            checkpsf = fitpsf(starcat.getPos(i));
            BVec checkpsf2(fitpsf2.getPsfOrder(),fitpsf2.getSigma());
            checkpsf2 = fitpsf2(starcat.getPos(i));
            double normsqdiff = (checkpsf2.vec()-checkpsf.vec()).TMV_normSq();
            rms += normsqdiff;
            ++count;
        }
        rms /= count;
        rms = sqrt(rms);
        dbg<<"fitpsf I/O rms = "<<rms<<std::endl;
        Test(rms < 1.e-4,"Fit PSF I/O");
    }
    params["fitpsf_ext"] = all_ext;

    // Measure shears
    ShearCatalog shearcat(incat,trans,fitpsf,params);
    ShearLog shearlog(params,"testshear.log");
    int nshear = shearcat.measureShears(im,weight_im.get(),shearlog);
    dbg<<shearlog<<std::endl;
    // There are 4557 galaxies in the file without error codes.
    // The code currently converges on more than 2800 of them,
    // although only about 1200 or so have no error flag.
    Test(nshear >= 2800,"Measure Shear");

    // Test I/O
    shearcat.write();
    ShearCatalog shearcat2(params);
    shearcat2.read();
    all_ext = params["shear_ext"];
    for(int k=0;k<int(all_ext.size());++k) {
        dbg<<"Test I/O with extension "<<all_ext[k]<<std::endl;
        params["shear_ext"] = all_ext[k];
        if (k > 0) shearcat2.read();
        Test(shearcat.size() == shearcat2.size(),"shearcat size I/O");
        for(int i=0;i<shearcat.size();++i) {
            Test(shearcat.getId(i) == shearcat2.getId(i),"shearcat id I/O");
            Test(std::abs(shearcat.getPos(i) - shearcat2.getPos(i)) <= 0.01,
                 "shearcat pos I/O");
            Test(std::abs(shearcat.getSky(i) - shearcat2.getSky(i)) <= 0.01,
                 "shearcat sky I/O");
            Test(std::abs(shearcat.getNoise(i) - shearcat2.getNoise(i)) <= 0.01,
                 "shearcat noise I/O");
            Test(shearcat.getFlags(i) == shearcat2.getFlags(i),
                 "shearcat flags I/O");
            dbg<<"skypos = "<<shearcat.getSkyPos(i)<<"  "<<shearcat2.getSkyPos(i)<<std::endl;
            Test(std::abs(shearcat.getSkyPos(i) - shearcat2.getSkyPos(i)) 
                 <= 0.01,
                 "shearcat skypos I/O");
            Test(std::abs(shearcat.getShear(i) - shearcat2.getShear(i)) <= 0.01,
                 "shearcat shear I/O");
            Test(std::abs(shearcat.getNu(i) - shearcat2.getNu(i)) <= 0.01,
                 "shearcat nu I/O");
#ifdef USE_TMV
            Test((DMatrix((shearcat.getCov(i) - shearcat2.getCov(i)))).norm()
                 <= 1.e-4*shearcat.getCov(i).norm(),
                 "shearcat cov I/O");
#else
            Test((shearcat.getCov(i) - shearcat2.getCov(i)).TMV_normInf()
                 <= 1.e-4*shearcat.getCov(i).norm(),
                 "shearcat cov I/O");
#endif
            Test(shearcat.getShape(i).size() == shearcat2.getShape(i).size(),
                 "shearcat shape.size I/O");
            Test(std::abs(shearcat.getShape(i).getSigma() -
                          shearcat2.getShape(i).getSigma()) <= 0.01,
                 "shearcat shape.getSigma I/O");
            Test((shearcat.getShape(i).vec() - shearcat2.getShape(i).vec()).TMV_normInf() 
                 <= 1.e-4*shearcat.getShape(i).vec().norm(),
                 "shearcat shape vector I/O");
        }
    }
    params["shear_ext"] = all_ext;

    std::cout<<"Passed full pipeline\n";
#endif

#ifdef TEST7
    {
        int order = 8;
        int orderSize = (order+1)*(order+2)/2;
        double sigma = 5.;

        // First lets make sure we understand Gary's LVector class:
        laguerre::LVector garyB(order);
        garyB.set(0,0,1.);
        garyB.set(1,0,std::complex<double>(0.4,0.7));
        garyB.set(2,0,std::complex<double>(0.1,0.2));
        garyB.set(1,1,0.3);

        BVec myB(order,sigma);
        myB(0) = 1.;
        myB(1) = 0.4;
        myB(2) = 0.7;
        myB(3) = 0.1;
        myB(4) = 0.2;
        myB(5) = 0.3;
        Test((myB.vec() - convertVector(garyB.rVector())).norm() <= 1.e-5,
             "compare my BVec with Gary's LVector");

        // Now check I(x,y)
        double x = -0.23;
        double y = 0.51;

        mv::DVector garyPsi = garyB.realPsi(order,x,y,laguerre::Ellipse());
        double garyIxy = garyPsi * garyB.rVector();
        // Gary normalizes his shapelets by 1/2Pi rather than 1/sqrt(Pi).
        garyIxy *= 2.*sqrtpi;

        DVector myPsi(orderSize);
        MakePsi(myPsi,std::complex<double>(x,y),order);
        double myIxy = EIGEN_ToScalar(EIGEN_Transpose(myPsi) * myB.vec());
        Test(std::abs(myIxy - garyIxy) <= 1.e-5,
             "compare MakePsi with Gary's realPsi");

        // Compare Gary's MakeLTransform to my CalculateZTransform
        std::complex<double> z(0.71,0.32);
        laguerre::Position<double> garyZ(std::real(z),std::imag(z));
        laguerre::LTransform garyZTransform = 
            laguerre::MakeLTransform(garyZ,order,order,true);
        mv::DMatrix garyZMatrix = garyZTransform.rMatrix();
        DMatrix myZMatrix(orderSize,orderSize);
        myZMatrix.setZero();
        CalculateZTransform(z,order,myZMatrix);
        Test((myZMatrix - convertMatrix(garyZMatrix)).norm() <= 1.e-5,
             "compare ZTransform with Gary's MakeLTransform");

        // Compare Gary's MakeLTransform to my CalculateMuTransform
        double mu = 0.37;
        laguerre::LTransform garyMuTransform = 
            laguerre::MakeLTransform(mu,order,order,true);
        mv::DMatrix garyMuMatrix = garyMuTransform.rMatrix();
        garyMuMatrix /= exp(2.*mu);
        DMatrix myMuMatrix(orderSize,orderSize);
        myMuMatrix.setZero();
        CalculateMuTransform(mu,order,myMuMatrix);
        Test((myMuMatrix - convertMatrix(garyMuMatrix)).norm() <= 1.e-5,
             "compare MuTransform with Gary's MakeLTransform");

        // Compare Gary's MakeLTransform to my CalculateGTransform
        std::complex<double> g(0.22,0.53);
        // Gary's constructor takes distortions, not shears, so convert.
        double gamma = std::abs(g);
        double delta = tanh(2.*atanh(gamma));
        double scale = delta/gamma;
        laguerre::Shear garyG(std::real(g)*scale,std::imag(g)*scale);
        laguerre::LTransform garyGTransform = 
            laguerre::MakeLTransform(garyG,order,order,true);
        mv::DMatrix garyGMatrix = garyGTransform.rMatrix();
        DMatrix myGMatrix(orderSize,orderSize);
        myGMatrix.setZero();
        CalculateGTransform(g,order,myGMatrix);
        Test((myGMatrix - convertMatrix(garyGMatrix)).norm() <= 1.e-5,
             "compare GTransform with Gary's MakeLTransform");

        DVector b_exp(45);
        // This is the b vector for an exponential disk with scale size
        // r0 = 5/1.67839, sheared by eta = arctanh(0.2), and measured
        // with sigma_GAL = 5:
        b_exp << 
            17.43818299 ,
            0, 0, 
            0.9998995137, 0, -4.241349805, 
            0, 0, 0, 0,
            0.07587021059, 0, -0.2836021410, 0, 3.483311775,
            0, 0, 0, 0, 0, 0, 
            0.006287543586, 0, -0.02080316557, 0, 0.3276615409, 0, -1.963729233,
            0, 0, 0, 0, 0, 0, 0, 0,
            0.0005450606239, 0, -0.001621332707, 0, 0.02997036976, 0, -0.1823596723, 0, 1.607970421;

        // Compare Gary's MakeLTransform with my CalculatePsfConvolve
        // First do PSF's with b00 = 1 and a single element = (0.1,0.2) or 0.1.
        for(int pplusq = 1, k=1; pplusq <= order; ++pplusq) {
            for(int p=pplusq, q=0; p>=q; --p,++q,++k) {
                dbg<<"p,q,k = "<<p<<','<<q<<','<<k<<std::endl;
                laguerre::LVector garyBPsf(order);
                garyBPsf.set(0,0,1.);
                BVec myBPsf(order,sigma_psf);
                myBPsf(0) = 1.;
                if (p == q) {
                    garyBPsf.set(p,q,std::complex<double>(0.1,0.0));
                    myBPsf(k) = 0.1;
                } else {
                    garyBPsf.set(p,q,std::complex<double>(0.1,0.2));
                    myBPsf(k++) = 0.1;
                    myBPsf(k) = 0.2;
                }
                double D = 1. / (1.+pow(sigma_psf/sigma,2));
                laguerre::LTransform garyPsfTransform = 
                    laguerre::MakeLTransform(garyBPsf,D,order,order,order);
                mv::DMatrix garyPsfMatrix = garyPsfTransform.rMatrix();
                DMatrix myPsfMatrix(orderSize,orderSize);
                myPsfMatrix.setZero();
                CalculatePsfConvolve(myBPsf,order,sigma,myPsfMatrix);
                garyPsfMatrix /= sqrtpi; // To Match my normalization
                dbg<<"Gary's convolution matrix/sqrt(pi) = "<<
                    convertMatrix(garyPsfMatrix)<<std::endl;
                dbg<<"My convolution matrix = "<<
                    myPsfMatrix<<std::endl;
                double normdiff = (myPsfMatrix-convertMatrix(garyPsfMatrix)).norm();
                dbg<<"Norm(diff) = "<<normdiff<<std::endl;
                Test(normdiff <= 1.e-5,
                     "compare PsfConvolve with Gary's MakeLTransform");
                dbg<<"C*b_exp = "<<myPsfMatrix * b_exp<<std::endl;
            }
        }
        // Next do some PSF's with lots of values != 0:
        // Also has mismatched orders (psf = 4, gal = order)
        for(int ip = 0; ip < NPSF; ++ip) {
            dbg<<"Start ip = "<<ip<<std::endl;
            BVec myBPsf(4,sigma_psf,bpsf_vecs[ip]);
            dbg<<"bpsf = "<<myBPsf.vec()<<std::endl;
            laguerre::LVector garyBPsf(4);
            for(int pplusq = 0, k=0; pplusq <= 4; ++pplusq) {
                for(int p=pplusq, q=0; p>=q; --p,++q,++k) {
                    if (p == q) {
                        garyBPsf.set(p,q,std::complex<double>(
                                bpsf_vecs[ip][k],0.));
                    } else {
                        garyBPsf.set(p,q,std::complex<double>(
                                bpsf_vecs[ip][k],bpsf_vecs[ip][k+1]));
                        ++k;
                    }
                }
            }
            dbg<<"Gary's bpsf = "<<
                convertVector(garyBPsf.rVector())<<std::endl;
            double D = 1. / (1.+pow(sigma_psf/sigma,2));
            dbg<<"D = "<<D<<std::endl;
            laguerre::LTransform garyPsfTransform = 
                laguerre::MakeLTransform(garyBPsf,D,order,order,4);
            mv::DMatrix garyPsfMatrix = garyPsfTransform.rMatrix();
            DMatrix myPsfMatrix(orderSize,orderSize);
            myPsfMatrix.setZero();
            CalculatePsfConvolve(myBPsf,order,sigma,myPsfMatrix);
            garyPsfMatrix /= sqrtpi; // To Match my normalization
            dbg<<"Gary's convolution matrix/sqrt(pi) = "<<
                convertMatrix(garyPsfMatrix)<<std::endl;
            dbg<<"My convolution matrix/2sqrt(pi) = "<<
                myPsfMatrix<<std::endl;
            double normdiff = (myPsfMatrix-convertMatrix(garyPsfMatrix)).norm();
            dbg<<"Norm(diff) = "<<normdiff<<std::endl;
            Test(normdiff <= 1.e-5,
                 "compare PsfConvolve with Gary's MakeLTransform");
        }
    }
    std::cout<<"Passed tests against Gary's code.\n";
#endif

    if (dbgout && dbgout != &std::cout) {delete dbgout; dbgout=0;}
    return 0;
}
#if 0
// Change to 1 to let gdb see where the program bombed out.
catch(int) {}
#else
#ifdef USE_TMV
catch (tmv::Error& e) {
    std::cerr<<"Caught "<<e<<std::endl;
    return 1;
}
#endif
catch (std::exception& e) {
    std::cerr<<"Caught std::exception: "<<e.what()<<std::endl;
}
catch (...) {
    std::cerr<<"Caught Unknown error\n";
}
#endif
