
#include "Ellipse.h"
#include "EllipseSolver.h"
#include "TMV.h"
#include <cmath>
#include "dbg.h"
#include <fstream>
#include "TimeVars.h"
#include "PsiHelper.h"
#include "Params.h"

#define N_FLUX_ATTEMPTS 0
#define MAXITER 4

static std::complex<double> addShears(
    const std::complex<double> g1, const std::complex<double> g2)
{
    double absg1 = std::abs(g1);
    if (absg1 == 0.) return g2;
    double absd1 = tanh(atanh(absg1)*2.);
    std::complex<double> d1 = absd1 / absg1 * g1;
    double absg2 = std::abs(g2);
    if (absg2 == 0.) return g1;
    double absd2 = tanh(atanh(absg2)*2.);
    std::complex<double> d2 = absd2 / absg2 * g2;

    double d2sq = absd2*absd2;
    double x = (1.-std::sqrt(1.-d2sq))/d2sq;
    std::complex<double> d3 = d1+d2*(1.+x*(std::real(d1)*d2-std::real(d2)*d1));
    d3 /= 1. + std::real(d1*std::conj(d2));
    double absd3 = std::abs(d3);
    double absg3 = tanh(atanh(absd3)/2.);
    std::complex<double> g3 = absg3/absd3*d3;
    xdbg<<"GammaAdd: "<<g1<<" + "<<g2<<" = "<<g3<<std::endl;
    return g3;
}

bool Ellipse::doMeasure(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>* psf,
    int order, double sigma, bool shouldUseInteg, long& flag, 
    tmv::Matrix<double>* cov, BVec* bRet, tmv::Matrix<double>* bCov)
{
    timeval tp;
    double t1=0.,t2=0.;

    // The non-linear solver is pretty sensitive to having a good
    // initial estimate.  So start with a simple estimate.
    if (order > 3) {
        doMeasure(pix,psf,order-2,sigma,shouldUseInteg,flag);
    }

    xdbg<<"Start DoMeasure: order = "<<order<<", psf = "<<bool(psf)<<std::endl;
    xdbg<<"fix = "<<_isFixedCen<<"  "<<_isFixedGamma<<"  "<<_isFixedMu<<std::endl;
    for(size_t i=0;i<pix.size();++i) xdbg<<"npix["<<i<<"] = "<<pix[i].size()<<std::endl;

    std::auto_ptr<BaseEllipseSolver> solver;

    tmv::Vector<double> x(5,0.);
    x(0) = std::real(_cen);  x(1) = std::imag(_cen); 
    x(2) = std::real(_gamma); x(3) = std::imag(_gamma);
    x(4) = _mu;
    tmv::Vector<double> xinit = x;
    tmv::Vector<double> f(5);

    if (shouldUseInteg && order <= 3) {
        if (_shouldDoTimings) {
            gettimeofday(&tp,0);
            t1 = tp.tv_sec + tp.tv_usec/1.e6;
        }
        // First we make the approximations the the galaxy is well-sampled
        // and has uniform variance.
        // Also, we start by just fitting the centroid.
        if (!_isFixedCen && (!_isFixedGamma || !_isFixedMu)) {
            if (psf) {
                // pixscale doesn't really matter unless we want an accurate B in
                // the end, so just use 1 here.
                solver.reset(new EllipseSolver2(pix,*psf,_fPsf,order,sigma,1.,
                                                false,true,true));
            } else {
                solver.reset(new EllipseSolver2(pix,order,sigma,1.,
                                                false,true,true));
            }
            xdbg<<"xinit (int fixgam) = "<<x<<std::endl;
            solver->useDogleg();
#ifdef NOTHROW
            solver->noUseCholesky();
#endif
            solver->setTol(1.e-3,1.e-8);
            if (XDEBUG) solver->setOutput(*dbgout);
            solver->setDelta0(0.01);
            solver->setMinStep(1.e-12);
            solver->setMaxIter(60);
            xdbg<<"Integrating solver, centroid only:\n";
            //if (XDEBUG) if (!solver->testJ(x,f,dbgout,1.e-5)) exit(1);
            bool ret = solver->solve(x,f);
            if (!ret) {
                dbg<<"failed integrating solver, centroid - x = "<<x<<std::endl;
                dbg<<"f = "<<f<<"  Norm(f) = "<<Norm(f)<<std::endl;
                dbg<<"b = "<<solver->getB().vec()<<std::endl;
                return false;
            }
            xdbg<<"Done: x (Integ, centroid only) = "<<x<<std::endl;
            xdbg<<"f = "<<f<<std::endl;
            xdbg<<"b = "<<solver->getB().vec()<<std::endl;
        }

        // Next allow the shear and/or mu to be fit as well:
#ifdef N_FLUX_ATTEMPTS
#if N_FLUX_ATTEMPTS > 0
        if (psf) {
            // pixscale doesn't really matter unless we want an accurate B in
            // the end, so just use 1 here.
            solver.reset(
                new EllipseSolver2(pix,*psf,_fPsf,order,sigma,1.,
                                   _isFixedCen,_isFixedGamma,_isFixedMu,true));
        } else {
            solver.reset(new EllipseSolver2(pix,order,sigma,1.,
                                            _isFixedCen,_isFixedGamma,_isFixedMu,true));
        }
        xdbg<<"xinit (integ) = "<<x<<std::endl;
        solver->useHybrid();
#ifdef NOTHROW
        solver->noUseCholesky();
        solver->noUseDirectH();
#endif
        solver->setTol(3.e-3,1.e-5);
        solver->setTau(1.0);
        if (XDEBUG) solver->setOutput(*dbgout);
        solver->setMinStep(1.e-12);
        solver->setMaxIter(10);
        xdbg<<"Integrating solver:\n";
        for(int iter = 0; iter < N_FLUX_ATTEMPTS; ++iter) {
            xdbg<<"Attempt #"<<iter<<std::endl;
            //if (XDEBUG) if (!solver->testJ(x,f,dbgout,1.e-5)) exit(1);
            solver->solve(x,f);
            xdbg<<"x => "<<x<<std::endl;
            xdbg<<"f => "<<f<<"  Norm(f) = "<<Norm(f)<<std::endl;
            xdbg<<"b => "<<solver->getB().vec()<<std::endl;
            if (Norm(f) < 1.e-2) break;
        } 
#endif
#endif
        // Repeat, but allow flux to change.
        if (psf) {
            solver.reset(
                new EllipseSolver2(pix,*psf,_fPsf,order,sigma,1.,
                                   _isFixedCen,_isFixedGamma,_isFixedMu));
        } else {
            solver.reset(new EllipseSolver2(pix,order,sigma,1.,
                                            _isFixedCen,_isFixedGamma,_isFixedMu));
        }
        xdbg<<"xinit (integ) = "<<x<<std::endl;
        solver->useDogleg();
#ifdef NOTHROW
        solver->noUseCholesky();
#endif
        solver->setTol(1.e-3,1.e-8);
        if (XDEBUG) solver->setOutput(*dbgout);
        solver->setDelta0(0.01);
        solver->setMinStep(1.e-12);
        solver->setMaxIter(10);
        xdbg<<"Final Integrating solver:\n";
        //if (XDEBUG) if (!solver->testJ(x,f,dbgout,1.e-5)) exit(1);
        //solver->verbose = true;
        solver->solve(x,f);
        xdbg<<"Done: x (Integrating) = "<<x<<std::endl;
        xdbg<<"f (integ) = "<<f<<std::endl;
        xdbg<<"b = "<<solver->getB().vec()<<std::endl;
        if (_shouldDoTimings) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            _times._tInteg += t2-t1;
        }
    }

#if 0
    // We should be close enough for the exact solver to work.
    // But just to be safe, fit each of centroid, gamma, mu separately first:
    if (_shouldDoTimings) {
        gettimeofday(&tp,0);
        t1 = tp.tv_sec + tp.tv_usec/1.e6;
    }
    if (!_isFixedCen) {
        if (psf) {
            solver.reset(new EllipseSolver(pix,*psf,_fPsf,order,sigma,
                                           false,true,true));
        } else {
            solver.reset(new EllipseSolver(pix,order,sigma,
                                           false,true,true));
        }
        xdbg<<"xinit = "<<x<<std::endl;
        solver->useDogleg();
#ifdef NOTHROW
        solver->noUseCholesky();
#endif
        solver->setTol(1.e-3,1.e-8);
        if (XDEBUG) solver->setOutput(*dbgout);
        solver->setDelta0(0.01);
        solver->setMinStep(1.e-12);
        solver->setMaxIter(50);
        xdbg<<"ML solver centroid:\n";
        bool rete= solver->solve(x,f);
        if (!ret) {
            dbg<<"ML solver, centroid, failed - x = "<<x<<std::endl;
            dbg<<"f = "<<f<<std::endl;
            dbg<<"b = "<<solver->getB().vec()<<std::endl;
            if (XDEBUG) solver->testJ(x,f,dbgout,1.e-5);
            if (_shouldDoTimings) {
                gettimeofday(&tp,0);
                t2 = tp.tv_sec + tp.tv_usec/1.e6;
                _times._tCentroid += t2-t1;
            }
            return false;
        }
        xdbg<<"Done: x (ML centroid only) = "<<x<<std::endl;
        xdbg<<"f = "<<f<<std::endl;
        xdbg<<"b = "<<solver->getB().vec()<<std::endl;
    }
    if (_shouldDoTimings) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        _times._tCentroid += t2-t1;
    }

    if (_shouldDoTimings) {
        gettimeofday(&tp,0);
        t1 = tp.tv_sec + tp.tv_usec/1.e6;
    }
    if (!_isFixedGamma) {
        // start with an initial guess from the latest b vector:
        if (solver.get()) {
            BVec b1 = solver->getB();
            std::complex<double> g1(b1(3)/b1(0)*sqrt(2.),-b1(4)/b1(0)*sqrt(2.));
            dbg<<"g1 = "<<g1<<std::endl;
            if (std::abs(g1) < 1.) {
                std::complex<double> g0(x[2],x[3]);
                std::complex<double> g = addShears(g0,g1);
                dbg<<"Initial gamma = "<<g<<std::endl;
                x[2] = std::real(g);
                x[3] = std::imag(g);
            }
        }

        if (psf) {
            solver.reset(new EllipseSolver(pix,*psf,_fPsf,order,sigma,
                                           true,false,true));
        } else {
            solver.reset(new EllipseSolver(pix,order,sigma,
                                           true,false,true));
        }
        xdbg<<"xinit = "<<x<<std::endl;
        solver->useDogleg();
#ifdef NOTHROW
        solver->noUseCholesky();
#endif
        solver->setTol(1.e-3,1.e-8);
        if (XDEBUG) solver->setOutput(*dbgout);
        solver->setDelta0(0.01);
        solver->setMinStep(1.e-12);
        solver->setMaxIter(50);
        xdbg<<"ML solver gamma:\n";
        bool ret = solver->solve(x,f);
        if (!ret) {
            dbg<<"ML solver, gamma, failed - x = "<<x<<std::endl;
            dbg<<"f = "<<f<<std::endl;
            dbg<<"b = "<<solver->getB().vec()<<std::endl;
            //if (XDEBUG) solver->testJ(x,f,dbgout,1.e-5);
            if (_shouldDoTimings) {
                gettimeofday(&tp,0);
                t2 = tp.tv_sec + tp.tv_usec/1.e6;
                _times._tGamma += t2-t1;
            }
            return false;
        }
        xdbg<<"Done: x (ML gamma only) = "<<x<<std::endl;
        xdbg<<"f = "<<f<<std::endl;
        xdbg<<"b = "<<solver->getB().vec()<<std::endl;
    }
    if (_shouldDoTimings) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        _times._tGamma += t2-t1;
    }

    if (_shouldDoTimings) {
        gettimeofday(&tp,0);
        t1 = tp.tv_sec + tp.tv_usec/1.e6;
    }
    if (!_isFixedMu) {
        if (psf) {
            solver.reset(new EllipseSolver(pix,*psf,_fPsf,order,sigma,
                                           true,true,false));
        } else {
            solver.reset(new EllipseSolver(pix,order,sigma,
                                           true,true,false));
        }
        xdbg<<"xinit = "<<x<<std::endl;
        solver->useDogleg();
#ifdef NOTHROW
        solver->noUseCholesky();
#endif
        solver->setTol(1.e-3,1.e-8);
        if (XDEBUG) solver->setOutput(*dbgout);
        solver->setDelta0(0.01);
        solver->setMinStep(1.e-12);
        solver->setMaxIter(50);
        xdbg<<"ML solver mu:\n";
        bool ret = solver->solve(x,f);
        if (!ret) {
            dbg<<"ML solver, mu, failed - x = "<<x<<std::endl;
            dbg<<"f = "<<f<<std::endl;
            dbg<<"b = "<<solver->getB().vec()<<std::endl;
            //if (XDEBUG) solver->testJ(x,f,dbgout,1.e-5);
            if (_shouldDoTimings) {
                gettimeofday(&tp,0);
                t2 = tp.tv_sec + tp.tv_usec/1.e6;
                _times._tGamma += t2-t1;
            }
            return false;
        }
        xdbg<<"Done: x (ML mu only) = "<<x<<std::endl;
        xdbg<<"f = "<<f<<std::endl;
        xdbg<<"b = "<<solver->getB().vec()<<std::endl;
    }
    if (_shouldDoTimings) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        _times._tGamma += t2-t1;
    }
#endif

    // Now fit everything, but first time, try to maintain the flux level
#ifdef N_FLUX_ATTEMPTS
#if N_FLUX_ATTEMPTS > 0
    if (_shouldDoTimings) {
        gettimeofday(&tp,0);
        t1 = tp.tv_sec + tp.tv_usec/1.e6;
    }
    if (psf) {
        solver.reset(new EllipseSolver(pix,*psf,_fPsf,order,sigma,
                                       _isFixedCen,_isFixedGamma,_isFixedMu,true));
    } else {
        solver.reset(new EllipseSolver(pix,order,sigma,
                                       _isFixedCen,_isFixedGamma,_isFixedMu,true));
    }
    xdbg<<"xinit = "<<x<<std::endl;
    solver->useHybrid();
#ifdef NOTHROW
    solver->noUseCholesky();
    solver->noUseDirectH();
#endif
    solver->setTol(3.e-3,1.e-5);
    solver->setTau(1.0);
    if (XDEBUG) solver->setOutput(*dbgout);
    solver->setMinStep(1.e-12);
    solver->setMaxIter(10);
    solver->setDelta0(0.01);
    xdbg<<"ML solver use flux:\n";
    for(int iter = 0; iter < N_FLUX_ATTEMPTS; ++iter) {
        xdbg<<"Attempt #"<<iter<<std::endl;
        //if (XDEBUG) if (!solver->testJ(x,f,dbgout,1.e-5)) exit(1);
        solver->solve(x,f);
        xdbg<<"x => "<<x<<std::endl;
        xdbg<<"f => "<<f<<"  Norm(f) = "<<Norm(f)<<std::endl;
        xdbg<<"b = "<<solver->getB().vec()<<std::endl;
        if (Norm(f) < 1.e-2) break;
    }
    xdbg<<"Done: x (ML fixed flux) = "<<x<<std::endl;
    xdbg<<"f = "<<f<<std::endl;
    xdbg<<"b = "<<solver->getB().vec()<<std::endl;
    if (_shouldDoTimings) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        _times._tFixFlux += t2-t1;
    }
#endif
#endif

    if (_shouldDoTimings) {
        gettimeofday(&tp,0);
        t1 = tp.tv_sec + tp.tv_usec/1.e6;
    }
    // Finally allow the flux change as needed.
    if (psf) {
        solver.reset(new EllipseSolver(pix,*psf,_fPsf,order,sigma,
                                       _isFixedCen,_isFixedGamma,_isFixedMu));
    } else {
        solver.reset(new EllipseSolver(pix,order,sigma,
                                       _isFixedCen,_isFixedGamma,_isFixedMu));
    }
    xdbg<<"xinit = "<<x<<std::endl;
    solver->useDogleg();
#ifdef NOTHROW
    solver->noUseCholesky();
#endif
    solver->setTol(1.e-8,1.e-25);
    if (XDEBUG) solver->setOutput(*dbgout);
    solver->setDelta0(0.01);
    solver->setMinStep(1.e-15);
    solver->setMaxIter(50);
    xdbg<<"Final solver:\n";
    for(int iter = 1;iter<=MAXITER;++iter) {
        //if (XDEBUG) if (!solver->testJ(x,f,dbgout,1.e-6)) exit(1);
        bool ret = solver->solve(x,f);
        if (ret) break;
        else if (iter == MAXITER) {
            if (_shouldDoTimings) {
                gettimeofday(&tp,0);
                t2 = tp.tv_sec + tp.tv_usec/1.e6;
                _times._tFinal += t2-t1;
            }
            return false;
        }
        dbg<<"Final solver pass failed - x = "<<x<<std::endl;
        dbg<<"b = "<<solver->getB().vec()<<std::endl;
        //if (XDEBUG) if (!solver->testJ(x,f,dbgout,1.e-5)) exit(1);
#if 1 
        // Try bumping it to get out of a local well and continue:
        BVec b1 = solver->getB();
        if (!_isFixedCen) {
            dbg<<"Current z = "<<std::complex<double>(x[0],x[1])<<std::endl;
            std::complex<double> z1(2.*b1(1),-2.*b1(2));
            z1 /= b1(0) - b1(5);
            if (std::abs(z1) < 1.) {
                dbg<<"Add "<<z1<<std::endl;
                x[0] += std::real(z1); x[1] += std::imag(z1); 
                dbg<<"New z = "<<std::complex<double>(x[0],x[1])<<std::endl;
            }
        }
#endif
#if 1
        if (!_isFixedGamma) {
            std::complex<double> g0(x[2],x[3]);
            dbg<<"Current gamma = "<<g0<<std::endl;
            std::complex<double> g1(std::sqrt(2.)*b1(3),-std::sqrt(2.)*b1(4));
            g1 /= b1(0) - (b1.getOrder() >= 4 ? b1(14) : 0.);
            if (std::abs(g1) < 1.) {
                dbg<<"Add "<<g1<<std::endl;
                g0 = addShears(g1,g0);
                x[2] = std::real(g0); x[3] = std::imag(g0); 
                dbg<<"New gamma = "<<std::complex<double>(x[2],x[3])<<std::endl;
            }
        }
#endif
#if 0
        if (!_isFixedMu) {
            dbg<<"Current mu = "<<x[4]<<std::endl;
            double m1 = -b1(5);
            m1 /= b1[0] - (b1.getOrder() >= 4 ? 2.*b1(14) : 0.);
            if (std::abs(m1) < 1.) {
                x[4] += m1; 
                dbg<<"New mu = "<<x[4]<<std::endl;
            }
        }
        dbg<<"New x = "<<x<<std::endl;
#endif

        if (psf) {
            solver.reset(new EllipseSolver(pix,*psf,_fPsf,order,sigma,
                                           _isFixedCen,_isFixedGamma,true));
        } else {
            solver.reset(new EllipseSolver(pix,order,sigma,
                                           _isFixedCen,_isFixedGamma,true));
        }

        solver->useHybrid();
#ifdef NOTHROW
        solver->noUseCholesky();
        solver->noUseDirectH();
#endif
        if (iter == 1) {
            solver->setTol(1.e-8,1.e-10);
            if (XDEBUG) solver->setOutput(*dbgout);
            solver->setDelta0(0.01);
            solver->setMinStep(1.e-15);
            solver->setMaxIter(50);
        } else {
            solver->setTol(10.*solver->getFTol(),10.*solver->getGTol());
        }
    }
    xdbg<<"Done: Final solver pass successful: x = "<<x<<std::endl;
    xdbg<<"f = "<<f<<std::endl;
    xdbg<<"b = "<<solver->getB().vec()<<std::endl;
    if (cov) {
        solver->useSVD();
        solver->getCovariance(*cov);
        xdbg<<"cov = "<<*cov<<std::endl;
    }

    // Check for flags
    solver->callF(x,f);
    if (Norm(f) > solver->getFTol()) {
        xdbg<<"Local minimum\n";
        flag |= SHEAR_LOCAL_MIN;
    }
    if (solver->getFTol() > 1.e-8) {
        xdbg<<"ftol was raised\n";
        flag |= SHEAR_POOR_FIT;
    }

    if (XDEBUG && psf) {
        for(size_t k=0;k<pix.size();++k) {
            double sig_psf = (*psf)[k].getSigma();
            double D = sigma*sigma / (sig_psf*sig_psf);
            xdbg<<"D = "<<D<<std::endl;
        }
    }
    if (_shouldDoTimings) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        _times._tFinal += t2-t1;
    }

    setCen(std::complex<double>(x(0),x(1)));
    setGamma(std::complex<double>(x(2),x(3))); 
    setMu(x(4));
    if (bRet) {
        *bRet = solver->getB();
        xdbg<<"bret = "<<bRet->vec()<<std::endl;
        if (bCov) solver->getBCov(*bCov);
    }
    return true;
}

