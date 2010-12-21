
#include "Ellipse.h"
#include "EllipseSolver.h"
#include <cmath>
#include "dbg.h"
#include <fstream>
#include "PsiHelper.h"
#include "Params.h"

#define MAX_ITER 3
#define MAX_START 0.1
#define MAX_SHIFT 0.3
//#define TRY_BUMP
//#define TRY_HYBRID

bool Ellipse::doMeasure(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>* psf,
    int galOrderInit, int galOrder2, double sigma, long& flag, double thresh,
    DMatrix* cov, BVec* bRet, DMatrix* bCov, 
    std::vector<PixelList>* pixels_model)
{
#ifdef TRY_BUMP
    const int NBUMPS_MU = 8;
    const int NBUMPS_GAMMA = 24;
    const int NBUMPS_BOTH = 16;
    const double x_bumps_mu[NBUMPS_MU][5] = {
        { 0.0, 0.0, 0.0, 0.0, 0.02 },
        { 0.0, 0.0, 0.0, 0.0, 0.05 },
        { 0.0, 0.0, 0.0, 0.0, 0.1 },
        { 0.0, 0.0, 0.0, 0.0, 0.2 },
        { 0.0, 0.0, 0.0, 0.0, -0.02 },
        { 0.0, 0.0, 0.0, 0.0, -0.05 },
        { 0.0, 0.0, 0.0, 0.0, -0.1 },
        { 0.0, 0.0, 0.0, 0.0, -0.2 } };
    const double x_bumps_gamma[NBUMPS_GAMMA][5] = {
        { 0.0, 0.0, 0.02, 0.0, 0.0 },
        { 0.0, 0.0, 0.05, 0.0, 0.0 },
        { 0.0, 0.0, 0.1, 0.0, 0.0 },
        { 0.0, 0.0, 0.2, 0.0, 0.0 },
        { 0.0, 0.0, -0.02, 0.0, 0.0 },
        { 0.0, 0.0, -0.05, 0.0, 0.0 },
        { 0.0, 0.0, -0.1, 0.0, 0.0 },
        { 0.0, 0.0, -0.2, 0.0, 0.0 },
        { 0.0, 0.0, 0.0, 0.02, 0.0 },
        { 0.0, 0.0, 0.0, 0.05, 0.0 },
        { 0.0, 0.0, 0.0, 0.1, 0.0 },
        { 0.0, 0.0, 0.0, 0.2, 0.0 },
        { 0.0, 0.0, 0.0, -0.02, 0.0 },
        { 0.0, 0.0, 0.0, -0.05, 0.0 },
        { 0.0, 0.0, 0.0, -0.1, 0.0 },
        { 0.0, 0.0, 0.0, -0.2, 0.0 },
        { 0.0, 0.0, 0.02, 0.02, 0.0 },
        { 0.0, 0.0, 0.02, -0.02, 0.0 },
        { 0.0, 0.0, -0.02, 0.02, 0.0 },
        { 0.0, 0.0, -0.02, -0.02, 0.0 },
        { 0.0, 0.0, 0.05, 0.05, 0.0 },
        { 0.0, 0.0, 0.05, -0.05, 0.0 },
        { 0.0, 0.0, -0.05, 0.05, 0.0 },
        { 0.0, 0.0, -0.05, -0.05, 0.0 } };
    const double x_bumps_both[NBUMPS_BOTH][5] = {
        { 0.0, 0.0, 0.02, 0.02, 0.02 },
        { 0.0, 0.0, 0.02, -0.02, 0.02 },
        { 0.0, 0.0, -0.02, 0.02, 0.02 },
        { 0.0, 0.0, -0.02, -0.02, 0.02 },
        { 0.0, 0.0, 0.02, 0.02, -0.02 },
        { 0.0, 0.0, 0.02, -0.02, -0.02 },
        { 0.0, 0.0, -0.02, 0.02, -0.02 },
        { 0.0, 0.0, -0.02, -0.02, -0.02 },
        { 0.0, 0.0, 0.05, 0.05, 0.05 },
        { 0.0, 0.0, 0.05, -0.05, 0.05 },
        { 0.0, 0.0, -0.05, 0.05, 0.05 },
        { 0.0, 0.0, -0.05, -0.05, 0.05 },
        { 0.0, 0.0, 0.05, 0.05, -0.05 },
        { 0.0, 0.0, 0.05, -0.05, -0.05 },
        { 0.0, 0.0, -0.05, 0.05, -0.05 },
        { 0.0, 0.0, -0.05, -0.05, -0.05 } };
#endif

    xdbg<<"Start DoMeasure: galOrder = "<<galOrderInit<<", psf = "<<bool(psf)<<std::endl;
    xdbg<<"fix = "<<_isFixedCen<<"  "<<_isFixedGamma<<"  "<<_isFixedMu<<std::endl;
    xdbg<<"Initial flag = "<<flag<<std::endl;
    int npix = 0;
    for(size_t i=0;i<pix.size();++i) {
        xdbg<<"npix["<<i<<"] = "<<pix[i].size()<<std::endl;
        npix += pix[i].size();
    }

    int galOrder = galOrderInit;
    int galSize = (galOrder+1)*(galOrder+2)/2;
    if (npix <= galSize) {
        while (npix <= galSize) { galSize -= galOrder+1; --galOrder; }
        dbg<<"Too few pixels ("<<npix<<") for given gal_order. \n";
        dbg<<"Reduced gal_order to "<<galOrder<<" gal_size = "<<galSize<<std::endl;
    }

    BVec b(galOrder,sigma);

    // iterate until the shape is round.
    for(int iter=1;iter<=MAX_ITER;++iter) {
        dbg<<"iter = "<<iter<<std::endl;

        DVector x(5,0.);
        if (!doMeasureShapelet(pix,psf,b,galOrder,galOrder2,bCov)) {
            xdbg<<"Could not measure a shapelet vector.\n";
            return false;
        }

        xdbg<<"b = "<<b<<std::endl;

        // Initial estimates from b:
        if (!_isFixedCen) {
            // b10/b00 ~= conj(zc)/2
            x(0) = b(1)/b(0)*2.;
            x(1) = -b(2)/b(0)*2.;
        }
        if (!_isFixedGamma) {
            // b20/b00 ~= conj(gamma)/sqrt(2)
            x(2) = b(3)/b(0)*sqrt(2.); 
            x(3) = -b(4)/b(0)*sqrt(2.); 
        }
        if (!_isFixedMu) {
            // b11/b00 ~= mu
            x(4) = b(5)/b(0);
        }

        // Don't start with more that 0.1 in any element.
        DVector xinit = x;
        for(int k=0;k<5;++k) {
            if (x[k] > MAX_START) {
                xdbg<<"Adjusting x["<<k<<"] from "<<x[k];
                x[k] = MAX_START;
                xdbg<<" to "<<x[k]<<std::endl;
            } else if (x[k] < -MAX_START) {
                xdbg<<"Adjusting x["<<k<<"] from "<<x[k];
                x[k] = -MAX_START;
                xdbg<<" to "<<x[k]<<std::endl;
            }
        }

        xdbg<<"x = "<<EIGEN_Transpose(x)<<std::endl;
        DVector f(5);

        EllipseSolver3 solver(b,galOrder2,_isFixedCen,_isFixedGamma,_isFixedMu);

        //solver.useDogleg();
        solver.useHybrid();
#ifdef NOTHROW
        solver.noUseCholesky();
#endif
        double ftol = thresh*thresh;
        double gtol = thresh*ftol;
        solver.setTol(ftol,gtol);
        if (XDEBUG) solver.setOutput(*dbgout);
        //solver.useVerboseOutput();
        solver.setDelta0(0.01);
        solver.setMinStep(1.e-6*gtol);
        solver.setMaxIter(200);
        solver.useSVD();
        if (solver.solve(x,f)) {
            dbg<<"Found good round frame:\n";
            dbg<<"x = "<<EIGEN_Transpose(x)<<std::endl;
            dbg<<"b = "<<EIGEN_Transpose(solver.getB().vec())<<std::endl;
        } else {
            dbg<<"findRoundFrame solver failed:\n";
            dbg<<"x = "<<EIGEN_Transpose(x)<<std::endl;
            dbg<<"b = "<<EIGEN_Transpose(solver.getB().vec())<<std::endl;
            //if (XDEBUG) if (!solver.testJ(x,f,dbgout,1.e-5)) exit(1);
#ifdef TRY_HYBRID
            dbg<<"Try Hybrid method.\n";
            solver.useHybrid();
            if (solver.solve(x,f)) {
                dbg<<"Found good round frame using Hybrid:\n";
                dbg<<"x = "<<EIGEN_Transpose(x)<<std::endl;
                dbg<<"b = "<<EIGEN_Transpose(solver.getB().vec())<<std::endl;
            } else {
                dbg<<"findRoundFrame solver failed again with Hybrid:\n";
                dbg<<"x = "<<EIGEN_Transpose(x)<<std::endl;
                dbg<<"b = "<<EIGEN_Transpose(solver.getB().vec())<<std::endl;
#endif
#ifdef TRY_BUMP
                bool found_better = false;
                if (iter < MAX_ITER) {
                    dbg<<"Try bumping:\n";
                    DVector x0 = x;
                    double best_f_norminf = f.normInf();
                    if (!_isFixedMu) for(int k=0;k<NBUMPS_MU;++k) {
                        DVector x1 = x0 + tmv::VectorViewOf(x_bumps_mu[k],5);
                        dbg<<"Try "<<x1<<std::endl;
                        solver.callF(x1,f);
                        double f_norminf = f.normInf();
                        dbg<<"f_norminf = "<<f_norminf<<std::endl;
                        if (f_norminf < best_f_norminf) {
                            dbg<<"Found a better x = "<<x1<<std::endl;
                            x = x1;  best_f_norminf = f_norminf; found_better=true;
                        }
                    }
                    if (!_isFixedGamma) for(int k=0;k<NBUMPS_GAMMA;++k) {
                        DVector x1 = x0 + tmv::VectorViewOf(x_bumps_gamma[k],5);
                        dbg<<"Try "<<x1<<std::endl;
                        solver.callF(x1,f);
                        double f_norminf = f.normInf();
                        dbg<<"f_norminf = "<<f_norminf<<std::endl;
                        if (f_norminf < best_f_norminf) {
                            dbg<<"Found a better x = "<<x1<<std::endl;
                            x = x1;  best_f_norminf = f_norminf; found_better=true;
                        }
                    }

                    if (!_isFixedMu && !_isFixedGamma) for(int k=0;k<NBUMPS_BOTH;++k) {
                        DVector x1 = x0 + tmv::VectorViewOf(x_bumps_both[k],5);
                        dbg<<"Try "<<x1<<std::endl;
                        solver.callF(x1,f);
                        double f_norminf = f.normInf();
                        dbg<<"f_norminf = "<<f_norminf<<std::endl;
                        if (f_norminf < best_f_norminf) {
                            dbg<<"Found a better x = "<<x1<<std::endl;
                            x = x1;  best_f_norminf = f_norminf; found_better=true;
                        }
                    }
                    if (!found_better) {
                        dbg<<"Couldn't find something better by bumping.\n";
                        return false;
                    }
                } else {
                    dbg<<"Last iteration, so don't try bumping.\n";
                    return false; 
                }
#else
                return false;
#endif
#ifdef TRY_HYBRID
            }
#endif
        }
        dbg<<"x = "<<x<<std::endl;
        if (iter < MAX_ITER) {
            // Don't shift by more that 0.3 in any element.
            for(int k=0;k<5;++k) {
                if (x[k] > MAX_SHIFT) {
                    xdbg<<"Adjusting x["<<k<<"] from "<<x[k];
                    x[k] = MAX_SHIFT;
                    xdbg<<" to "<<x[k]<<std::endl;
                } else if (x[k] < -MAX_SHIFT) {
                    xdbg<<"Adjusting x["<<k<<"] from "<<x[k];
                    x[k] = -MAX_SHIFT;
                    xdbg<<" to "<<x[k]<<std::endl;
                }
            }
        }

        shiftBy(std::complex<double>(x(0),x(1))*sigma,
                std::complex<double>(x(2),x(3)),
                x(4));

        dbg<<"ell => "<<*this<<std::endl;
        double normx = Norm(x);

        if (normx < thresh) {
            dbg<<"norm(x) = "<<normx<<" < "<<thresh<<
                ", so good enough to stop.\n";
            if (cov) {
                solver.useSVD();
                solver.getCovariance(*cov);
                xdbg<<"cov = "<<*cov<<std::endl;
            }
            // Check for flags
            xdbg<<"Checking for possible problems\n";
            if (!(b(0) > 0.)) {
                xdbg<<"getB = "<<solver.getB().vec()<<", negative flux.\n";
                xdbg<<"flag SHEAR_BAD_FLUX\n";
                flag |= SHEAR_BAD_FLUX;
            }
            double f_normInf = f.normInf();
            if (f_normInf < solver.getFTol()) {
                xdbg<<"Good solution\n";
            } else if (f_normInf < thresh) {
                xdbg<<"f is near solution, but above tolerance\n";
                xdbg<<"f.normInf = "<<f_normInf<<std::endl;
                xdbg<<"ftol = "<<solver.getFTol()<<std::endl;
                xdbg<<"flag SHEAR_POOR_FIT\n";
                flag |= SHEAR_POOR_FIT;
            } else {
                xdbg<<"Local minimum\n";
                xdbg<<"x = "<<x<<std::endl;
                xdbg<<"f = "<<f<<std::endl;
                xdbg<<"f.norm = "<<f.norm()<<std::endl;
                xdbg<<"f.normInf = "<<f_normInf<<std::endl;
                xdbg<<"ftol = "<<solver.getFTol()<<std::endl;
                xdbg<<"flag SHEAR_LOCAL_MIN\n";
                flag |= SHEAR_LOCAL_MIN;
            }
            break;
        } else if (iter == MAX_ITER) {
            dbg<<"norm(x) = "<<normx<<" > 1.e-4, but maxiter reached.\n";
        } else {
            dbg<<"norm(x) = "<<normx<<" > 1.e-4, so try again.\n";
        }
    }

    if (bRet) {
        if (galOrder < galOrderInit) {
            xdbg<<"flag SHAPE_REDUCED_ORDER\n";
            flag |= SHAPE_REDUCED_ORDER;
        }
        *bRet = b;
        xdbg<<"bret = "<<*bRet<<std::endl;
    }

    return true;
}


