
#include "Ellipse.h"
#include "EllipseSolver.h"
#include <cmath>
#include "dbg.h"
#include <fstream>
#include "PsiHelper.h"
#include "Params.h"

#define MAX_START 0.1
#define MAX_SHEAR 1.0
#define TRY_HYBRID
//#define TRY_FIXED_MU

bool Ellipse::doMeasure(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>* psf,
    int galOrderInit, int galOrder2, double sigma, long& flag, double thresh,
    DMatrix* cov, BVec* bRet, DMatrix* bCov, 
    std::vector<PixelList>* pixels_model)
{
    dbg<<"Start DoMeasure: galOrder = "<<galOrderInit<<", psf = "<<bool(psf)<<std::endl;
    dbg<<"fix = "<<_isFixedCen<<"  "<<_isFixedGamma<<"  "<<_isFixedMu<<std::endl;
    dbg<<"Initial flag = "<<flag<<std::endl;
    dbg<<"Thresh = "<<thresh<<std::endl;
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

    // The shapelet measurement becomes inaccurate for large shears
    // because the psf needs to be sheared, and the finite order of the 
    // decomosition means that we lose information off the end as the 
    // shear is larger.
    // The solution is to shear by at most 0.2 for the shapelet measurement,
    // even if the current shear is larger.
    double absg = std::abs(_gamma);
    if (absg > MAX_SHEAR) {
        _gamma /= (absg/0.1); 
        dbg<<"Scaled gamma back to "<<_gamma<<" for the shapelet measurement.\n";
    }

    if (!doMeasureShapelet(pix,psf,b,galOrder,galOrder2,bCov)) {
        xdbg<<"Could not measure a shapelet vector.\n";
        return false;
    }
    xdbg<<"b = "<<b<<std::endl;

    if (bRet) {
        if (galOrder < galOrderInit) {
            xdbg<<"flag SHAPE_REDUCED_ORDER\n";
            flag |= SHAPE_REDUCED_ORDER;
        }
        *bRet = b;
        xdbg<<"bret = "<<*bRet<<std::endl;
    }

    return findRoundFrame(b,psf,galOrder2,thresh,flag,cov);
}

bool Ellipse::findRoundFrame(
    const BVec& b, bool psf, int galOrder2, double thresh,
    long& flag, DMatrix* cov)
{
    DVector x(5,0.);

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
    if (!_isFixedMu && !psf) {
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

#ifdef NOTHROW
    solver.noUseCholesky();
#endif
    double ftol = thresh*thresh;
    double gtol = thresh*ftol;
    solver.setTol(ftol,gtol);
    solver.setMinStep(gtol*thresh);
    solver.setOutput(*dbgout);
    if (XDEBUG) solver.useExtraVerboseOutput();
    solver.setMinStep(1.e-6*gtol);
    solver.setMaxIter(200);
    solver.useSVD();
    do { // Just provide a place to break out to once we succeed.
        if (psf && !_isFixedMu) {
#ifdef TRY_FIXED_MU
            dbg<<"First keep Mu fixed\n.";
            // First try to find a solution with mu fixed.
            EllipseSolver3 solver2(
                b,galOrder2,_isFixedCen,_isFixedGamma,true);
            solver2.setTol(ftol,gtol);
            solver2.setOutput(*dbgout);
            if (XDEBUG) solver.useVerboseOutput();
            solver2.useHybrid();
            solver2.setMaxIter(200);
            if (solver2.solve(x,f) && NormInf(f) <= ftol) {
                dbg<<"Found good round frame without using mu:\n";
                dbg<<"x = "<<EIGEN_Transpose(x)<<std::endl;
                dbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
                break;
            }
            // Next let Mu help find a good round frame, but don't require 
            // b11 = 0.
            dbg<<"Now let mu vary, but don't try to zero b11.\n";
#endif
            solver.setDelta0(0.01);
            solver.useDogleg();
            solver.dontZeroB11();
        } else {
            dbg<<"!psf, so just use the Hybrid method.\n";
            // Otherwise, dogleg method isn't terribly effective.
            // Plus we're not going to care much if we only get a local 
            // min anyway, so just use the hybrid method which is more 
            // robust when we might hit a local minimum.
            solver.useHybrid();
        }
        if (solver.solve(x,f)) {
            dbg<<"Found good round frame:\n";
            dbg<<"x = "<<EIGEN_Transpose(x)<<std::endl;
            dbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
            break;
        } 
        dbg<<"findRoundFrame solver failed:\n";
        dbg<<"x = "<<EIGEN_Transpose(x)<<std::endl;
        dbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
        if (XDEBUG) if (!solver.testJ(x,f,dbgout,1.e-5)) exit(1);
#ifdef TRY_HYBRID
        if (psf && !_isFixedMu) {
            dbg<<"Try Hybrid method.\n";
            solver.useHybrid();
            if (solver.solve(x,f)) {
                dbg<<"Found good round frame using Hybrid:\n";
                dbg<<"x = "<<EIGEN_Transpose(x)<<std::endl;
                dbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
                break;
            } 
            dbg<<"findRoundFrame solver failed again with Hybrid:\n";
            dbg<<"x = "<<EIGEN_Transpose(x)<<std::endl;
            dbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
        }
#endif
        // If we get here, then nothing worked.
        // (Any success is followed by a break statement.)
        xdbg<<"flag SHEAR_DIDNT_CONVERGE\n";
        flag |= SHEAR_DIDNT_CONVERGE;
        return false;
    } while (false);
    dbg<<"x = "<<x<<std::endl;

    double sigma = b.getSigma();
    preShiftBy(std::complex<double>(x(0),x(1))*sigma,
               std::complex<double>(x(2),x(3)),
               x(4));

    dbg<<"ell => "<<*this<<std::endl;

    if (cov) {
        solver.useSVD();
        solver.getCovariance(*cov);
        xdbg<<"cov = "<<*cov<<std::endl;
    }
    // Check for flags
    xdbg<<"Checking for possible problems\n";
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
    return true;
}


