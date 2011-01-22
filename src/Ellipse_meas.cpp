
#include "Ellipse.h"
#include "EllipseSolver.h"
#include <cmath>
#include "dbg.h"
#include <fstream>
#include "PsiHelper.h"
#include "Params.h"

#define MAX_START 0.1

bool Ellipse::doMeasure(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>* psf,
    int galOrder, int galOrder2, int maxm,
    double sigma, long& flag, double thresh, DMatrix* cov)
{
    dbg<<"Start DoMeasure: galOrder = "<<galOrder<<", psf = "<<bool(psf)<<std::endl;
    dbg<<"fix = "<<_isFixedCen<<"  "<<_isFixedGamma<<"  "<<_isFixedMu<<std::endl;
    dbg<<"Thresh = "<<thresh<<std::endl;
    int npix = 0;
    for(size_t i=0;i<pix.size();++i) {
        xdbg<<"npix["<<i<<"] = "<<pix[i].size()<<std::endl;
        npix += pix[i].size();
    }

    int galSize = (galOrder+1)*(galOrder+2)/2;
    if (npix <= galSize) {
        dbg<<"Too few pixels ("<<npix<<") for given gal_order. \n";
        return false;
    }

    BVec b(galOrder,sigma);

    if (!doMeasureShapelet(pix,psf,b,galOrder,galOrder2,maxm)) {
        xdbg<<"Could not measure a shapelet vector.\n";
        return false;
    }
    if (!b(0) > 0) {
        xdbg<<"Bad flux in measured shapelet\n";
        return false;
    }
    xdbg<<"b = "<<b<<std::endl;

    return findRoundFrame(b,psf,galOrder2,thresh,flag,cov);
}

bool Ellipse::findRoundFrame(
    const BVec& b, bool psf, int galOrder2, double thresh,
    long& flag, DMatrix* cov)
{
    DVector x(5,0.);
    DVector f(5);

#if 0
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
#endif

    EllipseSolver3 solver(b,galOrder2,_isFixedCen,_isFixedGamma,_isFixedMu);

#ifdef NOTHROW
    solver.noUseCholesky();
#endif
    double ftol = thresh*thresh;
    double gtol = thresh*ftol;
    solver.setTol(ftol,gtol);
    solver.setMinStep(gtol*thresh);
    solver.setOutput(*dbgout);
    if (XDEBUG) solver.useVerboseOutput();
    solver.setMinStep(1.e-6*gtol);
    solver.setDelta0(0.01);
    solver.setMaxIter(200);
    if (psf && !_isFixedMu) {
        solver.setDelta0(0.01);
        solver.useDogleg();
        solver.dontZeroB11();
        solver.useSVD();
    } else {
        solver.useHybrid();
    }
    if (solver.solve(x,f)) {
        dbg<<"Found good round frame:\n";
        dbg<<"x = "<<EIGEN_Transpose(x)<<std::endl;
        dbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
        double f_normInf = f.normInf();
        if (psf && !_isFixedMu && !(f_normInf < solver.getFTol())) {
            xdbg<<"Oops, Local minimum, not real solution.\n";
            xdbg<<"f.norm = "<<f.norm()<<std::endl;
            xdbg<<"f.normInf = "<<f_normInf<<std::endl;
            xdbg<<"ftol = "<<solver.getFTol()<<std::endl;
            dbg<<"FLAG SHEAR_LOCAL_MIN\n";
            flag |= SHEAR_LOCAL_MIN;
            return false;
        }
    }  else {
        dbg<<"findRoundFrame solver failed:\n";
        dbg<<"x = "<<EIGEN_Transpose(x)<<std::endl;
        dbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
        //if (XDEBUG) if (!solver.testJ(x,f,dbgout,1.e-5)) exit(1);
        dbg<<"FLAG SHEAR_DIDNT_CONVERGE\n";
        flag |= SHEAR_DIDNT_CONVERGE;
        return false;
    }

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

    return true;
}


