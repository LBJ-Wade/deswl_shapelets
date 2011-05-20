
#include "Ellipse.h"
#include "EllipseSolver.h"
#include <cmath>
#include "dbg.h"
#include <fstream>
#include "PsiHelper.h"
#include "Params.h"

bool Ellipse::doMeasure(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>* psf,
    int galOrder, int galOrder2, int maxm,
    double sigma, long& flag, double thresh, DSmallMatrix22* cov)
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
    std::auto_ptr<DMatrix> bCov;
    bool doMeanLikelihood = psf && !_isFixedGamma;
    if (cov || doMeanLikelihood) 
        bCov.reset(new DMatrix(b.size(),b.size()));

    if (!doMeasureShapelet(pix,psf,b,galOrder,galOrder2,maxm,bCov.get())) {
        xdbg<<"Could not measure a shapelet vector.\n";
        return false;
    }
    if (!b(0) > 0) {
        xdbg<<"Bad flux in measured shapelet\n";
        return false;
    }
    dbg<<"b = "<<b<<std::endl;

    return findRoundFrame(b,psf,galOrder2,thresh,flag,bCov.get(),cov);
}

bool Ellipse::findRoundFrame(
    const BVec& b, bool psf, int galOrder2, double thresh,
    long& flag, const DMatrix* bCov, DSmallMatrix22* cov)
{
    double trCov = cov ? Trace(*cov) : 0.;
    bool doMeanLikelihood = psf && !_isFixedGamma;
    bool doMaxLikelihood = !doMeanLikelihood || trCov == 0.0;

    if (doMaxLikelihood) {
        DVector x(5);
        DVector f(5);

        x.setZero();

        EllipseSolver solver(b,galOrder2,_isFixedCen,_isFixedGamma,_isFixedMu);

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
            if (!doMeanLikelihood) {
                double f_normInf = f.TMV_normInf();
                if (psf && !_isFixedMu && !(f_normInf < solver.getFTol())) {
                    xdbg<<"Oops, Local minimum, not real solution.\n";
                    xdbg<<"f.norm = "<<f.norm()<<std::endl;
                    xdbg<<"f.normInf = "<<f_normInf<<std::endl;
                    xdbg<<"ftol = "<<solver.getFTol()<<std::endl;
                    dbg<<"FLAG SHEAR_LOCAL_MIN\n";
                    flag |= SHEAR_LOCAL_MIN;
                    return false;
                }
            }
        }  else {
            dbg<<"findRoundFrame solver failed:\n";
            dbg<<"x = "<<EIGEN_Transpose(x)<<std::endl;
            dbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
            if (!doMeanLikelihood) {
                //if (XDEBUG) if (!solver.testJ(x,f,dbgout,1.e-5)) exit(1);
                dbg<<"FLAG SHEAR_DIDNT_CONVERGE\n";
                flag |= SHEAR_DIDNT_CONVERGE;
                return false;
            }
        }

        double sigma = b.getSigma();
        preShiftBy(std::complex<double>(x(0),x(1))*sigma,
                   std::complex<double>(x(2),x(3)),
                   x(4));
        removeRotation();

        dbg<<"ell => "<<*this<<std::endl;

        if (cov) {
            Assert(bCov);
            solver.useSVD();
            DMatrix cov5(5,5);
            solver.getCovariance(*bCov,cov5);
            xdbg<<"cov5 = "<<cov5<<std::endl;
            *cov = cov5.TMV_subMatrix(2,4,2,4);
            if (!(cov->TMV_det() > 0.)) {
                dbg<<"cov has bad determinant: "<<
                    cov->TMV_det()<<std::endl;
                dbg<<"cov = "<<*cov<<std::endl;
                dbg<<"Full cov = "<<cov5<<std::endl;
                dbg<<"FLAG SHEAR_BAD_COVAR\n";
                flag |= SHEAR_BAD_COVAR;
            }
            // update trCov for next section
            trCov = cov ? Trace(*cov) : 0.;
        } else if (doMeanLikelihood) {
            // Then need good estimate of trCov for step size.
            Assert(bCov);
            solver.useSVD();
            DMatrix cov5(5,5);
            solver.getCovariance(*bCov,cov5);
            xdbg<<"cov5 = "<<cov5<<std::endl;
            trCov = Trace(cov5.TMV_subMatrix(2,4,2,4));
        }
    }

    if (doMeanLikelihood) {
        Assert(bCov);
        dbg<<"Starting likelihood-weighted mean calculation.\n";
        //dbg<<"b = "<<b<<std::endl;
        //dbg<<"bCov = "<<*bCov<<std::endl;
        std::complex<double> gamma = getGamma();
        dbg<<"gamma = "<<gamma<<std::endl;
        if (cov) dbg<<"cov = "<<*cov<<std::endl;
        dbg<<"trCov = "<<trCov<<std::endl;
        double g1c = real(gamma);
        double g2c = imag(gamma);
        double step = sqrt(trCov/2.);
        if (step > 0.1) step = 0.1;
        dbg<<"step = "<<step<<std::endl;
        const int Ngrid = 10;
        double suml=0.;
        std::complex<double> sumlg=0.;
        double sumlg1g1=0., sumlg1g2=0., sumlg2g2=0.;
        DMatrix like(2*Ngrid,2*Ngrid,0.);
        dbg<<"g1 = "<<(g1c+(-Ngrid+0.5)*step)<<" ... "<<(g1c+(Ngrid-0.5)*step)<<std::endl;
        dbg<<"g2 = "<<(g2c+(-Ngrid+0.5)*step)<<" ... "<<(g2c+(Ngrid-0.5)*step)<<std::endl;
        for(int i=-Ngrid;i<Ngrid;++i) for(int j=-Ngrid;j<Ngrid;++j) {
            double g1 = g1c + (i+0.5)*step;
            double g2 = g2c + (j+0.5)*step;
            std::complex<double> gtry(g1,g2);
            xdbg<<"gtry = "<<gtry<<std::endl;
            if (abs(gtry) >= 1.) {
                xdbg<<"|g| >= 1\n";
                continue;
            }
            Ellipse e1;
            e1.setGamma(gtry);
            e1.preShiftBy(0.,-gamma,0.);
            e1.removeRotation();
            std::complex<double> dg = e1.getGamma();
            xdbg<<"e1 = "<<e1<<std::endl;
            xdbg<<"dg = "<<dg<<std::endl;
            DMatrix S(6,b.size());
            calculateGTransform(dg,2,b.getOrder(),S);
            xdbg<<"S = "<<S<<std::endl;
            DVector b1 = S * b.vec();
            xdbg<<"Sheared b = "<<b1<<std::endl;
            DMatrix c1 = S * (*bCov) * S.transpose();
            xdbg<<"c1 = "<<c1<<std::endl;
            DVector b2 = b1.subVector(3,5);
            DMatrix c2 = c1.subMatrix(3,5,3,5);
            double chisq = b2 * (b2/c2);
            dbg<<"chisq("<<gtry<<") = "<<chisq<<std::endl;
            double l = exp(-chisq/2.);
            like(i+Ngrid,j+Ngrid) = l;
            sumlg += l * gtry;
            suml += l;
            sumlg1g1 += l * g1*g1;
            sumlg1g2 += l * g1*g2;
            sumlg2g2 += l * g2*g2;
        }
        sumlg /= suml;

        dbg<<"like = "<<like<<std::endl;
        dbg<<"likelihood-weighted mean = "<<sumlg<<std::endl;
        dbg<<"maximum likelihood value = "<<gamma<<std::endl;
        setGamma(sumlg);
        if (cov) {
            sumlg1g1 /= suml;  sumlg1g1 -= real(sumlg) * real(sumlg);
            sumlg1g2 /= suml;  sumlg1g2 -= real(sumlg) * imag(sumlg);
            sumlg2g2 /= suml;  sumlg2g2 -= imag(sumlg) * imag(sumlg);
            dbg<<"old cov = "<<(*cov)(0,0)<<"  "<<(*cov)(0,1)<<"  "<<(*cov)(1,1)<<std::endl;
            dbg<<"new cov = "<<sumlg1g1<<"  "<<sumlg1g2<<"  "<<sumlg2g2<<std::endl;
            (*cov)(0,0) = sumlg1g1;
            (*cov)(0,1) = sumlg1g2;
            (*cov)(1,0) = sumlg1g2;
            (*cov)(1,1) = sumlg2g2;
        }
    }

    return true;
}


