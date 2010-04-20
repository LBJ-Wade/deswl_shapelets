
#include "Ellipse.h"
#include "dbg.h"
#include "Params.h"
#include "BVec.h"
#include "Ellipse.h"
#include "Pixel.h"
#include "FittedPsf.h"
#include "Log.h"
#include "TimeVars.h"
#include "MeasureShearAlgo.h"

// TODO: This function is almost completely degenerage with the one in
//       MeasureShearAlgo.cpp.  I'm pretty sure these can be combined into
//       a single measureShear function.  This TODO is related to the
//       one below that the function is too long.  The similarities
//       will probably be more apparent when I break it up into smaller
//       chunks with a single driver function.

void measureMultiShear(
    Position& cen, 
    const std::vector<PixelList>& allpix,
    const std::vector<BVec>& psf,
    double galAperture, double maxAperture,
    int galOrder, int galOrder2,
    double fPsf, double minGalSize, bool fixCen,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    DSmallMatrix22& shearcov, BVec& shapelet,
    double& nu, long& flag)
{
    //TODO: This function is too long.  It should really be broken 
    //      up in to smaller units...
    //      In which case most of the overlap between this and 
    //      MeasureSingleShear1 would not have to be duplicated.

    dbg<<"Start MeasureMultiShear\n";
    dbg<<"cen = "<<cen<<std::endl;
    dbg<<"allpix.size = "<<allpix.size()<<std::endl;
    const int nExp = allpix.size();
    for(int i=0;i<nExp;++i)
        dbg<<"allpix["<<i<<"].size = "<<allpix[i].size()<<std::endl;
    dbg<<"psf.size = "<<psf.size()<<std::endl;
    Assert(psf.size() == allpix.size());

    try {
        // Find harmonic mean of psf sizes:
        // MJ: Is it correct to use the harmonic mean of sigma^2?
        xdbg<<"sigmaP values are: \n";
        double sigmaP = 0.;
        Assert(psf.size() > 0);
        for(int i=0;i<nExp;++i) {
            xdbg<<psf[i].getSigma()<<std::endl;
            sigmaP += 1./psf[i].getSigma();
        }
        sigmaP = double(psf.size()) / sigmaP;
        xdbg<<"Harmonic mean = "<<sigmaP<<std::endl;

        // We start with a native fit using sigma = 2*sigmaP:
        double sigmaObs = 2.*sigmaP;
        double galAp = galAperture * sigmaObs;
        xdbg<<"galAp = "<<galAperture<<" * "<<sigmaObs<<" = "<<galAp<<std::endl;
        if (maxAperture > 0. && galAp > maxAperture) {
            galAp = maxAperture;
            xdbg<<"      => "<<galAp<<std::endl;
        }
        xdbg<<"sigmaObs = "<<sigmaObs<<", sigmaP = "<<sigmaP<<std::endl;

        std::vector<PixelList> pix(allpix.size());
        int npix = 0;
        for(int i=0;i<nExp;++i) {
            getSubPixList(pix[i],allpix[i],galAp,flag);
            npix += pix[i].size();
        }
        xdbg<<"npix = "<<npix<<std::endl;

        // TODO: If we have single-epoch measurements, use those for starting
        // guess rather than crudeMeasure.
        Ellipse ell;
        ell.fixGam();
        if (fixCen) ell.fixCen();
        ell.crudeMeasure(pix,sigmaObs);
        xdbg<<"Crude Measure: centroid = "<<ell.getCen()<<", mu = "<<ell.getMu()<<std::endl;

        int go = 2;
        int galSize = (go+1)*(go+2)/2;

        if (times) ell.doTimings();
        long flag1=0;

        if (ell.measure(pix,go,sigmaObs,true,flag1)) {
            if (times) times->successNative(ell.getTimes());
            ++log._nsNative;
            dbg<<"Successful native fit:\n";
            dbg<<"Z = "<<ell.getCen()<<std::endl;
            dbg<<"Mu = "<<ell.getMu()<<std::endl;
        } else {
            if (times) times->failNative(ell.getTimes());
            ++log._nfNative;
            dbg<<"Native measurement failed\n";
            flag |= NATIVE_FAILED;
            if (!fixCen) cen += ell.getCen();
            return;
        }

        // Now we can do a deconvolving fit, but one that does not 
        // shear the coordinates.
        sigmaObs *= exp(ell.getMu());
        xdbg<<"sigmaObs -> "<<sigmaObs<<std::endl;
        if (sigmaObs < minGalSize*sigmaP) {
            dbg<<"skip: galaxy is too small -- "<<sigmaObs<<
                " psf size = "<<sigmaP<<std::endl;
            if (times) ++times->_nfSmall;
            ++log._nfSmall;
            flag |= TOO_SMALL;
            if (!fixCen) cen += ell.getCen();
            return;
        }
        galAp *= exp(ell.getMu());
        if (maxAperture > 0. && galAp > maxAperture) galAp = maxAperture;
        ell.setMu(0);

        // Check - mu should now be zero:
        // This one works
        //ell.measure(pix,go,sigmaObs,true,flag1);
        xxdbg<<"After native fit #2:\n";
        xxdbg<<"Mu = "<<ell.getMu()<<std::endl;

        npix = 0;
        for(int i=0;i<nExp;++i) {
            getSubPixList(pix[i],allpix[i],galAp,flag);
            npix += pix[i].size();
        }

        xdbg<<"npix = "<<npix<<std::endl;
        double sigpsq = pow(sigmaP,2);
        double sigma = pow(sigmaObs,2) - sigpsq;
        xdbg<<"sigmaP^2 = "<<sigpsq<<std::endl;
        xdbg<<"sigma^2 = sigmaObs^2 - sigmap^2 = "<<sigma<<std::endl;
        if (sigma < 0.1*sigpsq) sigma = 0.1*sigpsq;
        sigma = sqrt(sigma);
        xdbg<<"sigma = sqrt(sigma^2) "<<sigma<<std::endl;
        ell.setFP(fPsf);
        if (times) ell.resetTimes();
        flag1 = 0;
        if (ell.measure(pix,psf,go,sigma,false,flag1)) {
            if (times) times->successMu(ell.getTimes());
            ++log._nsMu;
            if (flag1) {
                Assert((flag1 & ~(SHEAR_LOCAL_MIN | SHEAR_POOR_FIT)) == false);
                // This is just bumps the Ellipse flags up two spots
                // from the SHEAR_* to SHAPE_* flags.
                flag1 <<= 2;
                Assert((flag1 & ~(SHAPE_LOCAL_MIN | SHAPE_POOR_FIT)) == false);
                flag |= flag1;
            }
            dbg<<"Successful deconvolving fit:\n";
            dbg<<"Mu = "<<ell.getMu()<<std::endl;
        } else {
            if (times) times->failMu(ell.getTimes());
            ++log._nfMu;
            dbg<<"Deconvolving measurement failed\n";
            flag |= DECONV_FAILED;
        }

#if 0
        // This doesn't work.  So for now, keep mu, sigma the same
        // I'd rather have mu = 0 and the right sigma.
        // Maybe I need to change to fit to solve for sigma directly
        // rather than solve for mu.
        sigma *= exp(ell.getMu());
        ell.setMu(0);

        // Check - mu should now be zero:
        b_gal.setSigma(sigma);
        dbg<<"Meausre with sigma = "<<sigma<<std::endl;
        ell.measure(pix,psf,go,sigma,false,flag1,&b_gal);
        dbg<<"After deconvolving fit #2:\n";
        dbg<<"Mu = "<<ell.getMu()<<std::endl;
        dbg<<"b_gal = "<<b_gal<<std::endl;
#endif

        // Measure the galaxy shape at the full order
        go = galOrder;
        galSize = (go+1)*(go+2)/2;
        if (npix <= galSize) {
            while (npix <= galSize) { galSize -= go+1; --go; }
            dbg<<"Too few pixels ("<<npix<<") for given galOrder. \n";
            dbg<<"Reduced galOrder to "<<go<<" galSize = "<<galSize<<std::endl;
            flag |= SHAPE_REDUCED_ORDER;
            if (go < 2) { 
                flag |= TOO_SMALL; 
                if (!fixCen) cen += ell.getCen();
                return; 
            }
        }
        shapelet.setSigma(sigma);
        std::complex<double> gale = 0.;
        ell.measureShapelet(pix,psf,shapelet,go);
        xdbg<<"Measured b_gal = "<<shapelet.vec()<<std::endl;

        // Also measure the isotropic significance
        BVec flux(0,sigma);
        DMatrix fluxCov(1,1);
        ell.measureShapelet(pix,psf,flux,0,&fluxCov);
        nu = flux(0) / sqrt(fluxCov(0,0));
        dbg<<"nu = "<<flux(0)<<" / sqrt("<<fluxCov(0,0)<<") = "<<nu<<std::endl;

        // If the b00 value in the shapelet doesn't match the direct flux
        // measurement, set a flag.
        if (!(flux(0) > 0.0  &&
              shapelet(0) >= flux(0)/3. &&
              shapelet(0) <= flux(0)*3.)) {
            dbg<<"Bad flux value: \n";
            dbg<<"flux = "<<flux(0)<<std::endl;
            dbg<<"shapelet = "<<shapelet.vec()<<std::endl;
            flag |= SHAPE_BAD_FLUX;
        }
        // If the above deconvolving fit failed, then return.
        if (flag & DECONV_FAILED) {
            if (!fixCen) cen += ell.getCen();
            return;
        }

        // Under normal circumstances, b20/b00 ~= conj(gamma)/sqrt(2)
        gale = std::complex<double>(shapelet(3),shapelet(4));
        gale /= shapelet(0);
        if (std::abs(gale) < 0.5) {
            ell.setGamma(conj(gale) * sqrt(2.));
            shear = ell.getGamma();
        }

        // Finally, we calculate the shear in the deconvolved galaxy.
        //ell.fixMu();
        ell.unfixGam();
        go = galOrder2;
        galSize = (go+1)*(go+2)/2;
        if (npix <= galSize) {
            while (npix <= galSize) { galSize -= go+1; --go; }
            dbg<<"Too few pixels ("<<npix<<") for given galOrder. \n";
            dbg<<"Reduced galOrder to "<<go<<" galSize = "<<galSize<<std::endl;
        }
        if (times) ell.resetTimes();
        DMatrix cov(5,5);
        if (ell.measure(pix,psf,go,sigma,false,flag,&cov)) {
            if (times) times->successGamma(ell.getTimes());
            ++log._nsGamma;
            dbg<<"Successful Gamma fit\n";
            dbg<<"Measured gamma = "<<ell.getGamma()<<std::endl;
        } else {
            if (times) times->failGamma(ell.getTimes());
            ++log._nfGamma;
            dbg<<"Gamma measurement failed\n";
            flag |= SHEAR_FAILED;
        }
        shear = ell.getGamma();
        DSmallMatrix22 cov1 = cov.TMV_subMatrix(2,4,2,4);
        if (!(cov1.TMV_det() > 0.)) flag |= SHEAR_BAD_COVAR;
        else shearcov = cov1;
        if (!fixCen) cen += ell.getCen();
#ifdef USE_TMV
    } catch (tmv::Error& e) {
        dbg<<"TMV Error thrown in MeasureSingleShear\n";
        dbg<<e<<std::endl;
        ++log._nfTmvError;
        flag |= TMV_EXCEPTION;
#endif
    } catch (std::exception& e) {
        dbg<<"std::exception thrown in MeasureSingleShear\n";
        dbg<<e.what()<<std::endl;
        ++log._nfOtherError;
        flag |= STD_EXCEPTION;
    } catch (...) {
        dbg<<"unkown exception in MeasureSingleShear\n";
        ++log._nfOtherError;
        flag |= UNKNOWN_EXCEPTION;
    } 
}


