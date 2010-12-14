
#include "dbg.h"
#include "Params.h"
#include "BVec.h"
#include "Ellipse.h"
#include "Pixel.h"
#include "Image.h"
#include "FittedPsf.h"
#include "Log.h"
#include "Params.h"
#include "MeasureShearAlgo.h"


void measureSingleShear1(
    Position& cen, const Image<double>& im, double sky,
    const Transformation& trans, const std::vector<BVec>& psf,
    double noise, double gain, const Image<double>* weightIm, 
    double galAperture, double maxAperture,
    int galOrder, int galOrder2,
    double fPsf, double minGalSize, bool fixCen,
    double xOffset, double yOffset,
    bool fixSigma, double fixSigmaValue,
    ShearLog& log, std::complex<double>& shear, 
    DSmallMatrix22& shearcov, BVec& shapelet,
    double& nu, long& flag)
{
    // Find harmonic mean of psf sizes:
    // MJ: Is it better to use the harmonic mean of sigma or sigma^2?
    //     (Or something else entirely?)
    double sigmaP = 0.;
    const int nPsf = psf.size();
    Assert(nPsf > 0);
    for(int i=0;i<nPsf;++i) {
        sigmaP += 1./psf[i].getSigma();
        xdbg<<"psf["<<i<<"] sigma = "<<psf[i].getSigma()<<std::endl;
    }
    sigmaP = double(nPsf) / sigmaP;
    xdbg<<"Harmonic mean = "<<sigmaP<<std::endl;

    //
    // Define initial sigma and aperture.
    // We start with a native fit using sigma = 2*sigmaP:
    // Or if we are fixed, then sigmaObs^2 = sigmaP^2 + sigma^2
    //
    double sigma = fixSigmaValue;
    double sigpsq = pow(sigmaP,2);
    double sigmaObs = 
        fixSigma ? 
        sqrt(sigpsq + sigma*sigma) : 
        2.*sigmaP;
    double galAp = galAperture * sigmaObs;
    xdbg<<"galap = "<<galAperture<<" * "<<sigmaObs<<" = "<<galAp<<std::endl;
    if (maxAperture > 0. && galAp > maxAperture) {
        galAp = maxAperture;
        xdbg<<"      => "<<galAp<<std::endl;
    }
    xdbg<<"sigma_obs = "<<sigmaObs<<", sigma_p = "<<sigmaP<<std::endl;
    long flag1=0;

    //
    // Load pixels from image
    //
    std::vector<PixelList> pix(1);
    getPixList(im,pix[0],cen,sky,noise,gain,weightIm,trans,
               galAp,xOffset,yOffset,flag);
    int npix = pix[0].size();
    xdbg<<"npix = "<<npix<<std::endl;
    if (npix < 10) {
        dbg<<"Too few pixels to continue: "<<npix<<std::endl;
        xdbg<<"flag LT10PIX\n";
        flag |= LT10PIX; // Should already be set by getPixList, but...
        xdbg<<"flag SHAPELET_NOT_DECONV\n";
        flag |= SHAPELET_NOT_DECONV;
        return;
    }
    int go = galOrder;
    int galSize = (go+1)*(go+2)/2;
    if (npix <= galSize) {
        while (npix <= galSize) { galSize -= go+1; --go; }
        dbg<<"Too few pixels ("<<npix<<") for given gal_order. \n";
        dbg<<"Reduced gal_order to "<<go<<" gal_size = "<<galSize<<std::endl;
        xdbg<<"flag SHAPE_REDUCED_ORDER\n";
        flag |= SHAPE_REDUCED_ORDER;
    }

    //
    // Do a crude measurement based on simple pixel sums.
    //
    Ellipse ell;
    ell.fixGam();
    if (fixCen) ell.fixCen();
    if (fixSigma) ell.fixMu();
    if (!fixCen || !fixSigma) {
        ell.crudeMeasure(pix[0],sigmaObs);
        xdbg<<"Crude Measure: centroid = "<<ell.getCen()<<
            ", mu = "<<ell.getMu()<<std::endl;
        sigmaObs *= exp(ell.getMu());
        ell.setMu(0.);
    }

    //
    // Correct the centroid first.
    //
    if (!fixCen) {
        flag1=0;
        ell.fixGam();
        ell.fixMu();
        shapelet.setSigma(sigmaObs);
        ell.measureShapelet(pix,shapelet,go,galOrder2);
        if (!ell.findRoundFrame(shapelet,galOrder2,flag1)) {
            ++log._nfCentroid;
            dbg<<"First centroid pass failed.\n";
            xdbg<<"flag CENTROID_FAILED\n";
            flag |= CENTROID_FAILED;
            xdbg<<"flag SHAPELET_NOT_DECONV\n";
            flag |= SHAPELET_NOT_DECONV;
            return;
        }
        dbg<<"After first centroid pass: cen = "<<ell.getCen()<<std::endl;

        // redo the pixlist:
        pix[0].clear();
        cen += ell.getCen();
        dbg<<"New center = "<<cen<<std::endl;
        getPixList(im,pix[0],cen,sky,noise,gain,weightIm,trans,
                   galAp,xOffset,yOffset,flag);
        npix = pix[0].size();
        if (npix < 10) {
            dbg<<"Too few pixels to continue: "<<npix<<std::endl;
            xdbg<<"flag LT10PIX\n";
            flag |= LT10PIX; // Should already be set by getPixList, but...
            xdbg<<"flag SHAPELET_NOT_DECONV\n";
            flag |= SHAPELET_NOT_DECONV;
            return;
        }
        go = galOrder;
        galSize = (go+1)*(go+2)/2;
        if (npix <= galSize) {
            while (npix <= galSize) { galSize -= go+1; --go; }
            dbg<<"Too few pixels ("<<npix<<") for given gal_order. \n";
            dbg<<"Reduced gal_order to "<<go<<" gal_size = "<<galSize<<std::endl;
            xdbg<<"flag SHAPE_REDUCED_ORDER\n";
            flag |= SHAPE_REDUCED_ORDER;
        }

// see Params.h
#ifdef RANDOMIZE_CENTER
        // The centroid doesn't necessarily get all the way to the perfect
        // centroid, so if the input guess is consistently biased in some 
        // direction (e.g. by a different convention for the (0,0) or (1,1) 
        // location), then this results in a bias in the shear values.  
        // To avoid that, we first get to what we think it the correct centroid.
        // Then we redo the aperture from this point.  Then we shift the 
        // centroid value by 0.1 pixel in a random direction, and resolve for
        // the centroid.  This is all pretty fast, so it's worth doing to
        // avoid the subtle bias that can otherwise result.

        std::complex<double> cen1 = ell.getCen();
        // get a offset of 0.1 pixels in a random direction.
        double x_offset, y_offset, rsq_offset;
        do {
            x_offset = double(rand())*2./double(RAND_MAX) - 1.;
            y_offset = double(rand())*2./double(RAND_MAX) - 1.;
            rsq_offset = x_offset*x_offset + y_offset*y_offset;
        } while (rsq_offset > 1. || rsq_offset == 0.);
        double r_offset = sqrt(rsq_offset);
        x_offset /= r_offset*10.;
        y_offset /= r_offset*10.;
        //std::cout<<"x_offset: "<<x_offset<<"   y_offset: "<<y_offset<<"\n";
        ell.setCen(std::complex<double>(x_offset,y_offset));
        dbg<<"After random offset "<<x_offset<<","<<y_offset<<
            ": cen = "<<cen1+ell.getCen()<<std::endl;

        // Now do it again.
        flag1=0;
        ell.measureShapelet(pix,shapelet,go,galOrder2);
        if (!ell.findRoundFrame(shapelet,galOrder2,flag1)) {
            ++log._nfCentroid;
            dbg<<"Second centroid pass failed.\n";
            xdbg<<"flag CENTROID_FAILED\n";
            flag |= CENTROID_FAILED;
            xdbg<<"flag SHAPELET_NOT_DECONV\n";
            flag |= SHAPELET_NOT_DECONV;
            return;
        }
        dbg<<"After second centroid pass: cen = "<<cen1+ell.getCen()<<std::endl;
        dbg<<"Which is now considered "<<ell.getCen()<<" relative to the new "
            "center value of "<<cen<<std::endl;
#else
        ell.setCen(std::complex<double>(0.,0.));
#endif
        ++log._nsCentroid;
    }

    //
    // Now adjust the sigma value 
    //
    flag1=0;
    if (!fixSigma) {
        ell.unfixMu();
        shapelet.setSigma(sigmaObs);
        ell.measureShapelet(pix,shapelet,go,galOrder2);
        if (ell.findRoundFrame(shapelet,galOrder2,flag1)) {
            ++log._nsNative;
            dbg<<"Successful native fit:\n";
            dbg<<"Z = "<<ell.getCen()<<std::endl;
            dbg<<"Mu = "<<ell.getMu()<<std::endl;
        } else {
            ++log._nfNative;
            dbg<<"Native measurement failed\n";
            xdbg<<"flag NATIVE_FAILED\n";
            flag |= NATIVE_FAILED;
            if (!fixCen) cen += ell.getCen();
            xdbg<<"flag SHAPELET_NOT_DECONV\n";
            flag |= SHAPELET_NOT_DECONV;
            return;
        }

        // Adjust the sigma value so we can reset mu = 0.
        sigmaObs *= exp(ell.getMu());
        if (sigmaObs < minGalSize*sigmaP) {
            dbg<<"skip: galaxy is too small -- "<<sigmaObs<<
                " psf size = "<<sigmaP<<std::endl;
            ++log._nfSmall;
            xdbg<<"flag TOO_SMALL\n";
            flag |= TOO_SMALL;
            if (!fixCen) cen += ell.getCen();
            xdbg<<"flag SHAPELET_NOT_DECONV\n";
            flag |= SHAPELET_NOT_DECONV;
            return;
        }
        galAp *= exp(ell.getMu());
        if (maxAperture > 0. && galAp > maxAperture) galAp = maxAperture;
        ell.setMu(0.);

    } 
    // Check - mu should now be zero:
    if (XDEBUG) {
        BVec b2 = shapelet;
        Ellipse e2 = ell;
        long flag2 = 0;
        b2.setSigma(sigmaObs);
        e2.measureShapelet(pix,b2,go,galOrder2);
        e2.findRoundFrame(b2,galOrder2,flag2);
        xdbg<<"After native fit #2:\n";
        xdbg<<"Mu = "<<e2.getMu()<<std::endl;
        xdbg<<"(Should be very close to zero.)\n";
#if 1
        double target_b00 = 7.411558529;
        xdbg<<"Rescaled to have b00 = "<<target_b00<<std::endl;
        xdbg<<b2.vec() * target_b00 / b2(0)<<std::endl;
#endif
    }

    //
    // Update deconvolved sigma
    //
    sigma = pow(sigmaObs,2) - sigpsq;
    xdbg<<"sigmaP^2 = "<<sigpsq<<std::endl;
    xdbg<<"sigma^2 = sigmaObs^2 - sigmap^2 = "<<sigma<<std::endl;
    if (sigma < 0.1*sigpsq) sigma = 0.1*sigpsq;
    sigma = sqrt(sigma);
    xdbg<<"sigma = "<<sigma<<std::endl;

    //
    // Update pixels
    //
    pix[0].clear();
    getPixList(im,pix[0],cen,sky,noise,gain,weightIm,trans,
               galAp,xOffset,yOffset,flag);
    npix = pix[0].size();
    xdbg<<"npix = "<<npix<<std::endl;
    if (npix < 10) {
        dbg<<"Too few pixels to continue: "<<npix<<std::endl;
        xdbg<<"flag LT10PIX\n";
        flag |= LT10PIX; // Should already be set by getPixList, but...
        xdbg<<"flag SHAPELET_NOT_DECONV\n";
        flag |= SHAPELET_NOT_DECONV;
        if (!fixCen) cen += ell.getCen();
        return;
    }
    go = galOrder;
    galSize = (go+1)*(go+2)/2;
    if (npix <= galSize) {
        while (npix <= galSize) { galSize -= go+1; --go; }
        dbg<<"Too few pixels ("<<npix<<") for given gal_order. \n";
        dbg<<"Reduced gal_order to "<<go<<" gal_size = "<<galSize<<std::endl;
        xdbg<<"flag SHAPE_REDUCED_ORDER\n";
        flag |= SHAPE_REDUCED_ORDER;
    }

    //
    // Adjust sigma based on deconvolved measurement.
    //
    flag1 = 0;
    shapelet.setSigma(sigma);
    ell.measureShapelet(pix,psf,shapelet,go,galOrder2);
#if 1
    double target_b00 = 12.39869805;
    xdbg<<"Rescaled to have b00 = "<<target_b00<<std::endl;
    xdbg<<shapelet.vec() * target_b00 / shapelet(0)<<std::endl;
#endif
    if (ell.findRoundFrame(shapelet,galOrder2,flag1)) {
        ++log._nsMu;
        if (flag1) {
            Assert(
                (flag1 & ~(SHEAR_LOCAL_MIN | SHEAR_POOR_FIT | SHEAR_BAD_FLUX)) 
                == false);
            // This is just bumps the Ellipse flags up two spots
            // from the SHEAR_* to SHAPE_* flags.
            flag1 <<= 2;
            Assert(
                (flag1 & ~(SHAPE_LOCAL_MIN | SHAPE_POOR_FIT | SHAPE_BAD_FLUX)) 
                == false);
            if (flag1 & SHAPE_LOCAL_MIN)
                xdbg<<"flag SHAPE_LOCAL_MIN\n";
            if (flag1 & SHAPE_POOR_FIT)
                xdbg<<"flag SHAPE_POOR_FIT\n";
            if (flag1 & SHAPE_BAD_FLUX)
                xdbg<<"flag SHAPE_BAD_FLUX\n";
            flag |= flag1;
        }
        dbg<<"Successful deconvolving fit:\n";
        xdbg<<"Mu = "<<ell.getMu()<<std::endl;
    } else {
        ++log._nfMu;
        dbg<<"Deconvolving measurement failed\n";
        xdbg<<"flag DECONV_FAILED\n";
        flag |= DECONV_FAILED;
        if (!fixCen) cen += ell.getCen();
        return;
    }
    sigma *= exp(ell.getMu());
    ell.setMu(0.);

    // Check - mu should now be close to zero:
    if (XDEBUG) {
        BVec b2 = shapelet;
        b2.setSigma(sigma);
        Ellipse e2 = ell;
        long flag2 = 0;
        e2.measureShapelet(pix,b2,go,galOrder2);
        e2.findRoundFrame(b2,galOrder2,flag2);
        xdbg<<"After native fit #2:\n";
        xdbg<<"Mu = "<<e2.getMu()<<std::endl;
        xdbg<<"(Should be close to zero, but probably not exact.)\n";
    }

    //
    // Remeasure the shape.
    //
    shapelet.setSigma(sigma);
    std::complex<double> gale = 0.;
    ell.measureShapelet(pix,psf,shapelet,go,galOrder2);
    dbg<<"ell.cen = "<<ell.getCen()<<std::endl;
    dbg<<"ell.gamma = "<<ell.getGamma()<<std::endl;
    dbg<<"ell.mu = "<<ell.getMu()<<std::endl;
    dbg<<"Measured deconvolved b_gal = "<<shapelet.vec()<<std::endl;

    //
    // Also measure the isotropic significance
    // TODO: Should this be in the reference frame where galaxy is round?
    //
    BVec flux(0,sigma);
    DMatrix fluxCov(1,1);
    ell.measureShapelet(pix,psf,flux,0,0,&fluxCov);
    nu = flux(0) / sqrt(fluxCov(0,0));
    dbg<<"nu = "<<flux(0)<<" / sqrt("<<fluxCov(0,0)<<") = "<<nu<<std::endl;
    // If the b00 value in the shapelet doesn't match the direct flux
    // measurement, set a flag.
    if (!(flux(0) > 0.0 &&
          shapelet(0) >= flux(0)/3. &&
          shapelet(0) <= flux(0)*3.)) {
        dbg<<"Bad flux value: \n";
        dbg<<"flux = "<<flux(0)<<std::endl;
        dbg<<"shapelet = "<<shapelet.vec()<<std::endl;
        xdbg<<"flag SHAPE_BAD_FLUX\n";
        flag |= SHAPE_BAD_FLUX;
    }


    //
    // Next, we calculate the shear for which this deconvolved 
    // measurement looks round.
    //
    ell.unfixGam();
    if (ell.findRoundFrame(shapelet,galOrder2,flag)) {
        dbg<<"Successful Gamma fit\n";
        xdbg<<"Measured gamma = "<<ell.getGamma()<<std::endl;
        xdbg<<"Mu  = "<<ell.getMu()<<std::endl;
        xdbg<<"Cen  = "<<ell.getCen()<<std::endl;
    } else {
        ++log._nfGamma;
        dbg<<"Measurement failed\n";
        xdbg<<"flag SHEAR_FAILED\n";
        flag |= SHEAR_FAILED;
    }

    //
    // Now we can update the pixels to be a circle in the frame where
    // the galaxy is round.
    //
    // TODO


    //
    // Then finally, we remeasure in the round frame and tweak gamma, et al
    // given the new measurement.
    //
    sigma *= exp(ell.getMu());
    ell.setMu(0.);
    shapelet.setSigma(sigma);
    ell.measureShapelet(pix,psf,shapelet,go,galOrder2);
    dbg<<"ell.cen = "<<ell.getCen()<<std::endl;
    dbg<<"ell.gamma = "<<ell.getGamma()<<std::endl;
    dbg<<"ell.mu = "<<ell.getMu()<<std::endl;
    dbg<<"Measured deconvolved b_gal = "<<shapelet.vec()<<std::endl;

    DMatrix cov(5,5);
    if (ell.findRoundFrame(shapelet,galOrder2,flag,&cov)) {
        ++log._nsGamma;
        dbg<<"Successful Gamma fit\n";
        xdbg<<"Measured gamma = "<<ell.getGamma()<<std::endl;
        xdbg<<"Mu  = "<<ell.getMu()<<std::endl;
        xdbg<<"Cen  = "<<ell.getCen()<<std::endl;
    } else {
        ++log._nfGamma;
        dbg<<"Measurement failed\n";
        xdbg<<"flag SHEAR_FAILED\n";
        flag |= SHEAR_FAILED;
    }

    //
    // Copy the shear and covariance to the output variables
    //
    shear = ell.getGamma();
    DSmallMatrix22 cov1 = cov.TMV_subMatrix(2,4,2,4);
    if (!(cov1.TMV_det() > 0.)) {
        dbg<<"cov1 has bad determinant: "<<cov1.TMV_det()<<std::endl;
        dbg<<"cov1 = "<<cov1<<std::endl;
        dbg<<"Full cov = "<<cov<<std::endl;
        xdbg<<"flag SHEAR_BAD_COVAR\n";
        flag |= SHEAR_BAD_COVAR;
    } else {
        shearcov = cov1;
    }
    if (!fixCen) cen += ell.getCen();
}

void measureSingleShear(
    Position& cen, const Image<double>& im, double sky,
    const Transformation& trans, const FittedPsf& fitpsf,
    double noise, double gain, const Image<double>* weightIm, 
    double galAperture, double maxAperture, 
    int galOrder, int galOrder2,
    double fPsf, double minGalSize, bool fixCen,
    double xOffset, double yOffset,
    bool fixSigma, double fixSigmaValue,
    ShearLog& log, std::complex<double>& shear, 
    DSmallMatrix22& shearcov, BVec& shapelet,
    double& nu, long& flag)
{
    // Get coordinates of the galaxy, and convert to sky coordinates
    try {
        // We don't need to save skypos.  We just want to catch the range
        // error here, so we don't need to worry about it for dudx, etc.
        Position skyPos;
        trans.transform(cen,skyPos);
        dbg<<"skypos = "<<skyPos<<std::endl;
    } catch (RangeException& e) {
        dbg<<"distortion range error: \n";
        xdbg<<"p = "<<cen<<", b = "<<e.getBounds()<<std::endl;
        ++log._nfRange1;
        xdbg<<"flag TRANSFORM_EXCEPTION\n";
        flag |= TRANSFORM_EXCEPTION;
        return;
    }

    // Calculate the psf from the fitted-psf formula:
    std::vector<BVec> psf(1, BVec(fitpsf.getPsfOrder(), fitpsf.getSigma()));
    try {
        dbg<<"for fittedpsf cen = "<<cen<<std::endl;
        psf[0] = fitpsf(cen);
    } catch (RangeException& e) {
        dbg<<"fittedpsf range error: \n";
        xdbg<<"p = "<<cen<<", b = "<<e.getBounds()<<std::endl;
        ++log._nfRange2;
        xdbg<<"flag FITTEDPSF_EXCEPTION\n";
        flag |= FITTEDPSF_EXCEPTION;
        return;
    }

    // Do the real meat of the calculation:
    dbg<<"measure single shear cen = "<<cen<<std::endl;
    try {
        measureSingleShear1(
            // Input data:
            cen, im, sky, trans, psf,
            // Noise variables:
            noise, gain, weightIm, 
            // Parameters:
            galAperture, maxAperture, galOrder, galOrder2,
            fPsf, minGalSize, fixCen, xOffset, yOffset,
            fixSigma, fixSigmaValue,
            // Log information
            log,
            // Ouput values:
            shear, shearcov, shapelet, nu, flag);
#ifdef USE_TMV
    } catch (tmv::Error& e) {
        dbg<<"TMV Error thrown in MeasureSingleShear\n";
        dbg<<e<<std::endl;
        ++log._nfTmvError;
        xdbg<<"flag TMV_EXCEPTION\n";
        flag |= TMV_EXCEPTION;
#endif
    } catch (std::exception& e) {
        dbg<<"std::exception thrown in MeasureSingleShear\n";
        dbg<<e.what()<<std::endl;
        ++log._nfOtherError;
        xdbg<<"flag STD_EXCEPTION\n";
        flag |= STD_EXCEPTION;
    } catch (...) {
        dbg<<"unkown exception in MeasureSingleShear\n";
        ++log._nfOtherError;
        xdbg<<"flag UNKNOWN_EXCEPTION\n";
        flag |= UNKNOWN_EXCEPTION;
    } 

}

