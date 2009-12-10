
#include "TMV.h"
#include "dbg.h"
#include "Params.h"
#include "BVec.h"
#include "Ellipse.h"
#include "Pixel.h"
#include "Image.h"
#include "FittedPsf.h"
#include "Log.h"
#include "Params.h"

void measureSingleShear1(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans, const std::vector<BVec>& psf,
    double noise, double gain, const Image<double>* weightIm, 
    double galAperture, double maxAperture,
    int galOrder, int galOrder2,
    double fPsf, double minGalSize,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    tmv::SmallMatrix<double,2,2>& shearcov, BVec& shapelet,
    double& nu, long& flag)
{
    // Find harmonic mean of psf sizes:
    // MJ: Is it correct to use the harmonic mean of sigma^2?
    double sigmaP = 0.;
    const int nPsf = psf.size();
    Assert(nPsf > 0);
    for(int i=0;i<nPsf;++i) sigmaP += 1./psf[i].getSigma();
    sigmaP = double(nPsf) / sigmaP;

    // We start with a native fit using sigma = 2*sigmaP:
    double sigmaObs = 2.*sigmaP;
    double galAp = galAperture * sigmaObs;
    xdbg<<"galap = "<<galAperture<<" * "<<sigmaObs<<" = "<<galAp<<std::endl;
    if (maxAperture > 0. && galAp > maxAperture) 
    {
        galAp = maxAperture;
        xdbg<<"      => "<<galAp<<std::endl;
    }
    xdbg<<"sigma_obs = "<<sigmaObs<<", sigma_p = "<<sigmaP<<std::endl;

    std::vector<PixelList> pix(1);
    getPixList(im,pix[0],cen,sky,noise,gain,weightIm,trans,galAp,flag);
    int npix = pix[0].size();
    xdbg<<"npix = "<<npix<<std::endl;

    Ellipse ell;
    ell.fixGam();
    ell.crudeMeasure(pix[0],sigmaObs);
    xdbg<<"Crude Measure: centroid = "<<ell.getCen()<<", mu = "<<ell.getMu()<<std::endl;

    int go = 2;
    int galSize = (go+1)*(go+2)/2;

    if (times) ell.doTimings();
    long flag1=0;
    if (ell.measure(pix,go,sigmaObs,true,flag1)) 
    {
        if (times) times->successNative(ell.getTimes());
        ++log._nsNative;
        dbg<<"Successful native fit:\n";
        dbg<<"Z = "<<ell.getCen()<<std::endl;
        dbg<<"Mu = "<<ell.getMu()<<std::endl;
    }
    else 
    {
        if (times) times->failNative(ell.getTimes());
        ++log._nfNative;
        dbg<<"Native measurement failed\n";
        flag |= NATIVE_FAILED;
        return;
    }

    // Now we can do a deconvolving fit, but one that does not 
    // shear the coordinates.
    sigmaObs *= exp(ell.getMu());
    if (sigmaObs < minGalSize*sigmaP) 
    {
        dbg<<"skip: galaxy is too small -- "<<sigmaObs<<
            " psf size = "<<sigmaP<<std::endl;
        ++log._nfSmall;
        flag |= TOO_SMALL;
        return;
    }
    galAp *= exp(ell.getMu());
    if (maxAperture > 0. && galAp > maxAperture) galAp = maxAperture;
    ell.setMu(0);

    // Check - mu should now be zero:
    // This one works
    //ell.measure(pix,go,sigmaObs,true,flag1);
    //xdbg<<"After native fit #2:\n";
    //xdbg<<"Mu = "<<ell.getMu()<<std::endl;

    pix[0].clear();
    getPixList(im,pix[0],cen,sky,noise,gain,weightIm,trans,galAp,flag);
    npix = pix[0].size();
    xdbg<<"npix = "<<npix<<std::endl;
    double sigpsq = pow(sigmaP,2);
    double sigma = pow(sigmaObs,2) - sigpsq;
    if (sigma < 0.1*sigpsq) sigma = 0.1*sigpsq;
    sigma = sqrt(sigma);
    xdbg<<"sigma = "<<sigma<<std::endl;
    ell.setFP(fPsf);
    if (times) ell.resetTimes();
    flag1 = 0;
    if (ell.measure(pix,psf,go,sigma,false,flag1)) 
    {
        if (times) times->successMu(ell.getTimes());
        ++log._nsMu;
        if (flag1) 
        {
            Assert((flag1 & ~(SHEAR_LOCAL_MIN | SHEAR_POOR_FIT)) == false);
            // This is just bumps the Ellipse flags up two spots
            // from the SHEAR_* to SHAPE_* flags.
            flag1 <<= 2;
            Assert((flag1 & ~(SHAPE_LOCAL_MIN | SHAPE_POOR_FIT)) == false);
            flag |= flag1;
        }
        dbg<<"Successful deconvolving fit:\n";
        xdbg<<"Mu = "<<ell.getMu()<<std::endl;
    }
    else 
    {
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
        dbg<<"Too few pixels ("<<npix<<") for given gal_order. \n";
        dbg<<"Reduced gal_order to "<<go<<" gal_size = "<<galSize<<std::endl;
        flag |= SHAPE_REDUCED_ORDER;
        if (go < 2) { flag |= TOO_SMALL; return; }
    }
    shapelet.setSigma(sigma);
    std::complex<double> gale = 0.;
    ell.measureShapelet(pix,psf,shapelet);
    xdbg<<"Measured b_gal = "<<shapelet.vec()<<std::endl;

    // Also measure the isotropic significance
    BVec flux(0,sigma);
    tmv::Matrix<double> fluxCov(1,1);
    ell.measureShapelet(pix,psf,flux,&fluxCov);
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
        flag |= SHAPE_BAD_FLUX;
    }

    // If the above deconvolving fit failed, then return.
    if (flag & DECONV_FAILED) return;

    // Under normal circumstances, b20/b00 ~= conj(gamma)/sqrt(2)
    gale = std::complex<double>(shapelet(3),shapelet(4));
    gale /= shapelet(0);
    if (std::abs(gale) < 0.5) 
    {
        ell.setGamma(conj(gale) * sqrt(2.));
        shear = ell.getGamma();
    }

    // Finally, we calculate the shear in the deconvolved galaxy.
    //ell.fixMu();
    ell.unfixGam();
    go = galOrder2;
    galSize = (go+1)*(go+2)/2;
    if (npix <= galSize) 
    {
        while (npix <= galSize) { galSize -= go+1; --go; }
        dbg<<"Too few pixels ("<<npix<<") for given gal_order. \n";
        dbg<<"Reduced gal_order to "<<go<<" gal_size = "<<galSize<<std::endl;
    }
    if (times) ell.resetTimes();
    tmv::Matrix<double> cov(5,5);
    if (ell.measure(pix,psf,go,sigma,false,flag,&cov)) 
    {
        if (times) times->successGamma(ell.getTimes());
        ++log._nsGamma;
        dbg<<"Successful Gamma fit\n";
        xdbg<<"Measured gamma = "<<ell.getGamma()<<std::endl;
        shear = ell.getGamma();
        tmv::SmallMatrix<double,2,2,tmv::RowMajor> cov1 = 
            cov.SubMatrix(2,4,2,4);
        if (!(cov1.Det() > 0.)) flag |= SHEAR_BAD_COVAR;
        else shearcov = cov1;
    }
    else 
    {
        if (times) times->failGamma(ell.getTimes());
        ++log._nfGamma;
        dbg<<"Measurement failed\n";
        flag |= SHEAR_FAILED;
        shear = ell.getGamma();
        tmv::SmallMatrix<double,2,2,tmv::RowMajor> cov1 = 
            cov.SubMatrix(2,4,2,4);
        if (!(cov1.Det() > 0.)) flag |= SHEAR_BAD_COVAR;
        else shearcov = cov1;
        return;
    }
}

void measureSingleShear(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans, const FittedPsf& fitpsf,
    double noise, double gain, const Image<double>* weightIm, 
    double galAperture, double maxAperture,
    int galOrder, int galOrder2,
    double fPsf, double minGalSize,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    tmv::SmallMatrix<double,2,2>& shearcov, BVec& shapelet,
    double& nu, long& flag)
{
    // Get coordinates of the galaxy, and convert to sky coordinates
    try 
    {
        // We don't need to save skypos.  We just want to catch the range
        // error here, so we don't need to worry about it for dudx, etc.
        Position skyPos;
        trans.transform(cen,skyPos);
        dbg<<"skypos = "<<skyPos<<std::endl;
    } 
    catch (RangeException& e) 
    {
        dbg<<"distortion range error: \n";
        xdbg<<"p = "<<cen<<", b = "<<e.getBounds()<<std::endl;
        if (times) ++times->_nfRange1;
        ++log._nfRange1;
        flag |= TRANSFORM_EXCEPTION;
        return;
    }

    // Calculate the psf from the fitted-psf formula:
    std::vector<BVec> psf(1,
                          BVec(fitpsf.getPsfOrder(), fitpsf.getSigma()));
    try 
    {
        dbg<<"for fittedpsf cen = "<<cen<<std::endl;
        psf[0] = fitpsf(cen);
    } 
    catch (RangeException& e) 
    {
        dbg<<"fittedpsf range error: \n";
        xdbg<<"p = "<<cen<<", b = "<<e.getBounds()<<std::endl;
        if (times) ++times->_nfRange2;
        ++log._nfRange2;
        flag |= FITTEDPSF_EXCEPTION;
        return;
    }

    // Do the real meat of the calculation:
    dbg<<"measure single shear cen = "<<cen<<std::endl;
    try 
    {
        measureSingleShear1(
            // Input data:
            cen, im, sky, trans, psf,
            // Noise variables:
            noise, gain, weightIm, 
            // Parameters:
            galAperture, maxAperture, galOrder, galOrder2, 
            fPsf, minGalSize, 
            // Time stats if desired:
            times,
            // Log information
            log,
            // Ouput values:
            shear, shearcov, shapelet, nu, flag);
    }
    catch (tmv::Error& e) 
    {
        dbg<<"TMV Error thrown in MeasureSingleShear\n";
        dbg<<e<<std::endl;
        ++log._nfTmvError;
        flag |= TMV_EXCEPTION;
    } 
    catch (...) 
    {
        dbg<<"unkown exception in MeasureSingleShear\n";
        ++log._nfOtherError;
        flag |= UNKNOWN_EXCEPTION;
    } 

}

