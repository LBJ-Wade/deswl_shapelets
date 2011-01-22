
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

#define MAX_ITER 3
#define MAX_DELTA_MU 0.2
#define MAX_DELTA_GAMMA 0.2

#if 0
static Position AddOffset(
    Position cen, std::complex<double> offset, const Transformation& trans)
{
    // ( u' ) = D ( x ) + ( u_off )
    // ( v' )     ( y )   ( v_off )
    //        = D ( x' )
    //            ( y' )
    // D ( x' ) = D ( x ) + ( u_off )
    //   ( y' )     ( y )   ( v_off )
    //   ( x' ) = ( x ) + D^-1 ( u_off )
    //   ( y' )   ( y )        ( v_off )
    DSmallMatrix22 localD;
    trans.getDistortion(cen,localD);
    tmv::SmallVector<double,2> xy;
    xy << real(offset) << imag(offset);
    xy = localD.inverse() * xy;
    Position ret(cen.getX() + xy(0), cen.getY() + xy(1));
    return ret;
}
#endif

void DoMeasureShear(
    const std::vector<PixelList>& allpix,
    const std::vector<BVec>& psf,
    double galAperture, double maxAperture,
    int galOrder, int galOrder2, int maxm,
    double minFPsf, double maxFPsf, double minGalSize, bool fixCen,
    bool fixSigma, double fixSigmaValue, bool nativeOnly,
    ShearLog& log, BVec& shapelet, 
    std::complex<double>& gamma, DSmallMatrix22& cov,
    double& nu, long& flag)
{
    try {
        dbg<<"Start MeasureShear\n";
        dbg<<"allpix.size = "<<allpix.size()<<std::endl;
        const int nExp = allpix.size();
        for(int i=0;i<nExp;++i)
            dbg<<"allpix["<<i<<"].size = "<<allpix[i].size()<<std::endl;
        dbg<<"psf.size = "<<psf.size()<<std::endl;
        Assert(psf.size() == allpix.size());

        // Find harmonic mean of psf sizes:
        // MJ: Is it better to use the harmonic mean of sigma or sigma^2?
        //     (Or something else entirely?)
        double sigmaP = 0.;
        const int nPsf = psf.size();
        Assert(nPsf > 0);
        for(int i=0;i<nPsf;++i) {
            sigmaP += 1./psf[i].getSigma();
            dbg<<"psf["<<i<<"] sigma = "<<psf[i].getSigma()<<std::endl;
        }
        sigmaP = double(nPsf) / sigmaP;
        dbg<<"Harmonic mean = "<<sigmaP<<std::endl;
        long flag1=0;
        std::complex<double> cen_offset = 0.;

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
        dbg<<"galap = "<<galAperture<<" * "<<sigmaObs<<" = "<<galAp<<std::endl;
        if (maxAperture > 0. && galAp > maxAperture) {
            galAp = maxAperture;
            dbg<<"      => "<<galAp<<std::endl;
        }
        dbg<<"sigma_obs = "<<sigmaObs<<", sigma_p = "<<sigmaP<<std::endl;

        //
        // Load pixels from main PixelLists.
        //
        std::vector<PixelList> pix(nExp);
        int npix = 0;
        for(int i=0;i<nExp;++i) {
            getSubPixList(pix[i],allpix[i],cen_offset,galAp,flag1);
            npix += pix[i].size();
        }
        dbg<<"npix = "<<npix<<std::endl;
        if (npix < 10) {
            dbg<<"Too few pixels to continue: "<<npix<<std::endl;
            dbg<<"FLAG LT10PIX\n";
            flag |= LT10PIX;
            dbg<<"FLAG SHAPELET_NOT_DECONV\n";
            flag |= SHAPELET_NOT_DECONV;
            return;
        }

        //
        // Do a crude measurement based on simple pixel sums.
        // TODO: If we have single-epoch measurements, use those for the
        //       starting guess rather than crudeMeasure.
        //
        Ellipse ell_init;
        if (fixCen) ell_init.fixCen();
        if (fixSigma) ell_init.fixMu();
        ell_init.fixGam();
        if (!fixCen || !fixSigma) {
            ell_init.crudeMeasure(pix[0],sigmaObs);
            xdbg<<"Crude Measure: centroid = "<<ell_init.getCen()<<
                ", mu = "<<ell_init.getMu()<<std::endl;
            sigmaObs *= exp(ell_init.getMu());
            ell_init.setMu(0.);
        }

        //
        // Correct the centroid first.
        //
        Ellipse ell_native = ell_init;
        if (fixCen) ell_native.fixCen();
        ell_native.fixMu();
        ell_native.fixGam();
        if (!fixCen) {
            flag1 = 0;
            if (!ell_native.measure(pix,2,2,2,sigmaObs,flag1,0.1)) {
                ++log._nfCentroid;
                dbg<<"First centroid pass failed.\n";
                flag |= flag1;
                dbg<<"FLAG CENTROID_FAILED\n";
                flag |= CENTROID_FAILED;
                dbg<<"FLAG SHAPELET_NOT_DECONV\n";
                flag |= SHAPELET_NOT_DECONV;
                return;
            }
            dbg<<"After first centroid pass: cen = "<<ell_native.getCen()<<std::endl;

            // redo the pixlist:
            cen_offset += ell_native.getCen();
            ell_native.setCen(0.);
            dbg<<"New center = "<<cen_offset<<std::endl;
            npix = 0;
            for(int i=0;i<nExp;++i) {
                getSubPixList(pix[i],allpix[i],cen_offset,galAp,flag1);
                npix += pix[i].size();
            }
            if (npix < 10) {
                dbg<<"Too few pixels to continue: "<<npix<<std::endl;
                dbg<<"FLAG LT10PIX\n";
                flag |= LT10PIX;
                dbg<<"FLAG SHAPELET_NOT_DECONV\n";
                flag |= SHAPELET_NOT_DECONV;
                return;
            }

            // The centroid doesn't necessarily get all the way to the perfect
            // centroid, so if the input guess is consistently biased in some 
            // direction (e.g. by a different convention for the (0,0) or (1,1) 
            // location), then this results in a bias in the shear values.  
            // To avoid that, we first get to what we think it the correct centroid.
            // Then we redo the aperture from this point.  Then we shift the 
            // centroid value by 0.1 arcsec in a random direction, and resolve for
            // the centroid.  This is all pretty fast, so it's worth doing to
            // avoid the subtle bias that can otherwise result.
            // get an offset of 0.1 arcsec in a random direction.
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
            ell_native.setCen(std::complex<double>(x_offset,y_offset));
            dbg<<"After random offset "<<x_offset<<","<<y_offset<<
                ": cen = "<<ell_native.getCen()<<std::endl;

            // Now do it again.
            flag1 = 0;
            if (!ell_native.measure(pix,2,2,2,sigmaObs,flag1,0.1)) {
                ++log._nfCentroid;
                dbg<<"Second centroid pass failed.\n";
                flag |= flag1;
                dbg<<"FLAG CENTROID_FAILED\n";
                flag |= CENTROID_FAILED;
                dbg<<"FLAG SHAPELET_NOT_DECONV\n";
                flag |= SHAPELET_NOT_DECONV;
                return;
            }
            dbg<<"After second centroid pass: cen = "<<
                ell_native.getCen()<<" relative to the new center: "<<
                cen_offset<<std::endl;

            ++log._nsCentroid;
        }

        //
        // Next find a good sigma value 
        //
        if (!fixSigma) {
            for(int iter=1;iter<=MAX_ITER;++iter) {
                dbg<<"Mu iter = "<<iter<<std::endl;
                flag1 = 0;
                ell_native.unfixMu();
                if (ell_native.measure(pix,2,2,2,sigmaObs,flag1,0.1)) {
                    // Don't ++nsNative yet.  Only success if also do round frame.
                    dbg<<"Successful native fit:\n";
                    dbg<<"Z = "<<ell_native.getCen()<<std::endl;
                    dbg<<"Mu = "<<ell_native.getMu()<<std::endl;
                } else {
                    ++log._nfNative;
                    dbg<<"Native measurement failed\n";
                    flag |= flag1;
                    dbg<<"FLAG NATIVE_FAILED\n";
                    flag |= NATIVE_FAILED;
                    dbg<<"FLAG SHAPELET_NOT_DECONV\n";
                    flag |= SHAPELET_NOT_DECONV;
                    return;
                }

                // Adjust the sigma value so we can reset mu = 0.
                double deltaMu = ell_native.getMu();
                cen_offset += ell_native.getCen();
                dbg<<"New center = "<<cen_offset<<std::endl;
                ell_native.setCen(0.);
                sigmaObs *= exp(ell_native.getMu());
                dbg<<"New sigmaObs = "<<sigmaObs<<std::endl;
                ell_native.setMu(0.);
                if (sigmaObs < minGalSize*sigmaP) {
                    dbg<<"skip: galaxy is too small -- "<<sigmaObs<<
                        " psf size = "<<sigmaP<<std::endl;
                    ++log._nfSmall;
                    flag |= flag1;
                    dbg<<"FLAG TOO_SMALL\n";
                    flag |= TOO_SMALL;
                    dbg<<"FLAG SHAPELET_NOT_DECONV\n";
                    flag |= SHAPELET_NOT_DECONV;
                    return;
                }

                // redo the pixlist:
                galAp = sigmaObs * galAperture;
                if (maxAperture > 0. && galAp > maxAperture) 
                    galAp = maxAperture;
                npix = 0;
                for(int i=0;i<nExp;++i) {
                    getSubPixList(pix[i],allpix[i],cen_offset,galAp,flag1);
                    npix += pix[i].size();
                }
                if (npix < 10) {
                    dbg<<"Too few pixels to continue: "<<npix<<std::endl;
                    dbg<<"FLAG LT10PIX\n";
                    flag |= LT10PIX;
                    dbg<<"FLAG SHAPELET_NOT_DECONV\n";
                    flag |= SHAPELET_NOT_DECONV;
                    return;
                }
                dbg<<"deltaMu = "<<deltaMu<<std::endl;
                if (std::abs(deltaMu) < MAX_DELTA_MU) {
                    dbg<<"deltaMu < "<<MAX_DELTA_MU<<std::endl;
                    break;
                } else if (iter < MAX_ITER-1) {
                    dbg<<"deltaMu >= "<<MAX_DELTA_MU<<std::endl;
                    continue;
                } else {
                    dbg<<"deltaMu >= "<<MAX_DELTA_MU<<std::endl;
                    dbg<<"But iter == "<<MAX_ITER<<", so stop.\n";
                }
            }
        }

        //
        // Next find the frame in which the native observation is round.
        //
        Ellipse ell_round = ell_native;
        if (fixCen) ell_round.fixCen();
        if (fixSigma) ell_round.fixMu();
        std::complex<double> gamma_prev = 0.;
        for(int iter=1;iter<=MAX_ITER;++iter) {
            dbg<<"Round iter = "<<iter<<std::endl;
            flag1 = 0;
            if (ell_round.measure(
                    pix,galOrder,galOrder2,maxm,sigmaObs,flag1,1.e-3)) {
                ++log._nsNative;
                dbg<<"Successful round fit:\n";
                dbg<<"Z = "<<ell_round.getCen()<<std::endl;
                dbg<<"Mu = "<<ell_round.getMu()<<std::endl;
                dbg<<"Gamma = "<<ell_round.getGamma()<<std::endl;
            } else {
                ++log._nfNative;
                dbg<<"Round measurement failed\n";
                flag |= flag1;
                dbg<<"FLAG NATIVE_FAILED\n";
                flag |= NATIVE_FAILED;
                dbg<<"FLAG SHAPELET_NOT_DECONV\n";
                flag |= SHAPELET_NOT_DECONV;
                return;
            }

            // Adjust the sigma value so we can reset mu = 0.
            if (!fixCen) cen_offset += ell_round.getCen();
            ell_round.setCen(0.);
            double deltaMu = ell_round.getMu();
            sigmaObs *= exp(ell_round.getMu());
            ell_round.setMu(0.);
            if (sigmaObs < minGalSize*sigmaP) {
                dbg<<"skip: galaxy is too small -- "<<sigmaObs<<
                    " psf size = "<<sigmaP<<std::endl;
                ++log._nfSmall;
                dbg<<"FLAG TOO_SMALL\n";
                flag |= TOO_SMALL;
                dbg<<"FLAG SHAPELET_NOT_DECONV\n";
                flag |= SHAPELET_NOT_DECONV;
                return;
            }
            ell_round.removeRotation();
            gamma = ell_round.getGamma();
            dbg<<"New center = "<<cen_offset<<std::endl;
            dbg<<"New sigmaObs = "<<sigmaObs<<std::endl;
            dbg<<"New gamma = "<<gamma<<std::endl;
            std::complex<double> deltaGamma = gamma - gamma_prev;
            gamma_prev = gamma;
            dbg<<"ell_round = "<<ell_round<<std::endl;

            // Get the pixels in an elliptical aperture based on the
            // observed shape.
            galAp = sigmaObs * galAperture;
            if (maxAperture > 0. && galAp > maxAperture) 
                galAp = maxAperture;
            npix = 0;
            for(int i=0;i<nExp;++i) {
                getSubPixList(pix[i],allpix[i],cen_offset,gamma,galAp,flag1);
                npix += pix[i].size();
            }
            dbg<<"npix = "<<npix<<std::endl;
            if (npix < 10) {
                dbg<<"Too few pixels to continue: "<<npix<<std::endl;
                dbg<<"FLAG LT10PIX\n";
                flag |= LT10PIX;
                dbg<<"FLAG SHAPELET_NOT_DECONV\n";
                flag |= SHAPELET_NOT_DECONV;
                return;
            }
            dbg<<"deltaMu = "<<deltaMu<<std::endl;
            dbg<<"deltaGamma = "<<deltaGamma<<std::endl;
            if (std::abs(deltaMu) < MAX_DELTA_MU &&
                std::abs(deltaGamma) < MAX_DELTA_GAMMA) {
                dbg<<"deltaMu < "<<MAX_DELTA_MU;
                dbg<<" and deltaGamma < "<<MAX_DELTA_GAMMA<<std::endl;
                break;
            } else if (iter < MAX_ITER-1) {
                dbg<<"deltaMu >= "<<MAX_DELTA_MU;
                dbg<<" or deltaGamma >= "<<MAX_DELTA_GAMMA<<std::endl;
                continue;
            } else {
                dbg<<"deltaMu >= "<<MAX_DELTA_MU;
                dbg<<" or deltaGamma >= "<<MAX_DELTA_GAMMA<<std::endl;
                dbg<<"But iter == "<<MAX_ITER<<", so stop.\n";
            }
        }

        if (nativeOnly) return;

        // Start with the specified fPsf, but allow it to increase up to
        // maxFPsf if there are any problems.
        long flag0 = flag;
        for(double fPsf=minFPsf; fPsf<=maxFPsf+1.e-3; fPsf+=0.5) {
            flag = flag0; // In case anything was set in a previous iteration.

            bool lastfpsf = (fPsf + 0.5 > maxFPsf+1.e-3);
            dbg<<"Start fPsf loop: fPsf = "<<fPsf<<std::endl;

            //
            // Set the sigma to use for the shapelet measurement:
            //
            // The sigma to use is 
            // sigma_s^2 = sigma_i^2 + fpsf sigma_p^2
            // And sigma_i^2 = sigma_o^2 - sigma_p^2
            // So sigma_s^2 = sigma_o^2 + (fpsf-1) sigma_p^2
            sigma = pow(sigmaObs,2) + (fPsf-1.) * sigpsq;
            xdbg<<"sigmaP^2 = "<<sigpsq<<std::endl;
            xdbg<<"sigmaObs^2 = "<<pow(sigmaObs,2)<<std::endl;
            xdbg<<"sigma^2 = sigmaObs^2 + (fPsf-1) * sigmap^2 = "<<sigma<<std::endl;
            if (sigma < 0.1*sigpsq) {
                xdbg<<"sigma too small, try larger fPsf\n";
                if (!lastfpsf) continue;
                dbg<<"FLAG TOO_SMALL\n";
                flag |= TOO_SMALL;
                return;
            }
            sigma = sqrt(sigma);
            dbg<<"sigma_s = "<<sigma<<std::endl;

            // 
            // Measure a deconvolving fit in the native frame.
            //
            shapelet.setSigma(sigma);
            if (ell_native.measureShapelet(
                    pix,psf,shapelet,galOrder,galOrder2,galOrder)) {
                dbg<<"Successful deconvolving fit:\n";
                ++log._nsMu;
            } else {
                dbg<<"Deconvolving measurement failed\n";
                if (!lastfpsf) continue;
                ++log._nfMu;
                dbg<<"FLAG DECONV_FAILED\n";
                flag |= DECONV_FAILED;
            }
            dbg<<"Measured deconvolved b_gal = "<<shapelet.vec()<<std::endl;

            //
            // Also measure the isotropic significance
            // TODO: Should this be in the frame where galaxy is round?
            //       Probably doesn't matter.
            //
            BVec flux(0,sigma);
            DMatrix fluxCov(1,1,0.);
            if (!ell_native.measureShapelet(pix,psf,flux,0,0,0,&fluxCov) ||
                !(flux(0) > 0) || !(fluxCov(0,0) > 0.) ||
                shapelet(0) >= flux(0)*3. || shapelet(0) <= flux(0)/3.) {
                // If the b00 value in the shapelet doesn't match the direct flux
                // measurement, set a flag.
                dbg<<"Bad flux value: \n";
                dbg<<"flux = "<<flux(0)<<std::endl;
                dbg<<"fluxCov = "<<fluxCov(0,0)<<std::endl;
                dbg<<"shapelet = "<<shapelet.vec()<<std::endl;
                dbg<<"flux ratio = "<<flux(0)/shapelet(0)<<std::endl;
                if (!lastfpsf) {
                    --log._nsMu;
                    continue;
                }
                dbg<<"FLAG SHAPE_BAD_FLUX\n";
                flag |= SHAPE_BAD_FLUX;
            }
            if (flux(0) > 0. && fluxCov(0,0) > 0.) {
                nu = flux(0) / sqrt(fluxCov(0,0));
                dbg<<"nu = "<<flux(0)<<" / sqrt("<<fluxCov(0,0)<<") = "<<nu<<std::endl;
            } else {
                nu = DEFVALNEG;
                dbg<<"nu set to error value = "<<nu<<std::endl;
            }

            //
            // Next, we find the frame where the deconvolved shape of the
            // galaxy is round.
            //
            for(int tryOrder=galOrder; tryOrder>=2; --tryOrder) {
                dbg<<"tryOrder = "<<tryOrder<<std::endl;
                if (tryOrder < galOrder) {
                    dbg<<"FLAG SHEAR_REDUCED_ORDER\n";
                    flag |= SHEAR_REDUCED_ORDER;
                }
                BVec shearedShape(tryOrder,sigma);
                if (maxm > tryOrder) maxm = tryOrder;

                if (XDEBUG) {
                    if (ell_round.measureShapelet(
                            pix,psf,shearedShape,tryOrder,galOrder2,tryOrder)) {
                        xdbg<<"Successful sheared deconvolving fit:\n";
                    } else {
                        xdbg<<"Sheared deconvolving measurement failed\n";
                    }
                    xdbg<<"Full m: shearedShape = "<<shearedShape<<std::endl;
                }

                if (ell_round.measureShapelet(
                        pix,psf,shearedShape,tryOrder,galOrder2,maxm)) {
                    dbg<<"Successful sheared deconvolving fit:\n";
                } else {
                    dbg<<"Sheared deconvolving measurement failed\n";
                    continue;
                }
                dbg<<"m <= "<<maxm<<": shearedShape = "<<shearedShape<<std::endl;

                Ellipse ell_shear = ell_round;
                if (fixCen) ell_shear.fixCen();
                if (fixSigma) ell_shear.fixMu();
                DMatrix cov5(5,5);
                long flag1=0;
                if (ell_shear.findRoundFrame(
                        shearedShape, true, galOrder2, 1.e-4, 
                        flag1, &cov5) && flag1==0) {
                    ++log._nsGamma;
                    ell_shear.removeRotation();
                    gamma = ell_shear.getGamma();
                    dbg<<"Successful shear measurement: "<<gamma<<std::endl;
                    cov = cov5.TMV_subMatrix(2,4,2,4);
                    if (!(cov.TMV_det() > 0.)) {
                        dbg<<"cov has bad determinant: "<<
                            cov.TMV_det()<<std::endl;
                        dbg<<"cov = "<<cov<<std::endl;
                        dbg<<"Full cov = "<<cov5<<std::endl;
                        dbg<<"FLAG SHEAR_BAD_COVAR\n";
                        flag |= SHEAR_BAD_COVAR;
                    }
                    return;
                } else {
                    ell_shear.removeRotation();
                    gamma = ell_shear.getGamma();
                    dbg<<"Unsuccessful shear measurement\n"; 
                    dbg<<"shear = "<<gamma<<std::endl;
                    if (tryOrder == 2) {
                        flag |= flag1;
                        dbg<<"flag = "<<flag<<std::endl;
                    }
                }
            }
            if (!lastfpsf) continue;
            ++log._nfGamma;
            dbg<<"FLAG SHEAR_FAILED\n";
            flag |= SHEAR_FAILED;
        }
#ifdef USE_TMV
    } catch (tmv::Error& e) {
        dbg<<"TMV Error thrown in MeasureShapes\n";
        dbg<<e<<std::endl;
        ++log._nfTmvError;
        dbg<<"FLAG TMV_EXCEPTION\n";
        flag |= TMV_EXCEPTION;
#endif
    } catch (std::exception& e) {
        dbg<<"std::exception thrown in MeasureShapes\n";
        dbg<<e.what()<<std::endl;
        ++log._nfOtherError;
        dbg<<"FLAG STD_EXCEPTION\n";
        flag |= STD_EXCEPTION;
    } catch (...) {
        dbg<<"unkown exception in MeasureShapes\n";
        ++log._nfOtherError;
        dbg<<"FLAG UNKNOWN_EXCEPTION\n";
        flag |= UNKNOWN_EXCEPTION;
    } 

}

