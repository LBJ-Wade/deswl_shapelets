
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

#define MAX_ITER 4
#define MAX_DELTA_MU 0.2
#define MAX_DELTA_GAMMA1 0.1
#define MAX_DELTA_GAMMA2 1.e-3
#define FAIL_DELTA_GAMMA2 1.e-1

void MeasureSingleShear(
    const std::vector<PixelList>& allpix,
    const std::vector<BVec>& psf,
    int& galorder, const ConfigFile& params,
    ShearLog& log, BVec& shapelet, 
    std::complex<double>& gamma, DSmallMatrix22& cov,
    double& nu, long& flag)
{
    double gal_aperture = params.read("shear_aperture",3.);
    double max_aperture = params.read("shear_max_aperture",0.);
    // Initial value, but also returned with actual final value.
    int galorder_init = params.read("shear_gal_order",6);
    galorder = galorder_init;
    int galorder2 = params.read("shear_gal_order2",20);
    int maxm = params.read("shear_maxm",galorder);
    int min_galorder = params.read("shear_min_gal_order",4);
    double min_fpsf = params.read("shear_f_psf",1.);
    double max_fpsf = params.read("shear_max_f_psf",min_fpsf);
    double min_galsize = params.read("shear_min_gal_size",0.);
    bool fixcen = params.read("shear_fix_centroid",false);
    bool fixsigma = params.keyExists("shear_force_sigma");
    double fixsigma_value = params.read("shear_force_sigma",0.);
    bool use_fake_pixels = params.read("shear_use_fake_pixels",false);
    double shear_inner_fake_aperture = params.read("shear_inner_fake_aperture",gal_aperture);
    double shear_outer_fake_aperture = params.read("shear_outer_fake_aperture",1.e100);
    double inner_fake_ap = 0.;
    double outer_fake_ap = 0.;

    try {
        dbg<<"Start MeasureSingleShear\n";
        dbg<<"allpix.size = "<<allpix.size()<<std::endl;
        const int nexp = allpix.size();
        for(int i=0;i<nexp;++i)
            dbg<<"allpix["<<i<<"].size = "<<allpix[i].size()<<std::endl;
        dbg<<"psf.size = "<<psf.size()<<std::endl;
        Assert(psf.size() == allpix.size());

        // Find harmonic mean of psf sizes:
        // MJ: Is it better to use the harmonic mean of sigma or sigma^2?
        //     (Or something else entirely?)
        double sigma_p = 0.;
        const int npsf = psf.size();
        Assert(npsf > 0);
        for(int i=0;i<npsf;++i) {
            sigma_p += 1./psf[i].getSigma();
            dbg<<"psf["<<i<<"] sigma = "<<psf[i].getSigma()<<std::endl;
        }
        sigma_p = double(npsf) / sigma_p;
        dbg<<"Harmonic mean = "<<sigma_p<<std::endl;
        long flag1=0;
        std::complex<double> cen_offset = 0.;

        //
        // Define initial sigma and aperture.
        // We start with a native fit using sigma = 2*sigma_p:
        // Or if we are fixed, then sigma_obs^2 = sigma_p^2 + sigma^2
        //
        double sigma = fixsigma_value;
        double sigpsq = pow(sigma_p,2);
        double sigma_obs = 
            fixsigma ? 
            sqrt(sigpsq + sigma*sigma) : 
            2.*sigma_p;
        double galap = gal_aperture * sigma_obs;
        dbg<<"galap = "<<gal_aperture<<" * "<<sigma_obs<<" = "<<galap<<std::endl;
        if (max_aperture > 0. && galap > max_aperture) {
            galap = max_aperture;
            dbg<<"      => "<<galap<<std::endl;
        }
        dbg<<"sigma_obs = "<<sigma_obs<<", sigma_p = "<<sigma_p<<std::endl;

        //
        // Load pixels from main PixelLists.
        //
        std::vector<PixelList> pix(nexp);
        int npix = 0;
        for(int i=0;i<nexp;++i) {
            if (use_fake_pixels) {
                inner_fake_ap = shear_inner_fake_aperture * sigma_obs;
                outer_fake_ap = shear_outer_fake_aperture * sigma_obs;
            }
            GetSubPixList(pix[i],allpix[i],cen_offset,0.,galap,
                          inner_fake_ap,outer_fake_ap,params,flag1);
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
        if (fixcen) ell_init.fixCen();
        if (fixsigma) ell_init.fixMu();
        ell_init.fixGam();
        if (!fixcen || !fixsigma) {
            ell_init.crudeMeasure(pix[0],sigma_obs);
            xdbg<<"Crude Measure: centroid = "<<ell_init.getCen()<<
                ", mu = "<<ell_init.getMu()<<std::endl;
            sigma_obs *= exp(ell_init.getMu());
            ell_init.setMu(0.);
        }

        //
        // Correct the centroid first.
        //
        Ellipse ell_native = ell_init;
        if (fixcen) ell_native.fixCen();
        ell_native.fixMu();
        ell_native.fixGam();
        if (!fixcen) {
            flag1 = 0;
            if (!ell_native.measure(pix,2,2,2,sigma_obs,flag1,0.1)) {
                ++log._nf_centroid;
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
            for(int i=0;i<nexp;++i) {
                if (use_fake_pixels) {
                    inner_fake_ap = shear_inner_fake_aperture * sigma_obs;
                    outer_fake_ap = shear_outer_fake_aperture * sigma_obs;
                }
                GetSubPixList(pix[i],allpix[i],cen_offset,0.,galap,
                              inner_fake_ap,outer_fake_ap,params,flag1);
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
            ell_native.setCen(std::complex<double>(x_offset,y_offset));
            dbg<<"After random offset "<<x_offset<<","<<y_offset<<
                ": cen = "<<ell_native.getCen()<<std::endl;

            // Now do it again.
            flag1 = 0;
            if (!ell_native.measure(pix,2,2,2,sigma_obs,flag1,0.1)) {
                ++log._nf_centroid;
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

            ++log._ns_centroid;
        }

        //
        // Next find a good sigma value 
        //
        if (!fixsigma) {
            for(int iter=1;iter<=MAX_ITER;++iter) {
                dbg<<"Mu iter = "<<iter<<std::endl;
                flag1 = 0;
                ell_native.unfixMu();
                if (ell_native.measure(pix,2,2,2,sigma_obs,flag1,0.1)) {
                    // Don't ++ns_native yet.  Only success if also do round frame.
                    dbg<<"Successful native fit:\n";
                    dbg<<"Z = "<<ell_native.getCen()<<std::endl;
                    dbg<<"Mu = "<<ell_native.getMu()<<std::endl;
                } else {
                    ++log._nf_native;
                    dbg<<"Native measurement failed\n";
                    flag |= flag1;
                    dbg<<"FLAG NATIVE_FAILED\n";
                    flag |= NATIVE_FAILED;
                    dbg<<"FLAG SHAPELET_NOT_DECONV\n";
                    flag |= SHAPELET_NOT_DECONV;
                    return;
                }

                // Adjust the sigma value so we can reset mu = 0.
                double delta_mu = ell_native.getMu();
                cen_offset += ell_native.getCen();
                dbg<<"New center = "<<cen_offset<<std::endl;
                ell_native.setCen(0.);
                sigma_obs *= exp(ell_native.getMu());
                dbg<<"New sigma_obs = "<<sigma_obs<<std::endl;
                ell_native.setMu(0.);
                if (sigma_obs < min_galsize*sigma_p) {
                    dbg<<"skip: galaxy is too small -- "<<sigma_obs<<
                        " psf size = "<<sigma_p<<std::endl;
                    ++log._nf_small;
                    flag |= flag1;
                    dbg<<"FLAG TOO_SMALL\n";
                    flag |= TOO_SMALL;
                    dbg<<"FLAG SHAPELET_NOT_DECONV\n";
                    flag |= SHAPELET_NOT_DECONV;
                    return;
                }

                // redo the pixlist:
                galap = sigma_obs * gal_aperture;
                if (max_aperture > 0. && galap > max_aperture) 
                    galap = max_aperture;
                npix = 0;
                for(int i=0;i<nexp;++i) {
                    if (use_fake_pixels) {
                        inner_fake_ap = shear_inner_fake_aperture * sigma_obs;
                        outer_fake_ap = shear_outer_fake_aperture * sigma_obs;
                    }
                    GetSubPixList(pix[i],allpix[i],cen_offset,0.,galap,
                                  inner_fake_ap,outer_fake_ap,params,flag1);
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
                dbg<<"delta_mu = "<<delta_mu<<std::endl;
                if (std::abs(delta_mu) < MAX_DELTA_MU) {
                    dbg<<"delta_mu < "<<MAX_DELTA_MU<<std::endl;
                    break;
                } else if (iter < MAX_ITER) {
                    dbg<<"delta_mu >= "<<MAX_DELTA_MU<<std::endl;
                    continue;
                } else {
                    dbg<<"delta_mu >= "<<MAX_DELTA_MU<<std::endl;
                    dbg<<"But iter == "<<MAX_ITER<<", so stop.\n";
                }
            }
        }
        

        //
        // Measure the isotropic significance
        //
        BVec flux(0,sigma);
        DMatrix flux_cov(1,1);
        if (!ell_native.measureShapelet(pix,psf,flux,0,0,0,&flux_cov) ||
            !(flux(0) > 0) || !(flux_cov(0,0) > 0.) ) {
            dbg<<"Failed flux measurement of bad flux value: \n";
            dbg<<"flux = "<<flux(0)<<std::endl;
            dbg<<"flux_cov = "<<flux_cov(0,0)<<std::endl;

            dbg<<"FLAG SHAPE_BAD_FLUX\n";
            flag |= SHAPE_BAD_FLUX;
        }
        dbg<<"flux = "<<flux<<std::endl;
        dbg<<"flux_cov = "<<flux_cov<<std::endl;
        if (flux(0) > 0. && flux_cov(0,0) > 0.) {
            nu = flux(0) / sqrt(flux_cov(0,0));
            dbg<<"nu = "<<flux(0)<<" / sqrt("<<flux_cov(0,0)<<") = "<<
                nu<<std::endl;
        } else {
            nu = DEFVALNEG;
            dbg<<"nu set to error value = "<<nu<<std::endl;
        }


        // 
        // Reduce order if necessary so that
        // (order+1)*(order+2)/2 < nu
        //
        galorder = galorder_init;
        if (params.read("shear_base_order_on_nu",true)) {
            int galsize;
            while (galorder > min_galorder) {
                if (maxm > galorder) maxm = galorder;
                galsize = (maxm+1)*(maxm+2)/2 + (2*maxm+1)*(galorder-maxm)/2;
                if (galsize <= nu) break;
                else --galorder;
            }
            if (galorder < galorder_init) {
                dbg<<"Reduced galorder to "<<galorder<<
                    " so that galsize = "<<galsize<<
                    " <= nu = "<<nu<<std::endl;
            }
        }


        //
        // Next find the frame in which the native observation is round.
        //
        Ellipse ell_round = ell_native;
        if (fixcen) ell_round.fixCen();
        if (fixsigma) ell_round.fixMu();
        for(int iter=1;iter<=MAX_ITER;++iter) {
            dbg<<"Round iter = "<<iter<<std::endl;
            flag1 = 0;
            std::complex<double> gamma_prev = ell_round.getGamma();
            if (ell_round.measure(
                    pix,galorder,galorder2,maxm,sigma_obs,flag1,1.e-2)) {
                ++log._ns_native;
                dbg<<"Successful round fit:\n";
                dbg<<"Z = "<<ell_round.getCen()<<std::endl;
                dbg<<"Mu = "<<ell_round.getMu()<<std::endl;
                dbg<<"Gamma = "<<ell_round.getGamma()<<std::endl;
            } else {
                ++log._nf_native;
                dbg<<"Round measurement failed\n";
                flag |= flag1;
                dbg<<"FLAG NATIVE_FAILED\n";
                flag |= NATIVE_FAILED;
                dbg<<"FLAG SHAPELET_NOT_DECONV\n";
                flag |= SHAPELET_NOT_DECONV;
                return;
            }

            // Adjust the sigma value so we can reset mu = 0.
            if (!fixcen) cen_offset += ell_round.getCen();
            ell_round.setCen(0.);
            double delta_mu = ell_round.getMu();
            sigma_obs *= exp(ell_round.getMu());
            ell_round.setMu(0.);
            if (sigma_obs < min_galsize*sigma_p) {
                dbg<<"skip: galaxy is too small -- "<<sigma_obs<<
                    " psf size = "<<sigma_p<<std::endl;
                ++log._nf_small;
                dbg<<"FLAG TOO_SMALL\n";
                flag |= TOO_SMALL;
                dbg<<"FLAG SHAPELET_NOT_DECONV\n";
                flag |= SHAPELET_NOT_DECONV;
                return;
            }
            gamma = ell_round.getGamma();
            dbg<<"New center = "<<cen_offset<<std::endl;
            if (std::abs(cen_offset) > 1.) {
                dbg<<"FLAG CENTROID_SHIFT\n";
                flag |= CENTROID_SHIFT;
            }
            dbg<<"New sigma_obs = "<<sigma_obs<<std::endl;
            dbg<<"New gamma = "<<gamma<<std::endl;
            std::complex<double> delta_gamma = gamma - gamma_prev;
            dbg<<"ell_round = "<<ell_round<<std::endl;

            // Get the pixels in an elliptical aperture based on the
            // observed shape.
            galap = sigma_obs * gal_aperture;
            if (max_aperture > 0. && galap > max_aperture) 
                galap = max_aperture;
            npix = 0;
            for(int i=0;i<nexp;++i) {
                if (use_fake_pixels) {
                    inner_fake_ap = shear_inner_fake_aperture * sigma_obs;
                    outer_fake_ap = shear_outer_fake_aperture * sigma_obs;
                }
                GetSubPixList(pix[i],allpix[i],cen_offset,gamma,galap,
                              inner_fake_ap,outer_fake_ap,params,flag1);
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
            dbg<<"delta_mu = "<<delta_mu<<std::endl;
            dbg<<"delta_gamma = "<<delta_gamma<<std::endl;
            if (std::abs(delta_mu) < MAX_DELTA_MU &&
                std::abs(delta_gamma) < MAX_DELTA_GAMMA1) {
                dbg<<"delta_mu < "<<MAX_DELTA_MU;
                dbg<<" and delta_gamma < "<<MAX_DELTA_GAMMA1<<std::endl;
                break;
            } else if (iter < MAX_ITER) {
                dbg<<"delta_mu >= "<<MAX_DELTA_MU;
                dbg<<" or delta_gamma >= "<<MAX_DELTA_GAMMA1<<std::endl;
                continue;
            } else {
                dbg<<"delta_mu >= "<<MAX_DELTA_MU;
                dbg<<" or delta_gamma >= "<<MAX_DELTA_GAMMA1<<std::endl;
                dbg<<"But iter == "<<MAX_ITER<<", so stop.\n";
            }
        }

        if (params.read("shear_native_only",false)) return;

        // Start with the specified fPsf, but allow it to increase up to
        // max_fpsf if there are any problems.
        long flag0 = flag;
        dbg<<"min_fpsf = "<<min_fpsf<<std::endl;
        dbg<<"max_fpsf = "<<max_fpsf<<std::endl;
        Assert(min_fpsf <= max_fpsf);
        for(double fPsf=min_fpsf; fPsf<=max_fpsf+1.e-3; fPsf+=0.5) {
            flag = flag0; // In case anything was set in a previous iteration.

            bool lastfpsf = (fPsf + 0.5 > max_fpsf+1.e-3);
            dbg<<"Start fPsf loop: fPsf = "<<fPsf<<std::endl;

            //
            // Set the sigma to use for the shapelet measurement:
            //
            // The sigma to use is 
            // sigma_s^2 = sigma_i^2 + fpsf sigma_p^2
            // And sigma_i^2 = sigma_o^2 - sigma_p^2
            // So sigma_s^2 = sigma_o^2 + (fpsf-1) sigma_p^2
            sigma = pow(sigma_obs,2) + (fPsf-1.) * sigpsq;
            xdbg<<"sigma_p^2 = "<<sigpsq<<std::endl;
            xdbg<<"sigma_obs^2 = "<<pow(sigma_obs,2)<<std::endl;
            xdbg<<"sigma^2 = sigma_obs^2 + (fPsf-1) * sigmap^2 = "<<sigma<<std::endl;
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
            DMatrix shapeCov(int(shapelet.size()),int(shapelet.size()));
            if (ell_native.measureShapelet(
                    pix,psf,shapelet,galorder,galorder2,galorder,&shapeCov)) {
                dbg<<"Successful deconvolving fit:\n";
                ++log._ns_mu;
            } else {
                dbg<<"Deconvolving measurement failed\n";
                if (!lastfpsf) continue;
                ++log._nf_mu;
                dbg<<"FLAG DECONV_FAILED\n";
                flag |= DECONV_FAILED;
            }
            dbg<<"Measured deconvolved b_gal = "<<shapelet.vec()<<std::endl;
            xdbg<<"shapeCov = "<<shapeCov<<std::endl;
            if (shapelet(0) >= flux(0)*3. || shapelet(0) <= flux(0)/3.) {
                // If the b00 value in the shapelet doesn't match the direct flux
                // measurement, set a flag.
                dbg<<"Bad flux value: \n";
                dbg<<"flux = "<<flux(0)<<std::endl;
                dbg<<"shapelet = "<<shapelet.vec()<<std::endl;
                dbg<<"flux ratio = "<<flux(0)/shapelet(0)<<std::endl;
            }

            //
            // Next, we find the frame where the deconvolved shape of the
            // galaxy is round.
            //
            for(int try_order=galorder; try_order>=min_galorder; --try_order) {
                dbg<<"try_order = "<<try_order<<std::endl;
                if (try_order < galorder) {
                    dbg<<"FLAG SHEAR_REDUCED_ORDER\n";
                    flag |= SHEAR_REDUCED_ORDER;
                }
                if (maxm > try_order) maxm = try_order;
                Ellipse ell_meas = ell_round;
                Ellipse ell_shear = ell_round;
                if (fixcen) ell_shear.fixCen();
                if (fixsigma) ell_shear.fixMu();
                //ell_shear.fixMu();
                double w = sqrt(sigma/sigma_p);
                bool success = false;
                cov.setZero();
                for(int iter=1;iter<=MAX_ITER;++iter) {
                    dbg<<"Shear iter = "<<iter<<std::endl;
                    flag1 = 0;
                    std::complex<double> gamma_prev = ell_shear.getGamma();
                    ell_meas.setGamma(
                        (w*ell_shear.getGamma() + ell_round.getGamma())/(w+1.));
                    if (ell_shear.measure(
                            pix,psf,try_order,galorder2,maxm,sigma,flag1,
                            1.e-2,&cov,&ell_meas)) {
                        dbg<<"Successful shear fit:\n";
                        dbg<<"Z = "<<ell_shear.getCen()<<std::endl;
                        dbg<<"Mu = "<<ell_shear.getMu()<<std::endl;
                        dbg<<"Gamma = "<<ell_shear.getGamma()<<std::endl;
                    } else {
                        dbg<<"Shear measurement failed\n";
                        success = false;
                        break;
                    }

                    gamma = ell_shear.getGamma();
                    dbg<<"New gamma = "<<gamma<<std::endl;
                    std::complex<double> delta_gamma = gamma - gamma_prev;
                    dbg<<"ell_shear = "<<ell_shear<<std::endl;

                    dbg<<"delta_gamma = "<<delta_gamma<<std::endl;
                    if (std::abs(delta_gamma) < MAX_DELTA_GAMMA2) {
                        dbg<<"delta_gamma < "<<MAX_DELTA_GAMMA2<<std::endl;
                        success = true;
                        flag1 = 0;
                        break;
                    } else if (iter == MAX_ITER) {
                        dbg<<"delta_gamma >= "<<MAX_DELTA_GAMMA2<<std::endl;
                        dbg<<"But iter == "<<MAX_ITER<<", so stop.\n";
                        success = std::abs(delta_gamma) < FAIL_DELTA_GAMMA2;
                    } else {
                        dbg<<"delta_gamma >= "<<MAX_DELTA_GAMMA2<<std::endl;
                    }
                }
                if (success) {
                    ++log._ns_gamma;
                    dbg<<"Successful shear measurement:\n";
                    dbg<<"shear = "<<gamma<<std::endl;
                    dbg<<"cov = "<<cov<<std::endl;
#if 0
                    ell_shear.correctForBias(
                        pix,psf,try_order,galorder2,maxm,sigma,&ell_meas);
                    gamma = ell_shear.getGamma();
                    dbg<<"after correct bias: shear => "<<gamma<<std::endl;
#endif
                    return;
                } else {
                    dbg<<"Unsuccessful shear measurement\n"; 
                    dbg<<"shear = "<<gamma<<std::endl;
                    if (try_order == min_galorder) {
                        flag |= flag1;
                        dbg<<"flag = "<<flag<<std::endl;
                    }
                }
            }
            if (!lastfpsf) continue;
            ++log._nf_gamma;
            dbg<<"FLAG SHEAR_FAILED\n";
            flag |= SHEAR_FAILED;
        }
#ifdef USE_TMV
    } catch (tmv::Error& e) {
        dbg<<"TMV Error thrown in MeasureShear\n";
        dbg<<e<<std::endl;
        ++log._nf_tmv_error;
        dbg<<"FLAG TMV_EXCEPTION\n";
        flag |= TMV_EXCEPTION;
#endif
    } catch (std::exception& e) {
        dbg<<"std::exception thrown in MeasureShear\n";
        dbg<<e.what()<<std::endl;
        ++log._nf_other_error;
        dbg<<"FLAG STD_EXCEPTION\n";
        flag |= STD_EXCEPTION;
    } catch (...) {
        dbg<<"unkown exception in MeasureShear\n";
        ++log._nf_other_error;
        dbg<<"FLAG UNKNOWN_EXCEPTION\n";
        flag |= UNKNOWN_EXCEPTION;
    } 

}

