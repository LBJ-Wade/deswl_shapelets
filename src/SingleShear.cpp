#include "types.h"
#include "Params.h"

#include "BVec.h"
#include "Ellipse.h"
#include "TMV.h"
#include "dbg.h"
#include "ConfigFile.h"
#include "Pixel.h"
#include "Image.h"
#include "FittedPSF.h"
#include "DoMeasure.h"
#include "TimeVars.h"
#include "Log.h"

#include <fstream>
#include <iostream>

void MeasureSingleShear(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans, const std::vector<BVec>& psf,
    double noise, double gain, const Image<double>* weight_im, 
    double gal_aperture, double max_aperture,
    int gal_order, int gal_order2,
    double f_psf, double min_gal_size,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    tmv::Matrix<double>& shearcov, BVec& shapelet,
    int32& flags)
{
  // Find harmonic mean of psf sizes:
  // MJ: Is it correct to use the harmonic mean of sigma^2?
  double sigma_p = 0.;
  Assert(psf.size() > 0);
  for(size_t i=0;i<psf.size();i++) sigma_p += 1./psf[i].GetSigma();
  sigma_p = double(psf.size()) / sigma_p;

  // We start with a native fit using sigma = 2*sigma_p:
  double sigma_obs = 2.*sigma_p;
  double galap = gal_aperture * sigma_obs;
  xdbg<<"galap = "<<gal_aperture<<" * "<<sigma_obs<<" = "<<galap<<std::endl;
  if (max_aperture > 0. && galap > max_aperture) {
    galap = max_aperture;
    xdbg<<"      => "<<galap<<std::endl;
  }
  xdbg<<"sigma_obs = "<<sigma_obs<<", sigma_p = "<<sigma_p<<std::endl;

  std::vector<std::vector<Pixel> > pix(1);
  int getpix_flag = 0;
  GetPixList(im,pix[0],cen,sky,noise,gain,weight_im,trans,galap,getpix_flag);
  if (getpix_flag) {
    dbg<<"skip: flag == "<<getpix_flag<<std::endl;
    if (times) times->nf_pixflag1++;
    if (getpix_flag & EDGE) {
      log.nf_edge1++;
      flags |= MSH_EDGE1;
    }
    if (flags & LT10PIX) {
      log.nf_npix1++;
      flags |= MSH_LT10PIX1;
    }
  }
  int npix = pix[0].size();
  xdbg<<"npix = "<<npix<<std::endl;

  Ellipse ell;
  ell.FixGam();
  ell.CrudeMeasure(pix[0],sigma_obs);
  xdbg<<"Crude Measure: centroid = "<<ell.GetCen()<<", mu = "<<ell.GetMu()<<std::endl;
  //double mu_1 = ell.GetMu();

  int go = 2;
  int gal_size = (go+1)*(go+2)/2;

  if (times) ell.DoTimings();
  if (ell.Measure(pix,go,sigma_obs,true)) {
    if (times) {
      times->ns_native++;
      times->ts_native_integ += ell.t_integ;
      times->ts_native_centroid += ell.t_centroid;
      times->ts_native_fixflux += ell.t_fixflux;
      times->ts_native_final += ell.t_final;
    }
    log.ns_native++;
    dbg<<"Successful native fit:\n";
    dbg<<"Z = "<<ell.GetCen()<<std::endl;
    dbg<<"Mu = "<<ell.GetMu()<<std::endl;
  }
  else {
    if (times) {
      times->nf_native++;
      times->tf_native_integ += ell.t_integ;
      times->tf_native_centroid += ell.t_centroid;
      times->tf_native_fixflux += ell.t_fixflux;
      times->tf_native_final += ell.t_final;
    }
    log.nf_native++;
    dbg<<"Native measurement failed\n";
    flags |= MSH_NATIVE_FAILED;
    shear = ell.GetGamma();
    ell.MeasureShapelet(pix,shapelet);
    return;
  }
  //double mu_2 = ell.GetMu();

  // Now we can do a deconvolving fit, but one that does not 
  // shear the coordinates.
  sigma_obs *= exp(ell.GetMu());
  if (sigma_obs < min_gal_size*sigma_p) {
    dbg<<"skip: galaxy is too small -- "<<sigma_obs<<" psf size = "<<sigma_p<<std::endl;
    if (times) times->nf_small++;
    log.nf_small++;
    flags |= MSH_TOOSMALL;
    shear = ell.GetGamma();
    ell.MeasureShapelet(pix,shapelet);
    return;
  }
  galap *= exp(ell.GetMu());
  if (max_aperture > 0. && galap > max_aperture) galap = max_aperture;
  ell.SetMu(0);

  // Check - mu should now be zero:
  // This one works
  //ell.Measure(pix,go,sigma_obs);
  //xdbg<<"After native fit #2:\n";
  //xdbg<<"Mu = "<<ell.GetMu()<<std::endl;

  pix[0].clear();
  getpix_flag=0;
  GetPixList(im,pix[0],cen,sky,noise,gain,weight_im,trans,galap,getpix_flag);
  if (getpix_flag) {
    dbg<<"skip: flag == "<<getpix_flag<<std::endl;
    if (times) times->nf_pixflag2++;
    if (flags & EDGE) {
      log.nf_edge2++;
      flags |= MSH_EDGE2;
    }
    if (flags & LT10PIX) {
      log.nf_npix2++;
      flags |= MSH_LT10PIX2;
    }
    shear = ell.GetGamma();
    ell.MeasureShapelet(pix,shapelet);
    return;
  }
  npix = pix[0].size();
  xdbg<<"npix = "<<npix<<std::endl;
  double sigpsq = pow(sigma_p,2);
  double sigma = pow(sigma_obs,2) - sigpsq;
  if (sigma < 0.1*sigpsq) sigma = 0.1*sigpsq;
  sigma = sqrt(sigma);
  xdbg<<"sigma = "<<sigma<<std::endl;
  ell.SetFP(f_psf);
  if (times) ell.ResetTimes();
  if (ell.Measure(pix,psf,go,sigma,false)) {
    if (times) {
      times->ns_mu++;
      times->ts_mu_integ += ell.t_integ;
      times->ts_mu_centroid += ell.t_centroid;
      times->ts_mu_fixflux += ell.t_fixflux;
      times->ts_mu_final += ell.t_final;
    }
    log.ns_mu++;
    dbg<<"Successful deconvolving fit:\n";
    xdbg<<"Mu = "<<ell.GetMu()<<std::endl;
  }
  else {
    if (times) {
      times->nf_mu++;
      times->tf_mu_integ += ell.t_integ;
      times->tf_mu_centroid += ell.t_centroid;
      times->tf_mu_fixflux += ell.t_fixflux;
      times->tf_mu_final += ell.t_final;
    }
    log.nf_mu++;
    dbg<<"Deconvolving measurement failed\n";
    flags |= MSH_DECONV_FAILED;
    shear = ell.GetGamma();
    ell.MeasureShapelet(pix,psf,shapelet);
    return;
  }
  //double mu_3 = ell.GetMu();

#if 0
  // This doesn't work.  So for now, keep mu, sigma the same
  // I'd rather have mu = 0 and the right sigma.
  // Maybe I need to change to fit to solve for sigma directly
  // rather than solve for mu.
  sigma *= exp(ell.GetMu());
  ell.SetMu(0);

  // Check - mu should now be zero:
  b_gal.SetSigma(sigma);
  dbg<<"Meausre with sigma = "<<sigma<<std::endl;
  ell.Measure(pix,psf,go,sigma,false,&b_gal);
  dbg<<"After deconvolving fit #2:\n";
  dbg<<"Mu = "<<ell.GetMu()<<std::endl;
  dbg<<"b_gal = "<<b_gal<<std::endl;
#endif

  // Measure the galaxy shape at the full order
  go = gal_order;
  gal_size = (go+1)*(go+2)/2;
  if (npix <= gal_size) {
    while (npix <= gal_size) { gal_size -= go+1; --go; }
    dbg<<"Too few pixels ("<<npix<<") for given gal_order. \n";
    dbg<<"Reduced gal_order to "<<go<<" gal_size = "<<gal_size<<std::endl;
  }
  //shapelet = new BVec(go,sigma);
  shapelet.SetSigma(sigma);
  std::complex<double> gale = 0.;
  ell.MeasureShapelet(pix,psf,shapelet);
  xdbg<<"Measured b_gal = "<<shapelet<<std::endl;

  // Under normal circumstances, b20/b00 ~= conj(gamma)/sqrt(2)
  gale = std::complex<double>(shapelet[3],shapelet[4]);
  gale /= shapelet[0];
  if (std::abs(gale) < 0.5) ell.SetGamma(conj(gale) * sqrt(2.));

  // Finally, we calculate the shear in the deconvolved galaxy.
  //ell.FixMu();
  ell.UnFixGam();
  go = gal_order2;
  gal_size = (go+1)*(go+2)/2;
  if (npix <= gal_size) {
    while (npix <= gal_size) { gal_size -= go+1; --go; }
    dbg<<"Too few pixels ("<<npix<<") for given gal_order. \n";
    dbg<<"Reduced gal_order to "<<go<<" gal_size = "<<gal_size<<std::endl;
  }
  if (times) ell.ResetTimes();
  tmv::Matrix<double> cov(5,5);
  if (ell.Measure(pix,psf,go,sigma,false,&cov)) {
    if (times) {
      times->ns_gamma++;
      times->ts_gamma_integ += ell.t_integ;
      times->ts_gamma_centroid += ell.t_centroid;
      times->ts_gamma_fixflux += ell.t_fixflux;
      times->ts_gamma_final += ell.t_final;
    }
    log.ns_gamma++;
    dbg<<"Successful Gamma fit\n";
    xdbg<<"Measured gamma = "<<ell.GetGamma()<<std::endl;
    shear = ell.GetGamma();
    shearcov = cov.SubMatrix(2,4,2,4);
    //shearcov.SetToIdentity(0.001);
  }
  else {
    if (times) {
      times->nf_gamma++;
      times->tf_gamma_integ += ell.t_integ;
      times->tf_gamma_centroid += ell.t_centroid;
      times->tf_gamma_fixflux += ell.t_fixflux;
      times->tf_gamma_final += ell.t_final;
    }
    log.nf_gamma++;
    dbg<<"Measurement failed\n";
    flags |= MSH_SHEAR_FAILED;
    shear = ell.GetGamma();
    shearcov = cov.SubMatrix(2,4,2,4);
    return;
  }
  //dbg<<"Stats: Mu: "<<mu_1<<"  "<<mu_2<<"  "<<mu_3<<std::endl;
  //dbg<<"       sigma = "<<sigma<<std::endl;
  //dbg<<"       gamma = "<<shear<<std::endl;

  // Finally measure the variance of the shear
  // TODO
  // (I'm not convinced that the above covariance matrix is a good estiamte.)
}

void MeasureSingleShear1(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans, const FittedPSF& fittedpsf,
    double noise, double gain, const Image<double>* weight_im, 
    double gal_aperture, double max_aperture,
    int gal_order, int gal_order2,
    double f_psf, double min_gal_size,
    OverallFitTimes* times, ShearLog& log,
    Position& skypos,
    std::complex<double>& shear, 
    tmv::Matrix<double>& shearcov, BVec& shapelet,
    int32& flags)
{
  // Get coordinates of the galaxy, and convert to sky coordinates
  try {
    trans.Transform(cen,skypos);
    dbg<<"skypos = "<<skypos<<std::endl;
  } catch (Range_error& e) {
    dbg<<"distortion range error: \n";
    xdbg<<"p = "<<cen<<", b = "<<e.b<<std::endl;
    if (times) times->nf_range1++;
    log.nf_range1++;
    flags |= MSH_TRANSFORM_EXCEPTION;
    return;
  }

  // Calculate the psf from the fitted-psf formula:
  std::vector<BVec> psf(1,
      BVec(fittedpsf.GetOrder(), fittedpsf.GetSigma()));
  try {
    dbg<<"for fittedpsf cen = "<<cen<<std::endl;
    psf[0] = fittedpsf(cen);
  } catch (Range_error& e) {
    dbg<<"fittedpsf range error: \n";
    xdbg<<"p = "<<cen<<", b = "<<e.b<<std::endl;
    if (times) times->nf_range2++;
    log.nf_range2++;
    flags |= MSH_FITTEDPSF_EXCEPTION;
    return;
  }

  // Do the real meat of the calculation:
  dbg<<"measure single shear cen = "<<cen<<std::endl;
  try { 
    MeasureSingleShear(
	// Input data:
	cen, im, sky, trans, psf,
	// Noise variables:
	noise, gain, weight_im, 
	// Parameters:
	gal_aperture, max_aperture, gal_order, gal_order2, 
	f_psf, min_gal_size, 
	// Time stats if desired:
	times,
	// Log information
	log,
	// Ouput values:
	shear, shearcov, shapelet, flags);
  } catch (tmv::Error& e) { 
    dbg<<"TMV Error thrown in MeasureSingleShear\n";
    dbg<<e<<std::endl;
    log.nf_tmverror++;
    flags |= MSH_TMV_EXCEPTION;
  } catch (...) {
    dbg<<"unkown exception in MeasureSingleShear\n";
    log.nf_othererror++;
    flags |= MSH_UNKNOWN_EXCEPTION;
  } 

}

