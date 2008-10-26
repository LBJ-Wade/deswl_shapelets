
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

#include <fstream>
#include <iostream>

//#define GAMMAOUT

// Define this to do the shapelet measurement.
//#define SHAPELET

void MeasureSingleShear(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans, const std::vector<BVec>& psf,
    double noise, double gain, const Image<double>* weight_im, 
    double gal_aperture, double max_aperture,
    int gal_order, int gal_order2,
    double f_psf, double min_gal_size,
    OverallFitTimes* times,
    bool& success, std::complex<double>& shear, 
    tmv::Matrix<double>& shearcov, BVec*& shapelet)
{
#ifdef GAMMAOUT
  static std::ofstream gammaout("gamma.out");
#endif
  success = false;

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
  int flag = 0;
  GetPixList(im,pix[0],cen,sky,noise,gain,weight_im,trans,galap,flag);
  if (flag) {
    dbg<<"skip: flag == "<<flag<<std::endl;
    if (times) times->nf_pixflag1++;
    return;
  }
  int npix = pix[0].size();
  xdbg<<"npix = "<<npix<<std::endl;
  if (npix < 10) {
    dbg<<"skip: npix == "<<npix<<std::endl;
    if (times) times->nf_npix1++;
    return;
  }

  Ellipse ell;
  ell.FixGam();
  ell.CrudeMeasure(pix[0],sigma_obs);
  xdbg<<"Crude Measure: centroid = "<<ell.GetCen()<<", mu = "<<ell.GetMu()<<std::endl;
  double mu_1 = ell.GetMu();

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
    dbg<<"Native measurement failed\n";
    return;
  }
  double mu_2 = ell.GetMu();

  // Now we can do a deconvolving fit, but one that does not 
  // shear the coordinates.
  sigma_obs *= exp(ell.GetMu());
  if (sigma_obs < min_gal_size*sigma_p) {
    dbg<<"skip: galaxy is too small -- "<<sigma_obs<<" psf size = "<<sigma_p<<std::endl;
    if (times) times->nf_small++;
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
  GetPixList(im,pix[0],cen,sky,noise,gain,weight_im,trans,galap,flag);
  if (flag) {
    dbg<<"skip: flag == "<<flag<<std::endl;
    if (times) times->nf_pixflag2++;
    return;
  }
  npix = pix[0].size();
  xdbg<<"npix = "<<npix<<std::endl;
  if (npix < 10) {
    dbg<<"skip: npix == "<<npix<<std::endl;
    if (times) times->nf_npix2++;
    return;
  }
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
    dbg<<"Deconvolving measurement failed\n";
    return;
  }
  double mu_3 = ell.GetMu();

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

#ifdef SHAPELET
  // Measure the galaxy shape at the full order
  go = gal_order;
  gal_size = (go+1)*(go+2)/2;
  if (npix <= gal_size) {
    while (npix <= gal_size) { gal_size -= go+1; --go; }
    dbg<<"Too few pixels ("<<npix<<") for given gal_order. \n";
    dbg<<"Reduced gal_order to "<<go<<" gal_size = "<<gal_size<<std::endl;
  }
  shapelet = new BVec(go,sigma);
  std::complex<double> gale = 0.;
  ell.MeasureShapelet(pix,psf,*shapelet);
  xdbg<<"b_gal = "<<*shapelet<<std::endl;
  gale = std::complex<double>((*shapelet)[3],(*shapelet)[4]);
  gale /= (*shapelet)[0];

  // Under normal circumstances, b20/b00 ~= conj(gamma)/sqrt(2)
  if (std::abs(gale) < 0.6) ell.SetGamma(conj(gale) * sqrt(2.));
#else
  shapelet = 0;
#endif

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
  //if (ell.Measure(pix,psf,go,sigma,false)) {
    if (times) {
      times->ns_gamma++;
      times->ts_gamma_integ += ell.t_integ;
      times->ts_gamma_centroid += ell.t_centroid;
      times->ts_gamma_fixflux += ell.t_fixflux;
      times->ts_gamma_final += ell.t_final;
    }
    dbg<<"Successful Gamma fit\n";
    xdbg<<"Measured gamma = "<<ell.GetGamma()<<std::endl;
    success = true;
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
    dbg<<"Measurement failed\n";
    return;
  }
  dbg<<"Stats: Mu: "<<mu_1<<"  "<<mu_2<<"  "<<mu_3<<std::endl;
  dbg<<"       sigma = "<<sigma<<std::endl;
#ifdef SHAPELET
  dbg<<"       b_20 = "<<gale<<"  gamma = "<<shear<<std::endl;
#else
  dbg<<"       gamma = "<<shear<<std::endl;
#endif
#ifdef GAMMAOUT
#ifdef _OPENMP
#pragma omp critical
#endif
  {
    gammaout<<real(gale)<<"  "<<imag(gale)<<"  "<<real(shear)<<"  "<<imag(shear)<<std::endl;
  }
#endif

  // Finally measure the variance of the shear
  // MJ
}

