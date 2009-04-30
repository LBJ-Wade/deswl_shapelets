
#include <fstream>
#include <iostream>

#include "ShearCatalog.h"
#include "ConfigFile.h"
#include "dbg.h"
#include "Params.h"
#include "BVec.h"
#include "Ellipse.h"
#include "TMV.h"
#include "Pixel.h"
#include "Image.h"
#include "FittedPSF.h"
#include "Log.h"
#include "Params.h"
#include "Form.h"
#include "WlVersion.h"
#include "TMV.h"

//#define SINGLEGAL 8146
//#define STARTAT 8000
//#define ENDAT 200

void MeasureSingleShear1(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans, const std::vector<BVec>& psf,
    double noise, double gain, const Image<double>* weight_im, 
    double gal_aperture, double max_aperture,
    int gal_order, int gal_order2,
    double f_psf, double min_gal_size, bool desqa,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    tmv::SmallMatrix<double,2,2>& shearcov, BVec& shapelet,
    long& flag)
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
  GetPixList(im,pix[0],cen,sky,noise,gain,weight_im,trans,galap,flag);
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
  long flag1=0;
  if (ell.Measure(pix,go,sigma_obs,true,flag1,desqa)) {
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
    flag |= NATIVE_FAILED;
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
    flag |= TOO_SMALL;
    return;
  }
  galap *= exp(ell.GetMu());
  if (max_aperture > 0. && galap > max_aperture) galap = max_aperture;
  ell.SetMu(0);

  // Check - mu should now be zero:
  // This one works
  //ell.Measure(pix,go,sigma_obs,true,flag1,desqa);
  //xdbg<<"After native fit #2:\n";
  //xdbg<<"Mu = "<<ell.GetMu()<<std::endl;

  pix[0].clear();
  GetPixList(im,pix[0],cen,sky,noise,gain,weight_im,trans,galap,flag);
  npix = pix[0].size();
  xdbg<<"npix = "<<npix<<std::endl;
  double sigpsq = pow(sigma_p,2);
  double sigma = pow(sigma_obs,2) - sigpsq;
  if (sigma < 0.1*sigpsq) sigma = 0.1*sigpsq;
  sigma = sqrt(sigma);
  xdbg<<"sigma = "<<sigma<<std::endl;
  ell.SetFP(f_psf);
  if (times) ell.ResetTimes();
  flag1 = 0;
  if (ell.Measure(pix,psf,go,sigma,false,flag1,desqa)) {
    if (times) {
      times->ns_mu++;
      times->ts_mu_integ += ell.t_integ;
      times->ts_mu_centroid += ell.t_centroid;
      times->ts_mu_fixflux += ell.t_fixflux;
      times->ts_mu_final += ell.t_final;
    }
    log.ns_mu++;
    if (flag1) {
      Assert((flag1 & ~(SHEAR_LOCAL_MIN | SHEAR_POOR_FIT)) == false);
      // This is just bumps the Ellipse flags up two spots
      // from the SHEAR_* to SHAPE_* flags.
      flag1 <<= 2;
      Assert((flag1 & ~(SHAPE_LOCAL_MIN | SHAPE_POOR_FIT)) == false);
      flag |= flag1;
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
    log.nf_mu++;
    dbg<<"Deconvolving measurement failed\n";
    flag |= DECONV_FAILED;
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
  ell.Measure(pix,psf,go,sigma,false,flag1,desqa,&b_gal);
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
    flag |= SHAPE_REDUCED_ORDER;
    if (go < 2) { flag |= TOO_SMALL; return; }
  }
  //shapelet = new BVec(go,sigma);
  shapelet.SetSigma(sigma);
  std::complex<double> gale = 0.;
  ell.MeasureShapelet(pix,psf,shapelet);
  xdbg<<"Measured b_gal = "<<shapelet<<std::endl;

  // Under normal circumstances, b20/b00 ~= conj(gamma)/sqrt(2)
  gale = std::complex<double>(shapelet[3],shapelet[4]);
  gale /= shapelet[0];
  if (std::abs(gale) < 0.5) {
    ell.SetGamma(conj(gale) * sqrt(2.));
    shear = ell.GetGamma();
  }

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
  if (ell.Measure(pix,psf,go,sigma,false,flag,desqa,&cov)) {
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
    tmv::SmallMatrix<double,2,2,tmv::RowMajor> cov1 = 
      cov.SubMatrix(2,4,2,4);
    if (!(cov1.Det() > 0.)) flag |= SHEAR_BAD_COVAR;
    else shearcov = cov1;
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
    flag |= SHEAR_FAILED;
    shear = ell.GetGamma();
    tmv::SmallMatrix<double,2,2,tmv::RowMajor> cov1 = 
      cov.SubMatrix(2,4,2,4);
    if (!(cov1.Det() > 0.)) flag |= SHEAR_BAD_COVAR;
    else shearcov = cov1;
    return;
  }
  //dbg<<"Stats: Mu: "<<mu_1<<"  "<<mu_2<<"  "<<mu_3<<std::endl;
  //dbg<<"       sigma = "<<sigma<<std::endl;
  //dbg<<"       gamma = "<<shear<<std::endl;

  // Finally measure the variance of the shear
  // TODO
  // (I'm not convinced that the above covariance matrix is a good estiamte.)
}

void MeasureSingleShear(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans, const FittedPSF& fitpsf,
    double noise, double gain, const Image<double>* weight_im, 
    double gal_aperture, double max_aperture,
    int gal_order, int gal_order2,
    double f_psf, double min_gal_size, bool desqa,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    tmv::SmallMatrix<double,2,2>& shearcov, BVec& shapelet,
    long& flag)
{
  // Get coordinates of the galaxy, and convert to sky coordinates
  try {
    // We don't need to save skypos.  We just want to catch the range
    // error here, so we don't need to worry about it for dudx, etc.
    Position skypos;
    trans.Transform(cen,skypos);
    dbg<<"skypos = "<<skypos<<std::endl;
  } catch (Range_error& e) {
    dbg<<"distortion range error: \n";
    xdbg<<"p = "<<cen<<", b = "<<e.b<<std::endl;
    if (times) times->nf_range1++;
    log.nf_range1++;
    flag |= TRANSFORM_EXCEPTION;
    return;
  }

  // Calculate the psf from the fitted-psf formula:
  std::vector<BVec> psf(1,
      BVec(fitpsf.GetPSFOrder(), fitpsf.GetSigma()));
  try {
    dbg<<"for fittedpsf cen = "<<cen<<std::endl;
    psf[0] = fitpsf(cen);
  } catch (Range_error& e) {
    dbg<<"fittedpsf range error: \n";
    xdbg<<"p = "<<cen<<", b = "<<e.b<<std::endl;
    if (times) times->nf_range2++;
    log.nf_range2++;
    flag |= FITTEDPSF_EXCEPTION;
    return;
  }

  // Do the real meat of the calculation:
  dbg<<"measure single shear cen = "<<cen<<std::endl;
  try { 
    MeasureSingleShear1(
	// Input data:
	cen, im, sky, trans, psf,
	// Noise variables:
	noise, gain, weight_im, 
	// Parameters:
	gal_aperture, max_aperture, gal_order, gal_order2, 
	f_psf, min_gal_size, desqa,
	// Time stats if desired:
	times,
	// Log information
	log,
	// Ouput values:
	shear, shearcov, shapelet, flag);
  } catch (tmv::Error& e) { 
    dbg<<"TMV Error thrown in MeasureSingleShear\n";
    dbg<<e<<std::endl;
    log.nf_tmverror++;
    flag |= TMV_EXCEPTION;
  } catch (...) {
    dbg<<"unkown exception in MeasureSingleShear\n";
    log.nf_othererror++;
    flag |= UNKNOWN_EXCEPTION;
  } 

}

int ShearCatalog::MeasureShears(const Image<double>& im,
    const Image<double>* weight_im, const Transformation& trans,
    const FittedPSF& fitpsf, ShearLog& log)
{
  int ngals = pos.size();
  dbg<<"ngals = "<<ngals<<std::endl;

  // Read some needed parameters
  double gal_aperture = params.get("shear_aperture");
  double max_aperture = params.read("shear_max_aperture",0.);
  int gal_order = params.get("shear_gal_order");
  int gal_order2 = params.read("shear_gal_order2",gal_order);
  double f_psf = params.get("shear_f_psf");
  double gain = params.read("image_gain",0.);
  double min_gal_size = params.get("shear_min_gal_size");
  bool desqa = params.read("des_qa",false);
  bool output_dots = params.read("output_dots",false);
  bool timing = params.read("timing",false);

  OverallFitTimes alltimes;

#ifdef ENDAT
  ngals = ENDAT;
#endif
  
  log.ngals = ngals;
#ifdef STARTAT
  log.ngals -= STARTAT;
#endif
#ifdef SINGLEGAL
  log.ngals = 1;
#endif
  log.ngoodin = std::count(flags.begin(),flags.end(),0);
  dbg<<log.ngoodin<<"/"<<log.ngals<<" galaxies with no input flags\n";

  // Main loop to measure shears
#ifdef _OPENMP
#pragma omp parallel 
  {
    try {
#endif
      OverallFitTimes times; // just for this thread
      ShearLog log1; // just for this thread
      log1.NoWriteLog();
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
      for(int i=0;i<ngals;i++) if (!flags[i]) {
#ifdef STARTAT
	if (i < STARTAT) continue;
#endif
#ifdef SINGLEGAL
	if (i < SINGLEGAL) continue;
	if (i > SINGLEGAL) break;
#endif
#ifdef _OPENMP
#pragma omp critical (output)
#endif
	{
	  if (output_dots) { std::cerr<<"."; std::cerr.flush(); }
	  dbg<<"galaxy "<<i<<":\n";
	  dbg<<"pos[i] = "<<pos[i]<<std::endl;
	}

	// Start with an error code of unknown failure, in case
	// something happens that I don't detect as an error.
	flags[i] = UNKNOWN_FAILURE;
	long flag1 = 0;
	MeasureSingleShear(
	    // Input data:
	    pos[i], im, sky[i], trans, 
	    // Fitted PSF
	    fitpsf,
	    // Noise variables:
	    noise[i], gain, weight_im, 
	    // Parameters:
	    gal_aperture, max_aperture, gal_order, gal_order2, 
	    f_psf, min_gal_size, desqa,
	    // Time stats if desired:
	    timing ? &times : 0, 
	    // Log information
	    log1,
	    // Ouput values:
	    shear[i], cov[i], shape[i], flag1);

	flags[i] = flag1;

	if (!flag1) {
	  ompdbg<<"Successful shear measurement: "<<shear[i]<<std::endl;
	}
	else {
	  ompdbg<<"Unsuccessful shear measurement\n"; 
	}

	if (timing) {
	  ompdbg<<"So far: ns = "<<times.ns_gamma<<",  nf = "<<times.nf_native;
	  ompdbg<<", "<<times.nf_mu<<", "<<times.nf_gamma<<std::endl;
	}

      }
#ifdef _OPENMP
#pragma omp critical (add_log)
#endif
      { 
	if (timing) alltimes += times;
	log += log1;
      }
#ifdef _OPENMP
    } 
    catch (...)
    {
      // This isn't supposed to happen.
      std::cerr<<"STATUS5BEG Caught some error in parallel region STATUS5END\n";
      exit(1);
    }
  }
#endif

  dbg<<log.ns_gamma<<" successful shear measurements, ";
  dbg<<ngals-log.ns_gamma<<" unsuccessful\n";
  log.ngood = std::count(flags.begin(),flags.end(),0);
  dbg<<log.ngood<<" with no flags\n";

  if (output_dots) {
    std::cerr
      <<std::endl
      <<"Success rate: "<<log.ns_gamma<<"/"<<log.ngoodin
      <<"  # with no flags: "<<log.ngood
      <<std::endl;
  }

  if (timing) {
    dbg<<"From timing structure:\n";
    dbg<<alltimes.ns_gamma<<" successful shear measurements, ";
    dbg<<alltimes.nf_native<<" + "<<alltimes.nf_mu;
    dbg<<" + "<<alltimes.nf_gamma<<" unsuccessful\n";
    std::cerr<<alltimes<<std::endl;
  }
  xdbg<<log<<std::endl;

  return log.ns_gamma;
}

ShearCatalog::ShearCatalog(const InputCatalog& incat,
    const Transformation& trans, const ConfigFile& _params) :
  id(incat.id), pos(incat.pos), sky(incat.sky), noise(incat.noise),
  flags(incat.flags), params(_params)
{
  dbg<<"Create ShearCatalog\n";

  // Fix flags to only have INPUT_FLAG set.
  // (I don't think this is necessary, but just in case.)
  for (size_t i=0; i<size(); i++) {
    flags[i] &= INPUT_FLAG;
  }

  // Calculate sky positions
  Position skypos_default(DEFVALNEG,DEFVALNEG);
  if (incat.ra.size() == 0 || incat.dec.size() == 0) {
    skypos.resize(size(),skypos_default);
    for(size_t i=0;i<size();i++) {
      try {
	trans.Transform(pos[i],skypos[i]);
	skypos[i] /= 3600.; // arcsec -> degrees
      } catch (Range_error& e) {
	xdbg<<"distortion range error\n";
	xdbg<<"p = "<<pos[i]<<", b = "<<e.b<<std::endl;
      }                       
    }
  } else {
    Assert(incat.ra.size() == size());
    Assert(incat.dec.size() == size());
    skypos.resize(size());
    for(size_t i=0;i<size();++i) {
      skypos[i] = Position(incat.ra[i],incat.dec[i]);
    }

    xdbg<<"Check transformation:\n";
    double rmserror = 0.;
    int count = 0;
    for(size_t i=0;i<size();i++) {
      try {
	Position temp;
	trans.Transform(pos[i],temp);
	temp /= 3600.; // arcsec -> degrees
	xdbg<<pos[i]<<"  "<<skypos[i]<<"  "<<temp;
	xdbg<<"  "<<temp-skypos[i]<<std::endl;
	rmserror += std::norm(temp-skypos[i]);
	count++;
      } catch (Range_error& e) {
	xdbg<<"distortion range error\n";
	xdbg<<"p = "<<pos[i]<<", b = "<<e.b<<std::endl;
      }
    }
    rmserror /= count;
    rmserror = std::sqrt(rmserror);
    xdbg<<"rms error = "<<rmserror*3600.<<" arcsec\n";
    if (rmserror > 0.1) {
      std::cout<<"STATUS3BEG Warning: Positions from WCS transformation have rms error of "<<rmserror<<" arcsec relative to ra, dec in catalog. STATUS3END"<<std::endl;
    }
  }

  shear.resize(id.size(),std::complex<double>(DEFVALPOS,DEFVALPOS));
  nu.resize(id.size(),DEFVALNEG);

  tmv::SmallMatrix<double,2,2> cov_default;
  cov_default = tmv::ListInit, DEFVALPOS, 0, 0, DEFVALPOS;
  cov.resize(id.size(),cov_default);

  int gal_order = params.get("shear_gal_order");
  BVec shape_default(gal_order,1.);
  shape_default.SetAllTo(DEFVALNEG);
  shape.resize(id.size(),shape_default);

  Assert(id.size() == size());
  Assert(pos.size() == size());
  Assert(sky.size() == size());
  Assert(noise.size() == size());
  Assert(flags.size() == size());
  Assert(skypos.size() == size());
  Assert(shear.size() == size());
  Assert(nu.size() == size());
  Assert(cov.size() == size());
  Assert(shape.size() == size());
}

ShearCatalog::ShearCatalog(const ConfigFile& _params) : params(_params)
{
  Read();

  Assert(id.size() == size());
  Assert(pos.size() == size());
  Assert(sky.size() == size());
  Assert(noise.size() == size());
  Assert(flags.size() == size());
  Assert(skypos.size() == size());
  Assert(shear.size() == size());
  Assert(nu.size() == size());
  Assert(cov.size() == size());
  Assert(shape.size() == size());
}

void ShearCatalog::WriteFits(std::string file) const
{
  Assert(id.size() == size());
  Assert(pos.size() == size());
  Assert(sky.size() == size());
  Assert(noise.size() == size());
  Assert(flags.size() == size());
  Assert(shear.size() == size());
  Assert(nu.size() == size());
  Assert(cov.size() == size());
  Assert(shape.size() == size());

  // ! means overwrite existing file
  CCfits::FITS fits("!"+file, CCfits::Write);


  const int nfields=17;
  std::vector<string> colnames(nfields);
  std::vector<string> colfmts(nfields);
  std::vector<string> colunits(nfields);

  colnames[0] = params.get("shear_id_col");
  colnames[1] = params.get("shear_x_col");
  colnames[2] = params.get("shear_y_col");
  colnames[3] = params.get("shear_sky_col");
  colnames[4] = params.get("shear_noise_col");
  colnames[5] = params.get("shear_flags_col");
  colnames[6] = params.get("shear_ra_col");
  colnames[7] = params.get("shear_dec_col");
  colnames[8] = params.get("shear_shear1_col");
  colnames[9] = params.get("shear_shear2_col");
  colnames[10] = params.get("shear_nu_col");
  colnames[11] = params.get("shear_cov00_col");
  colnames[12] = params.get("shear_cov01_col");
  colnames[13] = params.get("shear_cov11_col");
  colnames[14] = params.get("shear_order_col");
  colnames[15] = params.get("shear_sigma_col");
  colnames[16] = params.get("shear_coeffs_col");

  colfmts[0] = "1J"; // id
  colfmts[1] = "1D"; // x
  colfmts[2] = "1D"; // y
  colfmts[3] = "1D"; // sky
  colfmts[4] = "1D"; // noise
  colfmts[5] = "1J"; // flags
  colfmts[6] = "1D"; // ra
  colfmts[7] = "1D"; // dec
  colfmts[8] = "1D"; // shear1
  colfmts[9] = "1D"; // shear2
  colfmts[10] = "1D"; // nu
  colfmts[11] = "1D"; // cov00
  colfmts[12] = "1D"; // cov01
  colfmts[13] = "1D"; // cov11
  colfmts[14] = "1J"; // order
  colfmts[15] = "1D"; // sigma


  int ncoeff = shape[0].size();
  dbg<<"ncoeff = "<<ncoeff<<std::endl;
  std::stringstream coeff_form;
  coeff_form << ncoeff << "D";
  colfmts[16] = coeff_form.str(); // shapelet coeffs

  colunits[0] = "None";   // id
  colunits[1] = "Pixels"; // x
  colunits[2] = "Pixels"; // y
  colunits[3] = "ADU";    // sky
  colunits[4] = "ADU^2";  // noise
  colunits[5] = "None";   // flags
  colunits[6] = "Deg";    // ra
  colunits[7] = "Deg";    // dec
  colunits[8] = "None";   // shear1
  colunits[9] = "None";   // shear2
  colunits[10] = "None";  // nu
  colunits[11] = "None";  // cov00
  colunits[12] = "None";  // cov01
  colunits[13] = "None";  // cov11
  colunits[14] = "None";  // order
  colunits[15] = "Arcsec";// sigma
  colunits[16] = "None";  // coeffs

  dbg<<"Before Create table"<<std::endl;
  CCfits::Table* table;
  table = fits.addTable("shearcat",size(),colnames,colfmts,colunits);

  // Header keywords
  std::string tmvvers = tmv::TMV_Version();
  std::string wlvers = WlVersion();

  table->addKey("tmvvers", tmvvers, "version of TMV code");
  table->addKey("wlvers", wlvers, "version of weak lensing code");

  std::string str;
  double dbl;
  int intgr;
  //CCfitsWriteParKey(params, table, "version", str);
  CCfitsWriteParKey(params, table, "noise_method", str);
  CCfitsWriteParKey(params, table, "dist_method", str);

  CCfitsWriteParKey(params, table, "shear_aperture", dbl);
  CCfitsWriteParKey(params, table, "shear_max_aperture", dbl);
  CCfitsWriteParKey(params, table, "shear_gal_order", intgr);
  CCfitsWriteParKey(params, table, "shear_gal_order2", intgr);
  CCfitsWriteParKey(params, table, "shear_min_gal_size", dbl);
  CCfitsWriteParKey(params, table, "shear_f_psf", dbl);



  // data
  // make vector copies for writing
  std::vector<double> x(pos.size());
  std::vector<double> y(pos.size());
  std::vector<double> ra(size());
  std::vector<double> dec(size());
  std::vector<double> shear1(size());
  std::vector<double> shear2(size());
  std::vector<double> cov00(size());
  std::vector<double> cov01(size());

  for(size_t i=0;i<pos.size();i++) {
    x[i] = pos[i].GetX();
    y[i] = pos[i].GetY();
  }
  for(size_t i=0;i<size();i++) { 
    ra[i] = skypos[i].GetX();
    dec[i] = skypos[i].GetY();
  }
  for(size_t i=0;i<size();i++) { 
    shear1[i] = real(shear[i]);
    shear2[i] = imag(shear[i]);
  }
  std::vector<double> cov11(size());
  for(size_t i=0;i<size();i++) { 
    cov00[i] = cov[i](0,0);
    cov01[i] = cov[i](0,1);
    cov11[i] = cov[i](1,1);
  }


  int startrow=1;

  table->column(colnames[0]).write(id,startrow);
  table->column(colnames[1]).write(x,startrow);
  table->column(colnames[2]).write(y,startrow);
  table->column(colnames[3]).write(sky,startrow);
  table->column(colnames[4]).write(noise,startrow);
  table->column(colnames[5]).write(flags,startrow);
  table->column(colnames[6]).write(ra,startrow);
  table->column(colnames[7]).write(dec,startrow);
  table->column(colnames[8]).write(shear1,startrow);
  table->column(colnames[9]).write(shear2,startrow);
  table->column(colnames[10]).write(nu,startrow);
  table->column(colnames[11]).write(cov00,startrow);
  table->column(colnames[12]).write(cov01,startrow);
  table->column(colnames[13]).write(cov11,startrow);

  for (size_t i=0; i<size(); i++) {
    size_t row = i+1;
    long b_order = shape[i].GetOrder();
    double b_sigma = shape[i].GetSigma();

    table->column(colnames[14]).write(&b_order,1,row);
    table->column(colnames[15]).write(&b_sigma,1,row);
    double* cptr = (double *) shape[i].cptr();
    table->column(colnames[16]).write(cptr, ncoeff, 1, row);

  }



}




void ShearCatalog::WriteAscii(std::string file, std::string delim) const
{
  Assert(id.size() == size());
  Assert(pos.size() == size());
  Assert(sky.size() == size());
  Assert(noise.size() == size());
  Assert(flags.size() == size());
  Assert(skypos.size() == size());
  Assert(shear.size() == size());
  Assert(nu.size() == size());
  Assert(cov.size() == size());
  Assert(shape.size() == size());

  std::ofstream fout(file.c_str());
  if (!fout) {
    throw std::runtime_error("Error opening shear file");
  }

  Form hexform; hexform.hex().trail(0);

  for(size_t i=0;i<size();i++) {
    fout
      << id[i] << delim
      << pos[i].GetX() << delim
      << pos[i].GetY() << delim
      << sky[i] << delim
      << noise[i] << delim
      << hexform(flags[i]) << delim
      << skypos[i].GetX() << delim
      << skypos[i].GetY() << delim
      << real(shear[i]) << delim
      << imag(shear[i]) << delim
      << nu[i] << delim
      << cov[i](0,0) << delim
      << cov[i](0,1) << delim
      << cov[i](1,1) << delim
      << shape[i].GetOrder() << delim
      << shape[i].GetSigma();
    for(size_t j=0;j<shape[i].size();++j)
      fout << delim << shape[i][j];
    fout << std::endl;
  }
}

void ShearCatalog::Write() const
{
  std::vector<std::string> files = MultiName(params, "shear");  

  for(size_t i=0; i<files.size(); ++i) {
    const std::string& file = files[i];
    dbg<<"Writing shear catalog to file: "<<file<<std::endl;

    bool fitsio = false;
    if (params.keyExists("shear_io")) {
      std::vector<std::string> ios = params["shear_io"];
      Assert(ios.size() == files.size());
      fitsio = (ios[i] == "FITS");
    }
    else if (file.find("fits") != std::string::npos) 
      fitsio = true;

    if (fitsio) {
      WriteFits(file);
    } else {
      std::string delim = "  ";
      if (params.keyExists("shear_delim")) {
	std::vector<std::string> delims = params["shear_delim"];
	Assert(delims.size() == files.size());
	delim = delims[i];
      }
      else if (file.find("csv") != std::string::npos) delim = ",";
      WriteAscii(file,delim);
    }
  }
  dbg<<"Done Write ShearCatalog\n";
}

void ShearCatalog::ReadFits(std::string file)
{
  int hdu = params.read("shear_hdu",2);

  dbg<<"Opening FITS file at hdu "<<hdu<<std::endl;
  // true means read all as part of the construction
  CCfits::FITS fits(file, CCfits::Read, hdu-1, true);

  CCfits::ExtHDU& table=fits.extension(hdu-1);

  long nrows=table.rows();

  dbg<<"  nrows = "<<nrows<<std::endl;
  if (nrows <= 0) {
    throw std::runtime_error("nrows must be > 0");
  }

  std::string id_col=params.get("shear_id_col");
  std::string x_col=params.get("shear_x_col");
  std::string y_col=params.get("shear_y_col");
  std::string sky_col=params.get("shear_sky_col");
  std::string noise_col=params.get("shear_noise_col");
  std::string flags_col=params.get("shear_flags_col");
  std::string ra_col=params.get("shear_ra_col");
  std::string dec_col=params.get("shear_dec_col");
  std::string shear1_col=params.get("shear_shear1_col");
  std::string shear2_col=params.get("shear_shear2_col");
  std::string nu_col=params.get("shear_nu_col");
  std::string cov00_col=params.get("shear_cov00_col");
  std::string cov01_col=params.get("shear_cov01_col");
  std::string cov11_col=params.get("shear_cov11_col");
  std::string order_col=params.get("shear_order_col");
  std::string sigma_col=params.get("shear_sigma_col");
  std::string coeffs_col=params.get("shear_coeffs_col");

  long start=1;
  long end=nrows;

  dbg<<"Reading columns"<<std::endl;
  dbg<<"  "<<id_col<<std::endl;
  table.column(id_col).read(id, start, end);

  dbg<<"  "<<x_col<<"  "<<y_col<<std::endl;
  pos.resize(nrows);
  std::vector<double> x;
  std::vector<double> y;
  table.column(x_col).read(x, start, end);
  table.column(y_col).read(y, start, end);
  for(long i=0;i<nrows;++i) pos[i] = Position(x[i],y[i]);

  dbg<<"  "<<sky_col<<std::endl;
  table.column(sky_col).read(sky, start, end);

  dbg<<"  "<<noise_col<<std::endl;
  table.column(noise_col).read(noise, start, end);

  dbg<<"  "<<flags_col<<std::endl;
  table.column(flags_col).read(flags, start, end);

  dbg<<"  "<<ra_col<<"  "<<dec_col<<std::endl;
  skypos.resize(nrows);
  std::vector<double> ra;
  std::vector<double> dec;
  table.column(ra_col).read(ra, start, end);
  table.column(dec_col).read(dec, start, end);
  for(long i=0;i<nrows;++i) {
    skypos[i] = Position(ra[i],dec[i]);
  }

  dbg<<"  "<<shear1_col<<"  "<<shear2_col<<std::endl;
  shear.resize(nrows);
  std::vector<double> shear1;
  std::vector<double> shear2;
  table.column(shear1_col).read(shear1, start, end);
  table.column(shear2_col).read(shear2, start, end);
  for(long i=0;i<nrows;++i) {
    shear[i] = std::complex<double>(shear1[i],shear2[i]);
  }

  dbg<<"  "<<nu_col<<std::endl;
  table.column(nu_col).read(nu, start, end);

  dbg<<"  "<<cov00_col<<"  "<<cov01_col<<"  "<<cov11_col<<std::endl;
  cov.resize(nrows);
  std::vector<double> cov00;
  std::vector<double> cov01;
  std::vector<double> cov11;
  table.column(cov00_col).read(cov00, start, end);
  table.column(cov01_col).read(cov01, start, end);
  table.column(cov11_col).read(cov11, start, end);
  for(long i=0;i<nrows;++i) {
    cov[i] = tmv::ListInit, cov00[i], cov01[i], cov01[i], cov11[i];
  }

  // temporary
  std::vector<double> sigma;
  std::vector<int> order;
  table.column(sigma_col).read(sigma, start, end);
  table.column(order_col).read(order, start, end);

  shape.reserve(nrows);
  for (size_t i=0; i<size(); i++) {
    size_t row=i+1;

    shape.push_back(BVec(order[i],sigma[i]));
    int ncoeff=(order[i]+1)*(order[i]+2)/2;

    std::valarray<double> coeffs;
    table.column(coeffs_col).read(coeffs, row);

    double* ptri = (double* ) shape[i].cptr(); 
    for (int j=0; j<ncoeff; ++j) {
      ptri[i] = coeffs[i];
    }

  }
}

void ShearCatalog::ReadAscii(std::string file, std::string delim)
{
  std::ifstream fin(file.c_str());
  if (!fin) {
    throw std::runtime_error("Error opening stars file");
  }

  id.clear(); pos.clear(); sky.clear(); noise.clear(); flags.clear();
  skypos.clear(); shear.clear(); nu.clear(); cov.clear(); shape.clear();

  if (delim == "  ") {
    ConvertibleString flag;
    long id1,b_order;
    double x,y,sky1,noise1,ra,dec,s1,s2,nu1,c00,c01,c11,b_sigma;
    while ( fin >> id1 >> x >> y >> sky1 >> noise1 >> 
	flag >> ra >> dec >> s1 >> s2 >> nu1 >>
	c00 >> c01 >> c11 >> b_order >> b_sigma )
    {
      id.push_back(id1);
      pos.push_back(Position(x,y));
      sky.push_back(sky1);
      noise.push_back(noise1);
      flags.push_back(flag);
      skypos.push_back(Position(ra,dec));
      shear.push_back(std::complex<double>(s1,s2));
      nu.push_back(nu1);
      cov.push_back(tmv::SmallMatrix<double,2,2>());
      cov.back() = tmv::ListInit, c00, c01, c01, c11;
      shape.push_back(BVec(b_order,b_sigma));
      for(size_t j=0;j<shape.back().size();++j)
	fin >> shape.back()[j];
    }
  } else {
    if (delim.size() > 1) {
      // getline only works with single character delimiters.
      // Since I don't really expect a multicharacter delimiter to
      // be used ever, I'm just going to throw an exception here 
      // if we do need it, and I can write the workaround then.
      throw std::runtime_error(
	  "ReadAscii delimiter must be a single character");
    }
    char d = delim[0];
    long b_order;
    double x,y,ra,dec,s1,s2,b_sigma,c00,c01,c11;
    ConvertibleString temp;
    while (getline(fin,temp,d))
    {
      id.push_back(temp);
      getline(fin,temp,d); x = temp;
      getline(fin,temp,d); y = temp;
      pos.push_back(Position(x,y));
      getline(fin,temp,d); sky.push_back(temp);
      getline(fin,temp,d); noise.push_back(temp);
      getline(fin,temp,d); flags.push_back(temp);
      getline(fin,temp,d); ra = temp;
      getline(fin,temp,d); dec = temp;
      skypos.push_back(Position(ra,dec));
      getline(fin,temp,d); s1 = temp;
      getline(fin,temp,d); s2 = temp;
      shear.push_back(std::complex<double>(s1,s2));
      getline(fin,temp,d); nu.push_back(temp);
      getline(fin,temp,d); c00 = temp;
      getline(fin,temp,d); c01 = temp;
      getline(fin,temp,d); c11 = temp;
      cov.push_back(tmv::SmallMatrix<double,2,2>());
      cov.back() = tmv::ListInit, c00, c01, c01, c11;
      getline(fin,temp,d); b_order = temp;
      getline(fin,temp,d); b_sigma = temp;
      shape.push_back(BVec(b_order,b_sigma));
      for(size_t j=0;j<shape.back().size()-1;++j) {
	getline(fin,temp,d); shape.back()[j] = temp;
      }
      getline(fin,temp); shape.back()[shape.back().size()-1] = temp;
    }
  }
}


void ShearCatalog::Read()
{
  std::string file = Name(params,"shear",false,true);
  // false,true = input_prefix=false, mustexist=true.
  // It is an input here, but it is in the output_prefix directory.
  dbg<< "Reading Shear cat from file: " << file << std::endl;

  bool fitsio = false;
  if (params.keyExists("shear_io")) 
    fitsio = (params["shear_io"] == "FITS");
  else if (file.find("fits") != std::string::npos) 
    fitsio = true;

  if (fitsio) {
    ReadFits(file);
  } else {
    std::string delim = "  ";
    if (params.keyExists("shear_delim")) delim = params["shear_delim"];
    else if (file.find("csv") != std::string::npos) delim = ",";
    ReadAscii(file,delim);
  }
  dbg<<"Done Read ShearCatalog\n";
}

