
#include <fstream>
#include <iostream>

#include "ShearCatalog.h"
#include "ConfigFile.h"
#include "FitsFile.h"
#include "dbg.h"
#include "Params.h"
#include "BVec.h"
#include "Ellipse.h"
#include "TMV.h"
#include "Pixel.h"
#include "Image.h"
#include "FittedPSF.h"
#include "TimeVars.h"
#include "Log.h"

//#define SINGLEGAL 8146
//#define STARTAT 8000
//#define ENDAT 200

static void MeasureSingleShear1(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans, const std::vector<BVec>& psf,
    double noise, double gain, const Image<double>* weight_im, 
    double gal_aperture, double max_aperture,
    int gal_order, int gal_order2,
    double f_psf, double min_gal_size,
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
    flag |= NATIVE_FAILED;
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
    flag |= TOOSMALL;
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
    flag |= DECONV_FAILED;
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
    flag |= SHEAR_FAILED;
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

static void MeasureSingleShear(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans, const FittedPSF& fitpsf,
    double noise, double gain, const Image<double>* weight_im, 
    double gal_aperture, double max_aperture,
    int gal_order, int gal_order2,
    double f_psf, double min_gal_size,
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
      BVec(fitpsf.GetOrder(), fitpsf.GetSigma()));
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
	f_psf, min_gal_size, 
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
  // Read some needed parameters
  int ngals = pos.size();
  dbg<<"ngals = "<<ngals<<std::endl;

  Assert(params.keyExists("shear_aperture"));
  double gal_aperture = params["shear_aperture"];
  double max_aperture = 0.;
  if (params.keyExists("shear_max_aperture")) 
    max_aperture = params["shear_max_aperture"];

  Assert(params.keyExists("shear_gal_order"));
  int gal_order = params["shear_gal_order"];
  int gal_order2 = gal_order;
  if (params.keyExists("shear_gal_order2")) 
    gal_order2 = params["shear_gal_order2"];

  Assert(params.keyExists("shear_f_psf"));
  double f_psf = params["shear_f_psf"];

  double gain = 0.;
  if (params.keyExists("image_gain")) gain = params["image_gain"];

  Assert(params.keyExists("shear_min_gal_size"));
  double min_gal_size = params["shear_min_gal_size"];

  bool output_dots=false;
  if (params.keyExists("output_dots")) output_dots=true;
  bool timing=false;
  if (params.keyExists("timing")) timing=true;

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
      for(int i=0;i<ngals;i++) {
#ifdef STARTAT
	if (i < STARTAT) continue;
#endif
#ifdef SINGLEGAL
	if (i < SINGLEGAL) continue;
	if (i > SINGLEGAL) break;
#endif
#ifdef _OPENMP
#pragma omp critical
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
	    f_psf, min_gal_size, 
	    // Time stats if desired:
	    timing ? &times : 0, 
	    // Log information
	    log1,
	    // Ouput values:
	    shear[i], cov[i], shape[i], flag1);

	flags[i] = flag1;

	if (!flag1) {
	  dbg<<"Successful shear measurement: "<<shear[i]<<std::endl;
	}
	else {
	  dbg<<"Unsuccessful shear measurement\n"; 
	}

	if (timing) {
	  dbg<<"So far: ns = "<<times.ns_gamma<<",  nf = "<<times.nf_native;
	  dbg<<", "<<times.nf_mu<<", "<<times.nf_gamma<<std::endl;
	}

      }
#ifdef _OPENMP
#pragma omp critical
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

  int nsuccess = log.ns_gamma;
  dbg<<nsuccess<<" successful shear measurements, ";
  dbg<<ngals-nsuccess<<" unsuccessful.\n";

  if (timing) {
    dbg<<alltimes.ns_gamma<<" successful shear measurements, ";
    dbg<<alltimes.nf_native<<" + "<<alltimes.nf_mu;
    dbg<<" + "<<alltimes.nf_gamma<<" unsuccessful\n";
    std::cerr<<alltimes<<std::endl;
  }
  xdbg<<log<<std::endl;

  if (output_dots) { 
    std::cerr
      <<std::endl
      <<"Success rate: "<<nsuccess<<"/"<<ngals
      <<std::endl; 
  }

  return nsuccess;
}

ShearCatalog::ShearCatalog(const InputCatalog& incat,
    const Transformation& trans,
    const ConfigFile& _params, std::string key_prefix) :
  id(incat.id), pos(incat.pos), sky(incat.sky), noise(incat.noise),
  flags(incat.flags),
  params(_params), prefix(key_prefix)
{
  dbg<<"Create ShearCatalog\n";

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

  shear.resize(id.size(),0.);
  nu.resize(id.size(),DEFVALNEG);

  tmv::SmallMatrix<double,2,2> cov_default;
  cov_default = tmv::ListInit, DEFVALPOS, 0, 0, DEFVALPOS;
  cov.resize(id.size(),cov_default);

  int gal_order = params.get(prefix + "gal_order");
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

ShearCatalog::ShearCatalog(const ConfigFile& _params, std::string key_prefix) :
  params(_params), prefix(key_prefix)
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

  FitsFile fits(file, READWRITE, true);

  int ncoeff = shape[0].size();
  dbg<<"ncoeff = "<<ncoeff<<std::endl;

  std::stringstream coeff_form;
  coeff_form << ncoeff << "d";

  std::string id_col=params.get(prefix + "id_col");
  std::string x_col=params.get(prefix + "x_col");
  std::string y_col=params.get(prefix + "y_col");
  std::string sky_col=params.get(prefix + "sky_col");
  std::string noise_col=params.get(prefix + "noise_col");
  std::string flags_col=params.get(prefix + "flags_col");
  std::string ra_col=params.get(prefix + "ra_col");
  std::string dec_col=params.get(prefix + "dec_col");
  std::string shear1_col=params.get(prefix + "shear1_col");
  std::string shear2_col=params.get(prefix + "shear2_col");
  std::string nu_col=params.get(prefix + "nu_col");
  std::string cov00_col=params.get(prefix + "cov00_col");
  std::string cov01_col=params.get(prefix + "cov01_col");
  std::string cov11_col=params.get(prefix + "cov11_col");
  std::string order_col=params.get(prefix + "order_col");
  std::string sigma_col=params.get(prefix + "sigma_col");
  std::string coeffs_col=params.get(prefix + "coeffs_col");

  const int nfields=17;
  std::string table_cols[nfields] = {
    id_col,
    x_col,
    y_col,
    sky_col,
    noise_col,
    flags_col,
    ra_col,
    dec_col,
    shear1_col,
    shear2_col,
    nu_col,
    cov00_col,
    cov01_col,
    cov11_col,
    order_col,
    sigma_col,
    coeffs_col,
  };
  std::string table_types[nfields] = {
    "1j", // id
    "1d", // x
    "1d", // y
    "1d", // sky
    "1d", // noise
    "1i", // flags
    "1d", // ra
    "1d", // dec
    "1d", // shear1
    "1d", // shear2
    "1d", // nu
    "1d", // cov00
    "1d", // cov01
    "1d", // cov11
    "1j", // order
    "1d", // sigma
    coeff_form.str()
  };
  std::string table_units[nfields] = {
    "None",
    "Pixels",
    "Pixels",
    "ADU",
    "ADU^2",
    "None",
    "Deg",
    "Deg",
    "None",
    "None",
    "None",
    "None",
    "None",
    "None",
    "None",
    "Arcsec",
    "None"
  };

  // Create a binary table
  fits.CreateBinaryTable(size(),nfields,table_cols,table_types,table_units);

  // Write the header keywords
  fits.WriteParKey(params, "version", TSTRING);
  fits.WriteParKey(params, "noise_method", TSTRING);
  fits.WriteParKey(params, "dist_method", TSTRING);

  fits.WriteParKey(params, prefix + "aperture", TDOUBLE);
  fits.WriteParKey(params, prefix + "max_aperture", TDOUBLE);
  fits.WriteParKey(params, prefix + "gal_order", TLONG);
  fits.WriteParKey(params, prefix + "gal_order2", TLONG);
  fits.WriteParKey(params, prefix + "min_gal_size", TDOUBLE);
  fits.WriteParKey(params, prefix + "f_psf", TDOUBLE);

  fits.WriteColumn(TLONG, 1, 1, 1, size(), &id[0]);

  std::vector<double> x(size());
  std::vector<double> y(size());
  for(size_t i=0;i<size();i++) { 
    x[i] = pos[i].GetX();
    y[i] = pos[i].GetY();
  }
  fits.WriteColumn(TDOUBLE, 2, 1, 1, size(), &x[0]);
  fits.WriteColumn(TDOUBLE, 3, 1, 1, size(), &y[0]);

  fits.WriteColumn(TDOUBLE, 4, 1, 1, size(), &sky[0]);
  fits.WriteColumn(TDOUBLE, 5, 1, 1, size(), &noise[0]);
  fits.WriteColumn(TLONG, 6, 1, 1, size(), &flags[0]);

  std::vector<double> ra(size());
  std::vector<double> dec(size());
  for(size_t i=0;i<size();i++) { 
    ra[i] = skypos[i].GetX();
    dec[i] = skypos[i].GetY();
  }
  fits.WriteColumn(TDOUBLE, 7, 1, 1, size(), &ra[0]);
  fits.WriteColumn(TDOUBLE, 8, 1, 1, size(), &dec[0]);

  std::vector<double> shear1(size());
  std::vector<double> shear2(size());
  for(size_t i=0;i<size();i++) { 
    shear1[i] = real(shear[i]);
    shear2[i] = imag(shear[i]);
  }
  fits.WriteColumn(TDOUBLE, 9, 1, 1, size(), &shear1[0]);
  fits.WriteColumn(TDOUBLE, 10, 1, 1, size(), &shear2[0]);
  fits.WriteColumn(TDOUBLE, 11, 1, 1, size(), &nu[0]);

  std::vector<double> cov00(size());
  std::vector<double> cov01(size());
  std::vector<double> cov11(size());
  for(size_t i=0;i<size();i++) { 
    cov00[i] = cov[i](0,0);
    cov01[i] = cov[i](0,1);
    cov11[i] = cov[i](1,1);
  }
  fits.WriteColumn(TDOUBLE, 12, 1, 1, size(), &cov00[0]);
  fits.WriteColumn(TDOUBLE, 13, 1, 1, size(), &cov01[0]);
  fits.WriteColumn(TDOUBLE, 14, 1, 1, size(), &cov11[0]);
  
  for (size_t i=0; i<size(); i++) {
    size_t row = i+1;
    long b_order = shape[i].GetOrder();
    double b_sigma = shape[i].GetSigma();
    fits.WriteColumn(TLONG, 15, row, 1, 1, &b_order);
    fits.WriteColumn(TDOUBLE, 16, row, 1, 1, &b_sigma);
    fits.WriteColumn(TDOUBLE, 17, row, 1, ncoeff, shape[i].cptr());
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

  for(size_t i=0;i<size();i++) {
    fout
      << id[i] << delim
      << pos[i].GetX() << delim
      << pos[i].GetY() << delim
      << sky[i] << delim
      << noise[i] << delim
      << flags[i] << delim
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
  std::string file = Name(params,"shear");
  dbg<<"Writing to shear file: "<<file<<std::endl;

  if (file.find("fits") != std::string::npos) {
    WriteFits(file);
  } else {
    std::string delim = "  ";
    if (params.keyExists(prefix + "delim")) delim = params[prefix + "delim"];
    WriteAscii(file,delim);
  }
  dbg<<"Done Write\n";
}

void ShearCatalog::ReadFits(std::string file)
{
  int hdu = 2;
  if (params.keyExists(prefix + "hdu")) hdu = params[prefix + "hdu"];

  FitsFile fits(file);

  dbg<<"Moving to HDU #"<<hdu<<std::endl;
  fits.GotoHDU(hdu);

  long nrows=0;
  fits.ReadKey("NAXIS2", TLONG, &nrows);

  dbg<<"  nrows = "<<nrows<<std::endl;
  if (nrows <= 0) {
    throw std::runtime_error("nrows must be > 0");
  }

  std::string id_col=params.get(prefix + "id_col");
  std::string x_col=params.get(prefix + "x_col");
  std::string y_col=params.get(prefix + "y_col");
  std::string sky_col=params.get(prefix + "sky_col");
  std::string noise_col=params.get(prefix + "noise_col");
  std::string flags_col=params.get(prefix + "flags_col");
  std::string ra_col=params.get(prefix + "ra_col");
  std::string dec_col=params.get(prefix + "dec_col");
  std::string shear1_col=params.get(prefix + "shear1_col");
  std::string shear2_col=params.get(prefix + "shear2_col");
  std::string nu_col=params.get(prefix + "nu_col");
  std::string cov00_col=params.get(prefix + "cov00_col");
  std::string cov01_col=params.get(prefix + "cov01_col");
  std::string cov11_col=params.get(prefix + "cov11_col");
  std::string order_col=params.get(prefix + "order_col");
  std::string sigma_col=params.get(prefix + "sigma_col");
  std::string coeffs_col=params.get(prefix + "coeffs_col");

  dbg<<"Reading columns"<<std::endl;
  dbg<<"  "<<id_col<<std::endl;
  id.resize(nrows);
  fits.ReadScalarCol(id_col,TLONG,&id[0], nrows);

  dbg<<"  "<<x_col<<"  "<<y_col<<std::endl;
  pos.resize(nrows);
  std::vector<double> x(nrows);
  std::vector<double> y(nrows);
  fits.ReadScalarCol(x_col,TDOUBLE,&x[0], nrows);
  fits.ReadScalarCol(y_col,TDOUBLE,&y[0], nrows);
  for(long i=0;i<nrows;++i) pos[i] = Position(x[i],y[i]);

  dbg<<"  "<<sky_col<<std::endl;
  fits.ReadScalarCol(sky_col,TDOUBLE,&sky[0], nrows);

  dbg<<"  "<<noise_col<<std::endl;
  fits.ReadScalarCol(noise_col,TDOUBLE,&noise[0], nrows);

  dbg<<"  "<<flags_col<<std::endl;
  flags.resize(nrows,0);
  fits.ReadScalarCol(flags_col,TLONG,&flags[0], nrows);

  dbg<<"  "<<ra_col<<"  "<<dec_col<<std::endl;
  skypos.resize(nrows);
  std::vector<double> ra(nrows);
  std::vector<double> dec(nrows);
  fits.ReadScalarCol(ra_col,TDOUBLE,&ra[0], nrows);
  fits.ReadScalarCol(dec_col,TDOUBLE,&dec[0], nrows);
  for(long i=0;i<nrows;++i) skypos[i] = Position(ra[i],dec[i]);

  dbg<<"  "<<shear1_col<<"  "<<shear2_col<<std::endl;
  shear.resize(nrows);
  std::vector<double> shear1(nrows);
  std::vector<double> shear2(nrows);
  fits.ReadScalarCol(shear1_col,TDOUBLE,&shear1[0], nrows);
  fits.ReadScalarCol(shear2_col,TDOUBLE,&shear2[0], nrows);
  for(long i=0;i<nrows;++i) 
    shear[i] = std::complex<double>(shear1[i],shear2[i]);

  dbg<<"  "<<nu_col<<std::endl;
  nu.resize(nrows,0);
  fits.ReadScalarCol(nu_col,TDOUBLE,&nu[0], nrows);

  dbg<<"  "<<cov00_col<<"  "<<cov01_col<<"  "<<cov11_col<<std::endl;
  cov.resize(nrows);
  std::vector<double> cov00(nrows);
  std::vector<double> cov01(nrows);
  std::vector<double> cov11(nrows);
  fits.ReadScalarCol(cov00_col,TDOUBLE,&cov00[0], nrows);
  fits.ReadScalarCol(cov01_col,TDOUBLE,&cov01[0], nrows);
  fits.ReadScalarCol(cov11_col,TDOUBLE,&cov11[0], nrows);
  for(long i=0;i<nrows;++i) 
    cov[i] = tmv::ListInit, cov00[i], cov01[i], cov01[i], cov11[i];

  shape.reserve(nrows);
  for (size_t i=0; i<size(); i++) {
    size_t row=i+1;
    double b_sigma;
    long b_order;
    fits.ReadCell(order_col,TLONG,&b_order,row,1);
    fits.ReadCell(sigma_col,TDOUBLE,&b_sigma,row,1);
    shape.push_back(BVec(b_order,b_sigma));
    int ncoeff=(b_order+1)*(b_order+2)/2;
    fits.ReadCell(coeffs_col,TDOUBLE,shape[i].ptr(),row,ncoeff);
  }
}

void ShearCatalog::ReadAscii(std::string file, std::string delim)
{
  std::ifstream fin(file.c_str());
  if (!fin) {
    throw std::runtime_error("Error opening stars file");
  }

  if (delim == "  ") {
    long id1,flag,b_order;
    double x,y,sky1,noise1,ra,dec,s1,s2,nu1,b_sigma,c00,c01,c11;
    while ( fin >> id1 >> x >> y >> sky1 >> noise1 >> flag >>
	ra >> dec >> s1 >> s2 >> nu1 >>
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
  dbg<< "Reading PSF cat from file: " << file << std::endl;

  if (file.find("fits") != std::string::npos) {
    ReadFits(file);
  } else {
    std::string delim = "  ";
    if (params.keyExists(prefix + "delim")) delim = params[prefix + "delim"];
    ReadAscii(file,delim);
  }
  dbg<<"Done Read ShearCatalog\n";
}

