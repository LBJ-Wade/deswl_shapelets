
#include "ConfigFile.h"
#include "dbg.h"
#include "Params.h"
#include "Log.h"
#include "TimeVars.h"
#include "MultiShearCatalog.h"
#include "Ellipse.h"

int MultiShearCatalog::MeasureMultiShears(const Bounds& b, ShearLog& log)
{
  dbg<<"Start MeasureMultiShears for b = "<<b<<std::endl;
  int ngals = skypos.size();
  dbg<<"ngals = "<<ngals<<std::endl;

  // Read some needed parameters
  double gal_aperture = params.read<double>("shear_aperture");
  double max_aperture = params.read("shear_max_aperture",0.);
  int gal_order = params.read<int>("shear_gal_order");
  int gal_order2 = params.read("shear_gal_order2",gal_order);
  double f_psf = params.read<double>("shear_f_psf");
  double min_gal_size = params.read<double>("shear_min_gal_size");
  bool desqa = params.read("des_qa",false);
  bool output_dots = params.read("output_dots",false);
  bool timing = params.read("timing",false);

  OverallFitTimes alltimes;

  int nsuccess = 0;

#ifdef ENDAT
  ngals = ENDAT;
#endif
  
  // Main loop to measure shears
#ifdef _OPENMP
#pragma omp parallel 
  {
    try {
#endif
      OverallFitTimes times; // just for this thread
      ShearLog log1(params); // just for this thread
      log1.NoWriteLog();
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
      for(int i=0;i<ngals;++i) 
      {
	if (!b.Includes(skypos[i])) 
	{
	  xxdbg<<"Skipping galaxy "<<i<<" because "<<skypos[i]<<"not in bounds\n";
	  continue;
	}
	if (flags[i]) 
	{
	  xxdbg<<"Skipping galaxy "<<i<<" because flag = "<<flags[i]<<std::endl;
	  continue;
	}
#ifdef STARTAT
	if (i < STARTAT) continue;
#endif
#ifdef SINGLEGAL
	if (i < SINGLEGAL) continue;
	if (i > SINGLEGAL) break;
#endif

	if (output_dots)
#ifdef _OPENMP
#pragma omp critical (output)
#endif
	{
	  std::cerr<<"."; std::cerr.flush(); 
	}

	if (pixlist[i].size() == 0) 
	{
	  dbg<<"no valid single epoch images.\n";
	  flags[i] = NO_SINGLE_EPOCH_IMAGES;
	  continue;
	}

	// Start with an error code of unknown failure, in case
	// something happens that I don't detect as an error.
	flags[i] = UNKNOWN_FAILURE;
	long flag1 = 0;

	MeasureMultiShear(
	    // Input data:
	    skypos[i], pixlist[i], psflist[i],
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
	  dbg<<"Successful shear measurement: "<<shear[i]<<std::endl;
	  nsuccess++;
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

  dbg<<nsuccess<<" successful shear measurements in this pass.\n";
  dbg<<log.ns_gamma<<" successful shear measurements so far.\n";

  if (timing) {
    dbg<<"From timing structure:\n";
    dbg<<alltimes.ns_gamma<<" successful shear measurements, ";
    dbg<<alltimes.nf_native<<" + "<<alltimes.nf_mu;
    dbg<<" + "<<alltimes.nf_gamma<<" unsuccessful\n";
    std::cerr<<alltimes<<std::endl;
  }
  xdbg<<log<<std::endl;

  return nsuccess;
}

void MeasureMultiShear(
    const Position& cen, 
    const std::vector<PixelList>& allpix,
    const std::vector<BVec>& psf,
    double gal_aperture, double max_aperture,
    int gal_order, int gal_order2,
    double f_psf, double min_gal_size, bool desqa,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    tmv::SmallMatrix<double,2,2>& shearcov, BVec& shapelet,
    long& flag)
{
  //TODO: This function is too long.  It should really be broken 
  //      up in to smaller units...
  //      In which case most of the overlap between this and 
  //      MeasureSingleShear1 would not have to be duplicated.

  dbg<<"Start MeasureMultiShear\n";
  dbg<<"cen = "<<cen<<std::endl;
  dbg<<"allpix.size = "<<allpix.size()<<std::endl;
  for(size_t i=0;i<allpix.size();++i)
    dbg<<"allpix["<<i<<"].size = "<<allpix[i].size()<<std::endl;
  dbg<<"psf.size = "<<psf.size()<<std::endl;
  Assert(psf.size() == allpix.size());

  try 
  {
    // Find harmonic mean of psf sizes:
    // MJ: Is it correct to use the harmonic mean of sigma^2?
    xdbg<<"sigma_p values are: \n";
    double sigma_p = 0.;
    Assert(psf.size() > 0);
    for(size_t i=0;i<psf.size();++i) 
    {
      xdbg<<psf[i].GetSigma()<<std::endl;
      sigma_p += 1./psf[i].GetSigma();
    }
    sigma_p = double(psf.size()) / sigma_p;
    xdbg<<"Harmonic mean = "<<sigma_p<<std::endl;

    // We start with a native fit using sigma = 2*sigma_p:
    double sigma_obs = 2.*sigma_p;
    double galap = gal_aperture * sigma_obs;
    xdbg<<"galap = "<<gal_aperture<<" * "<<sigma_obs<<" = "<<galap<<std::endl;
    if (max_aperture > 0. && galap > max_aperture) {
      galap = max_aperture;
      xdbg<<"      => "<<galap<<std::endl;
    }
    xdbg<<"sigma_obs = "<<sigma_obs<<", sigma_p = "<<sigma_p<<std::endl;

    std::vector<PixelList> pix(allpix.size());
    int npix = 0;
    for(size_t i=0;i<pix.size();++i) {
      GetSubPixList(pix[i],allpix[i],galap,flag);
      npix += pix[i].size();
    }
    xdbg<<"npix = "<<npix<<std::endl;

    // TODO: If we have single-epoch measurements, use those for starting
    // guess rather than CrudeMeasure.
    Ellipse ell;
    ell.FixGam();
    ell.CrudeMeasure(pix,sigma_obs);
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
    xdbg<<"sigma_obs -> "<<sigma_obs<<std::endl;
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
    xxdbg<<"After native fit #2:\n";
    xxdbg<<"Mu = "<<ell.GetMu()<<std::endl;

    npix = 0;
    for(size_t i=0;i<pix.size();++i) {
      GetSubPixList(pix[i],allpix[i],galap,flag);
      npix += pix[i].size();
    }

    xdbg<<"npix = "<<npix<<std::endl;
    double sigpsq = pow(sigma_p,2);
    double sigma = pow(sigma_obs,2) - sigpsq;
    xdbg<<"sigma_p^2 = "<<sigpsq<<std::endl;
    xdbg<<"sigma^2 = sigma_obs^2 - sigmap^2 = "<<sigma<<std::endl;
    if (sigma < 0.1*sigpsq) sigma = 0.1*sigpsq;
    sigma = sqrt(sigma);
    xdbg<<"sigma = sqrt(sigma^2) "<<sigma<<std::endl;
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
      dbg<<"Mu = "<<ell.GetMu()<<std::endl;
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
      dbg<<"Measured gamma = "<<ell.GetGamma()<<std::endl;
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
      dbg<<"Gamma measurement failed\n";
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
    // TODO: Write function to directly measure shear variance
    // And also the BJ02 sigma_0 value.
    // (I'm not convinced that the above covariance matrix is a good estiamte.)
  } 
  catch (tmv::Error& e) 
  {
    dbg<<"TMV Error thrown in MeasureSingleShear\n";
    dbg<<e<<std::endl;
    log.nf_tmverror++;
    flag |= TMV_EXCEPTION;
  } catch (std::exception) {
    throw;
  } catch (...) {
    dbg<<"unkown exception in MeasureSingleShear\n";
    log.nf_othererror++;
    flag |= UNKNOWN_EXCEPTION;
  } 
}
