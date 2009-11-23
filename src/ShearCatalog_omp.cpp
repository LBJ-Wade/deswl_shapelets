
#include <iostream>

#include "ShearCatalog.h"
#include "ConfigFile.h"
#include "dbg.h"
#include "Params.h"
#include "Image.h"
#include "FittedPSF.h"
#include "Log.h"

//#define SINGLEGAL 8146
//#define STARTAT 8000
//#define ENDAT 200

int ShearCatalog::MeasureShears(const Image<double>& im,
    const Image<double>* weight_im, const Transformation& trans,
    const FittedPSF& fitpsf, ShearLog& log)
{
  int ngals = pos.size();
  dbg<<"ngals = "<<ngals<<std::endl;

  // Read some needed parameters
  double gal_aperture = params.read<double>("shear_aperture");
  double max_aperture = params.read("shear_max_aperture",0.);
  int gal_order = params.read<int>("shear_gal_order");
  int gal_order2 = params.read("shear_gal_order2",gal_order);
  double f_psf = params.read<double>("shear_f_psf");
  double gain = params.read("image_gain",0.);
  double min_gal_size = params.read<double>("shear_min_gal_size");
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
    try 
    {
#endif
      OverallFitTimes times; // just for this thread
      ShearLog log1(params); // just for this thread
      log1.NoWriteLog();
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
      for(int i=0;i<ngals;i++) 
      {
	if (flags[i]) continue;
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
	dbg<<"galaxy "<<i<<":\n";
	dbg<<"pos[i] = "<<pos[i]<<std::endl;

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
	    shear[i], cov[i], shape[i], nu[i], flag1);

	flags[i] = flag1;

	if (!flag1) 
	{
	  dbg<<"Successful shear measurement: "<<shear[i]<<std::endl;
	}
	else 
	{
	  dbg<<"Unsuccessful shear measurement\n"; 
	}

	if (timing) 
	{
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

  dbg<<log.ns_gamma<<" successful shear measurements, ";
  dbg<<ngals-log.ns_gamma<<" unsuccessful\n";
  log.ngood = std::count(flags.begin(),flags.end(),0);
  dbg<<log.ngood<<" with no flags\n";

  if (output_dots) 
  {
    std::cerr
      <<std::endl
      <<"Success rate: "<<log.ns_gamma<<"/"<<log.ngoodin
      <<"  # with no flags: "<<log.ngood
      <<std::endl;
  }

  if (timing) 
  {
    dbg<<"From timing structure:\n";
    dbg<<alltimes.ns_gamma<<" successful shear measurements, ";
    dbg<<alltimes.nf_native<<" + "<<alltimes.nf_mu;
    dbg<<" + "<<alltimes.nf_gamma<<" unsuccessful\n";
    std::cerr<<alltimes<<std::endl;
  }
  xdbg<<log<<std::endl;

  return log.ns_gamma;
}
