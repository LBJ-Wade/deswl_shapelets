
#include <sstream>

#include "PSFCatalog.h"
#include "dbg.h"
#include "Params.h"

//#define SINGLESTAR 18
//#define NSTARS 10
//#define STARTAT 85
//#define ENDAT 95

double PSFCatalog::EstimateSigma(const Image<double>& im,
    const Image<double>* weight_im, const Transformation& trans)
{
  // Initial sigma_p for shapelet measurements
  double sigma_p = 1.;
  if (params.keyExists("psf_seeing_est")) 
  {
    double seeing = params["psf_seeing_est"];
    // seeing is given as FWHM
    // for a gaussian 0.5 = exp(-((FWHM/2)/sigma)^2/2)
    // FWHM/sigma = 2.*sqrt(2 ln(2)) = 2.35
    sigma_p = seeing / 2.35;
  }

  // Calculate a good value of sigma to use:
  // (For this calculation, psfap is psf_aperture * 1 arcsec.)
  double gain = params.read("gain",0.);
  double psfap = params.read<double>("psf_aperture");
  dbg<<"psfap = "<<psfap<<std::endl;

  int nstars = pos.size();
  double meanmu = 0.;
  int count = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(guided) reduction(+ : meanmu) reduction(+ : count)
#endif
  for (int i=0; i<nstars; i++) if (!flags[i]) 
  {
    dbg<<"use i = "<<i<<std::endl;

    double sigma = sigma_p;
    long flag1 = 0; // Ignore flags set by CalcSigma
    CalcSigma(
	sigma,
	im, pos[i], sky[i], noise[i], gain, weight_im, 
	trans, psfap, flag1);
    // Ignore errors -- just don't add to meanmu
    if (flag1) continue;
    meanmu += log(sigma);
    count++;
  } // End omp parallel for

  if (count < nstars/3) 
  {
    std::ostringstream msgout;
    msgout<<"Too many objects were rejected. \n";
    msgout<<"nstars = "<<nstars<<", but only "<<count<<" successful measurements.\n";
    std::string msg = msgout.str();
    dbg<<msg<<std::endl;
    throw ProcessingError(msg);
  }
  meanmu /= count;
  xdbg<<"meanmu = "<<meanmu<<std::endl;
  xdbg<<"input sigma_p = "<<sigma_p<<std::endl;
  sigma_p = exp(meanmu);
  dbg<<"sigma_p = "<<sigma_p<<std::endl;
  return sigma_p;
}

int PSFCatalog::MeasurePSF(const Image<double>& im,
    const Image<double>* weight_im,
    const Transformation& trans, double sigma_p, PSFLog& log)
{
  // Read some needed parameters
  int psforder = params.read<int>("psf_order");
  bool output_dots = params.read("output_dots",false);
  bool desqa = params.read("des_qa",false);
  double gain = params.read("image_gain",0.);
  double psfap = params.read<double>("psf_aperture");
  dbg<<"psfap = "<<psfap<<std::endl;
  psfap *= sigma_p;  // arcsec
  dbg<<"psfap => "<<psfap<<std::endl;

  int nstars = size();
  dbg<<"nstars = "<<nstars<<std::endl;

#ifdef ENDAT
  nstars = ENDAT;
#endif

  log.nstars = nstars;
#ifdef STARTAT
  log.nstars -= STARTAT;
#endif
#ifdef SINGLESTAR
  log.nstars = 1;
#endif
  log.ngoodin = std::count(flags.begin(),flags.end(),0);
  dbg<<log.ngoodin<<"/"<<log.nstars<<" stars with no input flags\n";

  Assert(nstars<=int(pos.size()));
  Assert(nstars<=int(sky.size()));
  Assert(nstars<=int(noise.size()));
  Assert(nstars<=int(psf.size()));
  Assert(nstars<=int(nu.size()));
  Assert(nstars<=int(flags.size()));
  // Main loop to measure psf shapelets:
#ifdef _OPENMP
#pragma omp parallel 
  {
    try 
    {
#endif
      PSFLog log1(params);  // Just for this thread
      log1.NoWriteLog();
#ifdef _OPENMP
#pragma omp for schedule(guided)
#endif
      for(int i=0;i<nstars;i++) if (!flags[i]) 
      {
#ifdef STARTAT
	if (i < STARTAT) continue;
#endif
#ifdef SINGLESTAR
	if (i < SINGLESTAR) continue;
	if (i > SINGLESTAR) break;
	XDEBUG = true;
#endif
	if (output_dots)
#ifdef _OPENMP
#pragma omp critical (output)
#endif
	{
	  std::cerr<<"."; std::cerr.flush(); 
	}
	dbg<<"star "<<i<<":\n";
	dbg<<"pos["<<i<<"] = "<<pos[i]<<std::endl;

	// Start with an error code of unknown failure, in case
	// something happens that I don't detect as an error.
	flags[i] = UNKNOWN_FAILURE;
	dbg<<"Before MeasureSinglePSF1"<<std::endl;
	long flag1 = 0;
	MeasureSinglePSF(
	    // Input data:
	    pos[i], im, sky[i], trans, 
	    // Noise values:
	    noise[i], gain, weight_im,
	    // Parameters:
	    sigma_p, psfap, psforder, desqa,
	    // Log information
	    log1,
	    // Ouput value:
	    psf[i], nu[i], flag1);
	dbg<<"After MeasureSinglePSF"<<std::endl;

	flags[i] = flag1;
	if (!flag1) 
	{
	  dbg<<"Successful psf measurement: "<<psf[i]<<std::endl;
	}
	else 
	{
	  dbg<<"Unsuccessful psf measurement\n"; 
	}
      }
#ifdef _OPENMP
#pragma omp critical (add_log)
#endif
      {
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

  dbg<<log.ns_psf<<" successful star measurements, ";
  dbg<<nstars-log.ns_psf<<" unsuccessful\n";
  log.ngood = std::count(flags.begin(),flags.end(),0);
  dbg<<log.ngood<<" with no flags\n";

  if (output_dots) 
  {
    std::cerr
      <<std::endl
      <<"Success rate: "<<log.ns_psf<<"/"<<log.ngoodin
      <<"  # with no flags: "<<log.ngood
      <<std::endl;
  }

  return log.ns_psf;
}

