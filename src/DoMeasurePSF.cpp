
#include "Params.h"
#include "types.h"

#include "BVec.h"
#include "Ellipse.h"
#include "dbg.h"
#include "ConfigFile.h"
#include "Pixel.h"
#include "Image.h"
#include "FittedPSF.h"
#include "Transformation.h"
#include "DoMeasure.h"
#include "PsiHelper.h"
#include "Name.h"
#include "Log.h"

#include <fstream>
#include <iostream>

#ifdef _OPENMP
#undef _OPENMP
//#include <omp.h>
#endif

//#define SINGLESTAR 91
//#define NSTARS 10
//#define STARTAT 85
//#define ENDAT 95

int DoMeasurePSF(ConfigFile& params, PSFLog& log) 
{
  xdbg<<"Start DoMeasurePSF\n";

  // Load image:
  int image_hdu = 1;
  if (params.keyExists("image_hdu")) image_hdu = params["image_hdu"];
  Image<double> im(Name(params,"image",true),image_hdu);
  xdbg<<"Opened image "<<Name(params,"image",true)<<std::endl;

  // Read catalog info
  // (Also calculates the noise or opens the noise image as appropriate)
  std::vector<Position> all_pos;
  std::vector<double> all_sky;
  std::vector<double> all_noise;
  double gain;
  Image<double>* weight_im = 0;
  ReadCatalog(params,"starcat",all_pos,all_sky,all_noise,gain,weight_im);
  xdbg<<"Done read catalog "<<Name(params,"starcat",false,true)<<std::endl;

  // Fix sky if necessary
  if (all_sky.size() == 0) {
    size_t glob_sky = 0.;
    if (params.keyExists("sky")) glob_sky = params["sky"];
    else glob_sky = im.Median();
    dbg<<"Set global value of sky to "<<glob_sky<<std::endl;
    all_sky.resize(all_pos.size());
    fill(all_sky.begin(),all_sky.end(),glob_sky);
  }

  // Read distortion function
  Transformation trans(params);

  // Read some needed parameters
  int nstars = all_pos.size();
  dbg<<"nstars = "<<nstars<<std::endl;
  Assert(params.keyExists("psf_aperture"));
  double psfap = double(params["psf_aperture"]); 
  dbg<<"psfap = "<<psfap<<std::endl;
  Assert(params.keyExists("psf_order"));
  int psforder = params["psf_order"];
  bool output_dots=false;
  if (params.keyExists("output_dots")) output_dots=true;

  // Initial sigma_p for shapelet measurements
  double sigma_p = 1.;
  if (params.keyExists("seeing_est")) {
    double seeing = params["seeing_est"];
    // seeing is given as FWHM
    // for a gaussian 0.5 = exp(-((FWHM/2)/sigma)^2/2)
    // FWHM/sigma = 2.*sqrt(2 ln(2)) = 2.35
    sigma_p = seeing / 2.35;
  }

  // Calculate a good value of sigma to use:
  // (For this calculation, psfap is psf_aperture * 1 arcsec.)
  EstimateSigma(sigma_p,
      im,all_pos,all_sky,all_noise,gain,weight_im,trans,psfap);
  dbg<<"sigma_p = "<<sigma_p<<std::endl;
  psfap *= sigma_p;  // arcsec
  dbg<<"psfap => "<<psfap<<std::endl;

  // Setup output vectors
  std::vector<BVec> psf(nstars,BVec(psforder,sigma_p));
  std::vector<double> nu(nstars,0.);

  // Set up a default psf vector for output when an object measurement
  // fails
  BVec psf_default(psforder,sigma_p);
  for (size_t i=0; i<psf_default.size(); i++) {
    psf_default[i] = DEFVALNEG;
  }
  double nu_default = DEFVALNEG;

  // Array of flag values
  std::vector<int32> flagvec(nstars,0);

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

  // Main loop to measure psf shapelets:
#ifdef _OPENMP
#pragma omp parallel 
  { 
    try {
#pragma omp for schedule(guided)
#endif
      for(int i=0;i<nstars;i++) {
#ifdef STARTAT
	if (i < STARTAT) continue;
#endif
#ifdef SINGLESTAR
	if (i < SINGLESTAR) continue;
	if (i > SINGLESTAR) break;
	XDEBUG = true;
#endif
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  if (output_dots) { std::cerr<<"."; std::cerr.flush(); }
	  dbg<<"star "<<i<<":\n";
	}

	BVec psf1 = psf_default;
	double nu1 = nu_default;
	int32 flag1=0;
	try {
	  MeasureSinglePSF(
	      // Input data:
	      all_pos[i], im, all_sky[i], trans, 
	      // Noise values:
	      all_noise[i], gain, weight_im,
	      // Parameters:
	      sigma_p, psfap, psforder,
	      // Log information
	      log,
	      // Ouput value:
	      psf1, nu1, flag1);
	} catch (tmv::Error& e) {
	  dbg<<"TMV Error thrown in MeasureSinglePSF\n";
	  log.nf_tmverror++;
	  flag1 |= MPSF_TMV_EXCEPTION;
	} catch (...) {
	  dbg<<"unkown exception in MeasureSinglePSF\n";
	  log.nf_othererror++;
	  flag1 |= MPSF_UNKNOWN_EXCEPTION;
	}
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  flagvec[i] = flag1;
	  psf[i] = psf1;
	  nu[i] = nu1;
	  if (!flag1) {
	    dbg<<"Successful psf measurement: "<<psf1<<std::endl;
	  }
	  else {
	    dbg<<"Unsuccessful psf measurement\n"; 
	  }
	}
#ifdef SINGLESTAR
	exit(1);
#endif
      }
#ifdef _OPENMP
    }
    catch (...)
    { std::cerr<<"Caught some error in parallel region\n"; exit(1); }
  }
#endif


  int nsuccess = log.ns_psf;
  //for(int i=0;i<nstars;i++) if (!flagvec[i]) nsuccess++;

  dbg<<nsuccess<<" successful star measurements, ";
  dbg<<nstars-nsuccess<<" unsuccessful\n";

  // Output psf information:
  std::string psffile = Name(params,"psf");
  std::string psfdelim = "  ";
  if (params.keyExists("psf_delim")) psfdelim = params["psf_delim"];
  std::ofstream catout(psffile.c_str());
  Assert(catout);
  //catout << psforder <<"  "<< sigma_p <<std::endl;
  for(int i=0;i<nstars;i++) {
    DoMeasurePSFPrint(
	catout, 
	all_pos[i].GetX(),
	all_pos[i].GetY(),
	flagvec[i],
	nu[i],
	psforder, sigma_p, psf[i],
	psfdelim);
  }
  dbg<<"Done writing output psf catalog\n";

  // Fit the PSF with a polynomial:
  FittedPSF fittedpsf(psf,flagvec,all_pos,nu,sigma_p,params);
  dbg<<"Done fitting PSF\n";

  // Output fitted psf
  std::string fitpsffile = Name(params,"fitpsf");
  std::ofstream fitout(fitpsffile.c_str());
  fitout << fittedpsf;
  fitout.close();
  dbg<<"Done writing fitted PSF file\n";

  if (XDEBUG) {
    // Check fit:

    std::ifstream readfitout(fitpsffile.c_str());
    FittedPSF readfittedpsf(readfitout);
    for(int i=0;i<nstars;i++) if (!flagvec[i]) {
      xdbg<<"psf[i] = "<<psf[i]<<std::endl;
      BVec checkpsf(readfittedpsf.GetOrder(),readfittedpsf.GetSigma());
      checkpsf = readfittedpsf(all_pos[i]);
      xdbg<<"fittedpsf = "<<checkpsf<<std::endl;
      xdbg<<"Norm(diff) = "<<Norm(psf[i]-checkpsf)<<std::endl;
    }
  }

  dbg<<log<<std::endl;

  if (output_dots) { 
	  std::cerr
		  <<std::endl
		  <<"Success rate: "<<nsuccess<<"/"<<nstars
		  <<std::endl; 
  }

  return nsuccess;
}


void DoMeasurePSFPrint(
    std::ofstream& ostream,
    double x, double y, 
    int32 flags, 
    double nu, 
    int psforder, double sigma_p, const BVec& psf,
    const std::string& delim)
{
  ostream
    << x        << delim
    << y        << delim
    << flags    << delim
    << nu       << delim
    << psforder << delim
    << sigma_p  << delim;
  for(size_t i=0;i<psf.size();i++) ostream << psf[i] << delim;
  ostream << std::endl;
}


