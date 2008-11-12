
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

#include <fstream>
#include <iostream>

#ifdef _OPENMP
//#include <omp.h>
#endif

//#define SINGLESTAR 18
//#define NSTARS 10

int DoMeasurePSF(ConfigFile& params) 
{
  // Load image:
  int image_hdu = 1;
  if (params.keyExists("image_hdu")) image_hdu = params["image_hdu"];
  Image<double> im(Name(params,"image"),image_hdu);

  // Read catalog info
  // (Also calculates the noise or opens the noise image as appropriate)
  std::vector<Position> all_pos;
  std::vector<double> all_sky;
  std::vector<double> all_noise;
  double gain;
  Image<double>* weight_im = 0;
  ReadCatalog(params,"starcat",all_pos,all_sky,all_noise,gain,weight_im);

  // Read distortion function
  Transformation trans;
  if (params.keyExists("dist_ext") || params.keyExists("dist_file")) {
    std::string distfile = Name(params,"dist");
    std::ifstream distin(distfile.c_str());
    Assert(distin);
    distin >> trans;
  } // else stay with identity transformation.


  // Read some needed parameters
  int nstars = all_pos.size();
  dbg<<"nstars = "<<nstars<<std::endl;
  Assert(params.keyExists("psf_aperture"));
  double psfap = double(params["psf_aperture"]); 
  dbg<<"psfap = "<<psfap<<std::endl;
  Assert(params.keyExists("psf_order"));
  int psforder = params["psf_order"];
  bool output_dots=false;
  if (params.keyExists("ouput_dots")) output_dots=true;

  // Calculate a good value of sigma to use:
  // (For this calculation, psfap is psf_aperture * 1 arcsec.)
  double sigma_p = EstimateSigma(
      im,all_pos,all_sky,all_noise,gain,weight_im,trans,psfap);
  dbg<<"sigma_p = "<<sigma_p<<std::endl;
  psfap *= sigma_p;  // arcsec
  dbg<<"psfap => "<<psfap<<std::endl;

  // Setup output vectors
  std::vector<BVec*> psf(nstars,(BVec*)(0));
  std::vector<double> nu(nstars,0.);

  // Set up a default psf vector for output when an object measurement
  // fails
  BVec* psf_default = new BVec(psforder,sigma_p);
  for (int i=0; i<(*psf_default).size(); i++) {
    (*psf_default)[i] = DEFVALNEG;
  }
  double nu_default = DEFVALNEG;

  // Array of flag values
  vector<int32> flagvec(nstars,0);

#ifdef NSTARS
  nstars = NSTARS;
#endif
  // Main loop to measure psf shapelets:
#ifdef _OPENMP
#pragma omp parallel 
  { 
    try {
#pragma omp for schedule(guided)
#endif
      for(int i=0;i<nstars;i++) {
#ifdef SINGLESTAR
	if (i < SINGLESTAR) continue;
	if (i > SINGLESTAR) break;
	XDEBUG = true;
#endif
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  if (output_dots) { std::cout<<"."; std::cout.flush(); }
	  dbg<<"star "<<i<<":\n";
	}

	BVec* psf1(0);
	double nu1 = 0.;
	int32 flags=0;
	try {
	  MeasureSinglePSF(
	      // Input data:
	      all_pos[i], im, all_sky[i], trans, 
	      // Noise values:
	      all_noise[i], gain, weight_im,
	      // Parameters:
	      sigma_p, psfap, psforder,
	      // Ouput value:
	      psf1, nu1, flags);
	} catch (...) {
	  dbg<<"unkown exception in MeasureSingleShear\n";
	  flags |= DMPSF_MSPSF_UNKOWN_EXCEPTION;
	}
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  flagvec[i] = flags;
	  // only copy on success
	  if (psf1) {
	    psf[i] = psf1;
	    nu[i] = nu1;
	    dbg<<"Successful measurement: "<<*psf[i]<<std::endl;
	  } else {
	    dbg<<"Unsuccessful measurement\n"; 
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


  int nsuccess = 0;
  for(int i=0;i<nstars;i++) if (psf[i]) nsuccess++;
  dbg<<nsuccess<<" successful star measurements, ";
  dbg<<nstars-nsuccess<<" unsuccessful\n";

  // Output psf information:
  std::string outcatfile = Name(params,"psf");
  std::ofstream catout(outcatfile.c_str());
  Assert(catout);
  catout << psforder <<"  "<< sigma_p <<std::endl;
  for(int i=0;i<nstars;i++) if (psf[i]) {
    if (flagvec[i] == 0) {
      DoMeasurePSFPrint(
	  catout, 
	  all_pos[i].GetX(),
	  all_pos[i].GetY(),
	  flagvec[i],
	  nu[i],
	  psf[i]);
    } else {
      DoMeasurePSFPrint(
	  catout,
	  all_pos[i].GetX(),
	  all_pos[i].GetY(),
	  flagvec[i],
	  nu_default,
	  psf_default);
    }
  }
  dbg<<"Done writing output catalog\n";

  // Fit the PSF with a polynomial:
  FittedPSF fittedpsf(psf,all_pos,nu,sigma_p,params);
  dbg<<"Done fitting PSF\n";

  // Output fitted psf
  std::string outfitfile = Name(params,"fitpsf");
  std::ofstream fitout(outfitfile.c_str());
  fitout << fittedpsf;
  fitout.close();
  dbg<<"Done writing fitted PSF file\n";

  if (XDEBUG) {
    // Check fit:

    std::ifstream readfitout(outfitfile.c_str());
    FittedPSF readfittedpsf(readfitout);
    for(int i=0;i<nstars;i++) if (psf[i]) {
      xdbg<<"psf[i] = "<<*psf[i]<<std::endl;
      BVec checkpsf(readfittedpsf.GetOrder(),readfittedpsf.GetSigma());
      checkpsf = readfittedpsf(all_pos[i]);
      xdbg<<"fittedpsf = "<<checkpsf<<std::endl;
      xdbg<<"Norm(diff) = "<<Norm(*psf[i]-checkpsf)<<std::endl;
    }
  }

  // Cleanup memory
  for(int i=0;i<nstars;i++) if (psf[i]) delete psf[i];
  delete psf_default;

  return nsuccess;
}


void DoMeasurePSFPrint(
    std::ofstream& ostream,
    double x, double y, 
    int32 flags, 
    double nu, 
    BVec* psf)
{
  ostream
    << x       <<"  "
    << y       <<"  "
    << flags   <<"  "
    << nu      <<"  "
    << *psf
    << std::endl;
}


