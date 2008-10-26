
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
  Assert(params.keyExists("fits_ext"));
  Image<double> im(params["root"] + params["fits_ext"]);

  // Read catalog info
  // (Also calculates the noise or opens the noise image as appropriate)
  std::vector<Position> all_pos;
  std::vector<double> all_sky;
  std::vector<double> all_noise;
  double gain;
  Image<double>* weight_im = 0;
  ReadCatalog(params,all_pos,all_sky,all_noise,gain,weight_im);

  // Read distortion function
  Assert(params.keyExists("root"));
  Transformation trans;
  if (params.keyExists("dist_ext")) {
    std::string distfile = params["root"] + params["dist_ext"];
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
	MeasureSinglePSF(
	    // Input data:
	    all_pos[i], im, all_sky[i], trans, 
	    // Noise values:
	    all_noise[i], gain, weight_im,
	    // Parameters:
	    sigma_p, psfap, psforder,
	    // Ouput value:
	    psf1, nu1);
#ifdef _OPENMP
#pragma omp critical
#endif
	{
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
  Assert(params.keyExists("outcat_ext"));
  std::string outcatfile = params["root"] + params["outcat_ext"];
  std::ofstream catout(outcatfile.c_str());
  Assert(catout);
  catout << psforder <<"  "<< sigma_p <<std::endl;
  for(int i=0;i<nstars;i++) if (psf[i]) {
    catout << all_pos[i].GetX()<<"  "<<all_pos[i].GetY();
    catout << "  "<< nu[i] <<"  "<<psf[i]<<std::endl;
  }
  dbg<<"Done writing output catalog\n";

  // Fit the PSF with a polynomial:
  FittedPSF fittedpsf(psf,all_pos,nu,sigma_p,params);
  dbg<<"Done fitting PSF\n";

  // Output fitted psf
  Assert(params.keyExists("fitpsf_ext"));
  std::string outfitfile = params["root"] + params["fitpsf_ext"];
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

  return nsuccess;
}

