
#include "BVec.h"
#include "Ellipse.h"
#include "dbg.h"
#include "ConfigFile.h"
#include "Pixel.h"
#include "Image.h"
#include "FittedPSF.h"
#include "DoMeasure.h"
#include "TimeVars.h"
#include "PsiHelper.h"
#include "Name.h"

#include <fstream>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif


// The INVTRANS sections aren't needed for anything.  They just
// test the inverse transformation routines.
//#define INVTRANS

//#define SINGLEGAL 8146
//#define STARTAT 8000
//#define ENDAT 200

int DoMeasureShear(ConfigFile& params) 
{
  // Load image:
  Image<double> im(Name(params,"fits"));

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
  if (params.keyExists("dist_ext") || params.keyExists("dist_file")) {
    std::string distfile = Name(params,"dist");
    std::ifstream distin(distfile.c_str());
    Assert(distin);
    distin >> trans;
  } // else stay with identity transformation.


#ifdef INVTRANS
  Transformation invtrans;
  Bounds b;
  for(size_t i=0;i<all_pos.size();i++) b += all_pos[i];
  int invtrans_order = 3;
  invtrans.MakeInverseOf(trans,b,invtrans_order);
#endif

  // Read the fitted psf file
  std::string psffile = Name(params,"fitpsf");
  xdbg<<"Read fitted psf file "<<psffile<<std::endl;
  std::ifstream psfin(psffile.c_str());
  Assert(psfin);
  FittedPSF fittedpsf(psfin);
  xdbg<<"Done reading fittedpsf\n";

  // Read some needed parameters
  int ngals = all_pos.size();
  dbg<<"ngals = "<<ngals<<std::endl;
  Assert(params.keyExists("gal_aperture"));
  double gal_aperture = params["gal_aperture"];
  double max_aperture = 0.;
  if (params.keyExists("max_aperture")) 
    max_aperture = params["max_aperture"];
  Assert(params.keyExists("gal_order"));
  int gal_order = params["gal_order"];
  int gal_order2 = gal_order;
  if (params.keyExists("gal_order2")) gal_order2 = params["gal_order2"];
  Assert(params.keyExists("f_psf"));
  double f_psf = params["f_psf"];
  Assert(params.keyExists("min_gal_size"));
  double min_gal_size = params["min_gal_size"];
  bool output_dots=false;
  if (params.keyExists("output_dots")) output_dots=true;
  bool timing=false;
  if (params.keyExists("timing")) timing=true;

  // Setup output vectors
  std::vector<Position> skypos(ngals);
  std::vector<std::complex<double> > shear(ngals,0.);
  std::vector<tmv::Matrix<double> > shearcov(ngals,tmv::Matrix<double>(2,2));
  std::vector<bool> success(ngals,false);
  std::vector<BVec*> shapelet(ngals,(BVec*)(0));
  OverallFitTimes alltimes;

#ifdef ENDAT
  ngals = ENDAT;
#endif
  // Main loop to measure shears
#ifdef _OPENMP
  omp_set_num_threads(2);
#pragma omp parallel 
  {
#pragma omp critical 
    {
      dbg<<"start parallel region: thread "<<omp_get_thread_num()<<" of  "<<omp_get_num_threads()<<" threads  (max = "<<omp_get_max_threads()<<")"<<std::endl;
    }
    try {
#endif
      OverallFitTimes times; // just for this thread
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
      for(int i=0;i<ngals;i++) {
	if (output_dots) { std::cout<<"."; std::cout.flush(); }
#ifdef STARTAT
	if (i < STARTAT) continue;
#endif
#ifdef SINGLEGAL
	if (i < SINGLEGAL) continue;
	if (i > SINGLEGAL) break;
#endif
	dbg<<"galaxy "<<i<<":\n";
	dbg<<"all_pos[i] = "<<all_pos[i]<<std::endl;

	// Get coordinates of the galaxy, and convert to sky coordinates
	try {
	  trans.Transform(all_pos[i],skypos[i]);
	  dbg<<"skypos = "<<skypos[i]<<std::endl;
#ifdef INVTRANS
	  Position p2;
	  invtrans.Transform(skypos[i],p2);
	  dbg<<"p2 = "<<p2<<std::endl;
	  Assert(std::abs(p2-all_pos[i]) < 1.);
	  trans.InverseTransform(skypos[i],p2);
	  dbg<<"p2 = "<<p2<<std::endl;
	  Assert(std::abs(p2-all_pos[i]) < 1.e-5);
#endif
	} catch (Range_error& e) {
	  dbg<<"skip: distortion range error: \n";
	  xdbg<<"p = "<<all_pos[i]<<", b = "<<e.b<<std::endl;
	  times.nf_range1++;
	  continue;
	}

	// Calculate the psf from the fitted-psf formula:
	std::vector<BVec> psf(1,BVec(fittedpsf.GetOrder(),fittedpsf.GetSigma()));
	try {
	  dbg<<"for fittedpsf all_pos[i] = "<<all_pos[i]<<std::endl;
	  psf[0] = fittedpsf(all_pos[i]);
	} catch (Range_error& e) {
	  dbg<<"skip: fittedpsf range error: \n";
	  xdbg<<"p = "<<all_pos[i]<<", b = "<<e.b<<std::endl;
	  times.nf_range2++;
	  continue;
	}

	// Do the real meat of the calculation:
	bool success1;
	std::complex<double> shear1 = 0.;
	tmv::Matrix<double> shearcov1(2,2);
	BVec* shapelet1=0;
	dbg<<"measure single shear all_pos[i] = "<<all_pos[i]<<std::endl;
	try {
	  MeasureSingleShear(
	      // Input data:
	      all_pos[i], im, all_sky[i], trans, psf,
	      // Noise variables:
	      all_noise[i], gain, weight_im, 
	      // Parameters:
	      gal_aperture, max_aperture, gal_order, gal_order2, 
	      f_psf, min_gal_size, 
	      // Time stats if desired:
	      timing ? &times : 0,
	      // Ouput values:
	      success1, shear1, shearcov1, shapelet1);
	} catch (tmv::Error& e) {
	  dbg<<"skip: TMV Error thrown in MeasureSingleShear\n";
	  dbg<<e<<std::endl;
	  continue;
	}
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  success[i] = success1;
	  shear[i] = shear1;
	  shearcov[i] = shearcov1;
	  shapelet[i] = shapelet1;
	}

	if (timing) {
	  dbg<<"So far: ns = "<<times.ns_gamma<<",  nf = "<<times.nf_native;
	  dbg<<", "<<times.nf_mu<<", "<<times.nf_gamma<<std::endl;
	}

      }
      if (timing)
#ifdef _OPENMP
#pragma omp critical
#endif
      { 
	alltimes += times;
      }
#ifdef _OPENMP
    } 
    catch (...)
    {
      std::cout<<"Caught some kind of exception in the parallel region.\n";
    }
  }
  // End openmp parallel section.
#endif

  int nsuccess=0;
  if (timing) {
    dbg<<alltimes.ns_gamma<<" successful shear measurements, ";
    dbg<<alltimes.nf_native<<" + "<<alltimes.nf_mu;
    dbg<<" + "<<alltimes.nf_gamma<<" unsuccessful\n";
    nsuccess = alltimes.ns_gamma;
  } else {
    for(int i=0;i<ngals;i++) if (success[i]) nsuccess++;
    dbg<<nsuccess<<" successful shear measurements, ";
  }
  if (output_dots) { std::cout<<nsuccess<<std::endl; }

  // Output shear information:
  std::string outcatfile = Name(params,"outcat");
  std::ofstream catout(outcatfile.c_str());
  Assert(catout);
  for(int i=0;i<ngals;i++) if (success[i]) {
    catout << skypos[i].GetX()<<"   "<<skypos[i].GetY()<<"   ";;
    catout <<"  "<<real(shear[i])<<"   "<<imag(shear[i])<<"   ";
    catout <<shearcov[i](0,0)<<"   "<<shearcov[i](1,1)<<"   ";
    catout <<shearcov[i](0,1)<<std::endl;
  }
  dbg<<"Done writing output catalog\n";
  // TODO: Also output shapelets...

  if (timing) std::cout<<alltimes<<std::endl;

  // Cleanup memory
  if (weight_im) delete weight_im;
  for(int i=0;i<ngals;i++) if (shapelet[i]) delete shapelet[i];

  return nsuccess;
}

