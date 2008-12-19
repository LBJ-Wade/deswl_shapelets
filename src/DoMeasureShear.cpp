#include "Params.h"

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
#include "Log.h"

#include <fstream>
#include <iostream>

#include "SXCat.h"

#ifdef _OPENMP
//#include <omp.h>
#endif

//#define SINGLEGAL 8146
//#define STARTAT 8000
//#define ENDAT 200

// simpler io tests
void TestShearIO(ConfigFile& params, SHEAR_STRUCT& shcat)
{

  SHEAR_STRUCT test;
  ReadShearCat(params,test);

  std::stringstream err;

  if (test.id.size() != shcat.id.size()) {
    throw std::runtime_error("id in catalog is wrong size");
  }
  for (size_t i=0; i<test.id.size(); i++) {
    if (test.id[i] != shcat.id[i]) {
      err<<"id["<<i<<"] differs in catalog";
      throw std::runtime_error(err.str());
    }
  }

  if (test.shear_flags.size() != shcat.shear_flags.size()) {
      throw std::runtime_error("field shear_flags in catalog is wrong size");
  }
  for (size_t i=0; i<test.shear_flags.size(); i++) {
    if (test.shear_flags[i] != shcat.shear_flags[i]) {
      err<<"shear_flags["<<i<<"] differs in catalog";
      throw std::runtime_error(err.str());
    }
  }

  if (test.gal_order.size() != shcat.gal_order.size()) {
      throw std::runtime_error("field size_flags in catalog is wrong size");
  }
  for (size_t i=0; i<test.gal_order.size(); i++) {
    if (test.gal_order[i] != shcat.gal_order[i]) {
      err<<"gal_order["<<i<<"] differs in catalog";
      throw std::runtime_error(err.str());
    }
  }

}


// visual tests
static void TestReadShearCat(ConfigFile& params, SHEAR_STRUCT& shcat)
{
  SHEAR_STRUCT shcat_test;
  ReadShearCat(params, shcat_test);

  int gal_order = shcat.gal_order[0];
  int ncoeff = (gal_order+1)*(gal_order+2)/2;
  for (size_t i=0; i<shcat.id.size(); i++) {
    dbg<<i<<"diffs:"
      <<"\tid: "<<shcat.id[i]-shcat_test.id[i]

#ifdef SHEXTRA_PARS
      <<"\tsize_flags: "<<shcat.size_flags[i]-shcat_test.size_flags[i]
      <<"\tstar_flag: "<<shcat.star_flag[i]-shcat_test.star_flag[i]
      <<"\tsigma0: "<<shcat.sigma0[i]-shcat_test.sigma0[i]
#endif

      <<"\tflags: "<<shcat.shear_flags[i]-shcat_test.shear_flags[i]
      <<"\tshear1: "<<shcat.shear1[i]-shcat_test.shear1[i]
      <<"\tshear2: "<<shcat.shear2[i]-shcat_test.shear2[i]

      <<"\tcov00: "<<shcat.shear_cov00[i]-shcat_test.shear_cov00[i]
      <<"\tcov01: "<<shcat.shear_cov01[i]-shcat_test.shear_cov01[i]
      <<"\tcov11: "<<shcat.shear_cov11[i]-shcat_test.shear_cov11[i]

      <<"\torder: "<<shcat.gal_order[i]-shcat_test.gal_order[i]

      <<"\tshapelet[0]: "
      <<shcat.shapelets_prepsf[i][0]-shcat_test.shapelets_prepsf[i][0]

      <<"\tshapelet["<<ncoeff-1<<"]: "
      <<shcat.shapelets_prepsf[i][ncoeff-1]-
         shcat_test.shapelets_prepsf[i][ncoeff-1]

      <<std::endl;
  }
}


// This is now very DES centric
int DoMeasureShear_DES(ConfigFile& params, ShearLog& log) 
{
  // Load image:
  int image_hdu = 1;
  if (params.keyExists("image_hdu")) image_hdu = params["image_hdu"];
  Image<double> im(Name(params,"image",true),image_hdu);

  // Load weight image
  int weight_hdu = 1;
  if (params.keyExists("weight_hdu")) weight_hdu = params["weight_hdu"];
  //Image<double> weight_im(Name(params,"weight",true),weight_hdu);
  Image<double>* weight_im = new Image<double>(Name(params,"weight",true),weight_hdu);
  dbg<<"Opened weight image.\n";


  // The sextractor catalog
  SXCAT_STRUCT sxcat;
  ReadSXCat(params, sxcat);

  // Read distortion function
  Transformation trans(params);

  // These are unused
  double gain=0.0;


  // Read the fitted psf file
  FittedPSF fittedpsf;
  fittedpsf.Read(params);
  xdbg<<"Done reading fittedpsf\n";


  // Read some needed parameters
  int ngals = sxcat.pos.size();
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

  Assert(params.keyExists("shear_min_gal_size"));
  double min_gal_size = params["shear_min_gal_size"];

  bool output_dots=false;
  if (params.keyExists("output_dots")) output_dots=true;
  bool timing=false;
  if (params.keyExists("timing")) timing=true;

  // Setup output vectors
  std::vector<std::complex<double> > shear(ngals,0.);
  std::vector<tmv::Matrix<double> > shearcov(ngals,tmv::Matrix<double>(2,2));
  std::vector<BVec> shapelet(ngals,BVec(gal_order,DEFVALPOS));
  OverallFitTimes alltimes;
  // Vector of flag values
  std::vector<long> flagvec(ngals,0);

  // Default values when we have failure
  std::complex<double> shear_default = DEFVALNEG;

  tmv::Matrix<double> shearcov_default(2,2);
  shearcov_default(0,0) = DEFVALPOS; shearcov_default(0,1) = 0;
  shearcov_default(1,0) = 0;         shearcov_default(1,1) = DEFVALPOS;

  BVec shapelet_default(gal_order,DEFVALPOS);
  for (size_t i=0; i<shapelet_default.size(); i++) {
    shapelet_default[i] = DEFVALNEG;
  }

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
	  dbg<<"pos[i] = "<<sxcat.pos[i]<<std::endl;
	}

	// Inputs to shear measurement code
	std::complex<double> shear1 = shear_default;
	tmv::Matrix<double> shearcov1 = shearcov_default;
	BVec shapelet1 = shapelet_default;
	long flag1 = 0;

	MeasureSingleShear1(
	    // Input data:
	    sxcat.pos[i], im, sxcat.local_sky[i], trans, 
	    // Fitted PSF
	    fittedpsf,
	    // Noise variables:
	    sxcat.noise[i], gain, weight_im, 
	    // Parameters:
	    gal_aperture, max_aperture, gal_order, gal_order2, 
	    f_psf, min_gal_size, 
	    // Time stats if desired:
	    timing ? &times : 0, 
	    // Log information
	    log1,
	    // Ouput values:
	    shear1, shearcov1, shapelet1, flag1);

#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  // We always write out the answer for each object
	  // These are the default values if a catastrophic error occurred.
	  // Otherwise they are the best estimate up to the point of 
	  // any error.
	  // If there is no error, these are the correct estimates.
	  shear[i] = shear1;
	  shearcov[i] = shearcov1;
	  shapelet[i] = shapelet1;
	  flagvec[i] = flag1;
	  if (!flag1) {
	    dbg<<"Successful shear measurement: "<<shear1<<std::endl;
	  }
	  else {
	    dbg<<"Unsuccessful shear measurement\n"; 
	  }
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
      std::cerr<<"STATUS5BEG Caught some error in parallel region STATUS5END\n";
      exit(1);
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
    //for(int i=0;i<ngals;i++) if (flagvec[i] == 0) nsuccess++;
    nsuccess = log.ns_gamma;
    dbg<<nsuccess<<" successful shear measurements, ";
  }
  if (output_dots) { 
	  std::cerr
		  <<std::endl
		  <<"Success rate: "<<nsuccess<<"/"<<ngals
		  <<std::endl; 
  }

  // Output shear information:
  std::string shearfile = Name(params,"shear");





  // Get some info from the findstars run to make up for Joe's stupidity
  FINDSTARS_STRUCT fscat;
  ReadFindStarsCat(params, fscat);

  SHEAR_STRUCT shcat;
  ResizeShearCat(shcat, ngals, gal_order);

  for (int i=0; i<ngals; i++) {
    shcat.id[i] = sxcat.id[i];
    shcat.shear_flags[i] = flagvec[i];

    shcat.shear1[i] = real(shear[i]);
    shcat.shear2[i] = imag(shear[i]);

    shcat.shear_cov00[i] = shearcov[i](0,0);
    shcat.shear_cov01[i] = shearcov[i](0,1);
    shcat.shear_cov11[i] = shearcov[i](1,1);

    shcat.shapelets_prepsf[i] = shapelet[i];

    // copy in
#ifdef SHEXTRA_PARS
    shcat.sigma0[i] = fscat.sigma0[i];
    shcat.size_flags[i] = fscat.size_flags[i];
    shcat.star_flag[i] = fscat.star_flag[i];
#endif


  }

  if (1) {
    WriteShearCat(params, shcat);
    TestShearIO(params, shcat);
    if (0) {
      TestReadShearCat(params, shcat);
    }
  } else {
    std::string sheardelim = "  ";
    if (params.keyExists("shear_delim")) sheardelim = params["shear_delim"];
    std::ofstream catout(shearfile.c_str());
    Assert(catout);
    for(int i=0;i<ngals;i++) {
      DoMeasureShearPrint(
	  catout,
	  sxcat.ra[i], sxcat.dec[i],
	  flagvec[i],
	  shear[i], shearcov[i],
	  sheardelim);
    }
  }
  dbg<<"Done writing output shear catalog\n";
  // TODO: Also output shapelets...

  if (timing) std::cerr<<alltimes<<std::endl;
  xdbg<<log<<std::endl;

  // Cleanup memory
  if (weight_im) delete weight_im;

  return nsuccess;
}



int DoMeasureShear(ConfigFile& params, ShearLog& log) 
{
  // Load image:
  int image_hdu = 1;
  if (params.keyExists("image_hdu")) image_hdu = params["image_hdu"];
  Image<double> im(Name(params,"image",true),image_hdu);

  // Read catalog info
  // (Also calculates the noise or opens the noise image as appropriate)
  std::vector<Position> all_pos;
  std::vector<double> all_sky;
  std::vector<double> all_noise;
  std::vector<Position> all_skypos;
  double gain;
  Image<double>* weight_im = 0;
  ReadCatalog(params,"stars",all_pos,all_sky,all_noise,gain,weight_im,
      all_skypos);
  dbg<<"Finished Read cat\n";

  // Fix sky if necessary
  if (all_sky.size() == 0) {
    double glob_sky = 0.;
    if (params.keyExists("sky")) glob_sky = params["sky"];
    else glob_sky = im.Median();
    dbg<<"Set global value of sky to "<<glob_sky<<std::endl;
    all_sky.resize(all_pos.size());
    fill(all_sky.begin(),all_sky.end(),glob_sky);
  }

  // Read distortion function
  Transformation trans(params);
  Position skypos_default(DEFVALNEG,DEFVALNEG);
  if (all_skypos.size() == 0) {
    all_skypos.resize(all_pos.size(),skypos_default);
    for(size_t i=0;i<all_pos.size();i++) {
      try {
	trans.Transform(all_pos[i],all_skypos[i]);
	all_skypos[i] /= 3600.; // arcsec -> degrees
      } catch (Range_error& e) {
	xdbg<<"distortion range error\n";
	xdbg<<"p = "<<all_pos[i]<<", b = "<<e.b<<std::endl;
      }                       
    }
  } else if (XDEBUG) {
    xdbg<<"Check transformation:\n";
    double rmserror = 0.;
    int count = 0;
    for(size_t i=0;i<all_pos.size();i++) {
      try {
	Position temp;
	trans.Transform(all_pos[i],temp);
	temp /= 3600.; // arcsec -> degrees
	xdbg<<all_pos[i]<<"  "<<all_skypos[i]<<"  "<<temp;
	xdbg<<"  "<<temp-all_skypos[i]<<std::endl;
	rmserror += std::norm(temp-all_skypos[i]);
	count++;
      } catch (Range_error& e) {
	xdbg<<"distortion range error\n";
	xdbg<<"p = "<<all_pos[i]<<", b = "<<e.b<<std::endl;
      }
    }
    rmserror /= count;
    rmserror = std::sqrt(rmserror);
    xdbg<<"rms error = "<<rmserror*3600.<<" arcsec\n";
  }

  // Read the fitted psf file
  FittedPSF fittedpsf;
  fittedpsf.Read(params);
  xdbg<<"Done reading fittedpsf\n";

  // Read some needed parameters
  int ngals = all_pos.size();
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
  Assert(params.keyExists("shear_min_gal_size"));
  double min_gal_size = params["shear_min_gal_size"];
  bool output_dots=false;
  if (params.keyExists("output_dots")) output_dots=true;
  bool timing=false;
  if (params.keyExists("timing")) timing=true;

  // Setup output vectors
  std::vector<std::complex<double> > shear(ngals,0.);
  std::vector<tmv::Matrix<double> > shearcov(ngals,tmv::Matrix<double>(2,2));
  std::vector<BVec> shapelet(ngals,BVec(gal_order,DEFVALPOS));
  OverallFitTimes alltimes;
  // Vector of flag values
  std::vector<long> flagvec(ngals,0);

  // Default values when we have failure
  std::complex<double> shear_default = DEFVALNEG;

  tmv::Matrix<double> shearcov_default(2,2);
  shearcov_default(0,0) = DEFVALPOS; shearcov_default(0,1) = 0;
  shearcov_default(1,0) = 0;         shearcov_default(1,1) = DEFVALPOS;

  BVec shapelet_default(gal_order,DEFVALPOS);
  for (size_t i=0; i<shapelet_default.size(); i++) {
    shapelet_default[i] = DEFVALNEG;
  }

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
	  dbg<<"all_pos[i] = "<<all_pos[i]<<std::endl;
	}

	// Inputs to shear measurement code
	std::complex<double> shear1 = shear_default;
	tmv::Matrix<double> shearcov1 = shearcov_default;
	BVec shapelet1 = shapelet_default;
	long flag1 = 0;

	MeasureSingleShear1(
	    // Input data:
	    all_pos[i], im, all_sky[i], trans, 
	    // Fitted PSF
	    fittedpsf,
	    // Noise variables:
	    all_noise[i], gain, weight_im, 
	    // Parameters:
	    gal_aperture, max_aperture, gal_order, gal_order2, 
	    f_psf, min_gal_size, 
	    // Time stats if desired:
	    timing ? &times : 0, 
	    // Log information
	    log1,
	    // Ouput values:
	    shear1, shearcov1, shapelet1, flag1);

#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  // We always write out the answer for each object
	  // These are the default values if a catastrophic error occurred.
	  // Otherwise they are the best estimate up to the point of 
	  // any error.
	  // If there is no error, these are the correct estimates.
	  shear[i] = shear1;
	  shearcov[i] = shearcov1;
	  shapelet[i] = shapelet1;
	  flagvec[i] = flag1;
	  if (!flag1) {
	    dbg<<"Successful shear measurement: "<<shear1<<std::endl;
	  }
	  else {
	    dbg<<"Unsuccessful shear measurement\n"; 
	  }
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
      std::cerr<<"STATUS5BEG Caught some error in parallel region STATUS5END\n";
      exit(1);
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
    //for(int i=0;i<ngals;i++) if (flagvec[i] == 0) nsuccess++;
    nsuccess = log.ns_gamma;
    dbg<<nsuccess<<" successful shear measurements, ";
  }
  if (output_dots) { 
	  std::cerr
		  <<std::endl
		  <<"Success rate: "<<nsuccess<<"/"<<ngals
		  <<std::endl; 
  }

  // Output shear information:
  std::string shearfile = Name(params,"shear");
  std::string sheardelim = "  ";
  if (params.keyExists("shear_delim")) sheardelim = params["shear_delim"];
  std::ofstream catout(shearfile.c_str());
  Assert(catout);
  for(int i=0;i<ngals;i++) {
    DoMeasureShearPrint(
	catout,
	all_skypos[i].GetX(), all_skypos[i].GetY(),
	flagvec[i],
	shear[i], shearcov[i],
	sheardelim);
  }
  dbg<<"Done writing output shear catalog\n";
  // TODO: Also output shapelets...

  if (timing) std::cerr<<alltimes<<std::endl;
  dbg<<log<<std::endl;

  // Cleanup memory
  if (weight_im) delete weight_im;

  return nsuccess;
}

void DoMeasureShearPrint(
    std::ofstream& ostream,
    double x, double y, 
    long flags, 
    const std::complex<double>& shear, 
    const tmv::Matrix<double>& shearcov,
    const std::string& delim)
{
  ostream
    << x             << delim
    << y             << delim
    << flags         << delim
    << real(shear)   << delim
    << imag(shear)   << delim
    << shearcov(0,0) << delim
    << shearcov(0,1) << delim
    << shearcov(1,1)
    << std::endl;
}

