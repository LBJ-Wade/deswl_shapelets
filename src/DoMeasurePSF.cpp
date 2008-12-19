
#include "Params.h"

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

#include "SXCat.h"

#ifdef _OPENMP
#undef _OPENMP
//#include <omp.h>
#endif

//#define SINGLESTAR 18
//#define NSTARS 10
//#define STARTAT 85
//#define ENDAT 95


/*
void WritePSFFits(
    vector<long>& id, 
    vector<int32>& psf_flags,
    vector<double>& nu,
    vector<int32>& psf_order,
    vector<double>& sigmap,
    const BVec& psf)
{

  std::string testfile = "measurepsf.fits";
  FitsFile fits(testfile.c_str(), READWRITE, true);
  fitsfile* fptr = fits.get_fptr();

  int nfields=4;
  char *table_names[] =  
      {(char*)"id",
	(char*)"sigma0",
	(char*)"size_flags",
	(char*)"star_flag"};
  char *table_types[] =  
      {(char*)"1I",
	(char*)"1D",
	(char*)"1I",
	(char*)"1I"};
  char *table_units[] =  
      {(char*)"None",
	(char*)"pixels",
	(char*)"None",
	(char*)"None"};


  // Create a binary table
  int fits_status=0;
  int tbl_type = BINARY_TBL;
  fits_create_tbl(fptr, tbl_type, sxcat.id.size(), nfields, 
      table_names, table_types, table_units, NULL, &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error creating fits table: ";
    throw FitsException(serr);
  }


  // need a more general way to write columns
  fits_status=0;
  fits_write_col(
      fptr, TLONG, 1, 1, 1, sxcat.id.size(), &sxcat.id[0], &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing 'id' to fits column";
    throw FitsException(serr);
  }

  fits_status=0;
  fits_write_col(
      fptr, TDOUBLE, 2, 1, 1, sxcat.sigma.size(), &sxcat.sigma[0], 
      &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing 'sigma' to fits column: ";
    throw FitsException(serr);
  }

  fits_status=0;
  fits_write_col(
      fptr, TINT, 3, 1, 1, sxcat.id.size(), &sxcat.size_flags[0], 
      &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing 'size_flags' to fits column";
    throw FitsException(serr);
  }

  fits_status=0;
  fits_write_col(
      fptr, TINT, 4, 1, 1, sxcat.id.size(), &sxcat.star_flag[0], 
      &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing 'star_flag' to fits column";
    throw FitsException(serr);
  }




  fits.Close();

}

*/


static void TestPSFCatIO(ConfigFile& params, PSF_STRUCT& psfcat)
{

  PSF_STRUCT test;
  ReadPSFCat(params,test);

  std::stringstream err;

  if (test.id.size() != psfcat.id.size()) {
    throw std::runtime_error("id in catalog is wrong size");
  }
  for (size_t i=0; i<test.id.size(); i++) {
    if (test.id[i] != psfcat.id[i]) {
      err<<"id["<<i<<"] differs in catalog";
      throw std::runtime_error(err.str());
    }
  }

  if (test.psf_flags.size() != psfcat.psf_flags.size()) {
      throw std::runtime_error("field psf_flags in catalog is wrong size");
  }
  for (size_t i=0; i<test.psf_flags.size(); i++) {
    if (test.psf_flags[i] != psfcat.psf_flags[i]) {
      err<<"psf_flags["<<i<<"] differs in catalog";
      throw std::runtime_error(err.str());
    }
  }

  if (test.psf_order.size() != psfcat.psf_order.size()) {
      throw std::runtime_error("field size_flags in catalog is wrong size");
  }
  for (size_t i=0; i<test.psf_order.size(); i++) {
    if (test.psf_order[i] != psfcat.psf_order[i]) {
      err<<"psf_order["<<i<<"] differs in catalog";
      throw std::runtime_error(err.str());
    }
  }

}



static void ExtractStars(
    SXCAT_STRUCT& sxcat, 
    FINDSTARS_STRUCT& fscat, 
    FINDSTARS_STRUCT& starcat) 
{

  ResizeFindStarsCat(starcat, 0);

  dbg<<"Extract Stars:\n";
  dbg<<"ntot = "<<fscat.id.size()<<std::endl;
  Assert(fscat.star_flag.size() == fscat.id.size());
  Assert(fscat.sigma0.size() == fscat.id.size());
  Assert(fscat.size_flags.size() == fscat.id.size());
  Assert(sxcat.local_sky.size() == fscat.id.size());
  Assert(sxcat.noise.size() == fscat.id.size());
  Assert(sxcat.pos.size() == fscat.id.size());
  dbg<<"count stars = "<<std::count(fscat.star_flag.begin(),fscat.star_flag.end(),1)<<std::endl;
  for (size_t i=0; i<fscat.id.size(); i++) {
    if (fscat.star_flag[i] != 0) {
      starcat.id.push_back( fscat.id[i] );
      starcat.sigma0.push_back( fscat.sigma0[i] );
      starcat.size_flags.push_back( fscat.size_flags[i] );
      starcat.star_flag.push_back( fscat.star_flag[i] );

      starcat.local_sky.push_back( sxcat.local_sky[i] );
      starcat.noise.push_back( sxcat.noise[i] );
      starcat.pos.push_back( sxcat.pos[i] );
    }
  }
  dbg<<"nstars = "<<starcat.id.size()<<std::endl;
}


// This is now very DES specific in order to quickly meet Joe's demands
int DoMeasurePSF_DES(ConfigFile& params, PSFLog& log) 
{
  XDEBUG = true;

  xdbg<<"Start DoMeasurePSF\n";

  // Load image:
  int image_hdu = 1;
  if (params.keyExists("image_hdu")) image_hdu = params["image_hdu"];
  Image<double> im(Name(params,"image",true),image_hdu);
  xdbg<<"Opened image "<<Name(params,"image",true)<<std::endl;

  // Load weight image
  int weight_hdu = 1;
  if (params.keyExists("weight_hdu")) weight_hdu = params["weight_hdu"];
  //Image<double> weight_im(Name(params,"weight",true),weight_hdu);
  Image<double>* weight_im = new Image<double>(Name(params,"weight",true),weight_hdu);
  dbg<<"Opened weight image.\n";


  // Read some needed parameters
  Assert(params.keyExists("psf_aperture"));
  double psfap = double(params["psf_aperture"]); 
  dbg<<"psfap = "<<psfap<<std::endl;
  Assert(params.keyExists("psf_order"));
  int psforder = params["psf_order"];
  bool output_dots=false;
  if (params.keyExists("output_dots")) output_dots=true;



  // Load SExtractor catalog info
  SXCAT_STRUCT sxcat;
  ReadSXCat(params, sxcat);

  // The results of findstars
  FINDSTARS_STRUCT fscat;
  ReadFindStarsCat(params, fscat);

  // the subset of stars
  FINDSTARS_STRUCT starcat;
  ExtractStars(sxcat, fscat, starcat);


  // These are unused
  double gain=0.0;

  // Read distortion function
  Transformation trans(params);


  int nstars = starcat.id.size();
  dbg<<"nstars = "<<nstars<<std::endl;


  // Initial sigma_p for shapelet measurements
  double sigma_p = 1.;
  if (params.keyExists("psf_seeing_est")) {
    double seeing = params["psf_seeing_est"];
    // seeing is given as FWHM
    // for a gaussian 0.5 = exp(-((FWHM/2)/sigma)^2/2)
    // FWHM/sigma = 2.*sqrt(2 ln(2)) = 2.35
    sigma_p = seeing / 2.35;
  }


  // Calculate a good value of sigma to use:
  // (For this calculation, psfap is psf_aperture * 1 arcsec.)
  EstimateSigma(sigma_p,
      im,
      starcat.pos,starcat.local_sky,starcat.noise,
      gain,weight_im,trans,psfap);
  dbg<<"sigma_p = "<<sigma_p<<std::endl;
  psfap *= sigma_p;  // arcsec
  dbg<<"psfap => "<<psfap<<std::endl;


  // Output structure
  PSF_STRUCT psfcat;
  ResizePSFCat(psfcat, starcat.id.size(), psforder, sigma_p);
  // copy in ids of the stars
  for (size_t i=0; i<starcat.id.size(); i++) {
    psfcat.id[i] = starcat.id[i];
  }



  // Setup output vectors
  //std::vector<BVec> psf(nstars,BVec(psforder,sigma_p));
  //std::vector<double> nu(nstars,0.);

  // Set up a default psf vector for output when an object measurement
  // fails
  BVec psf_default(psforder,sigma_p);
  for (size_t i=0; i<psf_default.size(); i++) {
    psf_default[i] = DEFVALNEG;
  }
  double nu_default = DEFVALNEG;

  // Array of flag values
  //std::vector<int32> flagvec(nstars,0);

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
#endif
      PSFLog log1;  // Just for this thread
      log1.NoWriteLog();
#ifdef _OPENMP
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
	long flag1=0;
	try {
	  MeasureSinglePSF(
	      // Input data:
	      starcat.pos[i], im, starcat.local_sky[i], trans, 
	      // Noise values:
	      starcat.noise[i], gain, weight_im,
	      // Parameters:
	      sigma_p, psfap, psforder,
	      // Log information
	      log1,
	      // Ouput value:
	      psf1, nu1, flag1);
	} catch (tmv::Error& e) {
	  dbg<<"TMV Error thrown in MeasureSinglePSF\n";
	  dbg<<e<<std::endl;
	  log1.nf_tmverror++;
	  flag1 |= MPSF_TMV_EXCEPTION;
	} catch (...) {
	  dbg<<"unkown exception in MeasureSinglePSF\n";
	  log1.nf_othererror++;
	  flag1 |= MPSF_UNKNOWN_EXCEPTION;
	}
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  psfcat.psf_flags[i] = flag1;
	  psfcat.psf[i] = psf1;
	  psfcat.nu[i] = nu1;
	  if (!flag1) {
	    dbg<<"Successful psf measurement: "<<psf1<<std::endl;
	  }
	  else {
	    dbg<<"Unsuccessful psf measurement\n"; 
	  }
	}
#ifdef SINGLESTAR
	break;
#endif
      }
#ifdef _OPENMP
#pragma omp critical
#endif
      {
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
#endif


  int nsuccess = log.ns_psf;

  dbg<<nsuccess<<" successful star measurements, ";
  dbg<<nstars-nsuccess<<" unsuccessful\n";

/*
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
	starcat.pos[i].GetX(),
	starcat.pos[i].GetY(),
	psfcat.psf_flags[i],
	psfcat.nu[i],
	psforder, sigma_p, psfcat.psf[i],
	psfdelim);
  }
  */
  WritePSFCat(params, psfcat);
  TestPSFCatIO(params, psfcat);
  dbg<<"Done writing output psf catalog\n";


  // MJ -- Should we make a FittedPSFLog?  Are there parameters here that
  //       we want to ingest into a DB meta-table?
  //       If so, we would need to WritePSFCat after doing the
  //       FittedPSF calculation, since the fittedpsf file is not ingested.

  // Fit the PSF with a polynomial:
  FittedPSF fittedpsf(psfcat.psf,psfcat.psf_flags,starcat.pos,psfcat.nu,sigma_p,params);
  dbg<<"Done fitting PSF\n";


  // Output fitted psf
  /*
  std::string fitpsffile = Name(params,"fitpsf");
  std::ofstream fitout(fitpsffile.c_str());
  fitout << fittedpsf;
  fitout.close();
  */


  fittedpsf.Write(params);
  dbg<<"Done writing fitted PSF file\n";



  if (XDEBUG) {
    // Check fit:

    FittedPSF readfittedpsf;
    readfittedpsf.Read(params);
    for(int i=0;i<nstars;i++) if (!psfcat.psf_flags[i]) {
      xdbg<<"psf[i] = "<<psfcat.psf[i]<<std::endl;
      BVec checkpsf(readfittedpsf.GetOrder(),readfittedpsf.GetSigma());
      checkpsf = readfittedpsf(starcat.pos[i]);
      xdbg<<"fittedpsf = "<<checkpsf<<std::endl;
      xdbg<<"Norm(diff) = "<<Norm(psfcat.psf[i]-checkpsf)<<std::endl;
    }
  }
  xdbg<<log<<std::endl;

  if (output_dots) { 
    std::cerr
      <<std::endl
      <<"Success rate: "<<nsuccess<<"/"<<nstars
      <<std::endl; 
  }

  // testing
#if 0
  PSF_STRUCT testcat;
  ReadPSFCat(params, testcat);
  for (size_t i=0; i<psfcat.id.size(); i++) {
    if ( (psfcat.psf_flags[i] != testcat.psf_flags[i]) ||
	(psfcat.psf_order[i] != testcat.psf_order[i]) ) {
      std::cout
	<<i
	<<" psf_flags: "<<psfcat.psf_flags[i]<<"  "<<testcat.psf_flags[i]
	<<" psf_order: "<<psfcat.psf_order[i]<<"  "<<testcat.psf_order[i]
	<<std::endl;
    }
  }
#endif




  return nsuccess;
}

int DoMeasurePSF(ConfigFile& params, PSFLog& log) 
{
  XDEBUG = true;

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
  std::vector<Position> all_skypos;
  double gain;
  Image<double>* weight_im = 0;
  ReadCatalog(params,"starcat",all_pos,all_sky,all_noise,gain,weight_im,
      all_skypos);
  xdbg<<"Done read catalog "<<Name(params,"starcat",false,true)<<std::endl;

  // Fix sky if necessary
  if (all_sky.size() == 0) {
    double glob_sky = 0.;
    if (params.keyExists("image_sky")) glob_sky = params["image_sky"];
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
  if (params.keyExists("psf_seeing_est")) {
    double seeing = params["psf_seeing_est"];
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
  std::vector<long> flagvec(nstars,0);

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
#endif
      PSFLog log1;  // Just for this thread
      log1.NoWriteLog();
#ifdef _OPENMP
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
	long flag1=0;
	try {
	  MeasureSinglePSF(
	      // Input data:
	      all_pos[i], im, all_sky[i], trans, 
	      // Noise values:
	      all_noise[i], gain, weight_im,
	      // Parameters:
	      sigma_p, psfap, psforder,
	      // Log information
	      log1,
	      // Ouput value:
	      psf1, nu1, flag1);
	} catch (tmv::Error& e) {
	  dbg<<"TMV Error thrown in MeasureSinglePSF\n";
	  dbg<<e<<std::endl;
	  log1.nf_tmverror++;
	  flag1 |= MPSF_TMV_EXCEPTION;
	} catch (...) {
	  dbg<<"unkown exception in MeasureSinglePSF\n";
	  log1.nf_othererror++;
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
	break;
#endif
      }
#ifdef _OPENMP
#pragma omp critical
#endif
      {
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
  fittedpsf.Write(params);
  dbg<<"Done writing fitted PSF file\n";

  if (XDEBUG) {
    // Check fit:

    FittedPSF readfittedpsf;
    readfittedpsf.Read(params);
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
    long flags, 
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
    << sigma_p;
  for(size_t i=0;i<psf.size();i++) {
    ostream << delim << psf[i];
  }
  ostream << std::endl;
}


