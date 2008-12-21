
#include "ConfigFile.h"
#include "DoMeasure.h"
#include "TMV.h"
#include "dbg.h"
#include "Name.h"

#include "fitsio.h"
#include "FitsFile.h"

#include "SXCat.h"
#include "StarFinder.h"

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

std::ostream* dbgout = 0;
bool XDEBUG = false;


static void TestFindStarsIO(ConfigFile& params, FINDSTARS_STRUCT& fscat)
{

  FINDSTARS_STRUCT test;
  ReadFindStarsCat(params,test);

  std::stringstream err;

  if (test.id.size() != fscat.id.size()) {
    throw std::runtime_error("id in catalog is wrong size");
  }
  for (size_t i=0; i<test.id.size(); i++) {
    if (test.id[i] != fscat.id[i]) {
      err<<"id["<<i<<"] differs in catalog";
      throw std::runtime_error(err.str());
    }
  }

  if (test.size_flags.size() != fscat.size_flags.size()) {
      throw std::runtime_error("field size_flags in catalog is wrong size");
  }
  for (size_t i=0; i<test.size_flags.size(); i++) {
    if (test.size_flags[i] != fscat.size_flags[i]) {
      err<<"size_flags["<<i<<"] differs in catalog";
      throw std::runtime_error(err.str());
    }
  }

  if (test.star_flag.size() != fscat.star_flag.size()) {
      throw std::runtime_error("field size_flags in catalog is wrong size");
  }
  for (size_t i=0; i<test.star_flag.size(); i++) {
    if (test.star_flag[i] != fscat.star_flag[i]) {
      err<<"star_flag["<<i<<"] differs in catalog";
      throw std::runtime_error(err.str());
    }
  }


/*
  Assert(test.id.size() == fscat.id.size());
  Assert(test.sigma0.size() == fscat.id.size());
  Assert(test.size_flags.size() == fscat.id.size());
  Assert(test.star_flag.size() == fscat.id.size());
  for(size_t i=0;i<fscat.id.size();i++) {
    if (fscat.id[i] != test.id[i]) {
      dbg<<"id["<<i<<"] = "<<fscat.id[i]<<"  "<<test.id[i]<<std::endl;
      Assert(0);
    }
    if (fscat.sigma0[i] != test.sigma0[i]) {
      dbg<<"sigma0["<<i<<"] = "<<fscat.sigma0[i]<<"  "<<test.sigma0[i]<<std::endl;
      Assert(0);
    }
    if (fscat.size_flags[i] != test.size_flags[i]) {
      dbg<<"size_flags["<<i<<"] = "<<fscat.size_flags[i]<<"  "<<test.size_flags[i]<<std::endl;
      Assert(0);
    }
    if (fscat.star_flag[i] != test.star_flag[i]) {
      dbg<<"star_flag["<<i<<"] = "<<fscat.star_flag[i]<<"  "<<test.star_flag[i]<<std::endl;
      Assert(0);
    }
  }
*/
}


static void DoFindStars(ConfigFile& params, FindStarsLog& log)
{
  // Load image:
  int image_hdu = 1;
  if (params.keyExists("image_hdu")) image_hdu = params["image_hdu"];
  Image<double> im(Name(params,"image",true),image_hdu);
  dbg<<"Loaded image\n";

  // load weight image
  int weight_hdu = 1;
  if (params.keyExists("weight_hdu")) weight_hdu = params["weight_hdu"];
  Image<double> weight_im(Name(params,"weight",true),weight_hdu);
  dbg<<"Loaded weight image\n";

  SXCAT_STRUCT sxcat;
  ReadSXCat(params, sxcat);
  FINDSTARS_STRUCT fscat;
  ResizeFindStarsCat(fscat, sxcat.id.size());
  // copy over id information
  for (size_t i=0; i<sxcat.id.size(); i++) {
    fscat.id[i] = sxcat.id[i];
  }
  dbg<<"Read in SExtractor catalog\n";

  double psfap = double(params["psf_aperture"]); 

  // Read distortion function
  Transformation trans(params);
  dbg<<"Done read transformation\n";
  xdbg<<"Check transformation:\n";
  double rmserror = 0.;
  int count = 0;
  for(size_t i=0;i<sxcat.pos.size();i++) {
    try {
      Position skypos(sxcat.ra[i],sxcat.dec[i]);
      Position temp;
      trans.Transform(sxcat.pos[i],temp);
      temp /= 3600.; // arcsec -> degrees
      xdbg<<sxcat.pos[i]<<"  "<<skypos<<"  "<<temp;
      xdbg<<"  "<<(temp-skypos)*3600.<<std::endl;
      rmserror += std::norm(temp-skypos);
      count++;
    } catch (Range_error& e) {
      xdbg<<"distortion range error\n";
      xdbg<<"p = "<<sxcat.pos[i]<<", b = "<<e.b<<std::endl;
    }
  }
  rmserror /= count;
  rmserror = std::sqrt(rmserror);
  rmserror *= 3600.;  // degrees -> arcsec.
  xdbg<<"rms error = "<<rmserror<<" arcsec\n";
  // There seems to be a random error on the RA value of +- 0.05, which
  // gives an rms error of about 0.03.
  // I would have thought that with these being read from a binary FITS
  // file that there wouldn't be rounding issues.  But there seems to have
  // been rounding somewhere to the nearest 0.1 arcsec.
  if (rmserror > 0.1) { 
    std::cout<<"STATUS3BEG Warning: Positions from WCS transformation have rms error of "<<rmserror<<" arcsec relative to ra, dec in catalog. STATUS3END"<<std::endl;
  }

  // we are using weight image, so sk, noise, gain are dummy
  try {
    MeasureSigmas(
	im, 
	sxcat.pos, sxcat.local_sky, sxcat.noise, 0.0, 
	weight_im, 
	trans, 
	psfap, 
	fscat.sigma0,
	fscat.size_flags);
  } catch(...) {
    std::cerr<<"Caught signal from MeasureSigmas"<<std::endl;  
    throw;
  }
  dbg<<"Done MeasureSigmas\n";

  // Eventually separate this out
  dbg<<"Loading Star Finder"<<std::endl;
  StarFinder sf(params);
  dbg<<"Made StarFinder object\n";
  
  std::vector<int> starflags; 

  try {
    sf.RunFindStars(
	sxcat.flags,
	fscat.size_flags,
	sxcat.x,
	sxcat.y,
	fscat.sigma0,
	sxcat.mag,
	fscat.star_flag);
  } catch(StarFinderException& e) {
    dbg<<"Caught StarFinderException: "<<e.what()<<std::endl;  
    // Need to catch this here, so we can write the output file
    // with the sizes, even though we haven't figured out which 
    // objects are stars.
    std::string output_file = Name(params, "stars");
    dbg<<"Writing to stars file: "<<output_file<<std::endl;
    WriteFindStarsCat(params, fscat);
    throw;
  } catch (...) {
    throw;
  }
  dbg<<"After RunFindStars\n";

  if (params.keyExists("stars_maxoutmag")) {
    double maxoutmag = params["stars_maxoutmag"];
    for (size_t i=0; i<sxcat.pos.size(); i++) {
      if (sxcat.mag[i] > maxoutmag) fscat.star_flag[i] = 0;
    }
  }
  size_t star_count=0;
  for (size_t i=0; i<sxcat.pos.size(); i++) {
    if (fscat.star_flag[i] == 1) {
      star_count++;
    }
  }
  log.nobj = sxcat.pos.size();
  log.nstars = star_count;
  dbg<<"Starcount = "<<star_count<<std::endl;

  std::string output_file = Name(params, "stars");
  dbg<<"Writing to stars file: "<<output_file<<std::endl;

  if (star_count < 100) {
    std::cout<<"STATUS3BEG Warning: Only "<<star_count<<" stars found for Name="<<output_file<<". STATUS3END"<<std::endl;
  }

  std::string text_output_file = "findstars.dat";
  std::ofstream output(text_output_file.c_str());
  if (!output.is_open()) {
    throw std::runtime_error("Error opening stars file");
  }

  //std::string delim = params["stars_delim"];
  std::string delim = "  ";
  for (size_t i=0; i<sxcat.pos.size(); i++) {
    output
      <<fscat.id[i]<<delim
      <<fscat.sigma0[i]<<delim
      <<sxcat.mag[i]<<delim
      <<fscat.size_flags[i]<<delim
      <<fscat.star_flag[i]<<delim
      <<std::endl;
  }
  output.close();

  WriteFindStarsCat(params, fscat);
  TestFindStarsIO(params, fscat);

  dbg<<"Done Write\n";
}


int main(int argc, char **argv) try 
{
  // Read parameters
  if (argc < 2) {
    std::cout<<"STATUS5BEG Invalid command line for findstars. STATUS5END"<<std::endl;
    std::cerr<<"Usage: findstars configfile [param=value ...]\n";
    std::cerr<<"\tThe first parameter is the configuration file that has \n";
    std::cerr<<"\tall the parameters for this run. \n";
    std::cerr<<"\tThese values may be modified on the command line by \n";
    std::cerr<<"\tentering param/value pais as param=value. \n";
    std::cerr<<"Note: root is not usuallly given in the parameter file, \n";
    std::cerr<<"\tso the normal command line would be something like:\n";
    std::cerr<<"\tmeasurepsf measurepsf.config root=img123\n";
    return EXIT_FAILURE;
  }
  ConfigFile params(argv[1]);
  for(int k=2;k<argc;k++) params.Append(argv[k]);

  std::string logfile = "";
  if (params.keyExists("log_file") || params.keyExists("log_ext")) 
    logfile = Name(params,"log");

  std::string logdelim = "  ";
  if (params.keyExists("log_delim")) logdelim = params["log_delim"];

  std::string stars_file=Name(params,"stars");
  FindStarsLog log(logfile,stars_file); 
  // This automatically writes its output when it goes out of scope
  // whether that is naturally in after an exception is caught.
  // Log output is:  (all on one line)
  //    exitcode  nobjects  nstars_found

  try {
    // Setup debugging
    if (params.keyExists("verbose") && int(params["verbose"]) > 0) {
      if (params.keyExists("debug_file") || params.keyExists("debug_ext")) {
	dbgout = new std::ofstream(Name(params,"debug").c_str());
      }
      else dbgout = &std::cout;
      if (int(params["verbose"]) > 1) XDEBUG = true;
    }

#ifdef _OPENMP
    if (params.keyExists("omp_num_threads")) {
      int num_threads = params["omp_num_threads"];
      omp_set_num_threads(num_threads);
    }
#endif

    dbg<<"Config params = \n"<<params<<std::endl;

    DoFindStars(params,log);

    if (dbgout && dbgout != &std::cout) {delete dbgout; dbgout=0;}
    return EXIT_SUCCESS;
  }
#if 0
  // Change to 1 to let gdb see where the program bombed out.
  catch(int) {}
#else
  catch (file_not_found& e)
  {
    dbg<<"Caught \n"<<e.what()<<std::endl;
    std::cerr<<"Caught \n"<<e.what()<<std::endl;
    log.exitcode = FAILURE_FILE_NOT_FOUND;
    log.extraexitinfo = e.what();
    //if (Status(log.exitcode)==5) 
    return EXIT_FAILURE;
    //else return EXIT_SUCCESS;
  }
  catch (ConfigFile::file_not_found& e)
  {
    dbg<<"Caught \n"<<e.what()<<std::endl;
    std::cerr<<"Caught \n"<<e.what()<<std::endl;
    log.exitcode = FAILURE_CONFIGFILE_ERROR;
    log.extraexitinfo = e.what();
    //if (Status(log.exitcode)==5) 
    return EXIT_FAILURE;
    //else return EXIT_SUCCESS;
  }
  catch (ConfigFile::key_not_found& e)
  {
    dbg<<"Caught \n"<<e.what()<<std::endl;
    std::cerr<<"Caught \n"<<e.what()<<std::endl;
    log.exitcode = FAILURE_CONFIGFILE_ERROR;
    log.extraexitinfo = e.what();
    //if (Status(log.exitcode)==5) 
    return EXIT_FAILURE;
    //else return EXIT_SUCCESS;
  }
  catch (tmv::Error& e)
  {
    dbg<<"Caught \n"<<e<<std::endl;
    std::cerr<<"Caught \n"<<e<<std::endl;
    log.exitcode = FAILURE_TMV_ERROR;
    log.extraexitinfo = e.what();
    //if (Status(log.exitcode)==5) 
    return EXIT_FAILURE;
    //else return EXIT_SUCCESS;
  }
  catch(StarFinderException& e ) {
    dbg<<"Caught \n"<<e.what()<<std::endl;
    log.exitcode = FAILURE_STARFINDER_ERROR;
    log.extraexitinfo = e.what();
    //if (Status(log.exitcode)==5) 
    return EXIT_FAILURE;
    //else return EXIT_SUCCESS;
  }
  catch (std::exception& e)
  {
    dbg<<"Caught \n"<<e.what()<<std::endl;
    std::cerr<<"Caught \n"<<e.what()<<std::endl;
    log.exitcode = FAILURE_STD_EXCEPTION;
    log.extraexitinfo = e.what();
    //if (Status(log.exitcode)==5) 
    return EXIT_FAILURE;
    //else return EXIT_SUCCESS;
  }
  catch (ExitCode e)
  {
    dbg<<"Caught ExitCode "<<e<<std::endl;
    std::cerr<<"Caught ExitCode "<<e<<std::endl;
    log.exitcode = e;
    //if (Status(log.exitcode)==5) 
    return EXIT_FAILURE;
    //else return EXIT_SUCCESS;
  }
  catch (...)
  {
    log.exitcode = FAILURE;
    dbg<<"Caught Unknown error\n";
    std::cerr<<"Caught Unknown error\n";
    log.extraexitinfo="Caught unknown exception";
    //if (Status(log.exitcode)==5) 
    return EXIT_FAILURE;
    //else return EXIT_SUCCESS;
  }
#endif
}
catch (std::exception& e)
{
  dbg<<"Caught \n"<<e.what()<<std::endl;
  dbg<<"outside of the normal try..catch block.";
  dbg<<"Unable to write to the log file.\n";
  std::cerr<<"Caught \n"<<e.what()<<std::endl;
  std::cerr<<"outside of the normal try..catch block.";
  std::cout<<"STATUS5BEG Catastrophic error: "<<e.what()<<" STATUS5END\n";
  return EXIT_FAILURE;
}
catch (...)
{
  dbg<<"Cought an exception outside of the normal try..catch block.";
  dbg<<"Unable to write to the log file.\n";
  std::cerr<<"Cought an exception outside of the normal try..catch block.";
  std::cerr<<"Unable to write to the log file.\n";
  std::cout<<"STATUS5BEG Catastrophic error: Caught unknown exception STATUS5END\n";
  return EXIT_FAILURE;
}
