
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


void DoFindStars(ConfigFile& params, FindStarsLog& log)
{
  // Load image:
  int image_hdu = 1;
  if (params.keyExists("image_hdu")) image_hdu = params["image_hdu"];
  Image<double> im(Name(params,"image",true),image_hdu);

  // load weight image
  int weight_hdu = 1;
  if (params.keyExists("weight_hdu")) weight_hdu = params["weight_hdu"];
  Image<double> weight_im(Name(params,"weight",true),weight_hdu);

  SXCAT_STRUCT sxcat;
  ReadSXCat(params, sxcat);

  double psfap = double(params["psf_aperture"]); 

  // Read distortion function
  Transformation trans(params);

  // we are using weight image, so sk, noise, gain are dummy
  MeasureSigmas(
      im, 
      sxcat.pos, sxcat.local_sky, sxcat.noise, 0.0, 
      weight_im, 
      trans, 
      psfap, 
      sxcat.sigma,
      sxcat.size_flags);

  // Eventually separate this out
  std::cerr<<"Loading Star Finder"<<std::endl;
  StarFinder sf(params);

  
  std::vector<int> starflags; 
  sf.RunFindStars(
      sxcat.flags,
      sxcat.size_flags,
      sxcat.x,
      sxcat.y,
      sxcat.sigma,
      sxcat.mag,
      sxcat.star_flag);


  std::string output_file = Name(params, "allcat");

  std::cerr<<"Writing to allcat file: "<<output_file<<std::endl;
  std::ofstream output(output_file.c_str());

  std::string delim = params["allcat_delim"];
  for (size_t i=0; i<sxcat.pos.size(); i++) {
    output
      <<sxcat.id[i]<<delim
      <<sxcat.x[i]<<delim
      <<sxcat.y[i]<<delim
      <<sxcat.local_sky[i]<<delim
      <<sxcat.sigma[i]<<delim
      <<sxcat.mag[i]<<delim
      <<sxcat.size_flags[i]<<delim
      <<sxcat.star_flag[i]
      <<std::endl;
  }

  output_file = Name(params, "starcat");

  std::cerr<<"Writing to starcat file: "<<output_file<<std::endl;
  std::ofstream star_output(output_file.c_str());

  delim = params["starcat_delim"];
  int star_count=0;
  for (size_t i=0; i<sxcat.pos.size(); i++) {
    if (sxcat.star_flag[i] == 1) {
      star_count++;
      star_output
	<<sxcat.id[i]<<delim
	<<sxcat.x[i]<<delim
	<<sxcat.y[i]<<delim
	<<sxcat.local_sky[i]<<delim
	<<sxcat.sigma[i]<<delim
	<<sxcat.mag[i]<<delim
	<<sxcat.size_flags[i]<<delim
	<<sxcat.star_flag[i]
	<<std::endl;
    }
  }

  log.nobj = sxcat.pos.size();
  log.nstars = star_count;
}


int main(int argc, char **argv) try 
{
  // Read parameters
  if (argc < 2) {
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
  FindStarsLog log(logfile,logdelim); 
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
    if (params.keyExists("num_threads")) {
      int num_threads = params["num_threads"];
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
    return EXIT_FAILURE; // = 1 typically
  }
  catch (ConfigFile::file_not_found& e)
  {
    dbg<<"Caught \n"<<e.what()<<std::endl;
    std::cerr<<"Caught \n"<<e.what()<<std::endl;
    log.exitcode = FAILURE_CONFIGFILE_ERROR;
    return EXIT_FAILURE;
  }
  catch (ConfigFile::key_not_found& e)
  {
    dbg<<"Caught \n"<<e.what()<<std::endl;
    std::cerr<<"Caught \n"<<e.what()<<std::endl;
    log.exitcode = FAILURE_CONFIGFILE_ERROR;
    return EXIT_FAILURE;
  }
  catch (tmv::Error& e)
  {
    dbg<<"Caught \n"<<e<<std::endl;
    std::cerr<<"Caught \n"<<e<<std::endl;
    log.exitcode = FAILURE_TMV_ERROR;
    return EXIT_FAILURE;
  }
  catch (std::exception& e)
  {
    dbg<<"Caught \n"<<e.what()<<std::endl;
    std::cerr<<"Caught \n"<<e.what()<<std::endl;
    log.exitcode = FAILURE_STD_EXCEPTION;
    return EXIT_FAILURE;
  }
  catch (ExitCode e)
  {
    dbg<<"Caught ExitCode "<<e<<std::endl;
    std::cerr<<"Caught ExitCode "<<e<<std::endl;
    log.exitcode = e;
    return EXIT_FAILURE;
  }
  catch (...)
  {
    log.exitcode = FAILURE;
    dbg<<"Caught Unknown error\n";
    std::cerr<<"Caught Unknown error\n";
    return EXIT_FAILURE;
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
  std::cerr<<"Unable to write to the log file.\n";
  return EXIT_FAILURE;
}
catch (...)
{
  dbg<<"Cought an exception outside of the normal try..catch block.";
  dbg<<"Unable to write to the log file.\n";
  std::cerr<<"Cought an exception outside of the normal try..catch block.";
  std::cerr<<"Unable to write to the log file.\n";
  return EXIT_FAILURE;
}