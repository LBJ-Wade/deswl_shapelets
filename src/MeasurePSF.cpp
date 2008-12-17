
#include "ConfigFile.h"
#include "DoMeasure.h"
#include "TMV.h"
#include "dbg.h"
#include "Name.h"

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

std::ostream* dbgout = 0;
bool XDEBUG = false;

int main(int argc, char **argv) try 
{
  // Read parameters
  if (argc < 2) {
    std::cout<<"STATUS5BEG Invalid command line for measurepsf. STATUS5END"<<std::endl;
    std::cerr<<"Usage: measurepsf configfile [param=value ...]\n";
    std::cerr<<"The first parameter is the configuration file that has \n";
    std::cerr<<"all the parameters for this run. \n";
    std::cerr<<"These values may be modified on the command line by \n";
    std::cerr<<"entering param/value pais as param=value. \n";
    std::cerr<<"Note: root is not usuallly given in the parameter file, \n";
    std::cerr<<"so the normal command line would be something like:\n";
    std::cerr<<"measurepsf measurepsf.config root=img123\n";
    return EXIT_FAILURE;
  }
  ConfigFile params(argv[1]);
  for(int k=2;k<argc;k++) params.Append(argv[k]);

  std::string logfile = ""; // Default is to stdout
  if (params.keyExists("log_file") || params.keyExists("log_ext"))
    logfile = Name(params,"log");

  std::string logdelim = "  ";
  if (params.keyExists("log_delim")) logdelim = params["log_delim"];

  std::string psf_file=Name(params,"psf");
  PSFLog log(logfile,psf_file); 
  // This automatically writes its output when it goes out of scope
  // whether that is naturally in after an exception is caught.
  // Log output is:  (all on one line)
  //    exitcode  nstars  nsuccess 
  //    nfail_range_distortion
  //    nfail_edge  nfail_npix<10
  //    nfail_tmv_error  nfail_other_error
  //    nfail_psf_measurement


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

    DoMeasurePSF_DES(params,log);

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
    return EXIT_FAILURE; // = 1 typically
  }
  catch (ConfigFile::file_not_found& e)
  {
    dbg<<"Caught \n"<<e.what()<<std::endl;
    std::cerr<<"Caught \n"<<e.what()<<std::endl;
    log.exitcode = FAILURE_CONFIGFILE_ERROR;
    log.extraexitinfo = e.what();
    return EXIT_FAILURE;
  }
  catch (ConfigFile::key_not_found& e)
  {
    dbg<<"Caught \n"<<e.what()<<std::endl;
    std::cerr<<"Caught \n"<<e.what()<<std::endl;
    log.exitcode = FAILURE_CONFIGFILE_ERROR;
    log.extraexitinfo = e.what();
    return EXIT_FAILURE;
  }
  catch (tmv::Error& e)
  {
    dbg<<"Caught \n"<<e<<std::endl;
    std::cerr<<"Caught \n"<<e<<std::endl;
    log.exitcode = FAILURE_TMV_ERROR;
    log.extraexitinfo = e.what();
    return EXIT_FAILURE;
  }
  catch (std::exception& e)
  {
    dbg<<"Caught \n"<<e.what()<<std::endl;
    std::cerr<<"Caught \n"<<e.what()<<std::endl;
    log.exitcode = FAILURE_STD_EXCEPTION;
    log.extraexitinfo = e.what();
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
    dbg<<"Caught Unknown error\n";
    std::cerr<<"Caught Unknown error\n";
    log.exitcode = FAILURE;
    log.extraexitinfo = "Caught unknown exception";
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
