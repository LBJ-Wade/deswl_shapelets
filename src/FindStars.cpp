
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "TMV.h"
#include "ConfigFile.h"
#include "dbg.h"
#include "Name.h"
#include "StarCatalog.h"
#include "InputCatalog.h"
#include "StarFinder.h"

std::ostream* dbgout = 0;
bool XDEBUG = false;

static void DoFindStars(ConfigFile& params, FindStarsLog& log)
{
  // Read input catalog
  InputCatalog incat(params,"cat_");

  // Create StarCatalog from InputCatalog
  StarCatalog starcat(incat,params,"stars_");

  // Update the sizes to more robust values
  std::auto_ptr<Image<double> > weight_im;
  Image<double> im(params,weight_im);
  Transformation trans(params);
  starcat.CalcSizes(im,weight_im.get(),trans);

  try {
    starcat.FindStars(log);
  } catch(StarFinderException& e) {
    // Need to catch this here, so we can write the output file
    // with the sizes, even though we haven't figured out which 
    // objects are stars.
    dbg<<"Caught StarFinderException: "<<e.what()<<std::endl;  
    starcat.Write();
    throw;
  } catch (...) {
    dbg<<"Caught unknown exception\n";
    starcat.Write();
    throw;
  }
  dbg<<"After RunFindStars\n";

  // Write star catalog to file
  starcat.Write();

  xdbg<<"Log: \n"<<log<<std::endl;
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
    return EXIT_FAILURE;
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
  catch(StarFinderException& e ) {
    dbg<<"Caught \n"<<e.what()<<std::endl;
    log.exitcode = FAILURE_STARFINDER_ERROR;
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
    log.exitcode = FAILURE;
    dbg<<"Caught Unknown error\n";
    std::cerr<<"Caught Unknown error\n";
    log.extraexitinfo="Caught unknown exception";
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
