
#include "StarCatalog.h"
#include "InputCatalog.h"
#include "StarFinder.h"
#include "BasicSetup.h"

static void DoFindStars(ConfigFile& params, FindStarsLog& log)
{
  // Read image, transformation
  std::auto_ptr<Image<double> > weight_im;
  Image<double> im(params,weight_im);
  Transformation trans(params);

  // Read input catalog
  InputCatalog incat(params,&im);

  // Create StarCatalog from InputCatalog
  StarCatalog starcat(incat,params);

  // Update the sizes to more robust values
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
  ConfigFile params;
  if (BasicSetup(argc,argv,params,"findstars")) return EXIT_FAILURE;

  // Setup Log
  std::string logfile = ""; // Default is to stdout
  if (params.keyExists("log_file") || params.keyExists("log_ext")) 
    logfile = Name(params,"log");
  std::string stars_file=Name(params,"stars");
  FindStarsLog log(logfile,stars_file,des_qa); 

  try {
    DoFindStars(params,log);
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

  if (dbgout && dbgout != &std::cout) {delete dbgout; dbgout=0;}
  return EXIT_SUCCESS;
}
catch (std::exception& e)
{
  std::cerr<<"Fatal error: Caught \n"<<e.what()<<std::endl;
  if (des_qa)
    std::cout<<"STATUS5BEG Fatal error: "<<e.what()<<" STATUS5END\n";
  return EXIT_FAILURE;
}
catch (...)
{
  std::cerr<<"Fatal error: Cought an exception.\n";
  if (des_qa)
    std::cout<<"STATUS5BEG Fatal error: unknown exception STATUS5END\n";
  return EXIT_FAILURE;
}
