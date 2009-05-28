
#include "Image.h"
#include "FittedPSF.h"
#include "TimeVars.h"
#include "Log.h"
#include "CoaddCatalog.h"
#include "InputCatalog.h"
#include "BasicSetup.h"
#include <sys/time.h>

#include <iostream>
#include <fstream>

static void DoMeasureMultiShear(ConfigFile& params, ShearLog& log) 
{
  bool timing = params.read("timing",false);
  timeval tp;
  double t1=0.,t2=0.;

  CoaddCatalog coaddcat(params);
  dbg<<"Made coaddcat\n";

  coaddcat.ReadPixelLists();
  dbg<<"After ReadPixelLists\n";

  // coaddcat.MeasureShears();
  //dbg<<"After MeasureShears\n";


  return;
  /*
  if (timing) {
    gettimeofday(&tp,0);
    t1 = tp.tv_sec + tp.tv_usec/1.e6;
  }
  */

  /*
  // Load image:
  std::auto_ptr<Image<double> > weight_im;
  Image<double> im(params,weight_im);

  if (timing) {
    gettimeofday(&tp,0);
    t2 = tp.tv_sec + tp.tv_usec/1.e6;
    std::cout<<"Time: Open imgae = "<<t2-t1<<std::endl;
    t1 = t2;
  }

  // Read distortion function
  Transformation trans(params);

  if (timing) {
    gettimeofday(&tp,0);
    t2 = tp.tv_sec + tp.tv_usec/1.e6;
    std::cout<<"Time: Read Transformation = "<<t2-t1<<std::endl;
    t1 = t2;
  }

  // Read input catalog
  InputCatalog incat(params,&im);

  if (timing) {
    gettimeofday(&tp,0);
    t2 = tp.tv_sec + tp.tv_usec/1.e6;
    std::cout<<"Time: Read InputCatalog = "<<t2-t1<<std::endl;
    t1 = t2;
  }

  // Read the fitted psf file
  FittedPSF fitpsf(params);

  if (timing) {
    gettimeofday(&tp,0);
    t2 = tp.tv_sec + tp.tv_usec/1.e6;
    std::cout<<"Time: Read FittedPSF = "<<t2-t1<<std::endl;
    t1 = t2;
  }

  // Create shear catalog
  ShearCatalog shearcat(incat,trans,params);

  if (timing) {
    gettimeofday(&tp,0);
    t2 = tp.tv_sec + tp.tv_usec/1.e6;
    std::cout<<"Time: Create ShearCatalog = "<<t2-t1<<std::endl;
    t1 = t2;
  }

  // Measure shears and shapelet vectors
  int nshear = shearcat.MeasureShears(im,weight_im.get(),trans,fitpsf,log);

  if (timing) {
    gettimeofday(&tp,0);
    t2 = tp.tv_sec + tp.tv_usec/1.e6;
    std::cout<<"Time: Measure Shears = "<<t2-t1<<std::endl;
    t1 = t2;
  }

  // Write results to file
  shearcat.Write();

  if (timing) {
    gettimeofday(&tp,0);
    t2 = tp.tv_sec + tp.tv_usec/1.e6;
    std::cout<<"Time: Write ShearCatalog = "<<t2-t1<<std::endl;
    t1 = t2;
  }

  if (nshear == 0) throw std::runtime_error("No successful shear measurements");
  */

  xdbg<<"Log: \n"<<log<<std::endl;
}

int main(int argc, char **argv) try 
{
  ConfigFile params;
  if (BasicSetup(argc,argv,params,"multishear")) return EXIT_FAILURE;

  // Setup Log
  std::string logfile = ""; // Default is to stdout
  if (params.keyExists("log_file") || params.keyExists("log_ext")) 
    logfile = Name(params,"log");
  std::string multishear_file = Name(params,"multishear");
  ShearLog log(logfile,multishear_file,des_qa); 

  try {
    DoMeasureMultiShear(params,log);
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
    dbg<<"Caught Unknown error\n";
    std::cerr<<"Caught Unknown error\n";
    log.exitcode = FAILURE;
    log.extraexitinfo = "Caught unknown exception";
    //if (Status(log.exitcode)==5) 
    return EXIT_FAILURE;
    //else return EXIT_SUCCESS;
  }
#endif

  if (dbgout && dbgout != &std::cout) {delete dbgout; dbgout=0;}
  return EXIT_SUCCESS; // = 0 typically.  Defined in <cstdlib>
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
  std::cerr<<"Fatal error: Caught an exception.\n";
  if (des_qa) 
    std::cout<<"STATUS5BEG Fatal error: unknown exception STATUS5END\n";
  return EXIT_FAILURE;
}
