
#include <valarray>
#include "Image.h"
#include "FittedPsf.h"
#include "TimeVars.h"
#include "Log.h"
#include "ShearCatalog.h"
#include "InputCatalog.h"
#include "BasicSetup.h"
#include <sys/time.h>

static void doMeasureShear(ConfigFile& params, ShearLog& log) 
{
  bool timing = params.read("timing",false);
  timeval tp;
  double t1=0.,t2=0.;

  if (timing) {
    gettimeofday(&tp,0);
    t1 = tp.tv_sec + tp.tv_usec/1.e6;
  }

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
  FittedPsf fitpsf(params);

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
  int nshear = shearcat.measureShears(im,weight_im.get(),trans,fitpsf,log);

  if (timing) {
    gettimeofday(&tp,0);
    t2 = tp.tv_sec + tp.tv_usec/1.e6;
    std::cout<<"Time: Measure Shears = "<<t2-t1<<std::endl;
    t1 = t2;
  }

  // Write results to file
  shearcat.write();

  if (timing) {
    gettimeofday(&tp,0);
    t2 = tp.tv_sec + tp.tv_usec/1.e6;
    std::cout<<"Time: Write ShearCatalog = "<<t2-t1<<std::endl;
    t1 = t2;
  }

  if (nshear == 0) {
      throw ProcessingException(
          "No successful shear measurements");
  }

  xdbg<<"Log: \n"<<log<<std::endl;
}

int main(int argc, char **argv) try 
{
  ConfigFile params;
  if (BasicSetup(argc,argv,params,"measureshear")) return EXIT_FAILURE;

  // Setup Log
  std::string logfile = ""; // Default is to stdout
  if (params.keyExists("log_file") || params.keyExists("log_ext")) 
    logfile = makeName(params,"log",false,false);
  std::string shear_file = makeName(params,"shear",false,false);
  std::auto_ptr<ShearLog> log (
      new ShearLog(params,logfile,shear_file)); 

  try {
    doMeasureShear(params,*log);
  }
#if 0
  // Change to 1 to let gdb see where the program bombed out.
  catch(int) {}
#else
  CATCHALL
#endif

  if (dbgout && dbgout != &std::cout) {delete dbgout; dbgout=0;}
  return EXIT_SUCCESS; // = 0 typically.  Defined in <cstdlib>
}
catch (std::exception& e)
{
  std::cerr<<"Fatal error: Caught \n"<<e.what()<<std::endl;
  std::cout<<"STATUS5BEG Fatal error: "<<e.what()<<" STATUS5END\n";
  return EXIT_FAILURE;
}
catch (...)
{
  std::cerr<<"Fatal error: Cought an exception.\n";
  std::cout<<"STATUS5BEG Fatal error: unknown exception STATUS5END\n";
  return EXIT_FAILURE;
}
