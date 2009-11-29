
#include <valarray>
#include "StarCatalog.h"
#include "InputCatalog.h"
#include "StarFinder.h"
#include "BasicSetup.h"
#include <sys/time.h>

static void doFindStars(ConfigFile& params, FindStarsLog& log)
{
  bool timing = params.read("timing",false);
  timeval tp;
  double t1=0.,t2=0.;

  if (timing) {
    gettimeofday(&tp,0);
    t1 = tp.tv_sec + tp.tv_usec/1.e6;
  }

  // Read image, transformation
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

  // Create StarCatalog from InputCatalog
  StarCatalog starcat(incat,params);

  if (timing) {
    gettimeofday(&tp,0);
    t2 = tp.tv_sec + tp.tv_usec/1.e6;
    std::cout<<"Time: Make StarCatalog = "<<t2-t1<<std::endl;
    t1 = t2;
  }

  // Update the sizes to more robust values
  starcat.calculateSizes(im,weight_im.get(),trans);

  if (timing) {
    gettimeofday(&tp,0);
    t2 = tp.tv_sec + tp.tv_usec/1.e6;
    std::cout<<"Time: CalcSizes = "<<t2-t1<<std::endl;
    t1 = t2;
  }

  try {
    starcat.findStars(log);
  } catch(StarFinderException& e) {
    // Need to catch this here, so we can write the output file
    // with the sizes, even though we haven't figured out which 
    // objects are stars.
    dbg<<"Caught StarFinderException: "<<e.what()<<std::endl;  
    starcat.write();
    throw;
  } catch (...) {
    dbg<<"Caught unknown exception\n";
    starcat.write();
    throw;
  }
  dbg<<"After RunFindStars\n";

  if (timing) {
    gettimeofday(&tp,0);
    t2 = tp.tv_sec + tp.tv_usec/1.e6;
    std::cout<<"Time: FindStars = "<<t2-t1<<std::endl;
    t1 = t2;
  }

  // Write star catalog to file
  starcat.write();

  if (timing) {
    gettimeofday(&tp,0);
    t2 = tp.tv_sec + tp.tv_usec/1.e6;
    std::cout<<"Time: Write StarCatalog = "<<t2-t1<<std::endl;
    t1 = t2;
  }

  xdbg<<"Log: \n"<<log<<std::endl;
}

int main(int argc, char **argv) try 
{
  ConfigFile params;
  if (BasicSetup(argc,argv,params,"findstars")) return EXIT_FAILURE;

  // Setup Log
  std::string logfile = ""; // Default is to stdout
  if (params.keyExists("log_file") || params.keyExists("log_ext")) 
    logfile = makeName(params,"log",false,false);
  std::string stars_file=makeName(params,"stars",false,false);
  std::auto_ptr<FindStarsLog> log (
      new FindStarsLog(params,logfile,stars_file)); 

  try {
    doFindStars(params,*log);
  }
#if 0
  // Change to 1 to let gdb see where the program bombed out.
  catch(int) {}
#else
  CATCHALL
#endif

  if (dbgout && dbgout != &std::cout) {delete dbgout; dbgout=0;}
  return EXIT_SUCCESS;
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
