
#include <valarray>
#include "Image.h"
#include "Transformation.h"
#include "InputCatalog.h"
#include "StarCatalog.h"
#include "FittedPSF.h"
#include "PSFCatalog.h"
#include "StarFinder.h"
#include "ShearCatalog.h"
#include "TimeVars.h"
#include "Log.h"
#include "BasicSetup.h"

static void DoFullPipeline(ConfigFile& params,
    std::string logfile, std::auto_ptr<Log>& log)
{
  log.reset(
      new FindStarsLog(params,logfile,makeName(params,"stars",false,false)));

  // Load image:

  std::cerr<<"Running FindStars\n";
  std::auto_ptr<Image<double> > weight_im;
  Image<double> im(params,weight_im);

  // Read distortion function
  Transformation trans(params);

  // Read input catalog
  InputCatalog incat(params,&im);

  // Create StarCatalog from InputCatalog
  StarCatalog starcat(incat,params);

  // Update the sizes to more robust values
  starcat.calculateSizes(im,weight_im.get(),trans);

  try {
    starcat.findStars(static_cast<FindStarsLog&>(*log));
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

  // Write star catalog to file
  starcat.write();

  xdbg<<"FindStarsLog: \n"<<*log<<std::endl;

  log.reset(new PSFLog(params,logfile,makeName(params,"psf",false,false)));

  std::cerr<<"Determining PSF\n";
  // Create PSFCatalog from StarCatalog
  PSFCatalog psfcat(starcat,params);

  // Estimate the scale size to use for shapelet decompositions
  double sigma_p = psfcat.EstimateSigma(im,weight_im.get(),trans);

  // Do the actual PSF measurements
  psfcat.MeasurePSF(im,weight_im.get(),trans,sigma_p,
      static_cast<PSFLog&>(*log));

  // Write PSF catalog to file
  psfcat.Write();

  // MJ -- Should we make a FittedPSFLog?  Are there parameters here that
  //       we want to ingest into a DB meta-table?
  //       If so, we would need to WritePSFCat after doing the
  //       FittedPSF calculation, since the fittedpsf file is not ingested.

  // Fit the PSF with a polynomial:
  FittedPSF fitpsf(psfcat,params);

  // Write fitted psf to file
  fitpsf.Write();

  xdbg<<"PSFLog: \n"<<*log<<std::endl;

  log.reset(new ShearLog(params,logfile,makeName(params,"shear",false,false)));

  std::cerr<<"Measuring Shears\n";
  // Create shear catalog
  ShearCatalog shearcat(incat,trans,params);

  // Measure shears and shapelet vectors
  shearcat.MeasureShears(im,weight_im.get(),trans,fitpsf,
      static_cast<ShearLog&>(*log));

  // Write results to file
  shearcat.Write();

  xdbg<<"ShearLog: \n"<<*log<<std::endl;

}

int main(int argc, char **argv) try 
{
  ConfigFile params;
  if (BasicSetup(argc,argv,params,"fullpipe")) return EXIT_FAILURE;

  // Setup Log
  std::string logfile = ""; // Default is to stdout
  if (params.keyExists("log_file") || params.keyExists("log_ext")) 
    logfile = makeName(params,"log",false,false);
  std::auto_ptr<Log> log;

  try {
    DoFullPipeline(params,logfile,log);
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

