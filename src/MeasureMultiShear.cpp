
#include "Image.h"
#include "FittedPSF.h"
#include "TimeVars.h"
#include "Log.h"
#include "CoaddCatalog.h"
#include "MultiShearCatalog.h"
#include "BasicSetup.h"

#include <iostream>
#include <fstream>

static void DoMeasureMultiShear(ConfigFile& params, ShearLog& log) 
{
  CoaddCatalog coaddcat(params);
  dbg<<"Made coaddcat\n";

  MultiShearCatalog shearcat(coaddcat,params);
  dbg<<"Made multishearcat\n";

  int nshear = shearcat.MeasureMultiShears(log);
  dbg<<"After MeasureShears\n";

  shearcat.Write();
  dbg<<"After Write\n";

  if (nshear == 0) throw ProcessingError("No successful shear measurements");

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
  std::auto_ptr<ShearLog> log (
      new ShearLog(logfile,multishear_file,des_qa)); 

  try {
    DoMeasureMultiShear(params,*log);
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
