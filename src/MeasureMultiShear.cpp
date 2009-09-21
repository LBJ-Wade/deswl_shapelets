
#include <valarray>
#include "Image.h"
#include "FittedPSF.h"
#include "TimeVars.h"
#include "Log.h"
#include "CoaddCatalog.h"
#include "MultiShearCatalog.h"
#include "BasicSetup.h"

#include <iostream>
#include <fstream>
#include <valgrind/memcheck.h>

#ifdef __INTEL_COMPILER
// For VALGRIND.  
#pragma warning (disable : 593)
#pragma warning (disable : 1469)
#endif

static void DoMeasureMultiShear(ConfigFile& params, ShearLog& log) 
{
  bool output_dots = params.read("output_dots",false);

  CoaddCatalog coaddcat(params);
  dbg<<"Made coaddcat\n";

  MultiShearCatalog shearcat(coaddcat,params);
  dbg<<"Made multishearcat\n";
  dbg<<"Total bounds are "<<shearcat.skybounds<<std::endl;
  dbg<<"Area = "<<shearcat.skybounds.Area()/3600.<<" square arcmin\n";
  dbg<<" = "<<shearcat.skybounds.Area()/3600./3600.<<" square degrees\n";

  log.ngals = shearcat.size();
  log.ngoodin = std::count(shearcat.flags.begin(),shearcat.flags.end(),0);
  dbg<<log.ngoodin<<"/"<<log.ngals<<" galaxies with no input flags\n";

  long nshear = 0;
  double section_size = params.get("multishear_section_size");
  section_size *= 60.; // arcmin -> arcsec
  std::vector<Bounds> section_bounds = shearcat.SplitBounds(section_size);
  for(size_t i=0;i<section_bounds.size();++i) 
  {
    dbg<<"Starting section "<<(i+1)<<"/"<<section_bounds.size()<<std::endl;
    dbg<<section_bounds[i]<<std::endl;
    if (output_dots)
    {
      std::cerr<<"Starting section ";
      std::cerr<<(i+1)<<"/"<<section_bounds.size()<<std::endl;
    }
    int npix = shearcat.GetPixels(section_bounds[i]);
    if (output_dots)
      std::cerr<<npix<<" galaxies in this section.\n";

    long nshear1 = shearcat.MeasureMultiShears(section_bounds[i],log);

    nshear += nshear1;
    dbg<<"After MeasureShears: nshear = "<<nshear1<<"  "<<nshear<<std::endl;
    if (output_dots)
    {
      std::cerr<<std::endl;
      std::cerr<<nshear1<<" successful shear measurements ";
      std::cerr<<"(total "<<nshear<<")\n";
    }
    dbg<<"Valgrind Leak Check:\n";
    VALGRIND_DO_LEAK_CHECK;
    //if (VALGRIND_COUNT_ERRORS) break;
  }

  dbg<<"Done: "<<log.ns_gamma<<" successful shear measurements, ";
  dbg<<(log.ngoodin-log.ns_gamma)<<" unsuccessful, ";
  dbg<<(log.ngals-log.ngoodin)<<" with input flags\n";
  log.ngood = std::count(shearcat.flags.begin(),shearcat.flags.end(),0);
  dbg<<log.ngood<<" successful measurements with no measurement flags.\n";

  if (output_dots) 
  {
    std::cerr
      <<std::endl
      <<"Success rate: "<<log.ns_gamma<<"/"<<log.ngoodin
      <<"  # with no flags: "<<log.ngood
      <<std::endl;
#ifdef USE_POOL
    std::cerr<<"Memory used in PixelList pool = "<<
      pool_allocator<Pixel,PIXELLIST_BLOCK>::total_memory_used()/(1024.*1024.)<<
      " MB\n";
    std::cerr<<"Memory used in vector<PixelList> pool = "<<
      pool_allocator<PixelList,MULTI_BLOCK>::total_memory_used()/(1024.*1024.)<<
      " MB\n";
    std::cerr<<"Memory used in vector<vector<PixelList> > pool = "<<
      pool_allocator<std::vector<PixelList>,MULTI_BLOCK>::total_memory_used()/(1024.*1024.)<<
      " MB\n";
    std::cerr<<"Memory used in vector<int> pool = "<<
      pool_allocator<int,MULTI_BLOCK>::total_memory_used()/(1024.*1024.)<<
      " MB\n";
    std::cerr<<"Memory used in vector<vector<int> > pool = "<<
      pool_allocator<std::vector<int>,MULTI_BLOCK>::total_memory_used()/(1024.*1024.)<<
      " MB\n";
    std::cerr<<"Memory used in vector<Position> pool = "<<
      pool_allocator<Position,MULTI_BLOCK>::total_memory_used()/(1024.*1024.)<<
      " MB\n";
    std::cerr<<"Memory used in vector<vector<Position> > pool = "<<
      pool_allocator<std::vector<Position>,MULTI_BLOCK>::total_memory_used()/(1024.*1024.)<<
      " MB\n";
#endif
  }

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
      new ShearLog(params,logfile,multishear_file)); 

  try 
  {
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
  std::cout<<"STATUS5BEG Fatal error: "<<e.what()<<" STATUS5END\n";
  return EXIT_FAILURE;
}
catch (...)
{
  std::cerr<<"Fatal error: Caught an exception.\n";
  std::cout<<"STATUS5BEG Fatal error: unknown exception STATUS5END\n";
  return EXIT_FAILURE;
}
