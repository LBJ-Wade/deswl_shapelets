
#include "ConfigFile.h"
#include <iostream>
#include "DoMeasure.h"
#include "TMV.h"
#include "dbg.h"

#ifdef _OPENMP
#include <omp.h>
#endif

std::ostream* dbgout = 0;
bool XDEBUG = false;

std::string RootName(std::string fullname)
{
  size_t pos = fullname.rfind('.');
  Assert(pos != fullname.npos);
  std::string root;
  root = std::string(fullname,0,pos+1); // +1 to include the '.'
  return root;
}

int main(int argc, char **argv) 
#ifndef NOTHROW
  try 
#endif
{
  // Read parameters
  if (argc < 2) {
    std::cerr<<"Usage: measureshear fitsfile [param=value ...]\n";
    std::cerr<<"The first parameter is the input fits file name\n";
    std::cerr<<"The file measureshear.config is the configuration file.\n";
    std::cerr<<"Modifications to the parameters in the config file may be\n";
    std::cerr<<"made on the command line as param/value pairs: param=value. \n";
    return 1;
  }
  std::string infits = argv[1];
  ConfigFile params("measureshear.config");
  for(int k=2;k<argc;k++) params.Append(argv[k]);

  params["root"] = RootName(infits);

  // Setup debugging
  Assert(params.keyExists("verbose"));
  if (int(params["verbose"]) > 0) {
    if (params.keyExists("debug_file") && params["debug_file"] != "") 
      dbgout = new std::ofstream(params["debug_file"].c_str());
    else if (params.keyExists("debug_ext") && params["debug_ext"] != "") 
      dbgout = new std::ofstream((params["root"]+params["debug_ext"]).c_str());
    else dbgout = &std::cout;
    if (int(params["verbose"]) > 1) XDEBUG = true;
  }

  dbg<<"infits = "<<infits<<std::endl;
  dbg<<"params[root] = "<<params["root"]<<std::endl;
  dbg<<"params[fits_ext] = "<<params["fits_ext"]<<std::endl;
  dbg<<"sum = "<<params["root"] + params["fits_ext"]<<std::endl;
  Assert(infits == params["root"] + params["fits_ext"]);

#ifdef _OPENMP
  if (params.keyExists("num_threads")) {
    int num_threads = params["num_threads"];
    omp_set_num_threads(num_threads);
  }
#endif

  dbg<<"Config params = \n"<<params<<std::endl;

  DoMeasureShear(params);

  if (dbgout && dbgout != &std::cout) {delete dbgout; dbgout=0;}
  return 0;
}
#ifndef NOTHROW
#if 0
// Change to 1 to let gdb see where the program bombed out.
catch(int) {}
#else
catch (tmv::Error& e)
{
  std::cerr<<"Caught "<<e<<std::endl;
  return 1;
}
catch (std::exception& e)
{
  std::cerr<<"Caught std::exception: \n"<<e.what()<<std::endl;
}
catch (...)
{
  std::cerr<<"Caught Unknown error\n";
}
#endif
#endif
