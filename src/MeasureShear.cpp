
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
    std::cerr<<"Usage: measureshear configfile [param=value ...]\n";
    std::cerr<<"The first parameter is the configuration file that has \n";
    std::cerr<<"all the parameters for this run. \n";
    std::cerr<<"These values may be modified on the command line by \n";
    std::cerr<<"entering param/value pais as param=value. \n";
    std::cerr<<"Note: root is not usuallly given in the parameter file, \n";
    std::cerr<<"so the normal command line would be something like:\n";
    std::cerr<<"measureshear measureshear.config root=img123\n";
    return 1;
  }
  ConfigFile params(argv[1]);
  for(int k=2;k<argc;k++) params.Append(argv[k]);

  // Setup debugging
  if (params.keyExists("verbose") && int(params["verbose"]) > 0) {
    if (params.keyExists("debug_file") || params.keyExists("debug_ext")) {
      dbgout = new std::ofstream(Name(params,"debug").c_str());
    }
    else dbgout = &std::cout;
    if (int(params["verbose"]) > 1) XDEBUG = true;
  }


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
