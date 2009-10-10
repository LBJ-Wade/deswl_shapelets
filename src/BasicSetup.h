
#include <iostream>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "TMV.h"
#include "ConfigFile.h"
#include "dbg.h"
#include "Name.h"
#include "fp.h" // Generated with xxd -i fitsparams.config fp.h

std::ostream* dbgout = 0;
bool XDEBUG = false;

// Some things that are done at the beginning of each executable
inline int BasicSetup(int argc, char **argv,
    ConfigFile& params, std::string exec)
{
  // Check args:
  if (argc < 2) {
    std::cerr<<"Usage: "<<exec<<" configfile [param=value ...]\n";
    std::cerr<<"\tThe first parameter is the configuration file that has \n";
    std::cerr<<"\tall the parameters for this run. \n";
    std::cerr<<"\tThese values may be modified on the command line by \n";
    std::cerr<<"\tentering param/value pais as param=value. \n";
    std::cerr<<"Note: root is not usuallly given in the parameter file, \n";
    std::cerr<<"\tso the normal command line would be something like:\n";
    std::cerr<<"\t"<<exec<<" "<<exec<<".config root=img123\n";
    return 1;
  }

  // Read parameters
  // These are all the defaults, but might as well be explicit.
  params.setDelimiter("=");
  params.setInclude("+");
  params.setComment("#");
  params.Load(argv[1]);
  for(int k=2;k<argc;k++) params.Append(argv[k]);

  // Set number of openmp threads if necessary
#ifdef _OPENMP
  if (params.keyExists("omp_num_threads")) {
    int num_threads = params["omp_num_threads"];
    omp_set_num_threads(num_threads);
  }
#endif

  // Set root
  if ( !params.keyExists("root") && 
       !params.keyExists("image_file") &&
       params.keyExists("coaddcat_file") ) 
  {
    // Then use coaddcat_file rather than image_file to make root.
    params["image_file"] = params["coaddcat_file"];
  }
  SetRoot(params);

  // Setup debugging
  if (params.read("verbose",0) > 0) {
    if (params.read<int>("verbose") > 1) XDEBUG = true;
    if (params.keyExists("debug_file") || params.keyExists("debug_ext"))
    {
      std::string dbgfile = Name(params,"debug");
      dbgout = new std::ofstream(dbgfile.c_str());

#ifdef _OPENMP
      // For openmp runs, we use a cool feature known as threadprivate 
      // variables.  
      // In dbg.h, dbgout and XDEBUG are both set to be threadprivate.
      // This means that openmp sets up a separate value for each that
      // persists between threads.  
      // So here, we open a parallel block and initialize each thread's
      // copy of dbgout to be a different file.

      // To use this feature, dynamic threads must be off.  (Otherwise,
      // openmp doesn't know how many copies of each variable to make.)
      omp_set_dynamic(0); 

#pragma omp parallel copyin(dbgout, XDEBUG)
      {
        int threadnum = omp_get_thread_num();
        std::stringstream ss;
        ss << threadnum;
        std::string dbgfile2 = dbgfile + "_" + ss.str();
        if (threadnum > 0)
          // This is a memory leak, but a tiny one.
          dbgout = new std::ofstream(dbgfile2.c_str());
      }
#endif
    }
    else dbgout = &std::cout;
  }

  // Send TMV warnings to the debug output
  tmv::WriteWarningsTo(dbgout);

  // Read fits params
  dbg<<"Config params = \n"<<params<<std::endl;
  std::string fp((const char*)fitsparams_config,fitsparams_config_len);
  std::istringstream is(fp);
  params.Read(is);

  return 0;
}

#define CATCHALL \
catch (FileNotFound& e) \
{ \
  dbg<<"Caught \n"<<e.what()<<std::endl; \
  std::cerr<<"Caught \n"<<e.what()<<std::endl; \
  if (log.get()) { \
    log->exitcode = FAILURE_FILE_NOT_FOUND; \
    log->extraexitinfo = e.what(); \
  } \
  return EXIT_FAILURE; \
} \
catch (ParameterError& e) \
{ \
  dbg<<"Caught \n"<<e.what()<<std::endl; \
  std::cerr<<"Caught \n"<<e.what()<<std::endl; \
  if (log.get()) { \
    log->exitcode = FAILURE_PROCESSING_ERROR; \
    log->extraexitinfo = e.what(); \
  } \
  return EXIT_FAILURE; \
}  \
catch (ReadError& e) \
{ \
  dbg<<"Caught \n"<<e.what()<<std::endl; \
  std::cerr<<"Caught \n"<<e.what()<<std::endl; \
  if (log.get()) { \
    log->exitcode = FAILURE_READ_ERROR; \
    log->extraexitinfo = e.what(); \
  } \
  return EXIT_FAILURE; \
} \
catch (WriteError& e) \
{ \
  dbg<<"Caught \n"<<e.what()<<std::endl; \
  std::cerr<<"Caught \n"<<e.what()<<std::endl; \
  if (log.get()) { \
    log->exitcode = FAILURE_WRITE_ERROR; \
    log->extraexitinfo = e.what(); \
  } \
  return EXIT_FAILURE; \
} \
catch (ProcessingError& e) \
{ \
  dbg<<"Caught \n"<<e.what()<<std::endl; \
  std::cerr<<"Caught \n"<<e.what()<<std::endl; \
  if (log.get()) { \
    log->exitcode = FAILURE_PROCESSING_ERROR; \
    log->extraexitinfo = e.what(); \
  } \
  return EXIT_FAILURE; \
} \
catch (ConfigFile_FileNotFound& e) \
{ \
  dbg<<"Caught \n"<<e.what()<<std::endl; \
  std::cerr<<"Caught \n"<<e.what()<<std::endl; \
  if (log.get()) { \
    log->exitcode = FAILURE_FILE_NOT_FOUND; \
    log->extraexitinfo = e.what(); \
  } \
  return EXIT_FAILURE; \
} \
catch (ConfigFile_KeyNotFound& e) \
{ \
  dbg<<"Caught \n"<<e.what()<<std::endl; \
  std::cerr<<"Caught \n"<<e.what()<<std::endl; \
  if (log.get()) { \
    log->exitcode = FAILURE_PARAMETER_ERROR; \
    log->extraexitinfo = e.what(); \
  } \
  return EXIT_FAILURE; \
} \
catch (tmv::Error& e) \
{ \
  dbg<<"Caught \n"<<e<<std::endl; \
  std::cerr<<"Caught \n"<<e<<std::endl; \
  if (log.get()) { \
    log->exitcode = FAILURE_PROCESSING_ERROR; \
    log->extraexitinfo = e.what(); \
  } \
  return EXIT_FAILURE; \
} \
catch (std::exception& e) \
{ \
  dbg<<"Caught \n"<<e.what()<<std::endl; \
  std::cerr<<"Caught \n"<<e.what()<<std::endl; \
  if (log.get()) { \
    log->exitcode = FAILURE_PROCESSING_ERROR; \
    log->extraexitinfo = e.what(); \
  } \
  return EXIT_FAILURE; \
} \
catch (...) \
{ \
  dbg<<"Caught Unknown error\n"; \
  std::cerr<<"Caught Unknown error\n"; \
  if (log.get()) { \
    log->exitcode = FAILURE; \
    log->extraexitinfo = "Caught unknown exception"; \
  } \
  return EXIT_FAILURE; \
}

