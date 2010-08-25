
#include <iostream>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "MyMatrix.h"
#include "ConfigFile.h"
#include "dbg.h"
#include "Name.h"
#include "fp.h" // Generated with xxd -i fitsparams.config fp.h

std::ostream* dbgout = 0;
bool XDEBUG = false;

// Some things that are done at the beginning of each executable
inline int basicSetup(
    int argc, char **argv,
    ConfigFile& params, std::string exec)
{
    // Check args:
    if (argc < 2) {
        std::cerr<<
            "Usage: "<<exec<<" configfile [param=value ...]\n"
            "\tThe first parameter is the configuration file that has \n"
            "\tall the parameters for this run. \n"
            "\tThese values may be modified on the command line by \n"
            "\tentering param/value pais as param=value. \n"
            "Note: root is not usuallly given in the parameter file, \n"
            "\tso the normal command line would be something like:\n"
            "\t"<<exec<<" "<<exec<<".config root=img123\n";
        return 1;
    }

    // Read parameters
    // These are all the defaults, but might as well be explicit.
    params.setDelimiter("=");
    params.setInclude("+");
    params.setComment("#");
    params.load(argv[1]);
    for(int k=2;k<argc;k++) params.append(argv[k]);

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
    setRoot(params);

    // Setup debugging
    if (params.read("verbose",0) > 0) {
        if (params.read<int>("verbose") > 1) XDEBUG = true;
        if (params.keyExists("debug_file") || params.keyExists("debug_ext"))
        {
            std::string debugFile = makeName(params,"debug",false,false);
            dbgout = new std::ofstream(debugFile.c_str());

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
                int threadNum = omp_get_thread_num();
                std::stringstream ss;
                ss << threadNum;
                std::string debugFile2 = debugFile + "_" + ss.str();
                if (threadNum > 0) {
                    // This is a memory leak, but a tiny one.
                    dbgout = new std::ofstream(debugFile2.c_str());
                    dbgout->setf(std::ios_base::unitbuf);
                }
            }
#endif
        }
        else dbgout = &std::cout;
        dbgout->setf(std::ios_base::unitbuf);
    }

#ifdef USE_TMV
    // Send TMV warnings to the debug output
    tmv::WriteWarningsTo(dbgout);
#endif

    // Read fits params
    std::string fp((const char*)fitsparams_config,fitsparams_config_len);
    std::istringstream is(fp);
    params.read(is);
    dbg<<"Config params = \n"<<params<<std::endl;

    return 0;
}

#ifdef USE_TMV
#define CATCHTMV \
    catch (tmv::Error& e) { \
        dbg<<"Caught \n"<<e<<std::endl; \
        std::cerr<<"Caught \n"<<e<<std::endl; \
        if (log.get()) { \
            log->setExitCode(FAILURE_PROCESSING_ERROR, e.what()); \
        } \
        return EXIT_FAILURE; \
    }
#else
#define CATCHTMV
#endif

#define CATCHALL \
    CATCHTMV catch (FileNotFoundException& e) { \
        dbg<<"Caught \n"<<e.what()<<std::endl; \
        std::cerr<<"Caught \n"<<e.what()<<std::endl; \
        if (log.get()) { \
            log->setExitCode(FAILURE_FILE_NOT_FOUND, e.what()); \
        } \
        return EXIT_FAILURE; \
    } catch (ParameterException& e) { \
        dbg<<"Caught \n"<<e.what()<<std::endl; \
        std::cerr<<"Caught \n"<<e.what()<<std::endl; \
        if (log.get()) { \
            log->setExitCode(FAILURE_PROCESSING_ERROR, e.what()); \
        } \
        return EXIT_FAILURE; \
    } catch (ReadException& e) { \
        dbg<<"Caught \n"<<e.what()<<std::endl; \
        std::cerr<<"Caught \n"<<e.what()<<std::endl; \
        if (log.get()) { \
            log->setExitCode(FAILURE_READ_ERROR, e.what()); \
        } \
        return EXIT_FAILURE; \
    } catch (WriteException& e) { \
        dbg<<"Caught \n"<<e.what()<<std::endl; \
        std::cerr<<"Caught \n"<<e.what()<<std::endl; \
        if (log.get()) { \
            log->setExitCode(FAILURE_WRITE_ERROR, e.what()); \
        } \
        return EXIT_FAILURE; \
    } catch (ProcessingException& e) { \
        dbg<<"Caught \n"<<e.what()<<std::endl; \
        std::cerr<<"Caught \n"<<e.what()<<std::endl; \
        if (log.get()) { \
            log->setExitCode(FAILURE_PROCESSING_ERROR, e.what()); \
        } \
        return EXIT_FAILURE; \
    } catch (std::exception& e) { \
        dbg<<"Caught \n"<<e.what()<<std::endl; \
        std::cerr<<"Caught \n"<<e.what()<<std::endl; \
        if (log.get()) { \
            log->setExitCode(FAILURE_PROCESSING_ERROR, e.what()); \
        } \
        return EXIT_FAILURE; \
    } catch (...) { \
        dbg<<"Caught Unknown error\n"; \
        std::cerr<<"Caught Unknown error\n"; \
        if (log.get()) { \
            log->setExitCode(FAILURE, "Caught unknown exception"); \
        } \
        return EXIT_FAILURE; \
    }

