
//#include <valarray>
#include <sys/time.h>
#include "StarCatalog.h"
#include "InputCatalog.h"
#include "StarFinder.h"
#include "Scripts.h"
#include "BasicSetup.h"

#if 0
static void doFindStars(ConfigFile& params, FindStarsLog& log)
{
    bool isTiming = params.read("timing",false);
    timeval tp;
    double t1=0.,t2=0.;

    if (isTiming) {
        gettimeofday(&tp,0);
        t1 = tp.tv_sec + tp.tv_usec/1.e6;
    }

    // Read image, transformation
    std::auto_ptr<Image<double> > weightIm;
    Image<double> im(params,weightIm);

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Open imgae = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    // Read distortion function
    Transformation trans(params);

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Read Transformation = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    // Read input catalog
    InputCatalog inCat(params,&im);
    inCat.read();

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Read InputCatalog = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    // Create StarCatalog from InputCatalog
    StarCatalog starCat(inCat,params);

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Make StarCatalog = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    // Update the sizes to more robust values
    starCat.calculateSizes(im,weightIm.get(),trans);

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: CalcSizes = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    try {
        starCat.findStars(log);
    } catch(StarFinderException& e) {
        // Need to catch this here, so we can write the output file
        // with the sizes, even though we haven't figured out which 
        // objects are stars.
        dbg<<"Caught StarFinderException: "<<e.what()<<std::endl;  
        starCat.write();
        throw;
    } catch (...) {
        dbg<<"Caught unknown exception\n";
        starCat.write();
        throw;
    }
    dbg<<"After RunFindStars\n";

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: FindStars = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    // Write star catalog to file
    starCat.write();

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Write StarCatalog = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    xdbg<<"Log: \n"<<log<<std::endl;
}
#endif

int main(int argc, char **argv) try 
{
    ConfigFile params;
    if (basicSetup(argc,argv,params,"findstars")) return EXIT_FAILURE;

    // Setup Log
    std::string logFile = ""; // Default is to stdout
    if (params.keyExists("log_file") || params.keyExists("log_ext")) 
        logFile = makeName(params,"log",false,false);
    std::string starsFile=makeName(params,"stars",false,false);
    std::auto_ptr<FindStarsLog> log(
        new FindStarsLog(params,logFile,starsFile)); 

    try {
        bool isTiming = params.read("timing",false);
        timeval tp;
        double t1=0.,t2=0.;

        if (isTiming) {
            gettimeofday(&tp,0);
            t1 = tp.tv_sec + tp.tv_usec/1.e6;
        }

        // Read image, transformation
        std::auto_ptr<Image<double> > weightIm;
        Image<double> im(params,weightIm);

        if (isTiming) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: Open imgae = "<<t2-t1<<std::endl;
            t1 = t2;
        }

        // Read distortion function
        Transformation trans(params);

        if (isTiming) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: Read Transformation = "<<t2-t1<<std::endl;
            t1 = t2;
        }

        // Read input catalog
        InputCatalog inCat(params,&im);
        inCat.read();

        if (isTiming) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: Read InputCatalog = "<<t2-t1<<std::endl;
            t1 = t2;
        }

        std::auto_ptr<StarCatalog> starCat;
        doFindStars(params,*log,im,weightIm.get(),trans,inCat,starCat);
    }
#if 0
    // Change to 1 to let gdb see where the program bombed out.
    catch(int) {}
#else
    CATCHALL;
#endif

    if (dbgout && dbgout != &std::cout) {delete dbgout; dbgout=0;}
    return EXIT_SUCCESS;
} catch (std::exception& e) {
    std::cerr<<"Fatal error: Caught \n"<<e.what()<<std::endl;
    std::cout<<"STATUS5BEG Fatal error: "<<e.what()<<" STATUS5END\n";
    return EXIT_FAILURE;
} catch (...) {
    std::cerr<<"Fatal error: Cought an exception.\n";
    std::cout<<"STATUS5BEG Fatal error: unknown exception STATUS5END\n";
    return EXIT_FAILURE;
}
