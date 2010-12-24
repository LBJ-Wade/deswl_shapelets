
#include <sys/time.h>
#include "Image.h"
#include "Transformation.h"
#include "PsfCatalog.h"
#include "InputCatalog.h"
#include "FittedPsf.h"
#include "Log.h"
#include "BVec.h"
#include "Scripts.h"
#include "BasicSetup.h"

int main(int argc, char **argv) try 
{
    ConfigFile params;
    if (basicSetup(argc,argv,params,"measurepsf")) return EXIT_FAILURE;

    // Setup Log
    std::string logFile = ""; // Default is to stdout
    if (params.keyExists("log_file") || params.keyExists("log_ext"))
        logFile = makeName(params,"log",false,false);
    std::string psfFile=makeName(params,"psf",false,false);
    std::auto_ptr<PsfLog> log (
        new PsfLog(params,logFile,psfFile)); 

    try {
        bool isTiming = params.read("timing",false);
        timeval tp;
        double t1=0.,t2=0.;

        if (isTiming) {
            gettimeofday(&tp,0);
            t1 = tp.tv_sec + tp.tv_usec/1.e6;
        }

        // Load image:
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

        std::auto_ptr<StarCatalog> starCat;
        if (params.read("cat_all_stars",false)) {
            // Read input catalog
            InputCatalog inCat(params,&im);
            inCat.read();
            starCat.reset(new StarCatalog(inCat,params));
        } else {
            // Read star catalog info
            starCat.reset(new StarCatalog(params));
            starCat->read();
        }

        if (isTiming) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: Read StarCatalog = "<<t2-t1<<std::endl;
            t1 = t2;
        }

        std::auto_ptr<PsfCatalog> psfCat;
        std::auto_ptr<FittedPsf> fitPsf;
        double sigmaP = 0.;
        doMeasurePsf(
            params,*log,im,weightIm.get(),trans,*starCat,
            psfCat,fitPsf,sigmaP);
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
