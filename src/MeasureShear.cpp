
#include <sys/time.h>
#include "Image.h"
#include "InputCatalog.h"
#include "FittedPsf.h"
#include "ShearCatalog.h"
#include "Log.h"
#include "Scripts.h"
#include "BasicSetup.h"

int main(int argc, char **argv) try 
{
    ConfigFile params;
    if (basicSetup(argc,argv,params,"measureshear")) return EXIT_FAILURE;

    // Setup Log
    std::string logFile = ""; // Default is to stdout
    if (params.keyExists("log_file") || params.keyExists("log_ext")) 
        logFile = makeName(params,"log",false,false);
    std::string shearFile = makeName(params,"shear",false,false);
    std::auto_ptr<ShearLog> log (
        new ShearLog(params,logFile,shearFile)); 

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
            std::cout<<"Time: Open image = "<<t2-t1<<std::endl;
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

        bool nostars = params.read("cat_no_stars",false);
        if (!nostars) {
            // Read star catalog info
            StarCatalog starCat(params);
            starCat.read();

            if (isTiming) {
                gettimeofday(&tp,0);
                t2 = tp.tv_sec + tp.tv_usec/1.e6;
                std::cout<<"Time: Read StarCatalog = "<<t2-t1<<std::endl;
                t1 = t2;
            }

            // Flag known stars as too small to bother trying to measure 
            // the shear.
            inCat.flagStars(starCat);

            if (isTiming) {
                gettimeofday(&tp,0);
                t2 = tp.tv_sec + tp.tv_usec/1.e6;
                std::cout<<"Time: Flag stars = "<<t2-t1<<std::endl;
                t1 = t2;
            }
        }

        // Read the fitted psf file
        FittedPsf fitPsf(params);
        fitPsf.read();

        if (isTiming) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: Read FittedPSF = "<<t2-t1<<std::endl;
            t1 = t2;
        }

        std::auto_ptr<ShearCatalog> shearCat;
        doMeasureShear(
            params,*log,im,weightIm.get(),trans,inCat,fitPsf,
            shearCat);
    }
#if 0
    // Change to 1 to let gdb see where the program bombed out.
    catch(int) {}
#else
    CATCHALL;
#endif

    if (dbgout && dbgout != &std::cout) {delete dbgout; dbgout=0;}
    return EXIT_SUCCESS; // = 0 typically.  Defined in <cstdlib>
} catch (std::exception& e) {
    std::cerr<<"Fatal error: Caught \n"<<e.what()<<std::endl;
    std::cout<<"STATUS5BEG Fatal error: "<<e.what()<<" STATUS5END\n";
    return EXIT_FAILURE;
} catch (...) {
    std::cerr<<"Fatal error: Cought an exception.\n";
    std::cout<<"STATUS5BEG Fatal error: unknown exception STATUS5END\n";
    return EXIT_FAILURE;
}
