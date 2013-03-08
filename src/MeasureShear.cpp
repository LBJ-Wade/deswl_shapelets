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
    if (BasicSetup(argc,argv,params,"measureshear")) return EXIT_FAILURE;

    // Setup Log
    std::string log_file = ""; // Default is to stdout
    if (params.keyExists("log_file") || params.keyExists("log_ext")) 
        log_file = MakeName(params,"log",false,false);
    std::string shear_file = MakeName(params,"shear",false,false);
    std::auto_ptr<ShearLog> log (
        new ShearLog(params,log_file,shear_file)); 

    try {
        bool timing = params.read("timing",false);
        timeval tp;
        double t1=0.,t2=0.;

        if (timing) {
            gettimeofday(&tp,0);
            t1 = tp.tv_sec + tp.tv_usec/1.e6;
        }

        // Load image:
        std::auto_ptr<Image<double> > weight_image;
        Image<double> im(params,weight_image);

        if (timing) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: Open image = "<<t2-t1<<std::endl;
            t1 = t2;
        }

        // Read distortion function
        Transformation trans(params);

        if (timing) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: Read Transformation = "<<t2-t1<<std::endl;
            t1 = t2;
        }

        // Read input catalog
        InputCatalog incat(params,&im);
        incat.read();

        if (timing) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: Read InputCatalog = "<<t2-t1<<std::endl;
            t1 = t2;
        }

        bool nostars = params.read("cat_no_stars",false);
        if (!nostars) {
            // Read star catalog info
            StarCatalog starcat(params);
            starcat.read();

            if (timing) {
                gettimeofday(&tp,0);
                t2 = tp.tv_sec + tp.tv_usec/1.e6;
                std::cout<<"Time: Read StarCatalog = "<<t2-t1<<std::endl;
                t1 = t2;
            }

            // Flag known stars as too small to bother trying to measure 
            // the shear.
            incat.flagStars(starcat);

            if (timing) {
                gettimeofday(&tp,0);
                t2 = tp.tv_sec + tp.tv_usec/1.e6;
                std::cout<<"Time: Flag stars = "<<t2-t1<<std::endl;
                t1 = t2;
            }
        }

        // Read the fitted psf file
        FittedPsf fitpsf(params);
        fitpsf.read();

        if (timing) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: Read FittedPSF = "<<t2-t1<<std::endl;
            t1 = t2;
        }

        std::auto_ptr<ShearCatalog> shearcat;
        DoMeasureShear(
            params,*log,im,weight_image.get(),trans,incat,fitpsf,
            shearcat);
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
