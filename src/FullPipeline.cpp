
#include "Image.h"
#include "Transformation.h"
#include "InputCatalog.h"
#include "StarCatalog.h"
#include "FittedPsf.h"
#include "PsfCatalog.h"
#include "ShearCatalog.h"
#include "Log.h"
#include "Scripts.h"
#include "BasicSetup.h"

static void doFullPipeline(
    ConfigFile& params, std::string logFile, std::auto_ptr<Log>& log)
{
    // Load image:
    std::auto_ptr<Image<double> > weightIm;
    Image<double> im(params,weightIm);

    // Read distortion function
    Transformation trans(params);

    // Read input catalog
    InputCatalog inCat(params,&im);
    inCat.read();

    // Do FindStars script
    log.reset(
        new FindStarsLog(params,logFile,makeName(params,"stars",false,false)));
    std::auto_ptr<StarCatalog> starCat;
    doFindStars(
        params, static_cast<FindStarsLog&>(*log), im, weightIm.get(), trans,
        inCat, starCat);

    // Do MeasurePsf script
    log.reset(new PsfLog(params,logFile,makeName(params,"psf",false,false)));
    std::auto_ptr<PsfCatalog> psfCat;
    std::auto_ptr<FittedPsf> fitPsf;
    double sigmaP = 0.;
    doMeasurePsf(
        params, static_cast<PsfLog&>(*log), im, weightIm.get(), trans,
        *starCat, psfCat, fitPsf, sigmaP);

    // Flag stars, so don't try to measure shears for them.
    bool nostars = params.read("cat_no_stars",false);
    if (!nostars) inCat.flagStars(*starCat);

    // Do MeasusreShear script
    log.reset(new ShearLog(params,logFile,makeName(params,"shear",false,false)));
    std::auto_ptr<ShearCatalog> shearCat;
    doMeasureShear(
        params, static_cast<ShearLog&>(*log), im, weightIm.get(), trans,
        inCat, *fitPsf, shearCat);

    // Maybe do SplitStars script
    bool splitstars=params.read("splitstars",false);
    if (splitstars) {
        doSplitStars(
            params, logFile, log, im, weightIm.get(), trans,
            inCat, *starCat, sigmaP);
    }
}

int main(int argc, char **argv) try 
{
    ConfigFile params;
    if (basicSetup(argc,argv,params,"fullpipe")) return EXIT_FAILURE;

    // Setup Log
    std::string logFile = ""; // Default is to stdout
    if (params.keyExists("log_file") || params.keyExists("log_ext")) 
        logFile = makeName(params,"log",false,false);
    std::auto_ptr<Log> log;

    try {
        doFullPipeline(params,logFile,log);
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

