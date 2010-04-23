
#include <valarray>
#include "Image.h"
#include "Transformation.h"
#include "InputCatalog.h"
#include "StarCatalog.h"
#include "FittedPsf.h"
#include "PsfCatalog.h"
#include "StarFinder.h"
#include "ShearCatalog.h"
#include "TimeVars.h"
#include "Log.h"
#include "BasicSetup.h"

static void doFullPipeline(
    ConfigFile& params, std::string logFile, std::auto_ptr<Log>& log)
{
    log.reset(
        new FindStarsLog(params,logFile,makeName(params,"stars",false,false)));

    // Load image:

    std::cerr<<"Running FindStars\n";
    std::auto_ptr<Image<double> > weightIm;
    Image<double> im(params,weightIm);

    // Read distortion function
    Transformation trans(params);

    // Read input catalog
    InputCatalog inCat(params,&im);
    inCat.read();

    // Create StarCatalog from InputCatalog
    StarCatalog starCat(inCat,params);

    // Update the sizes to more robust values
    starCat.calculateSizes(im,weightIm.get(),trans);

    try {
        starCat.findStars(static_cast<FindStarsLog&>(*log));
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

    // Write star catalog to file
    starCat.write();

    xdbg<<"FindStarsLog: \n"<<*log<<std::endl;

    log.reset(new PsfLog(params,logFile,makeName(params,"psf",false,false)));

    std::cerr<<"Determining PSF\n";
    // Create PsfCatalog from StarCatalog
    PsfCatalog psfCat(starCat,params);

    // Estimate the scale size to use for shapelet decompositions
    double sigmaP = psfCat.estimateSigma(im,weightIm.get(),trans);

    // Do the actual PSF measurements
    psfCat.measurePsf(
        im,weightIm.get(),trans,sigmaP,static_cast<PsfLog&>(*log));

    // Write PSF catalog to file
    psfCat.write();

    // MJ -- Should we make a FittedPsfLog?  Are there parameters here that
    //       we want to ingest into a DB meta-table?
    //       If so, we would need to WritePsfCat after doing the
    //       FittedPsf calculation, since the fittedpsf file is not ingested.

    // Fit the PSF with a polynomial:
    FittedPsf fitPsf(psfCat,params);

    // Write fitted psf to file
    fitPsf.write();

    xdbg<<"PSFLog: \n"<<*log<<std::endl;

    log.reset(new ShearLog(params,logFile,makeName(params,"shear",false,false)));

    std::cerr<<"Measuring Shears\n";
    // Create shear catalog
    ShearCatalog shearCat(inCat,trans,fitPsf,params);

    // Measure shears and shapelet vectors
    shearCat.measureShears(im,weightIm.get(),static_cast<ShearLog&>(*log));

    // Write results to file
    shearCat.write();

    xdbg<<"ShearLog: \n"<<*log<<std::endl;

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

