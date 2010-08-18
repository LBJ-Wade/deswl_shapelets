
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
    FittedPsf fitPsf(psfCat,params,static_cast<PsfLog&>(*log));

    // Write fitted psf to file
    fitPsf.write();

    // Re-write the PSF catalog, since the interpolation may have changed
    // the flags.
    // TODO: It may be worth having a new routine that just updates the 
    // flags.  More efficient, since don't need to re-write everything.
    psfCat.write();

    xdbg<<"PSFLog: \n"<<*log<<std::endl;




    log.reset(new ShearLog(params,logFile,makeName(params,"shear",false,false)));

    std::cerr<<"Measuring Shears\n";
    // Create shear catalog
    ShearCatalog shearCat(inCat,trans,fitPsf,params);

    // Flag known stars as too small to bother trying to measure the shear.
    shearCat.flagStars(starCat);

    // Measure shears and shapelet vectors
    shearCat.measureShears(im,weightIm.get(),static_cast<ShearLog&>(*log));

    // Write results to file
    shearCat.write();

    xdbg<<"ShearLog: \n"<<*log<<std::endl;





    bool splitstars=params.read("splitstars",false);
    if (splitstars) {

      std::cerr<<"\n\nSplitting stars into two sets\n";

      // split the catalog and perform the psf analysis on both sets
      // separately.  We work in fits only to avoid dealing with all the
      // different naming cases.


      std::string star_file = makeFitsName(params, "stars");
      std::string psf_file = makeFitsName(params, "psf");
      std::string fitpsf_file = makeFitsName(params, "fitpsf");
      std::string shear_file = makeFitsName(params, "shear");

      // this will turn _stars.fits into _stars1.fits
      std::string star_file1 = addExtraToName(star_file, "1");
      std::string star_file2 = addExtraToName(star_file, "2");

      std::string psf_file1 = addExtraToName(psf_file, "1");
      std::string psf_file2 = addExtraToName(psf_file, "2");
      std::string fitpsf_file1 = addExtraToName(fitpsf_file, "1");
      std::string fitpsf_file2 = addExtraToName(fitpsf_file, "2");

      std::string shear_file1 = addExtraToName(shear_file, "1");
      std::string shear_file2 = addExtraToName(shear_file, "2");

      starCat.splitInTwo(star_file1, star_file2);

      //
      // set 1
      //
      std::cerr<<"\nDoing PSF/Shear for Set1\n";
      log.reset(new PsfLog(params,logFile,psf_file1));

      StarCatalog starCat1(params);
      starCat1.readFits(star_file1);

      PsfCatalog psfCat1(starCat1,params);
      psfCat1.measurePsf(
	  im,weightIm.get(),trans,sigmaP,static_cast<PsfLog&>(*log));

      psfCat1.writeFits(psf_file1);

      FittedPsf fitPsf1(psfCat1,params,static_cast<PsfLog&>(*log));

      fitPsf1.writeFits(fitpsf_file1);
      // rewrite since flags may have changed from outlier checking.
      psfCat1.writeFits(psf_file1);

      log.reset(new ShearLog(params,logFile,shear_file1));

      ShearCatalog shearCat1(inCat,trans,fitPsf1,params);

      shearCat1.measureShears(im,weightIm.get(),static_cast<ShearLog&>(*log));
      shearCat1.writeFits(shear_file1);

      //
      // set 2
      //
      std::cerr<<"\nDoing PSF/Shear for Set2\n";
      log.reset(new PsfLog(params,logFile,psf_file2));

      StarCatalog starCat2(params);
      starCat2.readFits(star_file2);

      PsfCatalog psfCat2(starCat2,params);
      psfCat2.measurePsf(
	  im,weightIm.get(),trans,sigmaP,static_cast<PsfLog&>(*log));

      psfCat2.writeFits(psf_file2);


      FittedPsf fitPsf2(psfCat2,params,static_cast<PsfLog&>(*log));

      fitPsf2.writeFits(fitpsf_file2);
      // rewrite since flags may have changed from outlier checking.
      psfCat2.writeFits(psf_file2);


      log.reset(new ShearLog(params,logFile,shear_file2));

      ShearCatalog shearCat2(inCat,trans,fitPsf2,params);

      shearCat2.measureShears(im,weightIm.get(),static_cast<ShearLog&>(*log));
      shearCat2.writeFits(shear_file2);



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

