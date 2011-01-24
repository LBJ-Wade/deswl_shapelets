
#include <sys/time.h>
#include "Image.h"
#include "Transformation.h"
#include "InputCatalog.h"
#include "StarCatalog.h"
#include "StarFinder.h"
#include "FittedPsf.h"
#include "PsfCatalog.h"
#include "ShearCatalog.h"
#include "Log.h"
#include "Scripts.h"

void doFindStars(
    ConfigFile& params, FindStarsLog& log,
    const Image<double>& im, const Image<double>* weightIm,
    const Transformation& trans, 
    const InputCatalog& inCat,
    std::auto_ptr<StarCatalog>& starCat)
{
    dbg<<"Starting FindStars script\n";

    bool isTiming = params.read("timing",false);
    timeval tp;
    double t1=0.,t2=0.;

    if (isTiming) {
        gettimeofday(&tp,0);
        t1 = tp.tv_sec + tp.tv_usec/1.e6;
    }

    // Create StarCatalog from InputCatalog
    starCat.reset(new StarCatalog(inCat,params));

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Make StarCatalog = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    // Update the sizes to more robust values
    starCat->calculateSizes(im,weightIm,trans);

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: CalcSizes = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    try {
        starCat->findStars(log);
    } catch(StarFinderException& e) {
        // Need to catch this here, so we can write the output file
        // with the sizes, even though we haven't figured out which 
        // objects are stars.
        dbg<<"Caught StarFinderException: "<<e.what()<<std::endl;  
        starCat->write();
        throw;
    } catch (...) {
        dbg<<"Caught unknown exception\n";
        starCat->write();
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
    starCat->write();

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Write StarCatalog = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    xdbg<<"FindStars Log: \n"<<log<<std::endl;
}

void doMeasurePsf(
    ConfigFile& params, PsfLog& log,
    const Image<double>& im, const Image<double>* weightIm,
    const Transformation& trans, 
    const StarCatalog& starCat,
    std::auto_ptr<PsfCatalog>& psfCat, std::auto_ptr<FittedPsf>& fitPsf,
    double& sigmaP)
{
    dbg<<"Starting MeasurePsf script\n";

    bool isTiming = params.read("timing",false);
    timeval tp;
    double t1=0.,t2=0.;

    if (isTiming) {
        gettimeofday(&tp,0);
        t1 = tp.tv_sec + tp.tv_usec/1.e6;
    }

    // Create PsfCatalog from StarCatalog
    psfCat.reset(new PsfCatalog(starCat,params));

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Create PSFCatalog = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    // Estimate the scale size to use for shapelet decompositions
    if (sigmaP == 0.)
        sigmaP = psfCat->estimateSigma(im,weightIm,trans);

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Estimate Sigma = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    // Do the actual PSF measurements
    int nPsf = psfCat->measurePsf(im,weightIm,trans,sigmaP,log);

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Measure PSF = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    // Write PSF catalog to file
    psfCat->write();

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Write PSFCatalog = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    if (nPsf == 0) {
        throw ProcessingException(
            "No successful PSF measurements");
    }

    // Fit the PSF with a polynomial:
    fitPsf.reset(new FittedPsf(*psfCat,params,log));

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Fit PSF = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    // Write fitted psf to file
    fitPsf->write();

    // Re-write the PSF catalog, since the interpolation may have changed
    // the flags.
    // TODO: It may be worth having a new routine that just updates the 
    // flags.  More efficient, since don't need to re-write everything.
    psfCat->write();

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Write FittedPSF = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    xdbg<<"PSF Log: \n"<<log<<std::endl;
}

void doMeasureShear(
    ConfigFile& params, ShearLog& log,
    const Image<double>& im, const Image<double>* weightIm,
    const Transformation& trans, 
    const InputCatalog& inCat, const FittedPsf& fitPsf,
    std::auto_ptr<ShearCatalog>& shearCat)
{
    dbg<<"Starting MeasureShear script\n";

    bool isTiming = params.read("timing",false);
    timeval tp;
    double t1=0.,t2=0.;

    if (isTiming) {
        gettimeofday(&tp,0);
        t1 = tp.tv_sec + tp.tv_usec/1.e6;
    }

    // Create shear catalog
    shearCat.reset(new ShearCatalog(inCat,trans,fitPsf,params));

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Create ShearCatalog = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    // Measure shears and shapelet vectors
    int nShear = shearCat->measureShears(im,weightIm,log);

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Measure Shears = "<<t2-t1<<std::endl;
        std::cout<<"Rate: "<<(t2-t1)/shearCat->size()<<" s / gal\n";
        t1 = t2;
    }

    // Write results to file
    shearCat->write();

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Write ShearCatalog = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    if (nShear == 0) {
        throw ProcessingException(
            "No successful shear measurements");
    }

    xdbg<<"Shear Log: \n"<<log<<std::endl;
}

void doSplitStars(
    ConfigFile& params, std::string logFile, std::auto_ptr<Log>& log,
    const Image<double>& im, const Image<double>* weightIm,
    const Transformation& trans, 
    const InputCatalog& inCat, const StarCatalog& starCat,
    double sigmaP)
{
    dbg<<"Starting SplitStars script\n";

    // split the catalog and perform the psf analysis on both sets
    // separately.  We work in fits only to avoid dealing with all the
    // different naming cases.

    std::string stars_file = makeFitsName(params, "stars");
    std::string psf_file = makeFitsName(params, "psf");
    std::string fitpsf_file = makeFitsName(params, "fitpsf");
    std::string shear_file = makeFitsName(params, "shear");

    std::string stars_file1 = addExtraToName(stars_file, "1");
    std::string stars_file2 = addExtraToName(stars_file, "2");
    starCat.splitInTwo(stars_file1, stars_file2);

    // set 1
    dbg<<"\nDoing PSF/Shear for Set1\n";
    ConfigFile params1 = params;
    params1["psf_file"] = addExtraToName(psf_file, "1");
    params1["fitpsf_file"] = addExtraToName(fitpsf_file, "1");
    params1["shear_file"] = addExtraToName(shear_file, "1");
    StarCatalog starCat1(params);
    starCat1.readFits(stars_file1);
    std::auto_ptr<PsfCatalog> psfCat1;
    std::auto_ptr<FittedPsf> fitPsf1;
    std::auto_ptr<ShearCatalog> shearCat1;
    log.reset(new PsfLog(params1,logFile,params1["psf_file"]));

    doMeasurePsf(
        params1, static_cast<PsfLog&>(*log),
        im, weightIm, trans, starCat1,
        psfCat1, fitPsf1, sigmaP);

    log.reset(new ShearLog(params,logFile,params1["shear_file"]));

    doMeasureShear(
        params1, static_cast<ShearLog&>(*log),
        im, weightIm, trans, inCat, *fitPsf1, shearCat1);

    // set 2
    dbg<<"\nDoing PSF/Shear for Set2\n";
    ConfigFile params2 = params;
    params2["psf_file"] = addExtraToName(psf_file, "2");
    params2["fitpsf_file"] = addExtraToName(fitpsf_file, "2");
    params2["shear_file"] = addExtraToName(shear_file, "2");
    StarCatalog starCat2(params2);
    starCat2.readFits(stars_file2);
    std::auto_ptr<PsfCatalog> psfCat2;
    std::auto_ptr<FittedPsf> fitPsf2;
    std::auto_ptr<ShearCatalog> shearCat2;
    log.reset(new PsfLog(params2,logFile,params2["psf_file"]));

    doMeasurePsf(
        params2, static_cast<PsfLog&>(*log),
        im, weightIm, trans, starCat2,
        psfCat2, fitPsf2, sigmaP);

    log.reset(new ShearLog(params,logFile,params2["shear_file"]));

    doMeasureShear(
        params2, static_cast<ShearLog&>(*log),
        im, weightIm, trans, inCat, *fitPsf2, shearCat2);
}

