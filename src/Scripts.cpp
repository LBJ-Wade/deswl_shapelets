
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

void DoFindStars(
    ConfigFile& params, FindStarsLog& log,
    const Image<double>& im, const Image<double>* weight_image,
    const Transformation& trans, 
    const InputCatalog& incat,
    std::auto_ptr<StarCatalog>& starcat)
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
    starcat.reset(new StarCatalog(incat,params));

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Make StarCatalog = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    // Update the sizes to more robust values
    starcat->calculateSizes(im,weight_image,trans);

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: CalcSizes = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    try {
        starcat->findStars(log);
    } catch(StarFinderException& e) {
        // Need to catch this here, so we can write the output file
        // with the sizes, even though we haven't figured out which 
        // objects are stars.
        dbg<<"Caught StarFinderException: "<<e.what()<<std::endl;  
        starcat->write();
        throw;
    } catch (...) {
        dbg<<"Caught unknown exception\n";
        starcat->write();
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
    starcat->write();

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Write StarCatalog = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    xdbg<<"FindStars Log: \n"<<log<<std::endl;
}

void DoMeasurePsf(
    ConfigFile& params, PsfLog& log,
    const Image<double>& im, const Image<double>* weight_image,
    const Transformation& trans, 
    const StarCatalog& starcat,
    std::auto_ptr<PsfCatalog>& psfcat, std::auto_ptr<FittedPsf>& fitpsf,
    double& sigma_p)
{
    dbg<<"Starting MeasurePsf script\n";

    bool isTiming = params.read("timing",false);
    timeval tp;
    double t1=0.,t2=0.;

    if (isTiming) {
        gettimeofday(&tp,0);
        t1 = tp.tv_sec + tp.tv_usec/1.e6;
    }

    if (params.read("psf_skip_measurements",false)) {
        // Option to read existing PsfCatalog rather than remeasure.
        // (Useful if you only want to redo the fitting step.)
        psfcat.reset(new PsfCatalog(params));
        psfcat->read();

        if (isTiming) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: Read PSFCatalog = "<<t2-t1<<std::endl;
            t1 = t2;
        }
    } else {
        // Create PsfCatalog from StarCatalog
        psfcat.reset(new PsfCatalog(starcat,params));

        if (isTiming) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: Create PSFCatalog = "<<t2-t1<<std::endl;
            t1 = t2;
        }

        // Estimate the scale size to use for shapelet decompositions
        if (sigma_p == 0.)
            sigma_p = psfcat->estimateSigma(im,weight_image,trans);

        if (isTiming) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: Estimate Sigma = "<<t2-t1<<std::endl;
            t1 = t2;
        }

        // Do the actual PSF measurements
        int npsf = psfcat->measurePsf(im,weight_image,trans,sigma_p,log);

        if (isTiming) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: Measure PSF = "<<t2-t1<<std::endl;
            t1 = t2;
        }

        // Write PSF catalog to file
        psfcat->write();

        if (isTiming) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            std::cout<<"Time: Write PSFCatalog = "<<t2-t1<<std::endl;
            t1 = t2;
        }

        if (npsf == 0) {
            throw ProcessingException(
                "No successful PSF measurements");
        }
    }


    // Fit the PSF with a polynomial:
    fitpsf.reset(new FittedPsf(*psfcat,params,log));

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Fit PSF = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    // Write fitted psf to file
    fitpsf->write();

    // Re-write the PSF catalog, since the interpolation may have changed
    // the flags.
    // TODO: It may be worth having a new routine that just updates the 
    // flags.  More efficient, since don't need to re-write everything.
    psfcat->write();

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Write FittedPSF = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    xdbg<<"PSF Log: \n"<<log<<std::endl;
}

void DoMeasureShear(
    ConfigFile& params, ShearLog& log,
    const Image<double>& im, const Image<double>* weight_image,
    const Transformation& trans, 
    const InputCatalog& incat, const FittedPsf& fitpsf,
    std::auto_ptr<ShearCatalog>& shearcat)
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
    shearcat.reset(new ShearCatalog(incat,trans,fitpsf,params));

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Create ShearCatalog = "<<t2-t1<<std::endl;
        t1 = t2;
    }

    // Measure shears and shapelet vectors
    int nShear = shearcat->measureShears(im,weight_image,log);

    if (isTiming) {
        gettimeofday(&tp,0);
        t2 = tp.tv_sec + tp.tv_usec/1.e6;
        std::cout<<"Time: Measure Shears = "<<t2-t1<<std::endl;
        std::cout<<"Rate: "<<(t2-t1)/shearcat->size()<<" s / gal\n";
        t1 = t2;
    }

    // Write results to file
    shearcat->write();

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

void DoSplitStars(
    ConfigFile& params, std::string log_file, std::auto_ptr<Log>& log,
    const Image<double>& im, const Image<double>* weight_image,
    const Transformation& trans, 
    const InputCatalog& incat, const StarCatalog& starcat,
    double sigma_p)
{
    dbg<<"Starting SplitStars script\n";

    // split the catalog and perform the psf analysis on both sets
    // separately.  We work in fits only to avoid dealing with all the
    // different naming cases.

    std::string stars_file = MakeFitsName(params, "stars");
    std::string psf_file = MakeFitsName(params, "psf");
    std::string fitpsf_file = MakeFitsName(params, "fitpsf");
    std::string shear_file = MakeFitsName(params, "shear");

    std::string stars_file1 = AddExtraToName(stars_file, "1");
    std::string stars_file2 = AddExtraToName(stars_file, "2");
    starcat.splitInTwo(stars_file1, stars_file2);

    // set 1
    dbg<<"\nDoing PSF/Shear for Set1\n";
    ConfigFile params1 = params;
    params1["psf_file"] = AddExtraToName(psf_file, "1");
    params1["fitpsf_file"] = AddExtraToName(fitpsf_file, "1");
    params1["shear_file"] = AddExtraToName(shear_file, "1");
    StarCatalog starcat1(params);
    starcat1.readFits(stars_file1);
    std::auto_ptr<PsfCatalog> psfcat1;
    std::auto_ptr<FittedPsf> fitpsf1;
    std::auto_ptr<ShearCatalog> shearcat1;
    log.reset(new PsfLog(params1,log_file,params1["psf_file"]));

    DoMeasurePsf(
        params1, static_cast<PsfLog&>(*log),
        im, weight_image, trans, starcat1,
        psfcat1, fitpsf1, sigma_p);

    log.reset(new ShearLog(params,log_file,params1["shear_file"]));

    DoMeasureShear(
        params1, static_cast<ShearLog&>(*log),
        im, weight_image, trans, incat, *fitpsf1, shearcat1);

    // set 2
    dbg<<"\nDoing PSF/Shear for Set2\n";
    ConfigFile params2 = params;
    params2["psf_file"] = AddExtraToName(psf_file, "2");
    params2["fitpsf_file"] = AddExtraToName(fitpsf_file, "2");
    params2["shear_file"] = AddExtraToName(shear_file, "2");
    StarCatalog starcat2(params2);
    starcat2.readFits(stars_file2);
    std::auto_ptr<PsfCatalog> psfcat2;
    std::auto_ptr<FittedPsf> fitpsf2;
    std::auto_ptr<ShearCatalog> shearcat2;
    log.reset(new PsfLog(params2,log_file,params2["psf_file"]));

    DoMeasurePsf(
        params2, static_cast<PsfLog&>(*log),
        im, weight_image, trans, starcat2,
        psfcat2, fitpsf2, sigma_p);

    log.reset(new ShearLog(params,log_file,params2["shear_file"]));

    DoMeasureShear(
        params2, static_cast<ShearLog&>(*log),
        im, weight_image, trans, incat, *fitpsf2, shearcat2);
}

