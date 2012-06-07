
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

static void DoFullPipeline(
    ConfigFile& params, std::string log_file, std::auto_ptr<Log>& log)
{
    bool output_info = params.read("output_info",true);

    // Load image:
    std::auto_ptr<Image<double> > weight_image;
    Image<double> im(params,weight_image);

    // Read distortion function
    Transformation trans(params);

    // Read input catalog
    InputCatalog incat(params,&im);
    incat.read();

    // Do FindStars script
    std::auto_ptr<StarCatalog> starcat;
    if ( (params.read("cat_all_stars",false) || 
          params.read("stars_trust_sg",false)) ) {
        starcat.reset(new StarCatalog(incat,params));
    } else {
        log.reset(
            new FindStarsLog(params,log_file,MakeName(params,"stars",false,false)));
        std::auto_ptr<StarCatalog> starcat;
        if (output_info) {
            std::cerr<<"Finding Stars"<<std::endl;
        }
        DoFindStars(
            params, static_cast<FindStarsLog&>(*log), im, weight_image.get(), trans,
            incat, starcat);
    }

    // Do MeasurePsf script
    log.reset(new PsfLog(params,log_file,MakeName(params,"psf",false,false)));
    std::auto_ptr<PsfCatalog> psfcat;
    std::auto_ptr<FittedPsf> fitpsf;
    double sigma_p = 0.;
    if (output_info) {
        std::cerr<<"Measuring PSF"<<std::endl;
    }
    DoMeasurePsf(
        params, static_cast<PsfLog&>(*log), im, weight_image.get(), trans,
        *starcat, psfcat, fitpsf, sigma_p);

    // Flag stars, so don't try to measure shears for them.
    bool nostars = params.read("cat_no_stars",false);
    if (!nostars) incat.flagStars(*starcat);

    // Do MeasusreShear script
    log.reset(new ShearLog(params,log_file,MakeName(params,"shear",false,false)));
    std::auto_ptr<ShearCatalog> shearcat;
    if (output_info) {
        std::cerr<<"Measuring Shear"<<std::endl;
    }
    DoMeasureShear(
        params, static_cast<ShearLog&>(*log), im, weight_image.get(), trans,
        incat, *fitpsf, shearcat);

    // Maybe do SplitStars script
    bool splitstars=params.read("splitstars",false);
    if (splitstars) {
        DoSplitStars(
            params, log_file, log, im, weight_image.get(), trans,
            incat, *starcat, sigma_p);
    }
}

int main(int argc, char **argv) try 
{
    ConfigFile params;
    if (BasicSetup(argc,argv,params,"fullpipe")) return EXIT_FAILURE;

    // Setup Log
    std::string log_file = ""; // Default is to stdout
    if (params.keyExists("log_file") || params.keyExists("log_ext")) 
        log_file = MakeName(params,"log",false,false);
    std::auto_ptr<Log> log;

    try {
        DoFullPipeline(params,log_file,log);
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

