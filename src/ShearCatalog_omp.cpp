
#include <iostream>

#include "ShearCatalog.h"
#include "ConfigFile.h"
#include "dbg.h"
#include "Params.h"
#include "Image.h"
#include "FittedPsf.h"
#include "Log.h"
#include "MeasureShearAlgo.h"

//#define SINGLEGAL 8146
//#define STARTAT 8000
//#define ENDAT 200

int ShearCatalog::measureShears(
    const Image<double>& im,
    const Image<double>* weightIm, ShearLog& log)
{

// should we ramdomize the trajectory of the
// centroid search?
// see Params.h
#ifdef RANDOMIZE_CENTER
    static bool first = true;
    if (first) {
        // initialize random seed:
        unsigned int seed=0;
        for (size_t i=0;i<this->size(); i++) {
          seed += i*_flags[i];
        }
        //srand ( time(NULL) );
	//std::cout<<"using seed: "<<seed<<"\n";
        srand (seed);
        first = false;
    }
#endif


    int nGals = size();
    dbg<<"ngals = "<<nGals<<std::endl;

    // Read some needed parameters
    double galAperture = _params.read<double>("shear_aperture");
    double maxAperture = _params.read("shear_max_aperture",0.);
    int galOrder = _params.read<int>("shear_gal_order");
    int galOrder2 = _params.read("shear_gal_order2",galOrder);
    double fPsf = _params.read<double>("shear_f_psf");
    double gain = _params.read("image_gain",0.);
    double minGalSize = _params.read<double>("shear_min_gal_size");
    bool galFixCen = _params.read("shear_fix_centroid",false);
    double xOffset = _params.read("cat_x_offset",0.);
    double yOffset = _params.read("cat_y_offset",0.);
    bool shouldOutputDots = _params.read("output_dots",false);
    bool isTiming = _params.read("timing",false);

    // This need to have been set.
    Assert(_trans);
    Assert(_fitPsf);

    OverallFitTimes allTimes;

#ifdef ENDAT
    nGals = ENDAT;
#endif

    log._nGals = nGals;
#ifdef STARTAT
    log._nGals -= STARTAT;
#endif
#ifdef SINGLEGAL
    log._nGals = 1;
#endif
    log._nGoodIn = std::count(_flags.begin(),_flags.end(),0);
    dbg<<log._nGoodIn<<"/"<<log._nGals<<" galaxies with no input flags\n";
    std::vector<Position> initPos = _pos;

    // Main loop to measure shears
#ifdef _OPENMP
#pragma omp parallel 
    {
        try {
#endif
            OverallFitTimes times; // just for this thread
            ShearLog log1(_params); // just for this thread
            log1.noWriteLog();
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for(int i=0;i<nGals;++i) {
                if (_flags[i]) {
                    xdbg<<i<<" skipped because has flag "<<_flags[i]<<std::endl;
                    continue;
                }
#ifdef STARTAT
                if (i < STARTAT) continue;
#endif
#ifdef SINGLEGAL
                if (i < SINGLEGAL) continue;
                if (i > SINGLEGAL) break;
#endif
                if (shouldOutputDots) {
#ifdef _OPENMP
#pragma omp critical (output)
#endif
                    {
                        std::cerr<<"."; std::cerr.flush(); 
                    }
                }
                dbg<<"galaxy "<<i<<":\n";
                dbg<<"pos = "<<_pos[i]<<std::endl;

                // Start with an error code of unknown failure, in case
                // something happens that I don't detect as an error.
                _flags[i] = UNKNOWN_FAILURE;
                long flag1 = 0;
                measureSingleShear(
                    // Input data:
                    _pos[i], im, _sky[i], *_trans, *_fitPsf,
                    // Noise variables:
                    _noise[i], gain, weightIm, 
                    // Parameters:
                    galAperture, maxAperture, galOrder, galOrder2, 
                    fPsf, minGalSize, galFixCen, xOffset, yOffset,
                    // Time stats if desired:
                    isTiming ? &times : 0, 
                    // Log information
                    log1,
                    // Ouput values:
                    _shear[i], _cov[i], _shape[i], _nu[i], flag1);

                dbg<<"After measureSingleShear, pos = "<<_pos[i]<<std::endl;

                _flags[i] = flag1;

                if (!flag1) {
                    dbg<<"Successful shear measurement: "<<_shear[i]<<std::endl;
                } else {
                    dbg<<"Unsuccessful shear measurement\n"; 
                }

                if (isTiming) {
                    dbg<<"So far: ns = "<<times._nsGamma<<
                        ",  nf = "<<times._nfNative;
                    dbg<<", "<<times._nfMu<<", "<<times._nfGamma<<std::endl;
                }

            }
#ifdef _OPENMP
#pragma omp critical (add_log)
#endif
            {
                if (isTiming) allTimes += times;
                log += log1;
            }
#ifdef _OPENMP
        } catch (std::exception& e) {
            // This isn't supposed to happen.
            std::cerr<<"Caught "<<e.what()<<std::endl;
            std::cerr<<"STATUS5BEG Caught some error in parallel region STATUS5END\n";
            exit(1);
        } catch (...) {
            std::cerr<<"STATUS5BEG Caught some error in parallel region STATUS5END\n";
            exit(1);
        }
    }
#endif

    dbg<<log._nsGamma<<" successful shear measurements, ";
    dbg<<nGals-log._nsGamma<<" unsuccessful\n";
    log._nGood = std::count(_flags.begin(),_flags.end(),0);
    dbg<<log._nGood<<" with no flags\n";

    if (shouldOutputDots) {
        std::cerr
            <<std::endl
            <<"Success rate: "<<log._nsGamma<<"/"<<log._nGoodIn
            <<"  # with no flags: "<<log._nGood
            <<std::endl;
    }

    if (isTiming) {
        dbg<<"From timing structure:\n";
        dbg<<allTimes._nsGamma<<" successful shear measurements, ";
        dbg<<allTimes._nfNative<<" + "<<allTimes._nfMu;
        dbg<<" + "<<allTimes._nfGamma<<" unsuccessful\n";
        std::cerr<<allTimes<<std::endl;
    }

    // Check if the input positions are biased with respect to the 
    // actual values.  If they are, then this can cause a bias in the
    // shears, so it's worth fixing.
    std::complex<double> meanOffset(0);
    int nOffset=0;
    for(int i=0;i<nGals;++i) if (!_flags[i]) {
        meanOffset += _pos[i] - initPos[i];
        ++nOffset;
    }
    meanOffset /= nOffset;
    if (std::abs(meanOffset) > 0.5) {
        std::cerr<<"STATUS3BEG Warning: A bias in the input positions found: "
            <<meanOffset<<".\nThis may bias the shear, so you should "
            <<"change cat_x_offset, cat_y_offset.  STATUS3END\n";
    }
    dbg<<"Found mean offset of "<<meanOffset<<" using "<<nOffset<<
        " galaxies with no error codes.\n";

    xdbg<<log<<std::endl;

    return log._nsGamma;
}
