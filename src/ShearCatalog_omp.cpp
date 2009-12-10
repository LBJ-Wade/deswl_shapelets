
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
    const Image<double>* weightIm, const Transformation& trans,
    const FittedPsf& fitPsf, ShearLog& log)
{
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
    bool shouldOutputDots = _params.read("output_dots",false);
    bool isTiming = _params.read("timing",false);

    OverallFitTimes allTimes;

#ifdef ENDAT
    nGals = ENDAT;
#endif

    log._nGals = nGals;
#ifdef STARTAT
    log._nGals -= STARTAT;
#endif
#ifdef SINGLEGAL
    log.nGals = 1;
#endif
    log._nGoodIn = std::count(_flags.begin(),_flags.end(),0);
    dbg<<log._nGoodIn<<"/"<<log._nGals<<" galaxies with no input flags\n";

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
                if (_flags[i]) continue;
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
                dbg<<"pos[i] = "<<_pos[i]<<std::endl;

                // Start with an error code of unknown failure, in case
                // something happens that I don't detect as an error.
                _flags[i] = UNKNOWN_FAILURE;
                long flag1 = 0;
                measureSingleShear(
                    // Input data:
                    _pos[i], im, _sky[i], trans, 
                    // Fitted PSF
                    fitPsf,
                    // Noise variables:
                    _noise[i], gain, weightIm, 
                    // Parameters:
                    galAperture, maxAperture, galOrder, galOrder2, 
                    fPsf, minGalSize, 
                    // Time stats if desired:
                    isTiming ? &times : 0, 
                    // Log information
                    log1,
                    // Ouput values:
                    _shear[i], _cov[i], _shape[i], _nu[i], flag1);

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
        } catch (...) {
            // This isn't supposed to happen.
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
    xdbg<<log<<std::endl;

    return log._nsGamma;
}
