
#include <sstream>

#include "PsfCatalog.h"
#include "dbg.h"
#include "Params.h"

//#define SINGLESTAR 18
//#define NSTARS 10
//#define STARTAT 85
//#define ENDAT 95

double PsfCatalog::estimateSigma(
    const Image<double>& im,
    const Image<double>* weightIm, const Transformation& trans)
{
    // Initial sigmaP for shapelet measurements
    double sigmaP = 1.;
    if (_params.keyExists("psf_seeing_est")) {
        double seeing = _params["psf_seeing_est"];
        // seeing is given as FWHM
        // for a gaussian 0.5 = exp(-((FWHM/2)/sigma)^2/2)
        // FWHM/sigma = 2.*sqrt(2 ln(2)) = 2.35
        sigmaP = seeing / 2.35;
    }

    // Calculate a good value of sigma to use:
    // (For this calculation, psfAp is psf_aperture * 1 arcsec.)
    double gain = _params.read("gain",0.);
    double psfAp = _params.read<double>("psf_aperture");
    dbg<<"psfap = "<<psfAp<<std::endl;

    const bool shouldUseShapeletSigma = true;

    int nStars = _pos.size();
    double meanMu = 0.;
    int count = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(guided) reduction(+ : meanMu) reduction(+ : count)
#endif
    for (int i=0; i<nStars; ++i) if (!_flags[i]) {
        dbg<<"use i = "<<i<<std::endl;

        double sigma = sigmaP;
        long flag1 = 0; // Ignore flags set by CalcSigma
        calculateSigma(
            sigma,
            im, _pos[i], _sky[i], _noise[i], gain, weightIm, 
            trans, psfAp, flag1, shouldUseShapeletSigma);
        // Ignore errors -- just don't add to meanMu
        if (flag1) continue;
        meanMu += log(sigma);
        ++count;
    } // End omp parallel for

    if (count < nStars/3) {
        std::ostringstream msgOut;
        msgOut<<"Too many objects were rejected. \n";
        msgOut<<"nstars = "<<nStars<<", but only "<<
            count<<" successful measurements.\n";
        std::string msg = msgOut.str();
        dbg<<msg<<std::endl;
        throw ProcessingException(msg);
    }
    meanMu /= count;
    xdbg<<"meanmu = "<<meanMu<<std::endl;
    xdbg<<"input sigma_p = "<<sigmaP<<std::endl;
    sigmaP = exp(meanMu);
    dbg<<"sigma_p = "<<sigmaP<<std::endl;
    return sigmaP;
}

int PsfCatalog::measurePsf(
    const Image<double>& im,
    const Image<double>* weightIm,
    const Transformation& trans, double sigmaP, PsfLog& log)
{
    // Read some needed parameters
    int psfOrder = _params.read<int>("psf_order");
    bool shouldOutputDots = _params.read("output_dots",false);
    double gain = _params.read("image_gain",0.);
    double psfAp = _params.read<double>("psf_aperture");
    dbg<<"psfap = "<<psfAp<<std::endl;
    psfAp *= sigmaP;  // arcsec
    dbg<<"psfap => "<<psfAp<<std::endl;

    int nStars = size();
    dbg<<"nstars = "<<nStars<<std::endl;

#ifdef ENDAT
    nStars = ENDAT;
#endif

    log._nStars = nStars;
#ifdef STARTAT
    log._nStars -= STARTAT;
#endif
#ifdef SINGLESTAR
    log._nStars = 1;
#endif
    log._nGoodIn = std::count(_flags.begin(),_flags.end(),0);
    dbg<<log._nGoodIn<<"/"<<log._nStars<<" stars with no input flags\n";

    Assert(nStars<=int(_pos.size()));
    Assert(nStars<=int(_sky.size()));
    Assert(nStars<=int(_noise.size()));
    Assert(nStars<=int(_psf.size()));
    Assert(nStars<=int(_nu.size()));
    Assert(nStars<=int(_flags.size()));
    // Main loop to measure psf shapelets:
#ifdef _OPENMP
#pragma omp parallel 
    {
        try {
#endif
            PsfLog log1(_params);  // Just for this thread
            log1.noWriteLog();
#ifdef _OPENMP
#pragma omp for schedule(guided)
#endif
            for(int i=0;i<nStars;++i) if (!_flags[i]) {
#ifdef STARTAT
                if (i < STARTAT) continue;
#endif
#ifdef SINGLESTAR
                if (i < SINGLESTAR) continue;
                if (i > SINGLESTAR) break;
                XDEBUG = true;
#endif
                if (shouldOutputDots) {
#ifdef _OPENMP
#pragma omp critical (output)
#endif
                    {
                        std::cerr<<"."; std::cerr.flush(); 
                    }
                }
                dbg<<"star "<<i<<":\n";
                dbg<<"pos["<<i<<"] = "<<_pos[i]<<std::endl;

                // Start with an error code of unknown failure, in case
                // something happens that I don't detect as an error.
                _flags[i] = UNKNOWN_FAILURE;
                dbg<<"Before MeasureSinglePSF1"<<std::endl;
                long flag1 = 0;
                measureSinglePsf(
                    // Input data:
                    _pos[i], im, _sky[i], trans, 
                    // Noise values:
                    _noise[i], gain, weightIm,
                    // Parameters:
                    sigmaP, psfAp, psfOrder, 
                    // Log information
                    log1,
                    // Ouput value:
                    _psf[i], _nu[i], flag1);
                dbg<<"After MeasureSinglePSF"<<std::endl;

                _flags[i] = flag1;
                if (!flag1) {
                    dbg<<"Successful psf measurement: "<<
                        _psf[i].vec()<<std::endl;
                } else {
                    dbg<<"Unsuccessful psf measurement\n"; 
                }
            }
#ifdef _OPENMP
#pragma omp critical (add_log)
#endif
            {
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

    dbg<<log._nsPsf<<" successful star measurements, ";
    dbg<<nStars-log._nsPsf<<" unsuccessful\n";
    log._nGood = std::count(_flags.begin(),_flags.end(),0);
    dbg<<log._nGood<<" with no flags\n";

    if (shouldOutputDots) {
        std::cerr
            <<std::endl
            <<"Success rate: "<<log._nsPsf<<"/"<<log._nGoodIn
            <<"  # with no flags: "<<log._nGood
            <<std::endl;
    }

    return log._nsPsf;
}

