
#include <sstream>

#include "PsfCatalog.h"
#include "Params.h"

//#define SINGLESTAR 18
//#define NSTARS 10
//#define STARTAT 85
//#define ENDAT 95

double PsfCatalog::estimateSigma(
    const Image<double>& im,
    const Image<double>* weight_image, const Transformation& trans)
{
    if (_params.keyExists("psf_force_sigma_p")) {
        double sigma_p = _params["psf_force_sigma_p"];
        dbg<<"Using forced value for sigma_p = "<<sigma_p<<std::endl;
        return sigma_p;
    }

    // Initial sigma_p for shapelet measurements
    double sigma_p = 1.;
    if (_params.keyExists("psf_seeing_est")) {
        double seeing = _params["psf_seeing_est"];
        // seeing is given as FWHM
        // for a gaussian 0.5 = exp(-((FWHM/2)/sigma)^2/2)
        // FWHM/sigma = 2.*sqrt(2 ln(2)) = 2.35
        sigma_p = seeing / 2.35;
    }

    // Calculate a good value of sigma to use:
    const bool use_shapelet_sigma = true;

    int nstars = _pos.size();
    double meanmu = 0.;
    int count = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(guided) reduction(+ : meanmu) reduction(+ : count)
#endif
    for (int i=0; i<nstars; ++i) if (!_flags[i]) {
        dbg<<"use i = "<<i<<std::endl;

        double sigma = sigma_p;
        long flag1 = 0; // Ignore flags set by CalcSigma
	double nu=-1; // use negative value so it will not calculate nu
        CalculateSigma(
	    sigma, nu,
            im, _pos[i], _sky[i], _noise[i], weight_image, 
            trans, _params, flag1, use_shapelet_sigma);
        // Ignore errors -- just don't add to meanmu
        if (flag1) continue;
        meanmu += log(sigma);
        ++count;
    } // End omp parallel for

    if (count < nstars/3) {
        std::ostringstream msgout;
        msgout<<"Too many objects were rejected. \n";
        msgout<<"nstars = "<<nstars<<", but only "<<
            count<<" successful measurements.\n";
        std::string msg = msgout.str();
        dbg<<msg<<std::endl;
        throw ProcessingException(msg);
    }
    meanmu /= count;
    xdbg<<"meanmu = "<<meanmu<<std::endl;
    xdbg<<"input sigma_p = "<<sigma_p<<std::endl;
    sigma_p = exp(meanmu);
    dbg<<"sigma_p = "<<sigma_p<<std::endl;
    return sigma_p;
}

int PsfCatalog::measurePsf(
    const Image<double>& im,
    const Image<double>* weight_image,
    const Transformation& trans, double sigma_p, PsfLog& log)
{
    // Read some needed parameters
    bool output_dots = _params.read("output_dots",false);
    double gain = _params.read("image_gain",0.);
    double psfAp = _params.read<double>("psf_aperture");
    bool psfFixCen = _params.read("psf_fix_centroid",false);
    double xOffset = _params.read("cat_x_offset",0.);
    double yOffset = _params.read("cat_y_offset",0.);
    dbg<<"psfap = "<<psfAp<<std::endl;

    int nstars = size();
    dbg<<"nstars = "<<nstars<<std::endl;

#ifdef ENDAT
    nstars = ENDAT;
#endif

    log._nstars = nstars;
#ifdef STARTAT
    log._nstars -= STARTAT;
#endif
#ifdef SINGLESTAR
    log._nstars = 1;
#endif
    log._ngoodin = std::count(_flags.begin(),_flags.end(),0);
    dbg<<log._ngoodin<<"/"<<log._nstars<<" stars with no input flags\n";

    Assert(nstars<=int(_pos.size()));
    Assert(nstars<=int(_sky.size()));
    Assert(nstars<=int(_noise.size()));
    Assert(nstars<=int(_psf.size()));
    Assert(nstars<=int(_nu.size()));
    Assert(nstars<=int(_flags.size()));
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
            for(int i=0;i<nstars;++i) if (!_flags[i]) {
#ifdef STARTAT
                if (i < STARTAT) continue;
#endif
#ifdef SINGLESTAR
                if (i < SINGLESTAR) continue;
                if (i > SINGLESTAR) break;
                XDEBUG = true;
#endif
                if (output_dots) {
#ifdef _OPENMP
#pragma omp critical (output)
#endif
                    {
                        std::cerr<<"."; std::cerr.flush(); 
                    }
                }
                dbg<<"star "<<i<<":\n";
                dbg<<"pos["<<i<<"] = "<<_pos[i]<<std::endl;

                dbg<<"Before MeasureSinglePSF1"<<std::endl;
                MeasureSinglePsf(
                    // Input data:
                    _pos[i], im, _sky[i], trans, 
                    // Noise values:
                    _noise[i], weight_image,
                    // Parameters:
                    sigma_p, _params,
                    // Log information
                    log1,
                    // Ouput value:
                    _psf[i], _nu[i], _flags[i]);
                dbg<<"After MeasureSinglePSF"<<std::endl;

                if (!_flags[i]) {
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

    dbg<<log._ns_psf<<" successful star measurements, ";
    dbg<<nstars-log._ns_psf<<" unsuccessful\n";
    log._ngood = std::count(_flags.begin(),_flags.end(),0);
    dbg<<log._ngood<<" with no flags\n";
    dbg<<"Breakdown of flags:\n";
    if (dbgout) PrintFlags(_flags,*dbgout);

    if (output_dots) {
        std::cerr
            <<std::endl
            <<"Success rate: "<<log._ns_psf<<"/"<<log._ngoodin
            <<"  # with no flags: "<<log._ngood
            <<std::endl;
        std::cerr<<"Breakdown of flags:\n";
        PrintFlags(_flags,std::cerr);
    }

    return log._ns_psf;
}

