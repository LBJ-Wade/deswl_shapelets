
#include <iostream>

#include "ShearCatalog.h"
#include "ConfigFile.h"
#include "dbg.h"
#include "Params.h"
#include "Image.h"
#include "FittedPsf.h"
#include "Log.h"
#include "MeasureShearAlgo.h"
#include "Ellipse.h"

//#define SINGLEGAL 13
//#define STARTAT 8000
//#define ENDAT 14

#ifdef SINGLEGAL
#undef _OPENMP
#endif

static bool GetPixPsf(
    const Image<double>& im, PixelList& pix,
    const Position pos, double sky, double noise, double gain,
    const Image<double>* weightIm, const Transformation& trans,
    double maxAperture, double xOffset, double yOffset, long& flags,
    const FittedPsf& fitPsf, BVec& psf, ShearLog& log)
{
    // Get the main PixelList for this galaxy:
    try {
        getPixList(
            im,pix,pos,sky,noise,gain,weightIm,
            trans,maxAperture,xOffset,yOffset,flags);
        int npix = pix.size();
        dbg<<"npix = "<<npix<<std::endl;
        if (npix < 10) {
            dbg<<"Too few pixels to continue: "<<npix<<std::endl;
            dbg<<"FLAG SHAPELET_NOT_DECONV\n";
            flags |= SHAPELET_NOT_DECONV;
            return false;
        }
    } catch (RangeException& e) {
        dbg<<"distortion range error: \n";
        xdbg<<"p = "<<pos<<", b = "<<e.getBounds()<<std::endl;
        ++log._nfRange1;
        dbg<<"FLAG TRANSFORM_EXCEPTION\n";
        flags |= TRANSFORM_EXCEPTION;
        return false;
    }

    // Get the interpolated psf for this galaxy:
    try {
        psf = fitPsf(pos);
    } catch (RangeException& e) {
        dbg<<"fittedpsf range error: \n";
        xdbg<<"p = "<<pos<<", b = "<<e.getBounds()<<std::endl;
        ++log._nfRange2;
        dbg<<"FLAG FITTEDPSF_EXCEPTION\n";
        flags |= FITTEDPSF_EXCEPTION;
        return false;
    }
    return true;
}
 
int ShearCatalog::measureShears(
    const Image<double>& im,
    const Image<double>* weightIm, ShearLog& log)
{
    static bool first = true;
#ifdef _OPENMP
#ifndef __INTEL_COMPILER
    // This should technically be guarded by a critical block, since 
    // it involves writing to global variables.  However, icpc seg faults
    // on the srand call if it is in a critical block.
    // I have no idea why, but it is clearly a bug in icpc.  v10.1 at least.
    // Anyway, it shouldn't be a problem, since this block should not 
    // be in a parallel region anyway, so there shouldn't be any problem
    // with multi-threading here anyway.
#pragma omp critical (srand)
#endif
#endif
    {
        if (first) {
            // initialize random seed:
            // NB: This will only stay deterministic if not using openmp.
            unsigned int seed=0;
            for (size_t i=0;i<this->size(); i++) {
                seed += i*_flags[i];
            }
            dbg<<"using seed: "<<seed<<"\n";
            srand (seed);
            first = false;
        }
    }

    int nGals = size();
    dbg<<"ngals = "<<nGals<<std::endl;

    // Read some needed parameters
    double galAperture = _params.read<double>("shear_aperture");
    double maxAperture = _params.read("shear_max_aperture",0.);
    int galOrder = _params.read<int>("shear_gal_order");
    int galOrder2 = _params.read<int>("shear_gal_order2");
    int maxm = _params.read("shear_maxm",galOrder);
    double minFPsf = _params.read("shear_f_psf",1.);
    double maxFPsf = _params.read("shear_max_f_psf",minFPsf);
    double gain = _params.read("image_gain",0.);
    double minGalSize = _params.read<double>("shear_min_gal_size");
    bool galFixCen = _params.read("shear_fix_centroid",false);
    bool galFixSigma = _params.keyExists("shear_force_sigma");
    double galFixSigmaValue = _params.read("shear_force_sigma",0.);
    double xOffset = _params.read("cat_x_offset",0.);
    double yOffset = _params.read("cat_y_offset",0.);
    bool shouldOutputDots = _params.read("output_dots",false);
    bool shouldOutputDesQa = _params.read("des_qa",false); 
    bool nativeOnly = _params.read("shear_native_only",false);

    // This need to have been set.
    Assert(_trans);
    Assert(_fitPsf);

#ifdef ENDAT
    for(int i=ENDAT; i<nGals; ++i) _flags[i] |= INPUT_FLAG;
    nGals = ENDAT;
#endif

    log._nGals = nGals;
#ifdef STARTAT
    for(int i=0; i<STARTAT; ++i) _flags[i] |= INPUT_FLAG;
    log._nGals -= STARTAT;
#endif
#ifdef SINGLEGAL
    for(int i=0; i<nGals; ++i) if (i != SINGLEGAL)  _flags[i] |= INPUT_FLAG;
    log._nGals = 1;
#endif
    log._nGoodIn = std::count(_flags.begin(),_flags.end(),0);
    dbg<<log._nGoodIn<<"/"<<log._nGals<<" galaxies with no input flags\n";
    std::vector<Position> initPos = _pos;

    // Main loop to measure shapes
#ifdef _OPENMP
#pragma omp parallel 
    {
        try {
#endif
            // Some variables to use just for this thread
            ShearLog log1(_params); 
            log1.noWriteLog();
            std::vector<PixelList> pix(1);
            std::vector<BVec> psf(
                1, BVec(_fitPsf->getPsfOrder(), _fitPsf->getSigma()));
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

                if (!GetPixPsf(
                        im,pix[0],_pos[i],_sky[i],_noise[i],gain,weightIm,
                        *_trans,maxAperture,xOffset,yOffset,_flags[i],
                        *_fitPsf,psf[0],log1)) {
                    continue;
                }

                // Now measure the shape and shear:
                DoMeasureShear(
                    // Input data:
                    pix, psf,
                    // Parameters:
                    galAperture, maxAperture, galOrder, galOrder2, maxm,
                    minFPsf, maxFPsf, minGalSize, galFixCen,
                    galFixSigma, galFixSigmaValue, nativeOnly,
                    // Log information
                    log1,
                    // Ouput values:
                    _shape[i], _shear[i], _cov[i], _nu[i], _flags[i]);

                if (!_flags[i]) {
                    dbg<<"Successful shear measurements: \n";
                    dbg<<"shape = "<<_shape[i]<<std::endl;
                    dbg<<"shear = "<<_shear[i]<<std::endl;
                    dbg<<"cov = "<<_cov[i]<<std::endl;
                    dbg<<"nu = "<<_nu[i]<<std::endl;
                } else {
                    dbg<<"Unsuccessful shear measurement\n"; 
                    dbg<<"flag = "<<_flags[i]<<std::endl;
                }
            }
            dbg<<"After loop"<<std::endl;
#ifdef _OPENMP
#pragma omp critical (add_log)
#endif
            {
                log += log1;
            }
#ifdef _OPENMP
        } catch (std::exception& e) {
            // This isn't supposed to happen.
            if (shouldOutputDesQa) {
                std::cerr<<"STATUS5BEG Caught error in parallel region STATUS5END\n";
            } 
            std::cerr<<"Caught "<<e.what()<<std::endl;
            std::cerr<<"Caught error in parallel region.  Aborting.\n";
            exit(1);
        } catch (...) {
            if (shouldOutputDesQa) {
                std::cerr<<"STATUS5BEG Caught error in parallel region STATUS5END\n";
            }
            std::cerr<<"Caught error in parallel region.  Aborting.\n";
            exit(1);
        }
    }
#endif
    dbg<<log._nsGamma<<" successful shape measurements, ";
    dbg<<nGals-log._nsGamma<<" unsuccessful\n";
    log._nGood = std::count(_flags.begin(),_flags.end(),0);
    dbg<<log._nGood<<" with no flags\n";
    dbg<<"Breakdown of flags:\n";
    if (dbgout) PrintFlags(_flags,*dbgout);

    if (shouldOutputDots) {
        std::cerr
            <<std::endl
            <<"Success rate: "<<log._nsGamma<<"/"<<log._nGoodIn
            <<"  # with no flags: "<<log._nGood
            <<std::endl;
        std::cerr<<"Breakdown of flags:\n";
        PrintFlags(_flags,std::cerr);
    }


    static const long ok_flags = (
        EDGE |
        SHEAR_REDUCED_ORDER |
        SHAPE_REDUCED_ORDER |
        SHAPE_POOR_FIT |
        SHAPE_LOCAL_MIN |
        SHAPE_BAD_FLUX );

    if (shouldOutputDots && !shouldOutputDesQa) {
        std::complex<double> meanShear = 0.;
        int nGoodShear = 0.;
        for(int i=0;i<nGals;++i) if (!(_flags[i] & ~ok_flags)) {
            meanShear += _shear[i];
            ++nGoodShear;
        }

        meanShear /= nGoodShear;
        dbg<<"nGoodShear = "<<nGoodShear<<"   meanShear = "<<meanShear<<std::endl;
        std::cerr
            <<"nGoodShear = "<<nGoodShear
            <<"   meanShear = "<<meanShear<<std::endl;
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
        if (shouldOutputDesQa) {
            std::cerr<<"STATUS3BEG Warning: A bias in the input positions found: "
                <<meanOffset<<".\nThis may bias the shear, so you should "
                <<"change cat_x_offset, cat_y_offset.  STATUS3END\n";
        }
    }
    dbg<<"Found mean offset of "<<meanOffset<<" using "<<nOffset<<
        " galaxies with no error codes.\n";

    xdbg<<log<<std::endl;

    return log._nsGamma;
}
