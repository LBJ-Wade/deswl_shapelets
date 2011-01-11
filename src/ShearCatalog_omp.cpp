
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

int ShearCatalog::measureShears(
    const Image<double>& im,
    const Image<double>* weightIm, ShearLog& log)
{
    static bool first = true;
#ifdef _OPENMP
#pragma omp critical (srand)
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

                BVec shearedShape(galOrder,1.);
                _flags[i] = 0;
                measureShapes(
                    // Input data:
                    _pos[i], im, _sky[i], *_trans, *_fitPsf,
                    // Noise variables:
                    _noise[i], gain, weightIm, 
                    // Parameters:
                    galAperture, maxAperture, galOrder, galOrder2,
                    minFPsf, maxFPsf, minGalSize, galFixCen, xOffset, yOffset,
                    galFixSigma, galFixSigmaValue,
                    // Log information
                    log1,
                    // Ouput values:
                    _shape[i], _shear[i], shearedShape, _nu[i],
                    _flags[i]);
                dbg<<"After measureShapes, pos = "<<_pos[i]<<std::endl;
                if (!_flags[i]) {
                    dbg<<"Successful shape measurements: \n";
                    dbg<<_shape[i]<<std::endl;
                    dbg<<shearedShape<<std::endl;
                } else {
                    dbg<<"Unsuccessful shape measurement\n"; 
                    dbg<<"flag = "<<_flags[i]<<std::endl;
                }

                // At this point, the shear value is just the shear of the
                // _observed_ galaxy.  We still need to convert this to
                // the shear of the frame in which the deconvolved
                // galaxy is round, using the observed shape in the 
                // sheared frame: shearedShape.
                Ellipse ell;
                ell.setGamma(_shear[i]);
                if (galFixCen) ell.fixCen();
                if (galFixSigma) ell.fixMu();
                DMatrix cov5(5,5);
                if (ell.findRoundFrame(
                        shearedShape, true, galOrder2, 1.e-4, 
                        _flags[i], &cov5)) {
                    ell.removeRotation();
                    _shear[i] = ell.getGamma();
                    dbg<<"Successful shear measurement: "<<_shear[i]<<std::endl;

                    DSmallMatrix22 cov2 = cov5.TMV_subMatrix(2,4,2,4);
                    if (!(cov2.TMV_det() > 0.)) {
                        dbg<<"cov2 has bad determinant: "<<
                            cov2.TMV_det()<<std::endl;
                        dbg<<"cov2 = "<<cov2<<std::endl;
                        dbg<<"Full cov = "<<cov5<<std::endl;
                        xdbg<<"flag SHEAR_BAD_COVAR\n";
                        _flags[i] |= SHEAR_BAD_COVAR;
                    } else {
                        _cov[i] = cov2;
                    }
                } else {
                    ell.removeRotation();
                    _shear[i] = ell.getGamma();
                    dbg<<"Unsuccessful shear measurement\n"; 
                    dbg<<"flag = "<<_flags[i]<<std::endl;
                    dbg<<"shear = "<<_shear[i]<<std::endl;
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

    if (shouldOutputDots && !shouldOutputDesQa) {
        std::complex<double> meanShear = 0.;
        int nGoodShear = 0.;
        for(int i=0;i<nGals;++i) if (!_flags[i]) {
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
