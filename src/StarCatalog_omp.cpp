
#include "StarCatalog.h"
#include "StarFinder.h"
#include "Bounds.h"
#include "ConfigFile.h"
#include "Image.h"
#include "Transformation.h"
#include "Params.h"
#include "Log.h"

//#define SINGLEGAL 111
//#define STARTAT 8000
//#define ENDAT 14

#ifdef SINGLEGAL
#undef _OPENMP
#endif

void StarCatalog::calculateSizes(
    const Image<double>& im, 
    const Image<double>*const weightIm, const Transformation& trans)
{
    Assert(int(_pos.size()) == size());
    Assert(int(_sky.size()) == size());
    Assert(int(_noise.size()) == size());
    Assert(int(_objSize.size()) == size());
    Assert(int(_flags.size()) == size());

    const int n = _pos.size();
    dbg<<"n = "<<n<<std::endl;
    double psfAp = _params.read<double>("psf_aperture"); 
    double gain = _params.read("image_gain",0.);
    double xOffset = _params.read("cat_x_offset",0.);
    double yOffset = _params.read("cat_y_offset",0.);

    const bool shouldUseShapeletSigma = _params.read(
        "stars_use_shapelet_sigma",true);
    const bool shouldOutputDots = _params.read("output_dots",false);

#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (int i=0; i<n; ++i) if (!_flags[i]) {
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

        dbg<<"use i = "<<i<<std::endl;

        // Negative value indicates not set yet.  Start with 1 then.
        if (_objSize[i] <= 0.) _objSize[i] = 1.;
        calculateSigma(
            _objSize[i],
            im, _pos[i], _sky[i], _noise[i], gain, weightIm, 
            trans, psfAp, xOffset, yOffset, 
            _flags[i], shouldUseShapeletSigma);
    }
    if (shouldOutputDots) std::cerr<<std::endl;
    dbg<<"Done MeasureSigmas\n";
}

