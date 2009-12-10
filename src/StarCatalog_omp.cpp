
#include "StarCatalog.h"
#include "StarFinder.h"
#include "Bounds.h"
#include "ConfigFile.h"
#include "Image.h"
#include "Transformation.h"
#include "Params.h"
#include "Log.h"

void StarCatalog::calculateSizes(
    const Image<double>& im, 
    const Image<double>*const weightIm, const Transformation& trans)
{
    Assert(_pos.size() == size());
    Assert(_sky.size() == size());
    Assert(_noise.size() == size());
    Assert(_objSize.size() == size());
    Assert(_flags.size() == size());

    const int n = _pos.size();
    dbg<<"n = "<<n<<std::endl;
    double psfAp = _params.read<double>("psf_aperture"); 
    double gain = _params.read("image_gain",0.);

    const bool shouldUseShapeletSigma = _params["stars_use_shapelet_sigma"];

#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (int i=0; i<n; ++i) if (!_flags[i]) {
        dbg<<"use i = "<<i<<std::endl;

        // Negative value indicates not set yet.  Start with 1 then.
        if (_objSize[i] <= 0.) _objSize[i] = 1.;
        calculateSigma(
            _objSize[i],
            im, _pos[i], _sky[i], _noise[i], gain, weightIm, 
            trans, psfAp, _flags[i], shouldUseShapeletSigma);
    }
    dbg<<"Done MeasureSigmas\n";
}

