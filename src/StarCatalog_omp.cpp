
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
    const Image<double>*const weight_image, const Transformation& trans)
{
    Assert(int(_pos.size()) == size());
    Assert(int(_sky.size()) == size());
    Assert(int(_noise.size()) == size());
    Assert(int(_objsize.size()) == size());
    Assert(int(_flags.size()) == size());

    const int n = _pos.size();
    dbg<<"n = "<<n<<std::endl;

    bool use_shapelet_sigma = _params.read("stars_use_shapelet_sigma",true);
    bool output_dots = _params.read("output_dots",false);

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
        if (output_dots) {
#ifdef _OPENMP
#pragma omp critical (output)
#endif
            {
                std::cerr<<"."; std::cerr.flush();
            }
        }

        dbg<<"use i = "<<i<<std::endl;

        // Negative value indicates not set yet.  Start with 1 then.
        if (_objsize[i] <= 0.) _objsize[i] = 1.;
        // set nu to positive so it will calculate it
        _nu[i]=1;
        CalculateSigma(
            _objsize[i],_nu[i],
            im, _pos[i], _sky[i], _noise[i], weight_image, 
            trans, _params, _flags[i], use_shapelet_sigma);
    }
    if (output_dots) std::cerr<<std::endl;
    dbg<<"Done MeasureSigmas\n";
}

