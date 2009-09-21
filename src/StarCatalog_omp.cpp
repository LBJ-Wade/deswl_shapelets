
#include "StarCatalog.h"
#include "StarFinder.h"
#include "Bounds.h"
#include "ConfigFile.h"
#include "Image.h"
#include "Transformation.h"
#include "Params.h"
#include "Log.h"

void StarCatalog::CalcSizes(const Image<double>& im, 
    const Image<double>*const weight_im, const Transformation& trans)
{
  Assert(pos.size() == size());
  Assert(sky.size() == size());
  Assert(noise.size() == size());
  Assert(objsize.size() == size());
  Assert(flags.size() == size());

  const int n = pos.size();
  dbg<<"n = "<<n<<std::endl;
  double psfap = params.read<double>("psf_aperture"); 
  double gain = params.read("image_gain",0.);

#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
  for (int i=0; i<n; i++) if (!flags[i]) {
    dbg<<"use i = "<<i<<std::endl;

    // Negative value indicates not set yet.  Start with 1 then.
    if (objsize[i] <= 0.) objsize[i] = 1.;
    CalcSigma(
	objsize[i],
	im, pos[i], sky[i], noise[i], gain, weight_im, 
	trans, psfap, flags[i]);
  }
  dbg<<"Done MeasureSigmas\n";
}

