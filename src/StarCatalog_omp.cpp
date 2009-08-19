
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
  std::ostream* dbgout_temp = dbgout;
  dbgout = 0;
#pragma omp parallel 
  {
#pragma omp for schedule(guided)
#endif
    for (int i=0; i<n; i++) if (!flags[i]) {
      ompdbg<<"use i = "<<i<<std::endl;

      try {
	// Negative value indicates not set yet.  Start with 1 then.
	if (objsize[i] <= 0.) objsize[i] = 1.;
	CalcSigma(
	    objsize[i],
	    im, pos[i], sky[i], noise[i], gain, weight_im, 
	    trans, psfap, flags[i]);
	ompdbg<<"objsize["<<i<<"]: "<<objsize[i]<<std::endl;
	ompdbg<<"flags["<<i<<"]: "<<flags[i]<<std::endl;
      } catch (tmv::Error& e) {
	ompdbg<<"Caught: "<<e<<std::endl;
	objsize[i] = DEFVALNEG;
	flags[i] |= TMV_EXCEPTION;
      } catch (std::exception& e) {
	ompdbg<<"Caught: "<<e.what()<<std::endl;
	objsize[i] = DEFVALNEG;
	flags[i] |= STD_EXCEPTION;
      } catch (...) {
	ompdbg<<"Caught unknown exception"<<std::endl;
	objsize[i] = DEFVALNEG;
	flags[i] |= UNKNOWN_EXCEPTION;
      }
    }
#ifdef _OPENMP
  } // End omp parallel 
  dbgout = dbgout_temp;
#endif
  dbg<<"Done MeasureSigmas\n";
}

