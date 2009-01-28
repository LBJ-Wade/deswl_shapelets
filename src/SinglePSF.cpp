
#include "DoMeasure.h"
#include "Pixel.h"
#include "Ellipse.h"
#include "Params.h"
#include <sstream>

#ifdef _OPENMP
#ifdef __PGI
#undef _OPENMP
#endif
#endif

double SingleSigma(
    const Image<double>& im, 
    const Position& pos,
    double sky,
    double noise,
    double gain,
    const Image<double>& weight_im, 
    const Transformation& trans, 
    double psfap,
    long& flag)
{
  double sigma = DEFVALNEG;

  std::vector<Pixel> pix;
  flag = 0;
  try {
    GetPixList(im, pix, pos, sky, noise, gain, &weight_im, trans, psfap, flag);
  } catch (Range_error& e) {
    xxdbg<<"transformation range error: \n";
    xxdbg<<"p = "<<pos<<", b = "<<e.b<<std::endl;
    flag = 1;
  }
  if (flag) return DEFVALNEG;
  xxdbg<<"npix = "<<pix.size()<<std::endl;

  Ellipse ell;
  ell.PeakCentroid(pix,psfap/3.);
  ell.CrudeMeasure(pix,1.); // sigma here = 1.
  xdbg<<"Crude Measure: centroid = "<<ell.GetCen();
  xdbg<<", mu = "<<ell.GetMu()<<std::endl;
  double mu = ell.GetMu();

  sigma = exp(mu);
  return sigma;
}

// estimate many sizes
void MeasureSigmas(
    const Image<double>& im,
    const std::vector<Position>& all_pos,
    const std::vector<double>& all_sky,
    const std::vector<double>& all_noise,
    double gain,
    const Image<double>& weight_im, 
    const Transformation& trans, 
    double psfap,
    std::vector<double>& sigmas,
    std::vector<long>& flags)
{
  sigmas.clear();
  sigmas.resize(all_pos.size(),0);
  flags.clear();
  flags.resize(all_pos.size(),0);

  int n = all_pos.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i=0; i<n; i++) {

    try {
      long flag1 = 0;
      double sigma1 = SingleSigma(
	  im, all_pos[i], all_sky[i], all_noise[i], gain, weight_im, 
	  trans, psfap, flag1);
#ifndef _OPENMP
      dbg<<"sigmas["<<i<<"]: "<<sigma1<<std::endl;
      dbg<<"flags["<<i<<"]: "<<flag1<<std::endl;
#endif
      sigmas[i] = sigma1;
      flags[i] = flag1;
#ifndef _OPENMP
    } catch (std::exception& e) {
      dbg<<"Caught: "<<e.what()<<std::endl;
      sigmas[i] = DEFVALNEG;
      flags[i] = 1;
#endif
    } catch (...) {
#ifndef _OPENMP
      dbg<<"Caught unknown exception"<<std::endl;
#endif
      sigmas[i] = DEFVALNEG;
      flags[i] = 1;
    }

  } // End omp parallel for
}

void EstimateSigma(
    double& sigma_p,
    const Image<double>& im,
    const std::vector<Position>& all_pos, const std::vector<double>& all_sky,
    const std::vector<double>& all_noise, double gain, 
    const Image<double>* weight_im, const Transformation& trans, double psfap)
{
  double meanmu=0.;
  int count=0;
  int nstars = all_pos.size();
  for(int i=0;i<nstars;i++) {
    xxdbg<<"star "<<i<<":\n";

    std::vector<Pixel> pix;
    long flag = 0;
    try {
      GetPixList(im,pix,all_pos[i],all_sky[i],all_noise[i],gain,
	  weight_im,trans,psfap,flag);
    } catch (Range_error& e) {
      xxdbg<<"skip: transformation range error: \n";
      xxdbg<<"p = "<<all_pos[i]<<", b = "<<e.b<<std::endl;
      continue;
    }
    if (flag) {
      xxdbg<<"skip: flag == "<<flag<<std::endl;
      continue;    
    }
    xxdbg<<"npix = "<<pix.size()<<std::endl;

    Ellipse ell;
    ell.PeakCentroid(pix,psfap/3.);
    ell.CrudeMeasure(pix,sigma_p); // sigma here = 1.
    xdbg<<"Crude Measure: centroid = "<<ell.GetCen();
    xdbg<<", mu = "<<ell.GetMu()<<std::endl;
    if (std::abs(ell.GetCen()) > 2.) {
      dbg<<"Large centroid excursion from CrudeMeasure for star "<<i<<"\n";
      dbg<<"Skipping this one for sigma_p calculation.\n";
      continue;
    }
    if (std::abs(ell.GetMu()) > 2.) {
      dbg<<"Large mu excursion from CrudeMeasure for star "<<i<<"\n";
      dbg<<"Skipping this one for sigma_p calculation.\n";
      continue;
    }
    meanmu += ell.GetMu();
    count++;
  }
  if (count < nstars/3) {
    std::ostringstream msgout;
    msgout<<"Too many objects were rejected. \n";
    msgout<<"nstars = "<<nstars<<", but only "<<count<<" successful measurements.\n";
    std::string msg = msgout.str();
    dbg<<msg<<std::endl;
    throw std::runtime_error(msg);
  }
  meanmu /= count;
  xdbg<<"meanmu = "<<meanmu<<std::endl;
  xdbg<<"input sigma_p = "<<sigma_p<<std::endl;
  sigma_p *= exp(meanmu);
  xdbg<<"sigma_p *= exp(mu) = "<<sigma_p<<std::endl;
}

void MeasureSinglePSF(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans,
    double noise, double gain, const Image<double>* weight_im,
    double sigma_p, double psfap, int psforder, PSFLog& log,
    BVec& psf, double& nu, long& flags)
{
  std::vector<std::vector<Pixel> > pix(1);
  long getpix_flag = 0;
  GetPixList(im,pix[0],cen,sky,noise,gain,weight_im,trans,psfap,getpix_flag);
  if (getpix_flag) {
    xdbg<<"skip: getpix flag == "<<getpix_flag<<std::endl;
    if (getpix_flag & EDGE) {
      flags |= MPSF_EDGE1;
      log.nf_edge++;
    }
    if (getpix_flag & LT10PIX) {
      flags |= MPSF_LT10PIX1;
      log.nf_npix++;
    }
  }
  int npix = pix[0].size();
  xdbg<<"npix = "<<npix<<std::endl;

  Ellipse ell;
  ell.FixGam();
  ell.FixMu();
  ell.PeakCentroid(pix[0],psfap/3.);
  ell.CrudeMeasure(pix[0],1.); // sigma here = 1.
  tmv::Matrix<double> cov(psf.size(),psf.size());
  if (ell.Measure(pix,psforder,sigma_p,false,0,&psf,&cov)) {
    xdbg<<"psf = "<<psf<<std::endl;
    nu = psf[0] / std::sqrt(cov(0,0));
    //nu = 1.;
    psf.Normalize();  // Divide by (0,0) element
    xdbg<<"Normalized psf: "<<psf<<std::endl;
    log.ns_psf++;
  }
  else {
    xdbg<<"Measurement failed\n";
    log.nf_psf++;
    nu = psf[0] / std::sqrt(cov(0,0));
    psf.Normalize();  // Divide by (0,0) element
    flags |= MPSF_MEASURE_FAILED;
  }
}

void MeasureSinglePSF1(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans,
    double noise, double gain, const Image<double>* weight_im,
    double sigma_p, double psfap, int psforder, PSFLog& log,
    BVec& psf, double& nu, long& flags)
{
  try {
    // We don't need to save skypos.  We just want to catch the range
    // error here, so we don't need to worry about it for dudx, etc.
    Position skypos;
    trans.Transform(cen,skypos);
    dbg<<"skypos = "<<skypos<<std::endl;
  } catch (Range_error& e) {
    xdbg<<"skip: transformation range error: \n";
    xdbg<<"p = "<<cen<<", b = "<<e.b<<std::endl;
    log.nf_range++;
    flags |= MPSF_TRANSFORM_EXCEPTION;
    return;
  }


  try {
    MeasureSinglePSF(
	// Input data:
	cen, im, sky, trans, 
	// Noise values:
	noise, gain, weight_im,
	// Parameters:
	sigma_p, psfap, psforder,
	// Log information
	log,
	// Ouput value:
	psf, nu, flags);
  } catch (tmv::Error& e) {
    dbg<<"TMV Error thrown in MeasureSinglePSF\n";
    dbg<<e<<std::endl;
    log.nf_tmverror++;
    flags |= MPSF_TMV_EXCEPTION;
  } catch (...) {
    dbg<<"unkown exception in MeasureSinglePSF\n";
    log.nf_othererror++;
    flags |= MPSF_UNKNOWN_EXCEPTION;
  }
}

