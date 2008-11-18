
#include "DoMeasure.h"
#include "Pixel.h"
#include "Ellipse.h"
#include "Params.h"

double EstimateSigma(
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
    int flag = 0;
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
    ell.CrudeMeasure(pix,1.); // sigma here = 1.
    xdbg<<"Crude Measure: centroid = "<<ell.GetCen();
    xdbg<<", mu = "<<ell.GetMu()<<std::endl;
    meanmu += ell.GetMu();
    count++;
  }
  meanmu /= count;
  xdbg<<"meanmu = "<<meanmu<<std::endl;
  double sigma_p = exp(meanmu);
  xdbg<<"sigma_p = "<<sigma_p<<std::endl;
  return sigma_p;
}

void MeasureSinglePSF(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans,
    double noise, double gain, const Image<double>* weight_im,
    double sigma_p, double psfap, int psforder, PSFLog& log,
    BVec& psf, double& nu, int& flags)
{
  flags = 0;
  std::vector<std::vector<Pixel> > pix(1);
  int getpix_flag = 0;
  try {
    GetPixList(im,pix[0],cen,sky,noise,gain,weight_im,trans,psfap,getpix_flag);
  }
  catch (Range_error& e) {
    xdbg<<"skip: transformation range error: \n";
    xdbg<<"p = "<<cen<<", b = "<<e.b<<std::endl;
    log.nf_range++;
    flags |= MPSF_TRANSFORM_EXCEPTION;
    return;
  }
  if (getpix_flag) {
    xdbg<<"skip: getpix flag == "<<getpix_flag<<std::endl;
    if (flags & EDGE) {
      flags |= MPSF_EDGE1;
      log.nf_edge++;
    }
    if (flags & LT10PIX) {
      flags |= MPSF_LT10PIX1;
      log.nf_npix++;
    }
    return;
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

