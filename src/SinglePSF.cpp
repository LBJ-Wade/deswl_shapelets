
#include "DoMeasure.h"
#include "Pixel.h"
#include "Ellipse.h"

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
    double sigma_p, double psfap, int psforder,
    BVec*& psf, double& nu)
{
  std::vector<std::vector<Pixel> > pix(1);
  int flag = 0;
  try {
    GetPixList(im,pix[0],cen,sky,noise,gain,weight_im,trans,psfap,flag);
  }
  catch (Range_error& e) {
    xdbg<<"skip: transformation range error: \n";
    xdbg<<"p = "<<cen<<", b = "<<e.b<<std::endl;
    psf = 0;
    return;
  }
  if (flag) {
    xdbg<<"skip: flag == "<<flag<<std::endl;
    psf = 0;
    return;
  }
  xdbg<<"npix = "<<pix[0].size()<<std::endl;
  xxdbg<<"pixels: \n";
  for(size_t k=0;k<pix[0].size();k++) {
    xxdbg<<pix[0][k].z<<"  "<<pix[0][k].I<<"  "<<pix[0][k].wt<<std::endl;
  }

  Ellipse ell;
  ell.FixGam();
  ell.FixMu();
  ell.PeakCentroid(pix[0],psfap/3.);
  ell.CrudeMeasure(pix[0],1.); // sigma here = 1.
  psf = new BVec(psforder,sigma_p);
  tmv::Matrix<double> cov(psf->size(),psf->size());
  if (ell.Measure(pix,psforder,sigma_p,false,0,psf,&cov)) {
  //if (ell.Measure(pix,psforder,sigma_p,false,0,psf)) {
    xdbg<<"psf = "<<*psf<<std::endl;
    nu = (*psf)[0] / std::sqrt(cov(0,0));
    //nu = 1.;
    psf->Normalize();  // Divide by (0,0) element
    xdbg<<"Normalized psf: "<<*psf<<std::endl;
  }
  else {
    xdbg<<"Measurement failed\n";
    delete psf;
    psf = 0;
  }
}

