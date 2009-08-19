
#include "Ellipse.h"
#include "EllipseSolver.h"
#include "TMV.h"
#include <cmath>
#include "dbg.h"
#include <fstream>
#include "TimeVars.h"
#include "PsiHelper.h"
#include "Params.h"

#define N_FLUX_ATTEMPTS 0
#define MAXITER 4

bool Ellipse::Measure(const std::vector<PixelList>& pix,
    const std::vector<BVec>& psf,
    int order, double sigma, bool use_integ, long& flag, bool desqa,
    tmv::Matrix<double>* cov, BVec* bret, tmv::Matrix<double>* bcov)
{ 
  return DoMeasure(pix,&psf,order,sigma,use_integ,flag,desqa,cov,bret,bcov); 
}

bool Ellipse::Measure(const std::vector<PixelList>& pix,
    int order, double sigma, bool use_integ, long& flag, bool desqa,
    tmv::Matrix<double>* cov, BVec* bret, tmv::Matrix<double>* bcov)
{
  return DoMeasure(pix,0,order,sigma,use_integ,flag,desqa,cov,bret,bcov);
}

void Ellipse::DoMeasureShapelet(const std::vector<PixelList>& pix,
    const std::vector<BVec>* psf, BVec& b, tmv::Matrix<double>* bcov) const
{
  xdbg<<"Start MeasureShapelet: order = "<<b.GetOrder()<<std::endl;
  xdbg<<"sigma = "<<b.GetSigma()<<std::endl;
  xdbg<<"el = "<<*this<<std::endl;
  
  // ( u )' = exp(-mu)/sqrt(1-gsq) ( 1-g1  -g2  ) ( u-uc )
  // ( v )                         ( -g2   1+g1 ) ( v-vc )
  // 
  // z' = u' + I v' 
  //    = exp(-mu)/sqrt(1-gsq)
  //             [ (1-g1)(u-uc)-g2(v-vc)-Ig2(u-uc)+I(1+g1)(v-vc) ]
  //    = exp(-mu)/sqrt(1-gsq) [ z-zc - g1(z-zc)* -Ig2(z-zc)* ]
  //    = exp(-mu)/sqrt(1-gsq) ( z-zc - g (z-zc)* )

  int order = b.GetOrder();
  double sigma = b.GetSigma();
  double gsq = std::norm(gamma);
  Assert(gsq < 1.);

  double m = exp(-mu)/sqrt(1.-gsq);

  size_t ntot = 0;
  for(size_t i=0;i<pix.size();i++) ntot += pix[i].size();
  dbg<<"ntot = "<<ntot<<" in "<<pix.size()<<" images\n";

  tmv::Vector<double> I(ntot);
  tmv::Vector<double> W(ntot);
  tmv::Vector<std::complex<double> > Z(ntot);

  for(size_t k=0,n=0;k<pix.size();k++) {
    double sigma_obs = psf ?
      sqrt(pow(sigma,2)+f_psf*pow((*psf)[k].GetSigma(),2)) :
      sigma;
    for(size_t i=0;i<pix[k].size();i++,n++) {
      I(n) = pix[k][i].I*pix[k][i].wt;
      W(n) = pix[k][i].wt;
      std::complex<double> z1 = pix[k][i].z;
      std::complex<double> z2 = m*(z1-cen) - m*gamma*conj(z1-cen);
      Z(n) = z2 / sigma_obs;
    }
  }

  size_t nsize = (order+1)*(order+2)/2;
  tmv::Matrix<double> A(ntot,nsize);
  MakePsi(A.View(),Z,order,&W);
  xxdbg<<"after makepsi\n";
  xxdbg<<"A.rows(0,10) = "<<A.Rows(0,10);

  if (psf) {
    for(size_t k=0,n=0,nx;k<pix.size();k++,n=nx) {
      size_t psfsize = (*psf)[k].size();
      BVec newpsf = (*psf)[k];
      tmv::Matrix<double> S(psfsize,psfsize);
      GTransform(gamma,newpsf.GetOrder(),S);
      newpsf = S * (*psf)[k];
      tmv::Matrix<double> D(psfsize,psfsize);
      MuTransform(mu,newpsf.GetOrder(),D);
      newpsf = D * newpsf;
      tmv::Matrix<double> C(nsize,nsize);
      PSFConvolve(newpsf,b.GetOrder(),b.GetSigma(),C);
      nx = n+pix[k].size();
      A.Rows(n,nx) *= C;
    }
    xxdbg<<"after psf correction\n";
  }
  b = I/A;
  if (bcov) {
    tmv::UpperTriMatrix<double> Rinv = A.QRD().GetR();
    Rinv.InvertSelf();
    *bcov = Rinv * Rinv.Transpose();
  }
}

void Ellipse::MeasureShapelet(const std::vector<PixelList>& pix,
    const std::vector<BVec>& psf, BVec& b) const
{ DoMeasureShapelet(pix,&psf,b); }

void Ellipse::MeasureShapelet(const std::vector<PixelList>& pix,
    BVec& b) const
{ DoMeasureShapelet(pix,0,b); }

