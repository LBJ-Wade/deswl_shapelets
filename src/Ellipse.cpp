
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

bool Ellipse::measure(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>& psf,
    int order, double sigma, bool shouldUseInteg, long& flag,
    tmv::Matrix<double>* cov, BVec* bRet, tmv::Matrix<double>* bCov)
{ 
    return doMeasure(pix,&psf,order,sigma,shouldUseInteg,flag,cov,bRet,bCov); 
}

bool Ellipse::measure(
    const std::vector<PixelList>& pix,
    int order, double sigma, bool shouldUseInteg, long& flag, 
    tmv::Matrix<double>* cov, BVec* bRet, tmv::Matrix<double>* bCov)
{
    return doMeasure(pix,0,order,sigma,shouldUseInteg,flag,cov,bRet,bCov);
}

void Ellipse::doMeasureShapelet(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>* psf, BVec& b, tmv::Matrix<double>* bCov) const
{
    xdbg<<"Start MeasureShapelet: order = "<<b.getOrder()<<std::endl;
    xdbg<<"sigma = "<<b.getSigma()<<std::endl;
    xdbg<<"el = "<<*this<<std::endl;

    // ( u )' = exp(-mu)/sqrt(1-gsq) ( 1-g1  -g2  ) ( u-uc )
    // ( v )                         ( -g2   1+g1 ) ( v-vc )
    // 
    // z' = u' + I v' 
    //    = exp(-mu)/sqrt(1-gsq)
    //             [ (1-g1)(u-uc)-g2(v-vc)-Ig2(u-uc)+I(1+g1)(v-vc) ]
    //    = exp(-mu)/sqrt(1-gsq) [ z-zc - g1(z-zc)* -Ig2(z-zc)* ]
    //    = exp(-mu)/sqrt(1-gsq) ( z-zc - g (z-zc)* )

    int order = b.getOrder();
    double sigma = b.getSigma();
    double gsq = std::norm(_gamma);
    Assert(gsq < 1.);

    double m = exp(-_mu)/sqrt(1.-gsq);

    int nTot = 0;
    const int nExp = pix.size();
    for(int i=0;i<nExp;++i) nTot += pix[i].size();
    dbg<<"ntot = "<<nTot<<" in "<<nExp<<" images\n";

    tmv::Vector<double> I(nTot);
    tmv::DiagMatrix<double> W(nTot);
    tmv::Vector<std::complex<double> > Z(nTot);

    for(int k=0,n=0;k<nExp;++k) {
        double sigma_obs = 
            psf ?
            sqrt(pow(sigma,2)+_fPsf*pow((*psf)[k].getSigma(),2)) :
            sigma;

        const int nPix = pix[k].size();
        for(int i=0;i<nPix;++i,++n) {
            I(n) = pix[k][i].getFlux()*pix[k][i].getInverseSigma();
            W(n) = pix[k][i].getInverseSigma();
            std::complex<double> z1 = pix[k][i].getPos();
            std::complex<double> z2 = m*(z1-_cen) - m*_gamma*conj(z1-_cen);
            Z(n) = z2 / sigma_obs;
        }
    }

    const int bSize = b.size();
    tmv::Matrix<double> A(nTot,bSize);
    makePsi(A.View(),Z,order,&W);

    if (psf) {
        for(int k=0,n=0,nx;k<nExp;++k,n=nx) {
            const int psfSize = (*psf)[k].size();
            BVec newPsf = (*psf)[k];
            tmv::Matrix<double> S(psfSize,psfSize);
            calculateGTransform(_gamma,newPsf.getOrder(),S);
            newPsf.vec() = S * (*psf)[k].vec();
            tmv::Matrix<double> D(psfSize,psfSize);
            calculateMuTransform(_mu,newPsf.getOrder(),D);
            newPsf.vec() = D * newPsf.vec();
            tmv::Matrix<double> C(bSize,bSize);
            calculatePsfConvolve(newPsf,b.getOrder(),b.getSigma(),C);
            const int nPix = pix[k].size();
            nx = n+nPix;
            A.Rows(n,nx) *= C;
        }
    }
    A.DivideUsing(tmv::SV);
    A.SaveDiv();
    const double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());
    A.SVD().Thresh(sqrtEps);
    dbg<<"For MeasureShapelet: svd = "<<A.SVD().GetS().diag()<<std::endl;
    dbg<<"Omitting last "<<A.rowsize()-A.SVD().GetKMax()<<" singular values\n";
    b.vec() = I/A;
    if (bCov) {
        A.InverseATA(*bCov);
    }
}

void Ellipse::measureShapelet(const std::vector<PixelList>& pix,
                              const std::vector<BVec>& psf, BVec& b, tmv::Matrix<double>* bCov) const
{ doMeasureShapelet(pix,&psf,b,bCov); }

void Ellipse::measureShapelet(const std::vector<PixelList>& pix,
                              BVec& b, tmv::Matrix<double>* bCov) const
{ doMeasureShapelet(pix,0,b,bCov); }

