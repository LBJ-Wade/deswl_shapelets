
#include "Ellipse.h"
#include "EllipseSolver.h"
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
    DMatrix* cov, BVec* bRet, DMatrix* bCov)
{ 
    return doMeasure(pix,&psf,order,sigma,shouldUseInteg,flag,cov,bRet,bCov); 
}

bool Ellipse::measure(
    const std::vector<PixelList>& pix,
    int order, double sigma, bool shouldUseInteg, long& flag, 
    DMatrix* cov, BVec* bRet, DMatrix* bCov)
{
    return doMeasure(pix,0,order,sigma,shouldUseInteg,flag,cov,bRet,bCov);
}

void Ellipse::doMeasureShapelet(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>* psf, BVec& b, int order, DMatrix* bCov) const
{
    xdbg<<"Start MeasureShapelet: order = "<<order<<std::endl;
    xdbg<<"b.order, sigma = "<<b.getOrder()<<", "<<b.getSigma()<<std::endl;
    xdbg<<"el = "<<*this<<std::endl;

    // ( u )' = exp(-mu)/sqrt(1-gsq) ( 1-g1  -g2  ) ( u-uc )
    // ( v )                         ( -g2   1+g1 ) ( v-vc )
    // 
    // z' = u' + I v' 
    //    = exp(-mu)/sqrt(1-gsq)
    //             [ (1-g1)(u-uc)-g2(v-vc)-Ig2(u-uc)+I(1+g1)(v-vc) ]
    //    = exp(-mu)/sqrt(1-gsq) [ z-zc - g1(z-zc)* -Ig2(z-zc)* ]
    //    = exp(-mu)/sqrt(1-gsq) ( z-zc - g (z-zc)* )

    double sigma = b.getSigma();
    double gsq = std::norm(_gamma);
    Assert(gsq < 1.);

    double m = exp(-_mu)/sqrt(1.-gsq);

    int nTot = 0;
    const int nExp = pix.size();
    for(int i=0;i<nExp;++i) nTot += pix[i].size();
    //dbg<<"ntot = "<<nTot<<" in "<<nExp<<" images\n";

    DVector I(nTot);
    DVector W(nTot);
    CDVector Z(nTot);

    for(int k=0,n=0;k<nExp;++k) {
        double sigma_obs = 
            psf ?
            sqrt(pow(sigma,2)+pow((*psf)[k].getSigma(),2)) :
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

    const int bSize = (order+1)*(order+2)/2;

    DMatrix A(nTot,bSize);
    Assert(nTot >= bSize); // Should have been addressed by calling routine.
    makePsi(A,TMV_vview(Z),order,&W);

    // TODO: Convert this to the CHIP_FRAME method.

    if (psf) {
        for(int k=0,n=0,nx;k<nExp;++k,n=nx) {
            const int psfSize = (*psf)[k].size();
            BVec newPsf = (*psf)[k];
            DMatrix S(psfSize,psfSize);
            calculateGTransform(_gamma,newPsf.getOrder(),S);
            newPsf.vec() = S * (*psf)[k].vec();
            DMatrix D(psfSize,psfSize);
            calculateMuTransform(_mu,newPsf.getOrder(),D);
            newPsf.vec() = D * newPsf.vec();
            DMatrix C(bSize,bSize);
            calculatePsfConvolve(newPsf,order,b.getSigma(),C);
            const int nPix = pix[k].size();
            nx = n+nPix;
            TMV_rowRange(A,n,nx) *= C;
        }
    }
#ifdef USE_TMV
    A.divideUsing(tmv::SV);
    A.saveDiv();
    const double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());
    A.svd().thresh(sqrtEps);
    dbg<<"For MeasureShapelet: svd = "<<A.svd().getS().diag()<<std::endl;
    if (A.svd().getKMax() < int(A.rowsize()))
        dbg<<"Omitting last "<<A.rowsize()-A.svd().getKMax()<<" singular values\n";
    b.vec().subVector(0,bSize) = I/A;
    b.vec().subVector(bSize,b.size()).setZero();
    if (bCov) {
        A.makeInverseATA(*bCov);
    }
#else
    const double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());
    Eigen::SVD<DMatrix> svd = A.svd().sort();
    const DMatrix& svd_u = svd.matrixU();
    const DVector& svd_s = svd.singularValues();
    const DMatrix& svd_v = svd.matrixV();
    double max = svd_s(0);
    int kmax = 0;
    double thresh = sqrtEps * max;
    while (kmax < svd_s.size() && svd_s(kmax) > thresh) ++kmax;
    dbg<<"For MeasureShapelet: svd = "<<svd_s.transpose()<<std::endl;
    if (kmax < svd_s.size())
        dbg<<"Omitting last "<<svd_s.size()-kmax<<" singular values\n";
    // USVtb = I
    // b = VS^-1UtI
    DVector temp = TMV_colRange(svd_u,0,kmax).transpose() * I;
    temp = svd_s.TMV_subVector(0,kmax).cwise().inverse().asDiagonal() * temp;
    b.vec().TMV_subVector(0,bSize) = TMV_colRange(svd_v,0,kmax) * temp;
    if (bCov) {
        // (AtA)^-1 = (VSUt USVt)^-1 = (V S^2 Vt)^-1 = V S^-2 Vt
        DMatrix temp2 = 
            svd_s.TMV_subVector(0,kmax).cwise().square().inverse().asDiagonal() *
            TMV_colRange(svd_v,0,kmax).transpose();
        *bCov = TMV_colRange(svd_v,0,kmax).transpose() * temp2;
    }
#endif
    if (order < b.getOrder()) {
        dbg<<"Need to zero the rest of b\n";
        dbg<<"bSize = "<<bSize<<"  b.size = "<<b.size()<<std::endl;
        // Zero out the rest of the shapelet vector:
        b.vec().TMV_subVector(bSize,b.size()).setZero();
    }
    xdbg<<"Done measure Shapelet\n";
    dbg<<"b = "<<b.vec()<<std::endl;
}

void Ellipse::doAltMeasureShapelet(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>* psf, BVec& b, int order,
    double pixScale, DMatrix* bCov) const
{
    xdbg<<"Start AltMeasureShapelet: order = "<<order<<std::endl;
    xdbg<<"b.order, sigma = "<<b.getOrder()<<", "<<b.getSigma()<<std::endl;
    xdbg<<"el = "<<*this<<std::endl;
    const int nExp = pix.size();

    const int bSize = (order+1)*(order+2)/2;
    b.vec().setZero();
    double sigma = b.getSigma();

    DDiagMatrix normD(bSize);
    double val = pixScale*pixScale/(sigma*sigma)/nExp;
    for(int n=0,k=0;n<=b.getOrder();++n) {
        for(int p=n,q=0;p>=q;--p,++q,++k) {
            if (p!=q) { normD(k) = val/2.; normD(++k) = val/2.; }
            else normD(k) = val;
        }
    }

    double gsq = std::norm(_gamma);
    Assert(gsq < 1.);
    double m = exp(-_mu)/sqrt(1.-gsq);

    for(int k=0;k<nExp;++k) {
        const int nPix = pix[k].size();
        DVector I(nPix);
        CDVector Z(nPix);

        double sigma_obs = 
            psf ?
            sqrt(pow(sigma,2)+pow((*psf)[k].getSigma(),2)) :
            sigma;

        for(int i=0;i<nPix;++i) {
            I(i) = pix[k][i].getFlux();
            std::complex<double> z1 = pix[k][i].getPos();
            std::complex<double> z2 = m*(z1-_cen) - m*_gamma*conj(z1-_cen);
            Z(i) = z2 / sigma_obs;
        }

        DMatrix A(nPix,bSize);
        makePsi(A,TMV_vview(Z),order);

        // b = C^-1 normD At I
        DVector b1 = A.transpose() * I;
        b1 = normD EIGEN_asDiag() * b1;

        if (psf) {
            const int psfSize = (*psf)[k].size();
            BVec newPsf = (*psf)[k];
            DMatrix S(psfSize,psfSize);
            calculateGTransform(_gamma,newPsf.getOrder(),S);
            newPsf.vec() = S * (*psf)[k].vec();
            DMatrix D(psfSize,psfSize);
            calculateMuTransform(_mu,newPsf.getOrder(),D);
            newPsf.vec() = D * newPsf.vec();
            DMatrix C(bSize,bSize);
            calculatePsfConvolve(newPsf,order,b.getSigma(),C);
#ifdef USE_TMV
            b1 /= C;
#else 
            C.lu().solve(b1,&b1);
#endif
        }
        b.vec() += b1;
    }
    if (order < b.getOrder()) {
        dbg<<"Need to zero the rest of b\n";
        dbg<<"bSize = "<<bSize<<"  b.size = "<<b.size()<<std::endl;
        // Zero out the rest of the shapelet vector:
        b.vec().TMV_subVector(bSize,b.size()).setZero();
    }
    dbg<<"Done measure Shapelet\n";
}

void Ellipse::measureShapelet(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>& psf, BVec& b, int order, DMatrix* bCov) const
{ doMeasureShapelet(pix,&psf,b,order,bCov); }

void Ellipse::measureShapelet(
    const std::vector<PixelList>& pix, BVec& b, int order, DMatrix* bCov) const
{ doMeasureShapelet(pix,0,b,order,bCov); }

void Ellipse::altMeasureShapelet(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>& psf, BVec& b, int order,
    double pixScale, DMatrix* bCov) const
{ doAltMeasureShapelet(pix,&psf,b,order,pixScale,bCov); }

void Ellipse::altMeasureShapelet(
    const std::vector<PixelList>& pix, BVec& b, int order,
    double pixScale, DMatrix* bCov) const
{ doAltMeasureShapelet(pix,0,b,order,pixScale,bCov); }

