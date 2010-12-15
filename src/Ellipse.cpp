
#include "Ellipse.h"
#include <cmath>
#include "dbg.h"
#include <fstream>
#include "PsiHelper.h"
#include "Params.h"

#define N_FLUX_ATTEMPTS 0
#define MAXITER 4

void Ellipse::shiftBy(
    const std::complex<double>& c2,
    const std::complex<double>& g2,
    const double m2)
{
    // Current values are c1, g1, m1.
    // These effect a transformation:
    //
    // z' = exp(-m1)/sqrt(1-|g1|^2) ( z-c1 - g1 (z-c1)* )
    //
    // Now execute a second transformation (c2,g2,m2):
    //
    // z'' = exp(-m2)/sqrt(1-|g2|^2) * 
    //    [ ( exp(-m1)/sqrt(1-|g1|^2) ( z-c1 - g1 (z-c1)* ) - c2 ) -
    //      g2 * ( exp(-m1)/sqrt(1-|g1|^2) ( z-c1 - g1 (z-c1)* ) - c2 )* ]
    //
    // We need to figure out what this corresponds to in terms of an 
    // effective:
    // z'' = exp(-m3)/sqrt(1-|g3|^2) ( z-c3 - g3 (z-c3)* )
    //
    // This is pretty complicated, but we can do it in steps.
    //
    // First, g3:
    // Taking just the g terms, we find that 
    // z'' =  [(1+g1* g2) z - (g1+g2) z*] / sqrt(1-|g1|^2) / sqrt(1-|g2|^2)
    //     =  (1+g1* g2) [z - (g1+g2)/(1+g1* g2) z*] 
    //                        / sqrt(1-|g1|^2) / sqrt(1-|g2|^2)
    // 
    // So we have an equation for g3:
    // g3/sqrt(1-|g3|^2) = (g1+g2)/(1+g1* g2) / sqrt(1-|g1|^2) / sqrt(1-|g2|^2)
    //
    // This isn't too hard to solve.  First do |g3|^2.  Then g3.
    // 
    // But look again at the above equation for z''. 
    // The first term is (1+g1* g2) which is complex.
    // This means it is both a dilation and a rotation.
    // But we aren't modelling rotations.  So we need to take it out.
    // We do this by multiplying z by (1+g1* g2)/|1+g1* g2| before
    // doing the rest of the stuff.
    //
    // So the real equation for g3 is:
    // g3/sqrt(1-|g3|^2) = (g1+g2)/|1+g1* g2| / sqrt(1-|g1|^2) / sqrt(1-|g2|^2)
    // 
    // Next m3:
    // We have a factor of |1+g1 g2*| from the shears.
    // And we also have exp(-m1) and exp(-m2). 
    // exp(-m3) = |1+g1 g2*| sqrt(1-g3^2) / sqrt(1-g1^2) / sqrt(1-g2^2)
    //
    // The m1 and m2 parts add an extra exp(-m1) * exp(-m2).
    //
    // Finally, the translation terms are a bit messy.
    // The translation terms in the equation for z'' are:
    //
    // exp(-m2)/sqrt(1-|g2|^2) * 
    //    [ ( exp(-m1)/sqrt(1-|g1|^2) ( -c1 - g1 (-c1)* ) - c2 ) -
    //      g2 * ( exp(-m1)/sqrt(1-|g1|^2) ( -c1 - g1 (-c1)* ) - c2 )* ]
    // 
    // This is some value, which we set equal to:
    //
    // exp(-m3)/sqrt(1-|g3|^2) (-c3 - g3(-c3)* )
    //
    // So solving for c3 - g3 c3* is straightforward.
    // c3 - g3 c3* = X
    // c3* - g3* c3 = X*
    // g3 c3* - |g3|^2 c3 = g3 X*
    // c3 - |g3|^2 c3 = X + g3 X*
    // c3 = (X + g3 X*) / (1-|g3|^2)
    //

    std::complex<double> c1 = _cen;
    std::complex<double> g1 = _gamma;
    double m1 = _mu;

    dbg<<"Start Ellipse::shiftBy\n";
    dbg<<"c1,g1,m1 = "<<c1<<"  "<<g1<<"  "<<m1<<std::endl;
    dbg<<"c2,g2,m2 = "<<c2<<"  "<<g2<<"  "<<m2<<std::endl;

    std::complex<double> g12 = 
        (g1+g2) / sqrt(1.-norm(g1)) / sqrt(1.-norm(g2));
    g12 /= abs(1.+conj(g1)*g2);
    // |g3|^2 / (1-|g3|^2) = |g12|^2
    // |g3|^2 = |g12|^2 - |g12|^2 |g3|^2
    // |g3|^2 = |g12|^2 / (1+|g12|^2)
    double g12sq = norm(g12);
    double g3sq = g12sq / (1.+g12sq);
    std::complex<double> g3 = g12 * sqrt(1.-g3sq);
    dbg<<"g3 = "<<g3<<", g3sq = "<<g3sq<<" = "<<norm(g3)<<std::endl;

    if (XDEBUG && g3sq > 0) {
        // Check with Miralda-Escude formula:
        std::complex<double> d1 = g1 * 2./(1.+norm(g1));
        std::complex<double> d2 = g2 * 2./(1.+norm(g2));
        std::complex<double> d3 = 
            d1 + d2 + 
            (1.-sqrt(1.-norm(d2))) * d2/norm(d2) * (conj(d1)*d2 - d1*conj(d2))/2.;
        d3 /= 1. + (real(d1)*real(d2) + imag(d1)*imag(d2));
        xdbg<<"Miralda-Escude formula gives d3 = "<<d3<<std::endl;
        xdbg<<"g3->d3 = "<<g3 * 2./(1.+norm(g3))<<std::endl;
    }

    double m3 = exp(-m1) * exp(-m2) * abs(1.+conj(g1)*g2);
    m3 = -log(m3);
    dbg<<"m3 = "<<m3<<std::endl;

    std::complex<double> temp = -c1 - g1 * conj(-c1);
    temp = exp(-m1)/sqrt(1.-norm(g1)) * temp - c2;
    temp = exp(-m2)/sqrt(1.-norm(g2)) * (temp - g2 * conj(temp));
    temp /= -exp(-m3)/sqrt(1.-norm(g3));
    std::complex<double> c3 = temp + g3*conj(temp) / (1.+g3sq);

    dbg<<"c3 = "<<c3<<std::endl;

    _cen = c3;
    _gamma = g3;
    _mu = m3;
}

bool Ellipse::measure(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>& psf,
    int order, double sigma, bool shouldUseInteg, long& flag,
    DMatrix* cov, BVec* bRet, DMatrix* bCov,
    std::vector<PixelList>* pixels_model)
{ 
    return doMeasure(pix,&psf,order,sigma,shouldUseInteg,flag,cov,bRet,bCov,
                     pixels_model); 
}

bool Ellipse::measure(
    const std::vector<PixelList>& pix,
    int order, double sigma, bool shouldUseInteg, long& flag, 
    DMatrix* cov, BVec* bRet, DMatrix* bCov,
    std::vector<PixelList>* pixels_model)
{
    return doMeasure(pix,0,order,sigma,shouldUseInteg,flag,cov,bRet,bCov,
                     pixels_model);
}

void Ellipse::doMeasureShapelet(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>* psf, BVec& b, int order, int order2,
    DMatrix* bCov) const
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
    //xdbg<<"Z = "<<Z<<std::endl;
    //xdbg<<"I = "<<I<<std::endl;
    //xdbg<<"W = "<<W<<std::endl;

    const int bsize = (order+1)*(order+2)/2;
    xdbg<<"bsize = "<<bsize<<std::endl;

    DMatrix A(nTot,bsize);
    Assert(nTot >= bsize); // Should have been addressed by calling routine.
    makePsi(A,TMV_vview(Z),order,&W);
    //xdbg<<"A = "<<A<<std::endl;

    if (psf) {
        for(int k=0,n=0,nx;k<nExp;++k,n=nx) {
            xdbg<<"psf = "<<(*psf)[k]<<std::endl;
            int psforder = (*psf)[k].getOrder();
            xdbg<<"psforder = "<<psforder<<std::endl;
            if (order2 > psforder) { psforder = order2; }
            xdbg<<"psforder => "<<psforder<<std::endl;
            const double psfsigma = (*psf)[k].getSigma();
            xdbg<<"sigma = "<<psfsigma<<std::endl;
            BVec newpsf(psforder,psfsigma);
            int psfsize = newpsf.size();
            xdbg<<"psfsize = "<<psfsize<<std::endl;
            DMatrix S(psfsize,psfsize);
            calculateGTransform(_gamma,psforder,S);
            newpsf.vec() = S.colRange(0,(*psf)[k].size()) * (*psf)[k].vec();
            xdbg<<"newpsf = "<<newpsf<<std::endl;
            DMatrix D(psfsize,psfsize);
            calculateMuTransform(_mu,psforder,D);
            newpsf.vec() = D * newpsf.vec();
            xdbg<<"newpsf => "<<newpsf<<std::endl;
            DMatrix C(bsize,bsize);
            calculatePsfConvolve(newpsf,order,b.getSigma(),C);
            const int nPix = pix[k].size();
            nx = n+nPix;
            TMV_rowRange(A,n,nx) *= C;
        }
    }
#ifdef USE_TMV
    A.saveDiv();
#if 1
    // QRP should be good enough, even if A is singular.
    // But also, A shouldn't really become singular here
    // (unlike in EllipseSolver where it sometimes does).
    A.divideUsing(tmv::QRP);
#else
    A.divideUsing(tmv::SV);
    const double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());
    A.svd().thresh(sqrtEps);
    xdbg<<"For MeasureShapelet: svd = "<<A.svd().getS().diag()<<std::endl;
    if (A.svd().getKMax() < int(A.rowsize()))
        xdbg<<"Omitting last "<<A.rowsize()-A.svd().getKMax()<<" singular values\n";
#endif
    b.vec().subVector(0,bsize) = I/A;
    xdbg<<"b = "<<b<<std::endl;
    xdbg<<"Norm(I) = "<<Norm(I)<<std::endl;
    xdbg<<"Norm(A*b) = "<<Norm(A*b.vec().subVector(0,bsize))<<std::endl;
    xdbg<<"Norm(I-A*b) = "<<Norm(I-A*b.vec().subVector(0,bsize))<<std::endl;
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
    xdbg<<"For MeasureShapelet: svd = "<<svd_s.transpose()<<std::endl;
    if (kmax < svd_s.size())
        xdbg<<"Omitting last "<<svd_s.size()-kmax<<" singular values\n";
    // USVtb = I
    // b = VS^-1UtI
    DVector temp = TMV_colRange(svd_u,0,kmax).transpose() * I;
    temp = svd_s.TMV_subVector(0,kmax).cwise().inverse().asDiagonal() * temp;
    b.vec().TMV_subVector(0,bsize) = TMV_colRange(svd_v,0,kmax) * temp;
    if (bCov) {
        // (AtA)^-1 = (VSUt USVt)^-1 = (V S^2 Vt)^-1 = V S^-2 Vt
        DMatrix temp2 = 
            svd_s.TMV_subVector(0,kmax).cwise().square().inverse().asDiagonal() *
            TMV_colRange(svd_v,0,kmax).transpose();
        *bCov = TMV_colRange(svd_v,0,kmax).transpose() * temp2;
    }
#endif
    if (order < b.getOrder()) {
        xdbg<<"Need to zero the rest of b\n";
        xdbg<<"bsize = "<<bsize<<"  b.size = "<<b.size()<<std::endl;
        // Zero out the rest of the shapelet vector:
        b.vec().TMV_subVector(bsize,b.size()).setZero();
    }
    xdbg<<"Done measure Shapelet\n";
    xdbg<<"b = "<<b.vec()<<std::endl;
}

void Ellipse::doAltMeasureShapelet(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>* psf, BVec& b, int order, int order2,
    double pixScale, DMatrix* bCov) const
{
    xdbg<<"Start AltMeasureShapelet: order = "<<order<<std::endl;
    xdbg<<"b.order, sigma = "<<b.getOrder()<<", "<<b.getSigma()<<std::endl;
    xdbg<<"el = "<<*this<<std::endl;
    const int nExp = pix.size();

    const int bsize = (order+1)*(order+2)/2;
    b.vec().setZero();
    double sigma = b.getSigma();

    DDiagMatrix normD(bsize);
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

        DMatrix A(nPix,bsize);
        makePsi(A,TMV_vview(Z),order);

        // b = C^-1 normD At I
        DVector b1 = A.transpose() * I;
        b1 = normD EIGEN_asDiag() * b1;
        if (psf) {
            int psforder = (*psf)[k].getOrder();
            if (order2 > psforder) { psforder = order2; }
            const double psfsigma = (*psf)[k].getSigma();
            BVec newpsf(psforder,psfsigma);
            const int psfsize = newpsf.size();
            DMatrix S(psfsize,psfsize);
            calculateGTransform(_gamma,psforder,S);
            newpsf.vec() = S.colRange(0,(*psf)[k].size()) * (*psf)[k].vec();
            DMatrix D(psfsize,psfsize);
            calculateMuTransform(_mu,psforder,D);
            newpsf.vec() = D * newpsf.vec();
            DMatrix C(bsize,bsize);
            calculatePsfConvolve(newpsf,order,b.getSigma(),C);
#ifdef USE_TMV
            b1 /= C;
#else 
            C.lu().solve(b1,&b1);
#endif
        }
        b.vec() += b1;
    }
    if (order < b.getOrder()) {
        xdbg<<"Need to zero the rest of b\n";
        xdbg<<"bsize = "<<bsize<<"  b.size = "<<b.size()<<std::endl;
        // Zero out the rest of the shapelet vector:
        b.vec().TMV_subVector(bsize,b.size()).setZero();
    }
    xdbg<<"Done measure Shapelet\n";
}

void Ellipse::measureShapelet(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>& psf, BVec& b, int order, int order2,
    DMatrix* bCov) const
{ doMeasureShapelet(pix,&psf,b,order,order2,bCov); }

void Ellipse::measureShapelet(
    const std::vector<PixelList>& pix, BVec& b, int order, int order2,
    DMatrix* bCov) const
{ doMeasureShapelet(pix,0,b,order,order2,bCov); }

void Ellipse::altMeasureShapelet(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>& psf, BVec& b, int order, int order2,
    double pixScale, DMatrix* bCov) const
{ doAltMeasureShapelet(pix,&psf,b,order,order2,pixScale,bCov); }

void Ellipse::altMeasureShapelet(
    const std::vector<PixelList>& pix, BVec& b, int order, int order2,
    double pixScale, DMatrix* bCov) const
{ doAltMeasureShapelet(pix,0,b,order,order2,pixScale,bCov); }

