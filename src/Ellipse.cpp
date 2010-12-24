
#include "Ellipse.h"
#include <cmath>
#include "dbg.h"
#include <fstream>
#include "PsiHelper.h"
#include "Params.h"

static std::complex<double> addShears(
    const std::complex<double> g1, const std::complex<double> g2)
{
    double absg1 = std::abs(g1);
    if (absg1 == 0.) return g2;
    double absd1 = tanh(atanh(absg1)*2.);
    std::complex<double> d1 = absd1 / absg1 * g1;
    double absg2 = std::abs(g2);
    if (absg2 == 0.) return g1;
    double absd2 = tanh(atanh(absg2)*2.);
    std::complex<double> d2 = absd2 / absg2 * g2;

    double d2sq = absd2*absd2;
    double x = (1.-std::sqrt(1.-d2sq))/d2sq;
    std::complex<double> d3 = d1+d2*(1.+x*(std::real(d1)*d2-std::real(d2)*d1));
    d3 /= 1. + std::real(d1*std::conj(d2));
    double absd3 = std::abs(d3);
    double absg3 = tanh(atanh(absd3)/2.);
    std::complex<double> g3 = absg3/absd3*d3;
    xdbg<<"GammaAdd: "<<g1<<" + "<<g2<<" = "<<g3<<std::endl;
    return g3;
}

void Ellipse::shiftBy(
    const std::complex<double>& c1,
    const std::complex<double>& g1,
    const double m1)
{
    // Current values are c2, g2, m2.
    // The model is that a round galaxy is first transformed by (c1,g1,m1).
    // And then by the current values of (c2,g2,m2).
    //
    // The first set (c1,g1,m1) effect a transformation:
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
    // doing the rest of the stuff.  This is allowed because we are 
    // modelling the galaxy to be initially round, so a rotation doesn't
    // change that.
    //
    // So the real equation for g3 is:
    // g3/sqrt(1-|g3|^2) = (g1+g2)/|1+g1* g2| / sqrt(1-|g1|^2) / sqrt(1-|g2|^2)
    // 
    // Next m3:
    // We have a factor of |1+g1* g2| from the shears.
    // And we also have exp(-m1) and exp(-m2). 
    // exp(-m3) = |1+g1* g2| exp(-m1 - m2)
    // m3 = m1 + m2 - ln(|1+g1* g2|)
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

    std::complex<double> c2 = _cen;
    std::complex<double> g2 = _gamma;
    double m2 = _mu;

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

    xdbg<<"Miralda-Escude formula gives g3 = "<<addShears(g1,g2)<<std::endl;
    xdbg<<"Reverse Miralda-Escude formula gives g3 = "<<addShears(g2,g1)<<std::endl;

    double m3 = m1 + m2 - log(abs(1.+conj(g1)*g2));
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
    int order, int order2, double sigma, long& flag, double thresh,
    DMatrix* cov, BVec* bRet, DMatrix* bCov,
    std::vector<PixelList>* pixels_model)
{ 
    return doMeasure(
        pix,&psf,order,order2,sigma,flag,thresh,cov,bRet,bCov,pixels_model); 
}

bool Ellipse::measure(
    const std::vector<PixelList>& pix,
    int order, int order2, double sigma, long& flag, double thresh,
    DMatrix* cov, BVec* bRet, DMatrix* bCov,
    std::vector<PixelList>* pixels_model)
{
    return doMeasure(
        pix,0,order,order2,sigma,flag,thresh,cov,bRet,bCov,pixels_model);
}

bool Ellipse::doMeasureShapelet(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>* psf, BVec& b, int& order, int order2,
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
    dbg<<"ntot = "<<nTot<<" in "<<nExp<<" images\n";

    DVector I(nTot);
    DVector W(nTot);
    CDVector Z(nTot);

    for(int k=0,n=0;k<nExp;++k) {
        double sigma_obs = 
            psf ?
            sqrt(pow(sigma,2)+pow((*psf)[k].getSigma(),2)) :
            sigma;
        dbg<<"sigma_obs["<<k<<"] = "<<sigma_obs<<std::endl;

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

    int bsize = (order+1)*(order+2)/2;
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
    A.divideUsing(tmv::SV);
    std::auto_ptr<tmv::MatrixView<double> > A_use(
        new tmv::MatrixView<double>(A.view()));
    A_use->svd();
    xdbg<<"svd = "<<A_use->svd().getS().diag()<<std::endl;
    while (A_use->svd().isSingular() || A_use->svd().condition() > 1.e3) {
        dbg<<"Poor condition in MeasureShapelet: \n";
        xdbg<<"svd = "<<A_use->svd().getS().diag()<<std::endl;
        xdbg<<"condition = "<<A_use->svd().condition()<<std::endl;
        if (order > 2) {
            bsize -= order+1; --order;
            dbg<<"Reducing order to "<<order<<std::endl;
            A_use.reset(new tmv::MatrixView<double>(A.colRange(0,bsize)));
            continue;
        } else {
            dbg<<"Order is already only 2.  Unable to measure shapelet.\n";
            return false;
        }
    }
    b.vec().subVector(0,bsize) = I/(*A_use);
    xdbg<<"b = "<<b<<std::endl;
    xdbg<<"Norm(I) = "<<Norm(I)<<std::endl;
    xdbg<<"Norm(A*b) = "<<Norm((*A_use)*b.vec().subVector(0,bsize))<<std::endl;
    xdbg<<"Norm(I-A*b) = "<<Norm(I-(*A_use)*b.vec().subVector(0,bsize))<<std::endl;
    while (!(b(0) > 0.)) {
        dbg<<"Calculated b vector has negative b(0):\n";
        dbg<<"b = "<<b<<std::endl;
        if (order > 2) {
            bsize -= order+1; --order;
            dbg<<"Reducing order to "<<order<<std::endl;
            A_use.reset(new tmv::MatrixView<double>(A.colRange(0,bsize)));
            b.vec().subVector(0,bsize) = I/(*A_use);
            continue;
        } else {
            dbg<<"Order is already only 2.  Unable to measure shapelet.\n";
            return false;
        }
    }
    if (bCov) {
        bCov->setZero();
        A_use->makeInverseATA(bCov->subMatrix(0,bsize,0,bsize));
    }
#else
    const double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());
    Eigen::SVD<DMatrix> svd = A.svd().sort();
    const DMatrix& svd_u = svd.matrixU();
    const DVector& svd_s = svd.singularValues();
    const DMatrix& svd_v = svd.matrixV();
    double max = svd_s(0);
    double min = svd_s(svd_s.size()-1);
    if (min < 1.e-6 * max) {
        xdbg<<"Poor condition in MeasureShapelet: \n";
        xdbg<<"svd = "<<svd_s.transpose()<<std::endl;
        xdbg<<"condition = "<<max/min<<std::endl;
        if (order > 2) {
            dbg<<"Reducing order to "<<order-1<<std::endl;
            return doMeasureShapelet(pix,psf,b,order-1,order2,bCov);
        } else {
            dbg<<"Order is already only 2.  Unable to measure shapelet.\n";
            return false;
        }
    }
    // USVtb = I
    // b = VS^-1UtI
    DVector temp = svd_u.transpose() * I;
    temp = svd_s.cwise().inverse().asDiagonal() * temp;
    b.vec().TMV_subVector(0,bsize) = svd_v * temp;
    if (!(b(0) > 0.)) {
        dbg<<"Calculated b vector has negative b(0):\n";
        dbg<<"b = "<<b<<std::endl;
        if (order > 2) {
            dbg<<"Reducing order to "<<order-1<<std::endl;
            return doMeasureShapelet(pix,psf,b,order-1,order2,bCov);
        } else {
            dbg<<"Order is already only 2.  Unable to measure shapelet.\n";
            return false;
        }
    }
    if (bCov) {
        bCov->setZero();
        // (AtA)^-1 = (VSUt USVt)^-1 = (V S^2 Vt)^-1 = V S^-2 Vt
        DMatrix temp2 = 
            svd_s.cwise().square().inverse().asDiagonal() * svd_v.transpose();
        bCov->TMV_subMatrix(0,bsize,0,bsize) = svd_v.transpose() * temp2;
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
    return true;
}

bool Ellipse::doAltMeasureShapelet(
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
    return true;
}

bool Ellipse::measureShapelet(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>& psf, BVec& b, int& order, int order2,
    DMatrix* bCov) const
{ return doMeasureShapelet(pix,&psf,b,order,order2,bCov); }

bool Ellipse::measureShapelet(
    const std::vector<PixelList>& pix, BVec& b, int& order, int order2,
    DMatrix* bCov) const
{ return doMeasureShapelet(pix,0,b,order,order2,bCov); }

bool Ellipse::altMeasureShapelet(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>& psf, BVec& b, int order, int order2,
    double pixScale, DMatrix* bCov) const
{ return doAltMeasureShapelet(pix,&psf,b,order,order2,pixScale,bCov); }

bool Ellipse::altMeasureShapelet(
    const std::vector<PixelList>& pix, BVec& b, int order, int order2,
    double pixScale, DMatrix* bCov) const
{ return doAltMeasureShapelet(pix,0,b,order,order2,pixScale,bCov); }

