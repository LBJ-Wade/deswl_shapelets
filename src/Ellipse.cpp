
#include "Ellipse.h"
#include <cmath>
#include "dbg.h"
#include <fstream>
#include "PsiHelper.h"
#include "Params.h"

// This is the Miralda-Escude formula for adding distortions.
// We convert to/from distortion to use their formula.
std::complex<double> addShears(
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
    //xdbg<<"GammaAdd: "<<g1<<" + "<<g2<<" = "<<g3<<std::endl;
    return g3;
}

void CalculateShift(
    std::complex<double> c1, std::complex<double> g1, std::complex<double> m1,
    std::complex<double> c2, std::complex<double> g2, std::complex<double> m2,
    std::complex<double>& c3, std::complex<double>& g3, 
    std::complex<double>& m3)
{
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
    // First, just take the terms with z or z*:
    //
    // z'' = exp(-m2)/sqrt(1-|g2|^2)/sqrt(1-|g1|^2) * 
    //       ( exp(-m1) (z - g1 z*) - g2 exp(-m1*) (z* - g1* z) )
    //     = exp(-m2)/sqrt(1-|g2|^2)/sqrt(1-|g1|^2) * 
    //       ( (exp(-m1) + exp(-m1*) g1* g2) z - 
    //         (exp(-m1) g1 + exp(-m1*) g2) z* )
    //     = exp(-m1-m2)/sqrt(1-|g2|^2)/sqrt(1-|g1|^2)
    //       ( (1 + R g1* g2) z - (g1 + R g2) z* )
    // where R == exp(2i imag(m1))
    //
    // So, 
    //
    // g3 = (g1 + R g2) / (1 + R g1* g2)
    // exp(-m3) = exp(-m1-m2) (1+R g1* g2) sqrt(1-|g3|^2)/
    //                   sqrt(1-|g1|^2) / sqrt(1-|g2|^2)
    // 
    // The g terms simplify (after a bit of algebra):
    // exp(-m3) = exp(-m1-m2) (1+R g1* g2)/|1+R g1* g2|
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

    dbg<<"Start CalculateShift\n";
    dbg<<"c1,g1,m1 = "<<c1<<"  "<<g1<<"  "<<m1<<std::endl;
    dbg<<"c2,g2,m2 = "<<c2<<"  "<<g2<<"  "<<m2<<std::endl;

    std::complex<double> Rg2 = std::polar(1.,2*imag(m1)) * g2;
    double normg1 = norm(g1);
    double normg2 = norm(g2);
    std::complex<double> denom = 1.+conj(g1)*Rg2;
    g3 = (g1+Rg2) / denom;
    dbg<<"g3 = "<<g3<<std::endl;

    //xdbg<<"Miralda-Escude formula gives g3 = "<<addShears(g1,g2)<<std::endl;

    m3 = m1 + m2 - std::complex<double>(0.,1.)*arg(denom);
    dbg<<"m3 = "<<m3<<std::endl;

    std::complex<double> X = -c1 - g1 * conj(-c1);
    X = exp(-m1)/sqrt(1.-normg1) * X - c2;
    X = exp(-m2)/sqrt(1.-normg2) * (X - g2 * conj(X));
    double normg3 = norm(g3);
    X /= -exp(-m3)/sqrt(1.-normg3);
    c3 = (X + g3*conj(X)) / (1.-normg3);
    dbg<<"c3 = "<<c3<<std::endl;
}

void Ellipse::preShiftBy(
    const std::complex<double>& c1,
    const std::complex<double>& g1,
    const std::complex<double>& m1)
{
    // Current values are c2, g2, m2.
    // The model is that a round galaxy is first transformed by (c1,g1,m1).
    // And then by the current values of (c2,g2,m2).
    CalculateShift(c1,g1,m1,_cen,_gamma,_mu,_cen,_gamma,_mu); 
}

void Ellipse::postShiftBy(
    const std::complex<double>& c2,
    const std::complex<double>& g2,
    const std::complex<double>& m2)
{
    // Current values are c1, g1, m1.
    // The model is that a round galaxy is first transformed by (c1,g1,m1).
    // And then by the provided values of (c2,g2,m2).
    CalculateShift(_cen,_gamma,_mu,c2,g2,m2,_cen,_gamma,_mu); 
}

void Ellipse::removeRotation()
{
    // This finds the rotation-free transformation that distorts a 
    // circle to the same final shape as the current transformation does.
    // In otherwords, we move the rotation part from the end of the 
    // transformation to the beginning, since a rotation of a circle 
    // is a null operation.
    // 
    // z' = exp(-m-it)/sqrt(1-|g|^2) ( z-c - g (z-c)* )
    //    = exp(-m)/sqrt(1-|g|^2) (exp(-it)(z-c) - g exp(-it) (z-c)*)
    //    = exp(-m)/sqrt(1-|g|^2) (exp(-it)(z-c) - g exp(-2it) (exp(-it)(z-c))*)
    //
    // So, g -> g exp(-2it)
    //     c -> c exp(-it)
    //     m -> m
    //     t -> 0

    if (imag(_mu) != 0.) { 
        std::complex<double> r = std::polar(1.,-imag(_mu));
        _gamma *= r*r;
        _cen *= r;
        _mu = real(_mu);
    }
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

    std::complex<double> m = exp(-_mu)/sqrt(1.-gsq);

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
            std::complex<double> z2 = m*((z1-_cen) - _gamma*conj(z1-_cen));
            Z(n) = z2 / sigma_obs;
        }
    }
    //xdbg<<"Z = "<<Z<<std::endl;
    //xdbg<<"I = "<<I<<std::endl;
    //xdbg<<"W = "<<W<<std::endl;

    int bsize = (order+1)*(order+2)/2;
    xdbg<<"bsize = "<<bsize<<std::endl;
    int bsize2 = (order2+1)*(order2+2)/2;
    xdbg<<"bsize2 = "<<bsize2<<std::endl;

    DMatrix A(nTot,bsize);
    Assert(nTot >= bsize); // Should have been addressed by calling routine.
    //xdbg<<"A = "<<A<<std::endl;

    if (psf) {
        for(int k=0,n=0,nx;k<nExp;++k,n=nx) {
            xdbg<<"psf = "<<(*psf)[k]<<std::endl;
            int psforder = (*psf)[k].getOrder();
            int newpsforder = std::max(psforder,order2);
            xdbg<<"psforder = "<<psforder<<std::endl;
            xdbg<<"newpsforder = "<<newpsforder<<std::endl;
            const double psfsigma = (*psf)[k].getSigma();
            xdbg<<"sigma = "<<psfsigma<<std::endl;
            BVec newpsf(newpsforder,psfsigma);
            int psfsize = (*psf)[k].size();
            int newpsfsize = newpsf.size();
            xdbg<<"psfsize = "<<psfsize<<std::endl;
            bool setnew = false;
            if (_gamma != 0.) {
                DMatrix S(newpsfsize,psfsize,0.);
                calculateGTransform(_gamma,newpsforder,psforder,S);
                newpsf.vec() = S * (*psf)[k].vec();
                setnew = true;
                xdbg<<"newpsf = "<<newpsf<<std::endl;
            }
            if (real(_mu) != 0.) {
                if (setnew) {
                    DMatrix D(newpsfsize,newpsfsize,0.);
                    calculateMuTransform(real(_mu),newpsforder,D);
                    newpsf.vec() = D * newpsf.vec();
                } else {
                    DMatrix D(newpsfsize,psfsize,0.);
                    calculateMuTransform(real(_mu),newpsforder,psforder,D);
                    newpsf.vec() = D * (*psf)[k].vec();
                    setnew = true;
                }
                newpsf.vec() *= exp(2.*real(_mu));
                xdbg<<"newpsf => "<<newpsf<<std::endl;
            }
            if (imag(_mu) != 0.) {
                if (setnew) {
                    DBandMatrix R(newpsfsize,newpsfsize,1,1,0.);
                    calculateThetaTransform(imag(_mu),newpsforder,R);
                    newpsf.vec() = R * newpsf.vec();
                } else {
                    DBandMatrix R(newpsfsize,psfsize,1,1,0.);
                    calculateThetaTransform(imag(_mu),newpsforder,psforder,R);
                    newpsf.vec() = R * (*psf)[k].vec();
                    setnew = true;
                }
                xdbg<<"newpsf => "<<newpsf<<std::endl;
            }
            DMatrix C(bsize2,bsize,0.);
            if (setnew) {
                calculatePsfConvolve(newpsf,order2,order,b.getSigma(),C);
            } else {
                calculatePsfConvolve((*psf)[k],order2,order,b.getSigma(),C);
            }

            const int nPix = pix[k].size();
            nx = n+nPix;
            DMatrix A1(nPix,bsize2);
            makePsi(A1,Z.TMV_subVector(n,nx),order2,&W);
            TMV_rowRange(A,n,nx) = A1 * C;
        }
    } else {
        makePsi(A,TMV_vview(Z),order,&W);
    }
    const double MAX_CONDITION = 1.e6;
#ifdef USE_TMV
    A.saveDiv();
    // I used to use SVD for this, rather than the QRP decomposition.
    // But QRP is significantly faster, and the condition estimate 
    // produced is good enough for what we need here.  
    // TODO: I need to add a real condition estimator to division methods
    // other than SVD in TMV.  LAPACK has this functionality...
    A.divideUsing(tmv::QRP);
    std::auto_ptr<tmv::MatrixView<double> > A_use(
        new tmv::MatrixView<double>(A.view()));
    A_use->qrpd();
    xdbg<<"R diag = "<<A_use->qrpd().getR().diag()<<std::endl;
    while (A_use->qrpd().isSingular() || 
           DiagMatrixViewOf(A_use->qrpd().getR().diag()).condition() > 
           MAX_CONDITION) {
        dbg<<"Poor condition in MeasureShapelet: \n";
        dbg<<"R diag = "<<A_use->qrpd().getR().diag()<<std::endl;
        xdbg<<"condition = "<<
            DiagMatrixViewOf(A_use->qrpd().getR().diag()).condition()<<
            std::endl;
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
    if (max > MAX_CONDITION * min) {
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
    const int bsize2 = (order2+1)*(order2+2)/2;
    b.vec().setZero();
    double sigma = b.getSigma();

    double gsq = std::norm(_gamma);
    Assert(gsq < 1.);
    std::complex<double> m = exp(-_mu)/sqrt(1.-gsq);

    for(int k=0;k<nExp;++k) {
        const int nPix = pix[k].size();
        DVector I(nPix);
        CDVector Z(nPix);

        double sigma_obs = 
            psf ?
            sqrt(pow(sigma,2)+pow((*psf)[k].getSigma(),2)) :
            sigma;
        dbg<<"sigma_obs = "<<sigma_obs<<std::endl;

        DDiagMatrix normD(bsize2);
        double val = pow(pixScale/sigma_obs,2) * exp(-2.*real(_mu)) / nExp;
        for(int n=0,i=0;n<=order2;++n) {
            for(int p=n,q=0;p>=q;--p,++q,++i) {
                if (p!=q) { normD(i) = val/2.; normD(++i) = val/2.; }
                else normD(i) = val;
            }
        }

        for(int i=0;i<nPix;++i) {
            I(i) = pix[k][i].getFlux();
            std::complex<double> z1 = pix[k][i].getPos();
            std::complex<double> z2 = m*((z1-_cen) - _gamma*conj(z1-_cen));
            Z(i) = z2 / sigma_obs;
        }

        DMatrix A(nPix,bsize2);
        makePsi(A,TMV_vview(Z),order2);

        // b = C^-1 normD At I
        DVector b1 = A.transpose() * I;
        dbg<<"b1 = "<<b1<<std::endl;
        b1 = normD EIGEN_asDiag() * b1;
        dbg<<"b1 => "<<b1<<std::endl;
        if (psf) {
            xdbg<<"psf = "<<(*psf)[k]<<std::endl;
            int psforder = (*psf)[k].getOrder();
            int newpsforder = std::max(psforder,order2);
            xdbg<<"psforder = "<<psforder<<std::endl;
            xdbg<<"newpsforder = "<<newpsforder<<std::endl;
            const double psfsigma = (*psf)[k].getSigma();
            xdbg<<"sigma = "<<psfsigma<<std::endl;
            BVec newpsf(newpsforder,psfsigma);
            int psfsize = (*psf)[k].size();
            int newpsfsize = newpsf.size();
            xdbg<<"psfsize = "<<psfsize<<std::endl;
            bool setnew = false;
            if (_gamma != 0.) {
                DMatrix S(newpsfsize,psfsize,0.);
                calculateGTransform(_gamma,newpsforder,psforder,S);
                newpsf.vec() = S * (*psf)[k].vec();
                setnew = true;
                xdbg<<"newpsf = "<<newpsf<<std::endl;
            }
            if (real(_mu) != 0.) {
                if (setnew) {
                    DMatrix D(newpsfsize,newpsfsize,0.);
                    calculateMuTransform(real(_mu),newpsforder,D);
                    newpsf.vec() = D * newpsf.vec();
                } else {
                    DMatrix D(newpsfsize,psfsize,0.);
                    calculateMuTransform(real(_mu),newpsforder,psforder,D);
                    newpsf.vec() = D * (*psf)[k].vec();
                    setnew = true;
                }
                newpsf.vec() *= exp(2.*real(_mu));
                xdbg<<"newpsf => "<<newpsf<<std::endl;
            }
            if (imag(_mu) != 0.) {
                if (setnew) {
                    DBandMatrix R(newpsfsize,newpsfsize,1,1,0.);
                    calculateThetaTransform(imag(_mu),newpsforder,R);
                    newpsf.vec() = R * newpsf.vec();
                } else {
                    DBandMatrix R(newpsfsize,psfsize,1,1,0.);
                    calculateThetaTransform(imag(_mu),newpsforder,psforder,R);
                    newpsf.vec() = R * (*psf)[k].vec();
                    setnew = true;
                }
                xdbg<<"newpsf => "<<newpsf<<std::endl;
            }

            DMatrix C(bsize2,bsize2,0.);
            calculatePsfConvolve(newpsf,order2,b.getSigma(),C);
#ifdef USE_TMV
            b1 /= C;
#else 
            C.lu().solve(b1,&b1);
#endif
            dbg<<"b1 /= C => "<<b1<<std::endl;
        }
        b.vec() += b1.TMV_subVector(0,bsize);
    }
    if (order < b.getOrder()) {
        xdbg<<"Need to zero the rest of b\n";
        xdbg<<"bsize = "<<bsize<<"  b.size = "<<b.size()<<std::endl;
        // Zero out the rest of the shapelet vector:
        b.vec().TMV_subVector(bsize,b.size()).setZero();
    }
    xdbg<<"Done (alt) measure Shapelet\n";
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

