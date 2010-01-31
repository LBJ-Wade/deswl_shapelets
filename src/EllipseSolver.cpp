
#include "EllipseSolver.h"
#include "TMV.h"
#include "dbg.h"
#include "PsiHelper.h"

//#define JTEST

#ifdef JTEST
#include "TestHelper.h"
#endif

// 
// EllipseSolver Implemetation
//

struct EllipseSolver::ESImpl 
{

    ESImpl(const std::vector<PixelList>& pix,
           int order, double sigma, 
           bool isFixedCen, bool isFixedGamma, bool isFixedMu, 
           bool shouldUseFlux);

    ESImpl(const std::vector<PixelList>& pix,
           const std::vector<BVec>& psf, 
           double fP, int order, double sigma, 
           bool isFixedCen, bool isFixedGamma, bool isFixedMu, 
           bool shouldUseFlux);

    void calculateF(const tmv::Vector<double>& x, tmv::Vector<double>& f) const;
    void calculateJ(
        const tmv::Vector<double>& x, const tmv::Vector<double>& f,
        tmv::Matrix<double>& J) const;

    void doF(const tmv::Vector<double>& x, tmv::Vector<double>& f) const;
    void doJ(
        const tmv::Vector<double>& x, const tmv::Vector<double>& f,
        tmv::Matrix<double>& J) const;

    size_t nsize, np2size;
    mutable BVec b;
    const std::vector<PixelList>& pix;
    const std::vector<BVec>* psf;
    double f_psf;
    std::vector<double> sigma_obs;
    tmv::Vector<double> I;
    tmv::Vector<std::complex<double> > Z;
    mutable tmv::Vector<std::complex<double> > Z1;
    tmv::DiagMatrix<double> W;
    mutable tmv::Matrix<double> psi;
    mutable tmv::Matrix<double> A_aux;
    mutable tmv::Matrix<double> AC;
    tmv::MatrixView<double> A;
    mutable std::vector<tmv::Matrix<double> > C;
    mutable std::vector<tmv::Matrix<double> > S;
    mutable std::vector<tmv::Matrix<double> > D;
    bool fixcen, fixgam, fixmu, useflux, numeric_j;
    tmv::Matrix<double> U;
    mutable tmv::Vector<double> xinit;
    mutable tmv::Vector<double> xx;
    mutable tmv::Vector<double> ff;
    mutable tmv::Matrix<double> jj;
    mutable tmv::Vector<double> x_short;
    mutable tmv::Vector<double> f_short;
    mutable tmv::Vector<double> Iresid;
    mutable tmv::Vector<double> Cb;
    mutable tmv::Matrix<double> Gb;
    mutable tmv::Matrix<double> dbdE;
    mutable tmv::Matrix<double> dAdEb;
    mutable tmv::Matrix<double> temp;
    mutable tmv::Vector<double> AtI;
    mutable tmv::Matrix<double> dCdE;
    mutable double fixuc,fixvc,fixg1,fixg2,fixm,flux;
    int maxpsforder;
    tmv::Matrix<double> _Gx;
    tmv::Matrix<double> _Gy;
    tmv::Matrix<double> _Gg1;
    tmv::Matrix<double> _Gg2;
    tmv::Matrix<double> _Gth;
    tmv::Matrix<double> _Gmu;
    mutable tmv::Matrix<double> _dTdE;
};

EllipseSolver::EllipseSolver(
    const std::vector<PixelList>& pix,
    int order, double sigma, 
    bool fixcen, bool fixgam, bool fixmu, bool useflux) :
    _pimpl(new ESImpl(pix,order,sigma,fixcen,fixgam,fixmu,useflux))
{}

EllipseSolver::EllipseSolver(
    const std::vector<PixelList>& pix,
    const std::vector<BVec>& psf, double fp,
    int order, double sigma, 
    bool fixcen, bool fixgam, bool fixmu, bool useflux) :
    _pimpl(new ESImpl(pix,psf,fp,order,sigma,fixcen,fixgam,fixmu,useflux))
{}

EllipseSolver::~EllipseSolver()
{ delete _pimpl; }

void EllipseSolver::calculateF(
    const tmv::Vector<double>& x, tmv::Vector<double>& f) const
{ _pimpl->calculateF(x,f); }

void EllipseSolver::calculateJ(
    const tmv::Vector<double>& x, const tmv::Vector<double>& f,
    tmv::Matrix<double>& J) const
{ 
    if (_pimpl->numeric_j) NLSolver::calculateJ(x,f,J);
    else _pimpl->calculateJ(x,f,J); 
}

static size_t calculateSumSize(const std::vector<PixelList>& v)
{
    size_t sum = 0;
    dbg<<"nPix = ";
    for(size_t i=0;i<v.size();++i) {
        dbg<< (i==0 ? "" : " + ") << v.size();
        sum += v[i].size();
    }
    dbg<<" = "<<sum<<std::endl;
    return sum;
}

EllipseSolver::ESImpl::ESImpl(
    const std::vector<PixelList>& _pix, 
    int order, double _sigma, 
    bool _fixcen, bool _fixgam, bool _fixmu, bool _useflux) :
    nsize((order+1)*(order+2)/2),
    np2size((order+3)*(order+4)/2), b(order,_sigma),
    pix(_pix), psf(0), f_psf(1.), sigma_obs(pix.size(),_sigma),
    I(calculateSumSize(pix)), Z(I.size()), Z1(I.size()), W(I.size()), 
    psi(I.size(),nsize,0.), A_aux(I.size(),np2size,0.), AC(0,0),
    A(A_aux.colRange(0,nsize)),
    fixcen(_fixcen), fixgam(_fixgam), fixmu(_fixmu), useflux(_useflux),
    numeric_j(false), U((fixcen?0:2)+(fixgam?0:2)+(fixmu?0:1),5,0.), 
    xinit(5,0.), xx(5), ff(6), jj(6,5),
    x_short(U.colsize()), f_short(U.colsize()+(useflux?1:0)),
    Iresid(I.size()), Cb(0), Gb(np2size,5),
    dbdE(nsize,5), dAdEb(I.size(),5), temp(nsize,5), AtI(np2size),
    dCdE(0,0), flux(0), maxpsforder(0),
    _Gx(np2size,nsize), _Gy(np2size,nsize), 
    _Gg1(np2size,nsize), _Gg2(np2size,nsize), 
    _Gth(np2size,nsize), _Gmu(np2size,nsize),
    _dTdE(0,0)
{
    size_t n=0;
    for(size_t k=0;k<pix.size();++k) {
        for(size_t i=0;i<pix[k].size();++i,++n) {
            I(n) = pix[k][i].getFlux()*pix[k][i].getInverseSigma();
            W(n) = pix[k][i].getInverseSigma();
            Z(n) = pix[k][i].getPos();
        }
    }
    Assert(n == I.size());

    A.divideUsing(tmv::QR);
    A.saveDiv();
    Assert(A.colsize() >= A.rowsize());
    // Otherwise there are too few pixels to constrain model.

    int k=0;
    if (!fixcen) { U(k++,0) = 1.; U(k++,1) = 1.; }
    if (!fixgam) { U(k++,2) = 1.; U(k++,3) = 1.; }
    if (!fixmu) { U(k,4) = 1.; }

    setupGx(_Gx,order+2,order);
    setupGy(_Gy,order+2,order);
    setupGg1(_Gg1,order+2,order);
    setupGg2(_Gg2,order+2,order);
    setupGth(_Gth,order+2,order);
    setupGmu(_Gmu,order+2,order);
}

static int calculateMaxPsfOrder(const std::vector<BVec>& psf)
{
    int maxpsforder = 0;
    for(size_t k=0;k<psf.size();++k) 
        if (psf[k].getOrder() > maxpsforder)
            maxpsforder = psf[k].getOrder();
    return maxpsforder;
}

#define ORDER_1 (std::max(maxpsforder,order))
#define ORDER_2 (std::max(maxpsforder,order+2))
#define ORDER_3 (std::max(maxpsforder+2,order))
#define SIZE_1 (ORDER_1+1)*(ORDER_1+2)/2
#define SIZE_1b (ORDER_1+3)*(ORDER_1+4)/2
#define SIZE_2 (ORDER_2+1)*(ORDER_2+2)/2
#define SIZE_3 (ORDER_3+1)*(ORDER_3+2)/2
#define SIZE_4 (maxpsforder+1)*(maxpsforder+2)/2
#define SIZE_5 (maxpsforder+3)*(maxpsforder+4)/2
EllipseSolver::ESImpl::ESImpl(
    const std::vector<PixelList>& _pix, 
    const std::vector<BVec>& _psf,
    double _fp, int order, double _sigma, 
    bool _fixcen, bool _fixgam, bool _fixmu, bool _useflux) :
    nsize((order+1)*(order+2)/2),
    np2size((order+3)*(order+4)/2), b(order,_sigma),
    pix(_pix), psf(&_psf), f_psf(_fp), sigma_obs(pix.size()),
    I(calculateSumSize(pix)), Z(I.size()), Z1(I.size()), W(I.size()), 
    psi(I.size(),nsize,0.), A_aux(I.size(),np2size), AC(I.size(),nsize),
    A(A_aux.colRange(0,nsize)), 
    fixcen(_fixcen), fixgam(_fixgam), fixmu(_fixmu), useflux(_useflux),
    numeric_j(false), U((fixcen?0:2)+(fixgam?0:2)+(fixmu?0:1),5,0.),
    xinit(5,0.), xx(5), ff(6), jj(6,5),
    x_short(U.colsize()), f_short(U.colsize()+(useflux?1:0)),
    Iresid(I.size()), Cb(nsize), Gb(np2size,5),
    dbdE(nsize,5), dAdEb(I.size(),5), temp(nsize,5), AtI(np2size),
    dCdE(nsize,nsize), flux(0), maxpsforder(calculateMaxPsfOrder(_psf)),
    _Gx(np2size,nsize), _Gy(np2size,nsize),
    _Gg1(SIZE_1b,SIZE_1), _Gg2(SIZE_1b,SIZE_1),
    _Gth(SIZE_1b,SIZE_1), _Gmu(SIZE_2,SIZE_3),
    _dTdE(SIZE_4,SIZE_4)
{
    for(size_t k=0,n=0;k<pix.size();++k) {
        for(size_t i=0;i<pix[k].size();++i,++n) {
            I(n) = pix[k][i].getFlux()*pix[k][i].getInverseSigma();
            W(n) = pix[k][i].getInverseSigma();
            Z(n) = pix[k][i].getPos();
        }
    }

    int maxpsforder = 0;
    C.reserve(pix.size());
    S.reserve(pix.size());
    D.reserve(pix.size());
    for(size_t k=0;k<pix.size();++k) {
        int psfsize = (*psf)[k].size();
        int psfsize2 = ((*psf)[k].getOrder()+3)*((*psf)[k].getOrder()+4)/2;
        C.push_back(tmv::Matrix<double>(nsize,nsize));
        S.push_back(tmv::Matrix<double>(psfsize,psfsize2));
        D.push_back(tmv::Matrix<double>(psfsize2,psfsize));
        if ((*psf)[k].getOrder() > maxpsforder)
            maxpsforder = (*psf)[k].getOrder();
        sigma_obs[k] = sqrt(pow(b.getSigma(),2)+f_psf*pow((*psf)[k].getSigma(),2));
    }

    AC.divideUsing(tmv::QR);
    AC.saveDiv();
    dbg<<"AC.colsize = nPixels = "<<I.size()<<" = "<<AC.colsize()<<std::endl;
    dbg<<"AC.rowsize = nModel("<<order<<") = "<<nsize<<" = "<<AC.rowsize()<<std::endl;
    Assert(AC.colsize() >= AC.rowsize());
    // Otherwise there are too few pixels to constrain model.

    int k=0;
    if (!fixcen) { U(k++,0) = 1.; U(k++,1) = 1.; }
    if (!fixgam) { U(k++,2) = 1.; U(k++,3) = 1.; }
    if (!fixmu) { U(k,4) = 1.; }

    setupGx(_Gx,order+2,order);
    setupGy(_Gy,order+2,order);
    setupGg1(_Gg1,ORDER_1+2,ORDER_1);
    setupGg2(_Gg2,ORDER_1+2,ORDER_1);
    setupGth(_Gth,ORDER_1+2,ORDER_1);
    setupGmu(_Gmu,ORDER_2,ORDER_3);
}
#undef ORDER_1
#undef ORDER_2
#undef ORDER_3
#undef SIZE_1
#undef SIZE_1b
#undef SIZE_2
#undef SIZE_3

void EllipseSolver::ESImpl::doF(
    const tmv::Vector<double>& x, tmv::Vector<double>& f) const
{
    Assert(x.size() == 5);
    Assert(f.size() == 6);

    // x(0),x(1) = u_cen, v_cen
    // x(2),x(3) = gamma_1, gamma_2
    // x(4) = mu
    //
    // ( u )' = exp(-mu)/sqrt(1-gsq) ( 1-g1  -g2  ) ( u-uc )
    // ( v )                         ( -g2   1+g1 ) ( v-vc )
    // 
    // z' = u' + I v' 
    //    = exp(-mu)/sqrt(1-gsq)
    //             [ (1-g1)(u-uc)-g2(v-vc)-Ig2(u-uc)+I(1+g1)(v-vc) ]
    //    = exp(-mu)/sqrt(1-gsq) [ z-zc - g1(z-zc)* -Ig2(z-zc)* ]
    //    = exp(-mu)/sqrt(1-gsq) ( z-zc - g (z-zc)* )

    // Guard against dangerous z,g,mu values that you can sometimes get
    // from the NLSolver.
    // Typically, these large values end up with overflows, which 
    // lead to nans and/or inf's in the resulting f.
    // Use 2 x the previous value;
    if (NormInf(x-xinit) > 4.) {
        f = b.vec().subVector(0,6);
        if (flux == 0.) flux = b(0);
        f /= flux;
        if (useflux) f(0) -= 1.;
        f *= 2.;
        return; 
    }

    std::complex<double> zc(x[0],x[1]);
    std::complex<double> g(x[2],x[3]);
    double gsq = std::norm(g);
    double mu = x[4];
    // Also guard against gsq > 1, since it leads to nan's with the sqrt(1-gsq)
    // factor below.
    if (gsq > 0.99 || (mu < -3. && Norm(x-xinit) > 0.3)) {
        f = b.vec().subVector(0,6);
        if (flux == 0.) flux = b(0);
        f /= flux;
        if (useflux) f(0) -= 1.;
        f *= 2.;
        return; 
    }

    double m0 = exp(-mu)/sqrt(1.-gsq);
    // z' = m*(z-zc) - m*g*conj(z-zc)
    //    = m*z - m*g*conj(z) - m*zc + m*g*conj(zc);
    Z1 = m0*Z;
    Z1 -= m0*g*Z.conjugate();
    Z1.addToAll(m0*(g*std::conj(zc)-zc));
    for(size_t k=0,n=0,nx;k<pix.size();++k,n=nx) {
        nx = n+pix[k].size();
        Z1.subVector(n,nx) /= sigma_obs[k];
    }
    makePsi(A,Z1,b.getOrder(),&W);

    if (psf) {
        for(size_t k=0,n=0,nx;k<pix.size();++k,n=nx) {
            nx = n+pix[k].size();
            BVec newpsf = (*psf)[k];
            const int psfsize = newpsf.size();
            tmv::MatrixView<double> S1 = S[k].colRange(0,psfsize);
            calculateGTransform(g,newpsf.getOrder(),S1);
            newpsf.vec() = S1 * (*psf)[k].vec();
            tmv::MatrixView<double> D1 = D[k].rowRange(0,psfsize);
            calculateMuTransform(mu,newpsf.getOrder(),D1);
            newpsf.vec() = D1 * newpsf.vec();
            calculatePsfConvolve(newpsf,b.getOrder(),b.getSigma(),C[k]);
            AC.rowRange(n,nx) = A.rowRange(n,nx) * C[k];
        }
        AC.resetDiv();
        if (AC.isSingular()) {
            // I used to keep going and just give an error. But these never 
            // eventually work, so just bail on the whole calculation.
            dbg<<"Found singular AC.  Bail.\n";
            throw tmv::Singular();
        }
        b.vec() = I/AC;
    } else {
        A.resetDiv();
        if (A.isSingular()) {
            dbg<<"Found singular A.  Bail.\n";
            throw tmv::Singular();
        }
        b.vec() = I/A;
    }

    f = b.vec().subVector(0,6);
    if (flux == 0.) flux = b(0);
    f /= flux;
    if (useflux) f(0) -= 1.;
}

void EllipseSolver::ESImpl::calculateF(
    const tmv::Vector<double>& x, tmv::Vector<double>& f) const
{
    Assert(x.size() == U.colsize());
    if (useflux) Assert(f.size() == U.colsize()+1);
    else Assert(f.size() == U.colsize());

    Assert(xx.size() == U.rowsize());
    xx = x*U;
    if (fixcen) { xx[0] = fixuc;  xx[1] = fixvc; }
    if (fixgam) { xx[2] = fixg1;  xx[3] = fixg2; }
    if (fixmu) { xx[4] = fixm; }

    doF(xx,ff);

    if (useflux) {
        f(0) = ff(0);
        f.subVector(1,f.size()) = U*ff.subVector(1,6);
    } else {
        f = U*ff.subVector(1,6);
    }
}

void EllipseSolver::ESImpl::doJ(
    const tmv::Vector<double>& x, const tmv::Vector<double>& ,
    tmv::Matrix<double>& J) const
{
    // A b = I
    // (A + dA/dE dE) (b + db) = I
    // A db + (dA/dE) dE b = I - A b
    // db/dE = -A^-1 (dA/dE) b
    //       = -A^-1 A_aux G b
    //
    // Well, it turns out that this is really an approximation,
    // since it relies on A b = I, which isn't exactly true.
    // b is solved in this equation in a maximum likelihood sense.
    // If Ab is not a perfect model of the intensity pattern,
    // then Ab != I.
    //
    // The correct equation is:
    //
    // At A b = At I
    //
    // Solving this for db/dE:
    //
    // (A + dA/dE dE)t (A + dA/dE dE) (b + db) = (A + dA/dE dE)t I
    // (dA/dE)t A b dE + At (dA/dE) b dE + At A db = (dA/dE)t I dE
    // (At A) db/dE = (dA/dE)t I - (dA/dE)t A b - At (dA/dE) b
    //
    // Now use A = QR and dA/dE = A_aux G to "simplify":
    //
    // Rt Qt Q R (db/dE) = (dA/dE)t (I-Ab) - Rt Qt (dA/dE) b
    // db/dE = R^-1 Rt^-1 (dA/dE)t (I-Ab) - R^-1 Qt (dA/dE) b
    //       = R^-1 (Rt^-1 Gt A_auxt (I-Ab) - Qt A_aux G b)
    //
    // When there is a psf to deal with, the basic equation
    // changes to:
    // A C b = I
    // Ct At A C b = Ct At I
    // 
    // Most of the math follows the same as before with the
    // exception of the calculation of dA/dE, which is now dAC/dE.
    // 
    // d(AC)/dE = dA/dE C + A dC/dE = Aaux G C + A dC/dE
    // 
    // C_pq^st = Sum_uv C_pq^stuv b_psf_uv
    // (dC/dE)_pq^st = Sum_uv C_pq^stuv G_uv^u'v' b_psf_u'v'
    // So dC/dE is obtained from the calculatePsfConvolve routine,
    // giving it the "psf" vector Gb_psf instead of b_psf.
    //
    // Alternatively, we can derive another version of dC/dE
    // by considering the convolution in two different frames:
    //
    // Viewed as transforming the psf by transformation T:
    // b_obs = C(T bpsf) b_i
    //
    // Or viewed as transformint the galaxy by T^-1
    // T^-1 b_obs = C(bpsf) T^-1 b_i
    // b_obs = T C(bpsf) T^-1 b_i
    //
    // Thus, C(T bpsf) = T C(bpsf) T^-1
    // C + dC = (T + dT/dE dE) C0 (T + dT/dE dE)^-1
    //  dC/dE = dT/dE C0 T^-1 - T C0 T^-1 dT/dE T^-1
    //        = dT/dE T^-1 [TC0T^-1] - [TC0T^-1] dT/dE T^-1
    //        = G T^-1 C - C G T^-1

    Assert(x.size() == 5);
    Assert(J.rowsize() == 5);
    Assert(J.colsize() == 6);

    std::complex<double> zc(x[0],x[1]);
    std::complex<double> g(x[2],x[3]);
    double fact = 1./(1.-std::norm(g)); // used below
    double mu = x[4];
    double m0 = exp(-mu)*sqrt(fact);

    if (psf) Iresid = I-AC*b.vec();
    else Iresid = I-A*b.vec();

    augmentPsi(A_aux,Z1,b.getOrder());

#ifdef JTEST
    tmv::Matrix<double> A_aux2(Z1.size(),np2size);
    makePsi(A_aux2.view(),Z1,b.getOrder()+2,&W);
    test(Norm(A_aux2-A_aux) < 1.e-8*Norm(A_aux),"A_aux");
#endif

    double g1 = std::real(g);
    double g2 = std::imag(g);

    tmv::ConstMatrixView<double> Gx = _Gx.subMatrix(0,np2size,0,nsize);
    tmv::ConstMatrixView<double> Gy = _Gy.subMatrix(0,np2size,0,nsize);
    tmv::ConstMatrixView<double> Gg1 = _Gg1.subMatrix(0,np2size,0,nsize);
    tmv::ConstMatrixView<double> Gg2 = _Gg2.subMatrix(0,np2size,0,nsize);
    tmv::ConstMatrixView<double> Gmu = _Gmu.subMatrix(0,np2size,0,nsize);
    tmv::ConstMatrixView<double> Gth = _Gth.subMatrix(0,np2size,0,nsize);

#ifdef JTEST
    // This section was only used for debugging purposes, but
    // I'm leaving it in, since the derivations of Gx[i] are
    // helpful in understanding some of the code following this section.
    tmv::Vector<double> f0(6);
    doF(x,f0);
    tmv::Matrix<double> J_num(6,5);
    tmv::Matrix<double> d2f_num(6,5);
    tmv::Vector<double> f1(6);
    tmv::Vector<double> f2(6);
    double dE = 1.e-8;
    for(int j=0;j<5;++j) {
        tmv::Vector<double> x1 = x;
        x1(j) -= dE;
        doF(x1,f1);
        tmv::Vector<double> x2 = x;
        x2(j) += dE;
        doF(x2,f2);
        J_num.col(j) = (f2-f1)/(2.*dE);
        d2f_num.col(j) = (f2+f1-2.*f0)/(dE*dE);
    }
    dbg<<"numerical df: "<<J_num<<std::endl;
    dbg<<"numerical d2f: "<<d2f_num<<std::endl;
    dbg<<"(f2-f0)/dE = "<<(f2-f0)/dE<<std::endl;
    dbg<<"(f0-f1)/dE = "<<(f0-f1)/dE<<std::endl;
    dbg<<"diff = "<<(f2-f0-f0+f1)/dE<<std::endl;
    doF(x,f0); // return A to original value

    tmv::Matrix<double> dAdEb_num(I.size(),5);
    tmv::Matrix<double> dAdEtI_num(nsize,5);
    for(int j=0;j<5;++j) {
        dbg<<"j = "<<j<<std::endl;
        tmv::Vector<double> x1 = x;
        x1[j] -= dE;
        std::complex<double> zc_1(x1[0],x1[1]);
        std::complex<double> g_1(x1[2],x1[3]);
        double mu_1 = x1[4];
        tmv::Vector<double> x2 = x;
        x2[j] += dE;
        std::complex<double> zc_2(x2[0],x2[1]);
        std::complex<double> g_2(x2[2],x2[3]);
        double mu_2 = x2[4];
        tmv::Matrix<double> A0(I.size(),nsize);
        tmv::Vector<std::complex<double> > z0(I.size());
        tmv::Matrix<double> A1(I.size(),nsize);
        tmv::Vector<std::complex<double> > z1(I.size());
        tmv::Matrix<double> A2(I.size(),nsize);
        tmv::Vector<std::complex<double> > z2(I.size());
        for(size_t k=0,n=0,nx;k<pix.size();++k,n=nx) {
            nx = n+pix[k].size();
            double m = m0/sigma_obs[k];
            z0.subVector(n,nx) = m*Z.subVector(n,nx);
            z0.subVector(n,nx) -= m*g*Z.conjugate().subVector(n,nx);
            z0.subVector(n,nx).addToAll(m*(g*std::conj(zc)-zc));
            double m_1 = exp(-mu_1)/sigma_obs[k]/sqrt(1.-std::norm(g_1));
            z1.subVector(n,nx) = m_1*Z.subVector(n,nx);
            z1.subVector(n,nx) -= m_1*g_1*Z.conjugate().subVector(n,nx);
            z1.subVector(n,nx).addToAll(m_1*(g_1*std::conj(zc_1)-zc_1));
            double m_2 = exp(-mu_2)/sigma_obs[k]/sqrt(1.-std::norm(g_2));
            z2.subVector(n,nx) = m_2*Z.subVector(n,nx);
            z2.subVector(n,nx) -= m_2*g_2*Z.conjugate().subVector(n,nx);
            z2.subVector(n,nx).addToAll(m_2*(g_2*std::conj(zc_2)-zc_2));
        }
        makePsi(A0.view(),z0,b.getOrder(),&W);
        makePsi(A1.view(),z1,b.getOrder(),&W);
        makePsi(A2.view(),z2,b.getOrder(),&W);
        tmv::Matrix<double> dAdE_num = (A2-A1)/(2.*dE);
        tmv::Matrix<double> d2A_num = (A2+A1-2.*A0)/(dE*dE);

        if (psf) {
            for(size_t k=0,n=0,nx;k<pix.size();++k,n=nx) {
                nx = n+pix[k].size();
                A0.rowRange(n,nx) *= C[k];

                tmv::Matrix<double> C1(nsize,nsize);
                BVec newpsf1 = (*psf)[k];
                ApplyG(g_1,newpsf1);
                ApplyMu(mu_1,newpsf1);
                calculatePsfConvolve(newpsf1,b.getOrder(),b.getSigma(),C1);
                A1.rowRange(n,nx) *= C1;

                tmv::Matrix<double> C2(nsize,nsize);
                BVec newpsf2 = (*psf)[k].vec();
                ApplyG(g_2,newpsf2);
                ApplyMu(mu_2,newpsf2);
                calculatePsfConvolve(newpsf2,b.getOrder(),b.getSigma(),C2);
                A2.rowRange(n,nx) *= C2;
            }
        }

        tmv::Matrix<double> dACdE_num = (A2-A1)/(2.*dE);
        tmv::Matrix<double> d2AC_num = (A2+A1-2.*A0)/(dE*dE);

        tmv::Matrix<double> dACdE(I.size(),nsize);
        tmv::Matrix<double> dAdE(I.size(),nsize);

        for(size_t k=0,n=0,nx;k<pix.size();++k,n=nx) {
            nx = n+pix[k].size();
            double m = m0/sigma_obs[k];

            std::vector<tmv::Matrix<double> > GG(
                5, tmv::Matrix<double>(np2size,nsize));

            // A_aux Gx[i] is the real dA/dE_i
            // given the fact that the u,v used to make the psi are given by:
            // u = exp(-mu)/sqrt(1-gsq) [ (1-g1) (x-xc) - g2 (y-yc) ]
            // v = exp(-mu)/sqrt(1-gsq) [ -g2 (x-xc) + (1+g1) (y-yc) ]

            // dA/dxc = -exp(-mu)/sqrt(1-gsq) [ dA/du (1-g1) + dA/dv (-g2) ]
            GG[0] = -m*((1.-g1)*Gx - g2*Gy);
            // dA/dyc = -exp(-mu)/sqrt(1-gsq) [ dA/du (-g2) + dA/dv (1+g1) ]
            GG[1] = -m*((1.+g1)*Gy - g2*Gx);
            // dA/dg1 = -exp(-mu)/sqrt(1-gsq) [ dA/du (x-xc) - dA/dv (y-yc) ]
            //             -1/2 exp(-mu) (1-gsq)^-3/2 (-2g1) [...]
            //   Use inverse transformation to simplify:
            //   exp(-mu)/sqrt(1-gsq) (x-xc) = ((1+g1) u + g2 v)/(1-gsq)
            //   exp(-mu)/sqrt(1-gsq) (y-yc) = (g2 u + (1-g1) v)/(1-gsq)
            // dA/dg1 = -1/(1-gsq) [ dA/du ((1+g1) u + g2 v - g1 u) - 
            //                         dA/dv (g2 u + (1-g1) v + g1 v ]
            //        = -1/(1-gsq) [ dA/du (u + g2 v) - dA/dv (g2 u + v) ]
            //        = -1/(1-gsq) [ (u dA/du - v dA/dv) - g2 (u dA/dv - v dA/du) ]
            GG[2] = -fact*(Gg1 - g2*Gth);
            // dA/dg2 = -exp(-mu)/sqrt(1-gsq) [ dA/du (y-yc) + dA/dv (x-xc) ]
            //             -1/2 exp(-mu) (1-gsq)^-3/2 (-2g2) [...]
            //        = -1/(1-gsq) [ dA/du (g2 u + (1-g1) v - g2 u) + 
            //                         dA/dv ((1+g1) u + g2 v - g2 v ]
            //        = -1/(1-gsq) [ dA/du (v - g1 v) + dA/dv (u + g1 u) ]
            //        = -1/(1-gsq) [ (u dA/dv + v dA/du) + g1 (u dA/dv - v dA/du) ]
            GG[3] = -fact*(Gg2 + g1*Gth);
            // dA/dmu = -(u dA/du + v dA/dv)
            GG[4] = -Gmu;

            dAdE.rowRange(n,nx) = A_aux.rowRange(n,nx)*GG[j];
            if (n == 0) {
                dbg<<"A_aux = "<<A_aux.rowRange(0,5)<<std::endl;
                dbg<<"GG["<<j<<"] = "<<GG[j]<<std::endl;
                dbg<<"dAdE = "<<dAdE.rowRange(0,5)<<std::endl;
                dbg<<"j = "<<j<<std::endl;
            }
            if (psf) {
                dACdE.rowRange(n,nx) = dAdE.rowRange(n,nx)*C[k];

                if (j>1) {
                    tmv::Matrix<double> C1(nsize,nsize);
                    BVec newpsf1 = (*psf)[k];
                    ApplyG(g_1,newpsf1);
                    ApplyMu(mu_1,newpsf1);
                    calculatePsfConvolve(newpsf1,b.getOrder(),b.getSigma(),C1);

                    tmv::Matrix<double> C2(nsize,nsize);
                    BVec newpsf2 = (*psf)[k];
                    ApplyG(g_2,newpsf2);
                    ApplyMu(mu_2,newpsf2);
                    calculatePsfConvolve(newpsf2,b.getOrder(),b.getSigma(),C2);

                    // For j=0,1 dC/dE = 0, so dA/dE is just AGC
                    // If not, we need to add A(dC/dE)
                    tmv::Matrix<double> dCdE_num = (C2-C1)/(2.*dE);
                    tmv::Matrix<double> d2C_num = (C2+C1-2.*C[k])/(dE*dE);

                    int n1 = ((*psf)[k].getOrder()+1)*((*psf)[k].getOrder()+2)/2;
                    int n2 = ((*psf)[k].getOrder()+3)*((*psf)[k].getOrder()+4)/2;
                    tmv::Matrix<double> S0(n1,n1);
                    tmv::Matrix<double> D0(n1,n1);
                    calculateGTransform(g,(*psf)[k].getOrder(),S0);
                    calculateMuTransform(mu,(*psf)[k].getOrder(),D0);
                    tmv::Matrix<double> T0 = D0*S0;
                    tmv::Matrix<double> S1(n1,n1);
                    tmv::Matrix<double> D1(n1,n1);
                    calculateGTransform(g_1,(*psf)[k].getOrder(),S1);
                    calculateMuTransform(mu_1,(*psf)[k].getOrder(),D1);
                    tmv::Matrix<double> T1 = D1*S1;
                    tmv::Matrix<double> D2(n1,n1);
                    tmv::Matrix<double> S2(n1,n1);
                    calculateGTransform(g_2,(*psf)[k].getOrder(),S2);
                    calculateMuTransform(mu_2,(*psf)[k].getOrder(),D2);
                    tmv::Matrix<double> T2 = D2*S2;
                    tmv::Matrix<double> dTdE_num = (T2-T1)/(2.*dE);
                    tmv::Matrix<double> d2T_num = (T2+T1-2.*T0)/(dE*dE);

                    tmv::Matrix<double> dTdE(n1,n1);
                    if (j==2) {
                        tmv::Matrix<double> Sx(n2,n2);
                        calculateGTransform(g,(*psf)[k].getOrder()+2,Sx);
                        tmv::Matrix<double> dSdE_num = (S2-S1)/(2.*dE);
                        tmv::Matrix<double> d2S_num = (S2+S1-2.*S0)/(dE*dE);
                        dbg<<"Gg1 size = "<<_Gg1.colsize()<<"  "<<
                            _Gg1.rowsize()<<std::endl;
                        dbg<<"Gth size = "<<_Gth.colsize()<<"  "<<
                            _Gth.rowsize()<<std::endl;
                        dbg<<"n1 = "<<n1<<", n2 = "<<n2<<std::endl;
                        tmv::Matrix<double> dSdE = 
                            (fact * Sx.rowRange(0,n1) * 
                             _Gg1.subMatrix(0,n2,0,n1)) +
                            (fact * g2 * Sx.rowRange(0,n1) * 
                             _Gth.subMatrix(0,n2,0,n1));
                        dbg<<"Norm(dS/dE - numeric dS/dE) = "<<
                            Norm(dSdE-dSdE_num)<<std::endl;
                        test(Norm(dSdE-dSdE_num) < 10.*dE*Norm(d2S_num),"dSdE");
                        dTdE = 
                            (fact * D[k].rowRange(0,n1) * Sx.rowRange(0,n1) * 
                             _Gg1.subMatrix(0,n2,0,n1)) +
                            (fact * g2 * D[k].rowRange(0,n1) * Sx.rowRange(0,n1) * 
                             _Gth.subMatrix(0,n2,0,n1));
                    } else if (j==3) {
                        tmv::Matrix<double> Sx(n2,n2);
                        calculateGTransform(g,(*psf)[k].getOrder()+2,Sx);
                        tmv::Matrix<double> dSdE_num = (S2-S1)/(2.*dE);
                        tmv::Matrix<double> d2S_num = (S2+S1-2.*S0)/(dE*dE);
                        tmv::Matrix<double> dSdE = 
                            (fact * Sx.rowRange(0,n1) * 
                             _Gg2.subMatrix(0,n2,0,n1)) -
                            (fact * g1 * Sx.rowRange(0,n1) * 
                             _Gth.subMatrix(0,n2,0,n1));
                        dbg<<"Norm(dS/dE - numeric dS/dE) = "<<
                            Norm(dSdE-dSdE_num)<<std::endl;
                        test(Norm(dSdE-dSdE_num) < 10.*dE*Norm(d2S_num),"dSdE");
                        dTdE = 
                            (fact * D[k].rowRange(0,n1) * Sx.rowRange(0,n1) * 
                             Gg2.subMatrix(0,n2,0,n1)) -
                            (fact * g1 * D[k].rowRange(0,n1) * Sx.rowRange(0,n1) * 
                             Gth.subMatrix(0,n2,0,n1));
                    } else if (j==4) {
                        tmv::Matrix<double> Dx(n2,n2);
                        calculateMuTransform(mu,(*psf)[k].getOrder()+2,Dx);
                        tmv::Matrix<double> dDdE_num = (D2-D1)/(2.*dE);
                        tmv::Matrix<double> d2D_num = (D2+D1-2.*D0)/(dE*dE);
                        tmv::Matrix<double> dDdE = 
                            _Gmu.subMatrix(0,n1,0,n2)*Dx.colRange(0,n1);
                        dbg<<"Norm(dD/dE - numeric dD/dE) = "<<
                            Norm(dDdE-dDdE_num)<<std::endl;
                        test(Norm(dDdE-dDdE_num) < 10.*dE*Norm(d2D_num),"dDdE");
                        dTdE = 
                            (_Gmu.subMatrix(0,n1,0,n2) *
                             Dx.colRange(0,n1) * S[k].colRange(0,n1));
                    }
                    dbg<<"Norm(dT/dE - numeric dT/dE) = "<<
                        Norm(dTdE-dTdE_num)<<std::endl;
                    test(Norm(dTdE-dTdE_num) < 10.*dE*Norm(dTdE_num),"dTdE");

                    BVec dbpsfdE((*psf)[k].getOrder(),(*psf)[k].getSigma());
                    dbpsfdE = dTdE*(*psf)[k].vec();
                    tmv::Matrix<double> dCdE(nsize,nsize);
                    calculatePsfConvolve(dbpsfdE,b.getOrder(),b.getSigma(),dCdE);
                    dbg<<"Norm(dC/dE - numeric dC/dE) = "<<
                        Norm(dCdE-dCdE_num)<<"\n";
                    test(Norm(dCdE-dCdE_num) < 10.*dE*Norm(d2C_num),"dCdE");
                    dACdE.rowRange(n,nx) += A.rowRange(n,nx)*dCdE;
                }
            }
        }
        dbg<<"Norm(dA/dE - numeric dA/dE) = "<<Norm(dAdE-dAdE_num)<<std::endl;
        dbg<<"Norm(dA/dE) = "<<Norm(dAdE)<<std::endl;
        dbg<<"Norm(d2A/dE2) = "<<Norm(d2A_num)<<std::endl;
        if (Norm(dAdE-dAdE_num) > 10.*dE*Norm(d2A_num)) {
            dbg<<"dAdE = "<<dAdE.rowRange(0,5)<<std::endl;
            dbg<<"dAdE_num = "<<dAdE_num.rowRange(0,5)<<std::endl;
            dbg<<"diff = "<<(dAdE-dAdE_num).rowRange(0,5)<<std::endl;
        }
        test(Norm(dAdE-dAdE_num) < 10.*dE*Norm(d2A_num),"dAdE");

        if (psf) {
            dbg<<"Norm(dAC/dE - numeric dAC/dE) = "<<
                Norm(dACdE-dACdE_num)<<std::endl;
            dbg<<"Norm(dAC/dE) = "<<Norm(dACdE)<<std::endl;
            dbg<<"Norm(d2A/dE2) = "<<Norm(d2A_num)<<std::endl;
            if (Norm(dACdE-dACdE_num) > 10.*dE*Norm(d2A_num)) {
                dbg<<"dACdE = "<<dACdE<<std::endl;
                dbg<<"dACdE_num = "<<dACdE_num<<std::endl;
                dbg<<"diff = "<<dACdE-dACdE_num<<std::endl;
            }
            test(Norm(dACdE-dACdE_num) < 10.*dE*Norm(d2AC_num),"dACdE");
            dAdE = dACdE;
        }

        dAdEb_num.col(j) = dAdE * b.vec();
        dAdEtI_num.col(j) = dAdE.transpose() * Iresid;

        // (At A) db/dE = (dA/dE)t (I - A b) - At (dA/dE) b
        tmv::Vector<double> dbdE_calc = dAdE.transpose() * Iresid;
        dbg<<"dAdEt I = "<<dbdE_calc<<std::endl;
        double kappa;
        if (psf) {
            dbdE_calc -= AC.transpose() * dAdE * b.vec();
            dbg<<"(ACt) dAdE b = "<<AC.transpose()*dAdE*b.vec()<<std::endl;
            dbg<<"dbdE_calc => "<<dbdE_calc<<std::endl;
            dbg<<"R.diag() = "<<AC.qrd().getR().diag()<<std::endl;
            dbdE_calc /= AC.qrd().getR().transpose();
            dbg<<"(/= Rt) dbdE_calc => "<<dbdE_calc<<std::endl;
            dbdE_calc /= AC.qrd().getR();
            dbg<<"(/= R) dbdE_calc => "<<dbdE_calc<<std::endl;
            kappa = AC.qrd().getR().DoCondition();
        } else {
            dbdE_calc -= A.transpose() * dAdE * b.vec();
            dbdE_calc /= A.qrd().getR().transpose();
            dbdE_calc /= A.qrd().getR();
            kappa = A.qrd().getR().DoCondition();
        }
        tmv::Vector<double> b0 = I/A0;
        tmv::Vector<double> b1 = I/A1;
        tmv::Vector<double> b2 = I/A2;
        tmv::Vector<double> dbdE_num = (b2-b1)/(2.*dE);
        tmv::Vector<double> d2b_num = (b2+b1-2.*b0)/(dE*dE);
        dbg<<"Norm(db/dE_calc - numeric db/dE) = "<<
            Norm(dbdE_calc-dbdE_num)<<std::endl;
        dbg<<"b0 = "<<b0<<std::endl;
        dbg<<"b1 = "<<b1<<std::endl;
        dbg<<"b2 = "<<b2<<std::endl;
        dbg<<"dbdE_calc = "<<dbdE_calc<<std::endl;
        dbg<<"dbdE_num = "<<dbdE_num<<std::endl;
        dbg<<"diff = "<<dbdE_calc-dbdE_num<<std::endl;
        dbg<<"d2b = "<<d2b_num<<std::endl;
        dbg<<"AtA dbdE_calc = "<<A0.Adjoint()*A0*dbdE_calc<<std::endl;
        dbg<<"AtA dbdE_num = "<<A0.Adjoint()*A0*dbdE_num<<std::endl;
        dbg<<"Norm(b0) = "<<Norm(b0)<<std::endl;
        dbg<<"Norm(d2b) = "<<Norm(d2b_num)<<std::endl;
        test(Norm(dbdE_calc-dbdE_num) < 10.*dE*kappa*kappa*Norm(d2b_num),
             "dbdE_calc");
        J.col(j) = dbdE_calc.subVector(0,6);
    }
    dbg<<"Semi-empirical: "<<J<<std::endl;
    dbg<<"Norm(diff) = "<<Norm(J-J_num)<<std::endl;
#endif

    // The dA/dE = AG model isn't completely accurate.
    // We need to expand out d/dE in terms of things that are
    // d/dx, d/dy, (x d/dx - y d/dy), etc.
    //
    // d/dE = dx/dE d/dx + dy/dE d/dy
    //
    // x = exp(-mu)/(1-gsq)^1/2 * ( (u-uc) - g1*(u-uc) - g2*(v-vc) )
    // y = exp(-mu)/(1-gsq)^1/2 * ( (v-vc) + g1*(v-vc) - g2*(u-uc) )
    //
    // The derivations of the GG matrices in the above JTEST
    // section are what we want.
    // Gb.col(i) = GG[i] * b

    dbdE.setZero();

    for(size_t k=0,n=0,nx;k<pix.size();++k,n=nx) {
        nx = n + pix[k].size();
        double m = m0/sigma_obs[k];
        if (psf) {
            Cb = C[k] * b.vec();

            // GG[0] = -m*((1.-g1)*Gx - g2*Gy);
            Gb.col(0) = -m * (1.-g1) * Gx * Cb;
            Gb.col(0) += m * g2 * Gy * Cb;

            // GG[1] = -m*((1.+g1)*Gy - g2*Gx);
            Gb.col(1) = -m * (1.+g1) * Gy * Cb;
            Gb.col(1) += m * g2 * Gx * Cb;

            // GG[2] = -fact*(Gg1 - g2*Gth);
            Gb.col(2) = -fact * Gg1 * Cb;
            Gb.col(2) += fact * g2 * Gth * Cb;

            // GG[3] = -fact*(Gg2 + g1*Gth);
            Gb.col(3) = -fact * Gg2 * Cb;
            Gb.col(3) -= fact * g1 * Gth * Cb;

            // GG[4] = -Gmu;
            Gb.col(4) = -Gmu * Cb;

            dAdEb.rowRange(n,nx) = A_aux.rowRange(n,nx) * Gb;
        } else {
            // GG[0] = -m*((1.-g1)*Gx - g2*Gy);
            Gb.col(0) = -m * (1.-g1) * Gx * b.vec();
            Gb.col(0) += m * g2 * Gy * b.vec();

            // GG[1] = -m*((1.+g1)*Gy - g2*Gx);
            Gb.col(1) = -m * (1.+g1) * Gy * b.vec();
            Gb.col(1) += m * g2 * Gx * b.vec();

            // GG[2] = -fact*(Gg1 - g2*Gth);
            Gb.col(2) = -fact * Gg1 * b.vec();
            Gb.col(2) += fact * g2 * Gth * b.vec();

            // GG[3] = -fact*(Gg2 + g1*Gth);
            Gb.col(3) = -fact * Gg2 * b.vec();
            Gb.col(3) -= fact * g1 * Gth * b.vec();

            // GG[4] = -Gmu;
            Gb.col(4) = -Gmu * b.vec();

            dAdEb.rowRange(n,nx) = A_aux.rowRange(n,nx) * Gb;
        }

        AtI = A_aux.rowRange(n,nx).transpose()*Iresid.subVector(n,nx);

        // Here, temp = (dAdE)t * I = GtAtI
        // GtAtI.col(i) = GG[i].transpose() * A_aux.transpose() * Iresid
        temp.col(0) = -m*(1.-g1) * Gx.transpose() * AtI;
        temp.col(0) += m*g2 * Gy.transpose() * AtI;

        temp.col(1) = -m*(1.+g1) * Gy.transpose() * AtI;
        temp.col(1) += m*g2 * Gx.transpose() * AtI;

        temp.col(2) = -fact * Gg1.transpose() * AtI;
        temp.col(2) += fact*g2 * Gth.transpose() * AtI;

        temp.col(3) = -fact * Gg2.transpose() * AtI;
        temp.col(3) -= fact*g1 * Gth.transpose() * AtI;

        temp.col(4) = -Gmu.transpose() * AtI;

        if (psf) {
            // with a psf, dAdE = AGC + A(dC/dE)
            // AGCb is already accounted for correctly by doing the Cb thing above.
            // for (dAdE)t, CtGtAtI is effected by the following line:
            size_t psfsize = (*psf)[k].size();

            temp = C[k].transpose() * temp;

            // The rest of this section does the A(dC/dE) terms, which
            // only exist for E = {g1,g2,mu}.  ie. not E = {x,y}.

            int n2 = ((*psf)[k].getOrder()+3)*((*psf)[k].getOrder()+4)/2;
            tmv::MatrixView<double> dTdE = _dTdE.subMatrix(0,psfsize,0,psfsize);
            tmv::ConstMatrixView<double> D1 = D[k].rowRange(0,psfsize);
            tmv::ConstMatrixView<double> S1 = S[k].colRange(0,psfsize);

            augmentGTransformCols(g,(*psf)[k].getOrder(),S[k].view());
            augmentMuTransformRows(mu,(*psf)[k].getOrder(),D[k].view());

#ifdef JTEST
            tmv::Matrix<double> S2(n2,n2);
            calculateGTransform(g,(*psf)[k].getOrder()+2,S2.view());
            dbg<<"S[k] = "<<S[k]<<std::endl;
            dbg<<"S2 = "<<S2.rowRange(0,psfsize)<<std::endl;
            dbg<<"diff = "<<S[k]-S2.rowRange(0,psfsize)<<std::endl;
            dbg<<"Norm(S[k]-S2) = "<<Norm(S[k]-S2.rowRange(0,psfsize))<<std::endl;
            dbg<<"Norm(S2) = "<<Norm(S2)<<std::endl;
            test(Norm(S[k]-S2.rowRange(0,psfsize)) < 1.e-8*Norm(S2),"AugmentG");
            tmv::Matrix<double> D2(n2,n2);
            calculateMuTransform(mu,(*psf)[k].getOrder()+2,D2.view());
            dbg<<"Norm(diff) = "<<Norm(D[k]-D2.colRange(0,psfsize))<<std::endl;
            dbg<<"Norm(D2) = "<<Norm(D2)<<std::endl;
            test(Norm(D[k]-D2.colRange(0,psfsize)) < 1.e-8*Norm(D2),"AugmentMu");
#endif

            // E = g1
            dTdE = fact*D1*S[k]*_Gg1.subMatrix(0,n2,0,psfsize);
            dTdE += fact*g2*D1*S[k]*_Gth.subMatrix(0,n2,0,psfsize);
            BVec dbpsfdE((*psf)[k].getOrder(),(*psf)[k].getSigma());
            dbpsfdE.vec() = dTdE*(*psf)[k].vec();
            calculatePsfConvolve(dbpsfdE,b.getOrder(),b.getSigma(),dCdE);
            dAdEb.col(2,n,nx) += A.rowRange(n,nx)*(dCdE*b.vec());
            temp.col(2) += dCdE.transpose() * AtI.subVector(0,nsize);

            // E = g2
            dTdE = fact*D1*S[k]*_Gg2.subMatrix(0,n2,0,psfsize);
            dTdE -= fact*g1*D1*S[k]*_Gth.subMatrix(0,n2,0,psfsize);
            dbpsfdE.vec() = dTdE*(*psf)[k].vec();
            calculatePsfConvolve(dbpsfdE,b.getOrder(),b.getSigma(),dCdE);
            dAdEb.col(3,n,nx) += A.rowRange(n,nx)*(dCdE*b.vec());
            temp.col(3) += dCdE.transpose() * AtI.subVector(0,nsize);

            // E = mu
            dTdE = _Gmu.subMatrix(0,psfsize,0,n2)*D[k]*S1;
            dbpsfdE.vec() = dTdE*(*psf)[k].vec();
            calculatePsfConvolve(dbpsfdE,b.getOrder(),b.getSigma(),dCdE);
            dAdEb.col(4,n,nx) += A.rowRange(n,nx)*(dCdE*b.vec());
            temp.col(4) += dCdE.transpose() * AtI.subVector(0,nsize);
        }
        dbdE += temp;
    }
#ifdef JTEST
    dbg<<"Norm(dAdEtI-dAdEtI_num) = "<<Norm(dbdE-dAdEtI_num)<<std::endl;
    dbg<<"Norm(dAdEb-dAdEb_num) = "<<Norm(dAdEb-dAdEb_num)<<std::endl;
#endif

    // db/dE = R^-1 (Rt^-1 dAdEt (I-Ab) - Qt dAdE b)
    // So far dbdE is just dAdEt (I-Ab)
    if (psf) {
        dbdE /= AC.qrd().getR().transpose();
        //dbdE -= AC.qrd().getQ().transpose() * dAdEb;
        temp = dAdEb / AC.qrd().getQ();
        dbdE -= temp;
        dbdE /= AC.qrd().getR();
    } else {
        dbdE /= A.qrd().getR().transpose();
        //dbdE -= A.qrd().getQ().transpose() * dAdEb;
        temp = dAdEb / A.qrd().getQ();
        dbdE -= temp;
        dbdE /= A.qrd().getR();
    }
    J = dbdE.rowRange(0,6);
    J /= flux;

#ifdef JTEST
    dbg<<"analytic: "<<J.Clip(1.e-5)<<std::endl;
    dbg<<"Norm(diff) = "<<Norm(J-J_num)<<std::endl;
#endif
}

void EllipseSolver::ESImpl::calculateJ(
    const tmv::Vector<double>& x, const tmv::Vector<double>& f,
    tmv::Matrix<double>& J) const
{
    Assert(x.size() == U.colsize());
    Assert(J.rowsize() == U.colsize());
    if (useflux) {
        Assert(J.colsize() == U.colsize()+1);
    } else {
        Assert(J.colsize() == U.colsize());
    }

    xx = x*U;
    if (fixcen) { xx[0] = fixuc;  xx[1] = fixvc; }
    if (fixgam) { xx[2] = fixg1;  xx[3] = fixg2; }
    if (fixmu) { xx[4] = fixm; }

    doJ(xx,f,jj);
    if (useflux) {
        J.row(0) = jj.row(0)*U.transpose();
        J.rowRange(1,J.colsize()) = U*jj.rowRange(1,6)*U.transpose();
    } else {
        J = U*jj.rowRange(1,6)*U.transpose();
    }
}

void EllipseSolver::useNumericJ() { _pimpl->numeric_j = true; }

const BVec& EllipseSolver::getB() const { return _pimpl->b; }

void EllipseSolver::getBCov(tmv::Matrix<double>& bcov) const 
{
    // It's not completely obvious that the first one is correct.
    // We find b as the solution of:
    // A C b = I
    // This means that 
    // cov(Cb) = (At A)^-1
    // But the propagation of covariance says that 
    // cov(Cb) = C cov(b) Ct
    // So we can solve for cov(b) as:
    // cov(b) = C^-1 cov(Cb) Ct^-1 =
    //        = C^-1 (At A)^-1 Ct^-1
    //        = (Ct At A C)^-1
    //        = ((AC)t (AC))^-1
    // This is what AC.makeInverseATA returns.
    if (_pimpl->psf) _pimpl->AC.makeInverseATA(bcov);
    else _pimpl->A.makeInverseATA(bcov);
}

void EllipseSolver::getCovariance(tmv::Matrix<double>& cov) const 
{
    tmv::Matrix<double> cov1(_pimpl->x_short.size(),_pimpl->x_short.size());
    NLSolver::getCovariance(cov1);
    cov = _pimpl->U.transpose()*cov1*_pimpl->U;
    dbg<<"getCovariance:\n";
    dbg<<"cov1 = "<<cov1<<std::endl;
    dbg<<"full cov = "<<cov<<std::endl;
}

void EllipseSolver::getInverseCovariance(tmv::Matrix<double>& invcov) const 
{
    tmv::Matrix<double> invcov1(_pimpl->x_short.size(),_pimpl->x_short.size());
    NLSolver::getInverseCovariance(invcov1);
    invcov = _pimpl->U.transpose()*invcov1*_pimpl->U;
    dbg<<"getCovariance:\n";
    dbg<<"invcov1 = "<<invcov1<<std::endl;
    dbg<<"full invcov = "<<invcov<<std::endl;
}

void EllipseSolver::callF(
    const tmv::Vector<double>& x, tmv::Vector<double>& f) const
{
    Assert(x.size() == 5);
    Assert(f.size() == 5);
    _pimpl->xinit = x;
    if (_pimpl->fixcen) { _pimpl->fixuc = x[0]; _pimpl->fixvc = x[1]; }
    if (_pimpl->fixgam) { _pimpl->fixg1 = x[2]; _pimpl->fixg2 = x[3]; }
    if (_pimpl->fixmu) { _pimpl->fixm = x[4]; }

    _pimpl->x_short = _pimpl->U * x;

    _pimpl->flux = 0;

    calculateF(_pimpl->x_short,_pimpl->f_short);

    if (_pimpl->useflux) 
        f = _pimpl->f_short.subVector(1,_pimpl->f_short.size()) * _pimpl->U;
    else f = _pimpl->f_short*_pimpl->U;
}

bool EllipseSolver::solve(tmv::Vector<double>& x, tmv::Vector<double>& f) const
{
    Assert(x.size() == 5);
    Assert(f.size() == 5);
    _pimpl->xinit = x;
    if (_pimpl->fixcen) { _pimpl->fixuc = x[0]; _pimpl->fixvc = x[1]; }
    if (_pimpl->fixgam) { _pimpl->fixg1 = x[2]; _pimpl->fixg2 = x[3]; }
    if (_pimpl->fixmu) { _pimpl->fixm = x[4]; }

    _pimpl->x_short = _pimpl->U * x;

    _pimpl->flux = 0;
    bool ret;
    try {
        ret = NLSolver::solve(_pimpl->x_short,_pimpl->f_short);
    } catch (...) {
        ret = false;
    }

    x = _pimpl->x_short*_pimpl->U;
    if (_pimpl->fixcen) { x[0] = _pimpl->fixuc; x[1] = _pimpl->fixvc; }
    if (_pimpl->fixgam) { x[2] = _pimpl->fixg1; x[3] = _pimpl->fixg2; }
    if (_pimpl->fixmu) { x[4] = _pimpl->fixm; }
    if (_pimpl->useflux) 
        f = _pimpl->f_short.subVector(1,_pimpl->f_short.size()) * _pimpl->U;
    else f = _pimpl->f_short*_pimpl->U;

    return ret;
}

void EllipseSolver::calculateNumericH(
    const tmv::Vector<double>& x, const tmv::Vector<double>& f,
    tmv::SymMatrix<double>& h) const
{
    Assert(x.size() == 5);
    Assert(f.size() == 5);
    Assert(h.size() == 5);
    _pimpl->xinit = x;
    if (_pimpl->fixcen) { _pimpl->fixuc = x[0]; _pimpl->fixvc = x[1]; }
    if (_pimpl->fixgam) { _pimpl->fixg1 = x[2]; _pimpl->fixg2 = x[3]; }
    if (_pimpl->fixmu) { _pimpl->fixm = x[4]; }

    _pimpl->x_short = _pimpl->U * x;

    tmv::SymMatrix<double> h_short(_pimpl->x_short.size(),0.);

    _pimpl->flux = 0;
    NLSolver::calculateNumericH(_pimpl->x_short,_pimpl->f_short,h_short);

    h = tmv::SymMatrix<double>(tmv::Matrix<double>(
            _pimpl->U.transpose() * h_short * _pimpl->U));
    if (_pimpl->fixcen) { h(0,0) = h(1,1) = 1.0; }
    if (_pimpl->fixgam) { h(2,2) = h(3,3) = 1.0; }
    if (_pimpl->fixmu) { h(4,4) = 1.0; }
}

bool EllipseSolver::testJ(
    const tmv::Vector<double>& x, tmv::Vector<double>& f,
    std::ostream* os, double relerr) const 
{
    Assert(x.size() == 5);
    Assert(f.size() == 5);
    _pimpl->xinit = x;
    if (_pimpl->fixcen) { _pimpl->fixuc = x[0]; _pimpl->fixvc = x[1]; }
    if (_pimpl->fixgam) { _pimpl->fixg1 = x[2]; _pimpl->fixg2 = x[3]; }
    if (_pimpl->fixmu) { _pimpl->fixm = x[4]; }

    _pimpl->x_short = _pimpl->U * x;

    _pimpl->flux = 0;
    return NLSolver::testJ(_pimpl->x_short,_pimpl->f_short,os,relerr);
}

