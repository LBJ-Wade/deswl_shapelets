
#include "EllipseSolver.h"
#include "dbg.h"
#include "PsiHelper.h"

//#define JTEST
//#define ALWAYS_NUMERIC_J

#ifdef JTEST
#include "TestHelper.h"
bool shouldShowTests = true;
bool shouldThrow = false;
std::string lastSuccess = "";
std::ostream* testout = &std::cout;
#endif

struct EllipseSolver3::ESImpl3
{

    ESImpl3(
        const BVec& b0, int order,
        bool isFixedCen, bool isFixedGamma, bool isFixedMu);

    void calculateF(const DVector& x, DVector& f) const;
    void calculateJ(const DVector& x, const DVector& f, DMatrix& J) const;

    void doF(const DVector& x, DVector& f) const;
    void doJ(const DVector& x, const DVector& f, DMatrix& J) const;

    int b0order, bxorder, bxorderp2;
    mutable BVec b0;
    mutable BVec bx;
    mutable DVector bxsave;
    int b0size, bxsize, bxsizep2;
    mutable DMatrix Daug;
    mutable DMatrix Saug;
    mutable DMatrix Taug;
    mutable DMatrixView D;
    mutable DMatrixView S;
    mutable DMatrixView T;
    mutable DVector dbdE;
    mutable DVector Db0;
    mutable DVector SDb0;
    mutable DVector GDb0;
    mutable DVector GthDb0;
    bool fixcen, fixgam, fixmu, numeric_j;
    DMatrix U;
    mutable DVector xinit;
    mutable DVector xx;
    mutable DVector ff;
    mutable DMatrix jj;
    mutable DVector x_short;
    mutable DVector f_short;
    mutable double fixuc,fixvc,fixg1,fixg2,fixm;
    DMatrix Gx;
    DMatrix Gy;
    DMatrix Gg1;
    DMatrix Gg2;
    DMatrix Gth;
    DMatrix Gmu;
};

EllipseSolver3::EllipseSolver3(
    const BVec& b0, int order,
    bool fixcen, bool fixgam, bool fixmu) :
    _pimpl(new ESImpl3(b0,order,fixcen,fixgam,fixmu))
{}

EllipseSolver3::~EllipseSolver3() 
{ delete _pimpl; }

void EllipseSolver3::calculateF(const DVector& x, DVector& f) const
{ _pimpl->calculateF(x,f); }

void EllipseSolver3::calculateJ(
    const DVector& x, const DVector& f, DMatrix& J) const
{ 
#ifdef ALWAYS_NUMERIC_J
    NLSolver::calculateJ(x,f,J);
#else
    if (_pimpl->numeric_j) NLSolver::calculateJ(x,f,J);
    else _pimpl->calculateJ(x,f,J); 
#endif
}

EllipseSolver3::ESImpl3::ESImpl3(
    const BVec& _b0, int _order,
    bool _fixcen, bool _fixgam, bool _fixmu) :
    b0order(_b0.getOrder()), bxorder(_order), bxorderp2(bxorder+2),
    b0(_b0), bx(bxorder,b0.getSigma()), bxsave(bx.size()),
    b0size(b0.size()), bxsize(bx.size()), 
    bxsizep2((bxorderp2+1)*(bxorderp2+2)/2),
    Daug(bxsize,bxsizep2), Saug(bxsize,bxsizep2), Taug(bxsize,bxsizep2),
    D(Daug.colRange(0,bxsize)),
    S(Saug.colRange(0,bxsize)),
    T(Taug.colRange(0,bxsize)),
    dbdE(6), Db0(bxsize), SDb0(bxsize), GDb0(bxsizep2), GthDb0(bxsizep2),
    fixcen(_fixcen), fixgam(_fixgam), fixmu(_fixmu),
    numeric_j(false), U((fixcen?0:2)+(fixgam?0:2)+(fixmu?0:1),5), 
    xinit(5), xx(5), ff(5), jj(5,5),
    x_short(U.TMV_colsize()), f_short(U.TMV_colsize()),
    Gx(bxsizep2,bxsize), Gy(bxsizep2,bxsize), 
    Gg1(bxsizep2,bxsize), Gg2(bxsizep2,bxsize), 
    Gth(bxsizep2,bxsize), Gmu(bxsize,bxsizep2)
{
    //xdbg<<"EllipseSolver3: \n";
    //xdbg<<"b0 = "<<b0.vec()<<std::endl;
    //xdbg<<"b0.order, sigma = "<<b0.getOrder()<<"  "<<b0.getSigma()<<std::endl;
    //xdbg<<"order = "<<order<<std::endl;
    //xdbg<<"bx.order, sigma = "<<bx.getOrder()<<"  "<<bx.getSigma()<<std::endl;
    //xdbg<<"fixcen,gam,mu = "<<fixcen<<"  "<<fixgam<<"  "<<fixmu<<std::endl;
    U.setZero();
    xinit.setZero();

    int k=0;
    if (!fixcen) { U(k++,0) = 1.; U(k++,1) = 1.; }
    if (!fixgam) { U(k++,2) = 1.; U(k++,3) = 1.; }
    if (!fixmu) { U(k,4) = 1.; }

    setupGx(Gx,bxorderp2,bxorder);
    setupGy(Gy,bxorderp2,bxorder);
    setupGg1(Gg1,bxorderp2,bxorder);
    setupGg2(Gg2,bxorderp2,bxorder);
    setupGth(Gth,bxorderp2,bxorder);
    setupGmu(Gmu,bxorder,bxorderp2);

    if (fixcen) { 
        T.setToIdentity();
        TMV_colRange(Taug,bxsize,bxsizep2).setZero(); 
    }
    if (fixgam) { 
        S.setToIdentity();
        TMV_colRange(Saug,bxsize,bxsizep2).setZero(); 
    }
    if (fixmu) { 
        D.setToIdentity();
        TMV_colRange(Daug,bxsize,bxsizep2).setZero(); 
    }
}

void EllipseSolver3::ESImpl3::doF(const DVector& x, DVector& f) const
{
    Assert(x.size() == 5);
    Assert(f.size() == 5);

    //xdbg<<"Start doF\n";
    //xdbg<<"x = "<<x<<std::endl;

    if ((x-xinit).TMV_normInf() > 4.) {
        dbg<<"bad x-xinit: "<<(x-xinit)<<std::endl;
        xdbg<<"x = "<<x<<std::endl;
        f = 2.* bx.vec().TMV_subVector(1,6) / bx(0);
        return;
    }

    std::complex<double> zc(x[0],x[1]);
    std::complex<double> g(x[2],x[3]);
    double mu = x[4];

    if (norm(g) > 0.99) {
        dbg<<"bad gsq = "<<norm(g)<<std::endl;
        xdbg<<"x = "<<x<<std::endl;
        f = 2.* bx.vec().TMV_subVector(1,6) / bx(0);
        return;
    }
    int order = bx.getOrder();

    bxsave = bx.vec();
    bx = b0;
    //xdbg<<"Start with b0 = "<<b0.vec()<<std::endl;
    //xdbg<<"bx = "<<bx.vec()<<std::endl;
    if (!fixmu) {
        calculateMuTransform(mu,order,Daug);
        bx.vec() = D * bx.vec();
        //xdbg<<"mu = "<<mu<<" bx => "<<bx.vec()<<std::endl;
    } else if (fixm != 0.) {
        bx.vec() = D * bx.vec();
    }
    if (!fixgam) {
        calculateGTransform(g,order,Saug);
        bx.vec() = S * bx.vec();
        //xdbg<<"g = "<<g<<" bx => "<<bx.vec()<<std::endl;
    } else if (fixg1 != 0. || fixg2 != 0.) {
        bx.vec() = S * bx.vec();
    }
    if (!fixcen) {
        calculateZTransform(zc,order,Taug);
        bx.vec() = T * bx.vec();
        //xdbg<<"zc = "<<zc<<" bx => "<<bx.vec()<<std::endl;
    } else if (fixuc != 0. || fixvc != 0.) {
        bx.vec() = T * bx.vec();
    }

    if (bx(0) <= 0.) {
        dbg<<"bad bx(0): "<<bx(0)<<std::endl;
        xdbg<<"x = "<<x<<std::endl;
        xdbg<<"bx = "<<x<<std::endl;
        f = 2.* bxsave.TMV_subVector(1,6) / bxsave(0);
        bx.vec() = bxsave;
        return;
    }

    f = bx.vec().TMV_subVector(1,6) / bx(0);
}

void EllipseSolver3::ESImpl3::doJ(
    const DVector& x, const DVector& f, DMatrix& J) const
{
    // bx = T S D b0
    //
    // dbx/dzc = (dT/dzc) S D b0
    // dbx/dg  = T (dS/dg) D b0
    // dbx/dmu = T S (dD/dmu) b0
    //
    // f = bx(1:6) / bx(0)
    // df/dE = dbx/dE(1:6) / bx(0) - (1/bx(0)^2) (dbx(0)/dE) bx(1:6)
    //       = ( dbx/dE(1:6) - dbx/dE(0) f ) / bx(0)

    //xdbg<<"Start doJ\n";
    //xdbg<<"x = "<<x<<std::endl;
    //xdbg<<"f = "<<f<<std::endl;
    //xdbg<<"J = "<<J<<std::endl;

    Assert(x.size() == 5);
    Assert(J.TMV_rowsize() == 5);
    Assert(J.TMV_colsize() == 5);

    std::complex<double> zc(x[0],x[1]);
    std::complex<double> g(x[2],x[3]);
    double mu = x[4];
    double fact = 1./(1.-norm(g));

    if (!fixcen || !fixgam) {
        // Both cen and gam use this:
        Db0 = D.colRange(0,b0size) * b0.vec();
    }

    // Leave off one factor of bx(0) for now.
    // We apply it at the end to the whole matrix.
    if (!fixcen) {
        // dT/dx = T Gx;
        // dT/dy = T Gy;
        // db/dx = T Gx S D b0
        // db/dy = T Gy S D b0
        augmentZTransformCols(zc,bxorder,Taug);
        SDb0 = S * Db0;
        GDb0 = Gx * SDb0;
        dbdE = Taug.rowRange(0,6) * GDb0;
        //dbg<<"dbdE = "<<dbdE<<std::endl;
        J.col(0) = dbdE.subVector(1,6) - dbdE(0) * f;
        //dbg<<"J(0) = "<<J.col(0)<<std::endl;
        GDb0 = Gy * SDb0;
        dbdE = Taug.rowRange(0,6) * GDb0;
        //dbg<<"dbdE = "<<dbdE<<std::endl;
        J.col(1) = dbdE.subVector(1,6) - dbdE(0) * f;
        //dbg<<"J(1) = "<<J.col(1)<<std::endl;
    } else {
        J.col(0).setZero();
        //dbg<<"J(0) = "<<J.col(0)<<std::endl;
        J.col(1).setZero();
        //dbg<<"J(1) = "<<J.col(1)<<std::endl;
    }

    if (!fixgam) {
        // dS/dg1 = fact * S * (Gg1 + g2 * Gth);
        // dS/dg2 = fact * S * (Gg2 - g1 * Gth);
        augmentGTransformCols(g,bxorder,Saug);
        double g1 = real(g);
        double g2 = imag(g);
        GDb0 = Gg1 * Db0;
        GthDb0 = Gth * Db0;
        GDb0 += g2*GthDb0;
        SDb0 = fact * Saug * GDb0;
        dbdE = T.rowRange(0,6) * SDb0;
        //dbg<<"dbdE = "<<dbdE<<std::endl;
        J.col(2) = dbdE.subVector(1,6) - dbdE(0) * f;
        //dbg<<"J(2) = "<<J.col(2)<<std::endl;
        GDb0 = Gg2 * Db0;
        GDb0 -= g1*GthDb0;
        SDb0 = fact * Saug * GDb0;
        dbdE = T.rowRange(0,6) * SDb0;
        //dbg<<"dbdE = "<<dbdE<<std::endl;
        J.col(3) = dbdE.subVector(1,6) - dbdE(0) * f;
        //dbg<<"J(3) = "<<J.col(3)<<std::endl;
    } else {
        J.col(2).setZero();
        //dbg<<"J(2) = "<<J.col(2)<<std::endl;
        J.col(3).setZero();
        //dbg<<"J(3) = "<<J.col(3)<<std::endl;
    }

    if (!fixmu) {
        // dD/dmu = -D GmuT
        augmentMuTransformCols(mu,bxorder,Daug);
        GDb0 = -Gmu.transpose().colRange(0,b0size) * b0.vec();
        Db0 = Daug * GDb0;
        SDb0 = S * Db0; 
        dbdE = T.rowRange(0,6) * SDb0;
        //dbg<<"dbdE = "<<dbdE<<std::endl;
        J.col(4) = dbdE.subVector(1,6) - dbdE(0) * f;
        //dbg<<"J(4) = "<<J.col(4)<<std::endl;
    } else {
        J.col(4).setZero();
        //dbg<<"J(4) = "<<J.col(4)<<std::endl;
    }
    J /= bx(0);
    //dbg<<"J = "<<J<<std::endl;

#ifdef JTEST
    int order = bx.getOrder();
    double dE = 1.e-6;
    dbg<<"dE = "<<dE<<std::endl;
    int orderp2 = order+2;
    testout = dbgout;
    
    // This section was only used for debugging purposes, but
    // I'm leaving it in, since the derivations of Gx[i] are
    // helpful in understanding some of the code following this section.
    
    //
    // dD/dmu:
    //
    DMatrix Dbig(bxsizep2,bxsizep2);
    DMatrix D0(bxsize,bxsize);
    DMatrix D1(bxsize,bxsize);
    DMatrix D2(bxsize,bxsize);
    calculateMuTransform(mu,orderp2,Dbig);
    calculateMuTransform(mu,order,D0);
    calculateMuTransform(mu-dE,order,D1);
    calculateMuTransform(mu+dE,order,D2);
    DMatrix Gmubig(bxsize,bxsizep2);
    setupGmu(Gmubig,order,orderp2);

    DMatrix dDdmu_num = (D2-D1)/(2.*dE);
    DMatrix d2D_num = (D2+D1-2.*D0)/(dE*dE);
    // The - and transpose are really to change the diagonal of Gmu
    // from -1 to +1.  dD/dmu is really D Gmu + 2D.
    // But because of the structure of Gmu (being anti-symmetric with -1's
    // on the diagonal), this is equivalent to -D Gmu.transpose().
    DMatrix dDdmu = -Dbig.rowRange(0,bxsize) * Gmubig.transpose();
    dbg<<"Gmu = "<<Gmubig.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"D = "<<Dbig.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"dDdmu = "<<dDdmu.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"dDdmu_num = "<<dDdmu_num.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"Norm(dD/dmu - numeric dD/dmu) = "<<Norm(dDdmu-dDdmu_num)<<std::endl;
    dbg<<"dE*Norm(d2D_num) = "<<dE*Norm(d2D_num)<<std::endl;
    test(Norm(dDdmu-dDdmu_num) < 30.*dE*Norm(d2D_num),"dDdmu");

    //
    // dS/dg1
    //
    DMatrix Sbig(bxsizep2,bxsizep2);
    DMatrix S0(bxsize,bxsize);
    DMatrix S1(bxsize,bxsize);
    DMatrix S2(bxsize,bxsize);
    calculateGTransform(g,orderp2,Sbig);
    calculateGTransform(g,order,S0);
    calculateGTransform(g-std::complex<double>(dE,0),order,S1);
    calculateGTransform(g+std::complex<double>(dE,0),order,S2);
    DMatrix Gg1big(bxsizep2,bxsize);
    DMatrix Gthbig(bxsizep2,bxsize);
    setupGg1(Gg1big,orderp2,order);
    setupGth(Gthbig,orderp2,order);

    DMatrix dSdg1_num = (S2-S1)/(2.*dE);
    DMatrix d2S_num = (S2+S1-2.*S0)/(dE*dE);
    DMatrix dSdg1 = fact * Sbig.rowRange(0,bxsize) * (Gg1big + imag(g) * Gthbig);
    dbg<<"Gg1 = "<<Gg1big.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"Gth = "<<Gthbig.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"S = "<<Sbig.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"dSdg1 = "<<dSdg1.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"dSdg1_num = "<<dSdg1_num.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"Norm(dS/dg1 - numeric dS/dg1) = "<<Norm(dSdg1-dSdg1_num)<<std::endl;
    dbg<<"dE*Norm(d2S_num) = "<<dE*Norm(d2S_num)<<std::endl;
    test(Norm(dSdg1-dSdg1_num) < 30.*dE*Norm(d2S_num),"dSdg1");

    //
    // dS/dg2
    //
    calculateGTransform(g-std::complex<double>(0,dE),order,S1);
    calculateGTransform(g+std::complex<double>(0,dE),order,S2);
    DMatrix Gg2big(bxsizep2,bxsize);
    setupGg2(Gg2big,orderp2,order);

    DMatrix dSdg2_num = (S2-S1)/(2.*dE);
    d2S_num = (S2+S1-2.*S0)/(dE*dE);
    DMatrix dSdg2 = fact * Sbig.rowRange(0,bxsize) * (Gg2big - real(g) * Gthbig);
    dbg<<"Gg2 = "<<Gg2big.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"Gth = "<<Gthbig.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"S = "<<Sbig.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"dSdg2 = "<<dSdg2.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"dSdg2_num = "<<dSdg2_num.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"Norm(dS/dg2 - numeric dS/dg2) = "<<Norm(dSdg2-dSdg2_num)<<std::endl;
    dbg<<"dE*Norm(d2S_num) = "<<dE*Norm(d2S_num)<<std::endl;
    test(Norm(dSdg2-dSdg2_num) < 30.*dE*Norm(d2S_num),"dSdg2");

    //
    // dT/dx
    //
    DMatrix Tbig(bxsizep2,bxsizep2);
    DMatrix T0(bxsize,bxsize);
    DMatrix T1(bxsize,bxsize);
    DMatrix T2(bxsize,bxsize);
    calculateZTransform(zc,orderp2,Tbig);
    calculateZTransform(zc,order,T0);
    calculateZTransform(zc-std::complex<double>(dE,0),order,T1);
    calculateZTransform(zc+std::complex<double>(dE,0),order,T2);
    DMatrix Gxbig(bxsizep2,bxsize);
    setupGx(Gxbig,orderp2,order);

    DMatrix dTdx_num = (T2-T1)/(2.*dE);
    DMatrix d2T_num = (T2+T1-2.*T0)/(dE*dE);
    DMatrix dTdx = Tbig.rowRange(0,bxsize) * Gxbig;
    dbg<<"Gx = "<<Gxbig.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"T = "<<Tbig.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"dTdx = "<<dTdx.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"dTdx_num = "<<dTdx_num.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"Norm(dT/dx - numeric dT/dx) = "<<Norm(dTdx-dTdx_num)<<std::endl;
    dbg<<"dE*Norm(d2T_num) = "<<dE*Norm(d2T_num)<<std::endl;
    test(Norm(dTdx-dTdx_num) < 30.*dE*Norm(d2T_num),"dTdx");

    //
    // dT/dy
    //
    calculateZTransform(zc-std::complex<double>(0,dE),order,T1);
    calculateZTransform(zc+std::complex<double>(0,dE),order,T2);
    DMatrix Gybig(bxsizep2,bxsize);
    setupGy(Gybig,orderp2,order);

    DMatrix dTdy_num = (T2-T1)/(2.*dE);
    d2T_num = (T2+T1-2.*T0)/(dE*dE);
    DMatrix dTdy = Tbig.rowRange(0,bxsize) * Gybig;
    dbg<<"Gy = "<<Gybig.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"T = "<<Tbig.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"dTdy = "<<dTdy.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"dTdy_num = "<<dTdy_num.subMatrix(0,6,0,6)<<std::endl;
    dbg<<"Norm(dT/dy - numeric dT/dy) = "<<Norm(dTdy-dTdy_num)<<std::endl;
    dbg<<"dE*Norm(d2T_num) = "<<dE*Norm(d2T_num)<<std::endl;
    test(Norm(dTdy-dTdy_num) < 30.*dE*Norm(d2T_num),"dTdy");

    //
    // J 
    //
    
    DMatrix J_num(5,5);
    DVector x0 = x;
    DVector f0(5);
    doF(x0,f0);
    DVector bx0 = bx.vec().subVector(0,6);
    dbg<<"x0 = "<<x0<<std::endl;
    dbg<<"b0 = "<<bx0<<std::endl;
    dbg<<"f0 = "<<f0<<std::endl;

    for(int k=0;k<5;k++) {
        dbg<<"k = "<<k<<std::endl;

        DVector f1(5);
        DVector x1 = x; x1(k) -= dE;
        doF(x1,f1);
        DVector b1 = bx.vec().subVector(0,6);
        dbg<<"x1 = "<<x1<<std::endl;
        dbg<<"b1 = "<<b1<<std::endl;
        dbg<<"f1 = "<<f1<<std::endl;

        DVector f2(5);
        DVector x2 = x; x2(k) += dE;
        doF(x2,f2);
        DVector b2 = bx.vec().subVector(0,6);
        dbg<<"x2 = "<<x2<<std::endl;
        dbg<<"b2 = "<<b2<<std::endl;
        dbg<<"f2 = "<<f2<<std::endl;

        DVector dbdE_num = (b2-b1)/(2.*dE);
        dbg<<"dbdE_num = "<<dbdE_num<<std::endl;
        DVector d2bdE = (b2+b1-2.*bx0)/(dE*dE);
        dbg<<"d2bdE = "<<d2bdE<<std::endl;

        DVector dbdE(6);
        double g1 = real(g);
        double g2 = imag(g);
        Db0 = D.colRange(0,b0size) * b0.vec();
        dbg<<"Db0 = "<<Db0<<std::endl;
        switch (k) {
          case 0: SDb0 = S * Db0;
                  dbg<<"SDb0 = "<<SDb0<<std::endl;
                  GDb0 = Gx * SDb0;
                  dbg<<"GDb0 = "<<GDb0<<std::endl;
                  dbdE = Tbig.rowRange(0,6) * GDb0;
                  dbg<<"dbdx = "<<dbdE<<std::endl;
                  dbg<<"dbdx = dTdx S D b0 = "<<
                      dTdx.rowRange(0,6) * S * Db0<<std::endl;
                  break;
          case 1: SDb0 = S * Db0;
                  dbg<<"SDb0 = "<<SDb0<<std::endl;
                  GDb0 = Gy * SDb0;
                  dbg<<"GDb0 = "<<GDb0<<std::endl;
                  dbdE = Tbig.rowRange(0,6) * GDb0;
                  dbg<<"dbdy = "<<dbdE<<std::endl;
                  dbg<<"dbdy = dTdy S D b0 = "<<
                      dTdy.rowRange(0,6) * S * Db0<<std::endl;
                  break;
          case 2: GDb0 = Gg1 * Db0;
                  dbg<<"GDb0 = "<<GDb0<<std::endl;
                  GthDb0 = Gth * Db0;
                  dbg<<"GthDb0 = "<<GthDb0<<std::endl;
                  GDb0 += g2*GthDb0;
                  dbg<<"GDb0 = "<<GDb0<<std::endl;
                  SDb0 = fact * Sbig.rowRange(0,bxsize) * GDb0;
                  dbg<<"SDb0 = "<<SDb0<<std::endl;
                  dbdE = T.rowRange(0,6) * SDb0;
                  dbg<<"dbdg1 = "<<dbdE<<std::endl;
                  dbg<<"dbdg1 = T dSdg1 D b0 = "<<
                      T.rowRange(0,6) * dSdg1 * Db0<<std::endl;
                  break;
          case 3: GDb0 = Gg2 * Db0;
                  dbg<<"GDb0 = "<<GDb0<<std::endl;
                  GthDb0 = Gth * Db0;
                  dbg<<"GthDb0 = "<<GthDb0<<std::endl;
                  GDb0 -= g1*GthDb0;
                  dbg<<"GDb0 = "<<GDb0<<std::endl;
                  SDb0 = fact * Sbig.rowRange(0,bxsize) * GDb0;
                  dbg<<"SDb0 = "<<SDb0<<std::endl;
                  dbdE = T.rowRange(0,6) * SDb0;
                  dbg<<"dbdg2 = "<<dbdE<<std::endl;
                  dbg<<"dbdg2 = T dSdg2 D b0 = "<<
                      T.rowRange(0,6) * dSdg2 * Db0<<std::endl;
                  break;
          case 4: GDb0 = -Gmu.transpose().colRange(0,b0size) * b0.vec();
                  dbg<<"GDb0 = "<<GDb0<<std::endl;
                  Db0 = Dbig.rowRange(0,bxsize) * GDb0;
                  dbg<<"Db0 = "<<Db0<<std::endl;
                  SDb0 = S * Db0; 
                  dbg<<"SDb0 = "<<SDb0<<std::endl;
                  dbdE = T.rowRange(0,6) * SDb0;
                  dbg<<"dbdmu = "<<dbdE<<std::endl;
                  dbg<<"dbdmu = T S dDdmu b0 = "<<
                      T.rowRange(0,6) * S * dDdmu.colRange(0,b0size) * b0.vec()
                      <<std::endl;
                  break;
        }
        dbg<<"dbdE = "<<dbdE<<std::endl;
        dbg<<"Norm(dbdE-dbdE_num) = "<<Norm(dbdE-dbdE_num)<<std::endl;
        dbg<<"dE*Norm(d2bdE) = "<<dE*Norm(d2bdE)<<std::endl;
        test(Norm(dbdE-dbdE_num) < 30.*dE*Norm(d2bdE),"dbdE");

        DVector dfdE = dbdE.subVector(1,6) - dbdE(0) * f0;
        dfdE /= bx(0);
        dbg<<"dfdE from dbdE = "<<dfdE<<std::endl;
        dbg<<"dfdE from dbdE_num = "<<
            (dbdE_num.subVector(1,6) - dbdE_num(0) * f0)/bx(0)<<std::endl;;
        dbg<<"dfdE in J = "<<J.col(k)<<std::endl;
        J_num.col(k) = (f2-f1)/(2.*dE);
        dbg<<"dfdE_num = "<<J_num.col(k)<<std::endl;
        DVector d2f_num = (f2+f1-2.*f0)/(dE*dE);
        dbg<<"d2f_num = "<<d2f_num<<std::endl;
        dbg<<"Norm(dfdE_num-dfdE) ="<<Norm(J_num.col(k)-dfdE)<<std::endl;
        dbg<<"dE*Norm(d2f_num) ="<<dE*Norm(d2f_num)<<std::endl;
        test(Norm(J_num.col(k)-dfdE) < 30.*dE*Norm(d2f_num),"dfdE");
    }

    dbg<<"J_num = "<<J_num<<std::endl;
    dbg<<"analytic: "<<J<<std::endl;
    if (fixcen) J_num.colRange(0,2).setZero(); 
    if (fixgam) J_num.colRange(2,4).setZero(); 
    if (fixmu) J_num.col(4).setZero(); 
    dbg<<"J_num => "<<J_num<<std::endl;
    dbg<<"Norm(diff) = "<<Norm(J-J_num)<<std::endl;
    dbg<<"dE*Norm(J) = "<<dE*Norm(J_num)<<std::endl;
    test(Norm(J_num-J) < dE*Norm(J_num),"dfdE");
    doF(x,f0); // return everything to original values
#endif
}

void EllipseSolver3::ESImpl3::calculateF(const DVector& x, DVector& f) const
{
    Assert(x.size() == U.TMV_colsize());
    Assert(f.size() == U.TMV_colsize());

    Assert(xx.size() == U.TMV_rowsize());
    EIGEN_Transpose(xx) = EIGEN_Transpose(x)*U;
    if (fixcen) { xx[0] = fixuc;  xx[1] = fixvc; }
    if (fixgam) { xx[2] = fixg1;  xx[3] = fixg2; }
    if (fixmu) { xx[4] = fixm; }

    doF(xx,ff);

    f = U*ff;
}

void EllipseSolver3::ESImpl3::calculateJ(
    const DVector& x, const DVector& f, DMatrix& j) const
{
    Assert(x.size() == U.TMV_colsize());
    Assert(f.size() == U.TMV_colsize());
    Assert(j.TMV_rowsize() == U.TMV_colsize());
    Assert(j.TMV_colsize() == U.TMV_colsize());

    EIGEN_Transpose(xx) = EIGEN_Transpose(x)*U;
    if (fixcen) { xx[0] = fixuc;  xx[1] = fixvc; }
    if (fixgam) { xx[2] = fixg1;  xx[3] = fixg2; }
    if (fixmu) { xx[4] = fixm; }

    doJ(xx,ff,jj);
    j = U*jj*U.transpose();
}

void EllipseSolver3::useNumericJ() { _pimpl->numeric_j = true; }

const BVec& EllipseSolver3::getB() const { return _pimpl->bx; }

void EllipseSolver3::getCovariance(DMatrix& cov) const 
{
    dbg<<"getCovariance:\n";
    DMatrix cov1(_pimpl->x_short.size(),_pimpl->x_short.size());
    dbg<<"cov1 = "<<cov1<<std::endl;
    NLSolver::getCovariance(cov1);
    dbg<<"cov1 => "<<cov1<<std::endl;
    cov = _pimpl->U.transpose()*cov1*_pimpl->U;
    dbg<<"full cov = "<<cov<<std::endl;
}

void EllipseSolver3::getInverseCovariance(DMatrix& invcov) const 
{
    DMatrix invcov1(_pimpl->x_short.size(),_pimpl->x_short.size());
    NLSolver::getInverseCovariance(invcov1);
    invcov = _pimpl->U.transpose()*invcov1*_pimpl->U;
    dbg<<"getInverseCovariance:\n";
    dbg<<"invcov1 = "<<invcov1<<std::endl;
    dbg<<"full invcov = "<<invcov<<std::endl;
}

void EllipseSolver3::callF(const DVector& x, DVector& f) const
{
    Assert(x.size() == 5);
    Assert(f.size() == 5);
    _pimpl->xinit = x;
    if (_pimpl->fixcen) { _pimpl->fixuc = x[0]; _pimpl->fixvc = x[1]; }
    if (_pimpl->fixgam) { _pimpl->fixg1 = x[2]; _pimpl->fixg2 = x[3]; }
    if (_pimpl->fixmu) { _pimpl->fixm = x[4]; }

    _pimpl->x_short = _pimpl->U * x;

    calculateF(_pimpl->x_short,_pimpl->f_short);

    EIGEN_Transpose(f) = EIGEN_Transpose(_pimpl->f_short)*_pimpl->U;
}

bool EllipseSolver3::solve(DVector& x, DVector& f) const
{
    Assert(x.size() == 5);
    Assert(f.size() == 5);
    _pimpl->xinit = x;
    if (_pimpl->fixcen) { _pimpl->fixuc = x[0]; _pimpl->fixvc = x[1]; }
    if (_pimpl->fixgam) { _pimpl->fixg1 = x[2]; _pimpl->fixg2 = x[3]; }
    if (_pimpl->fixmu) { _pimpl->fixm = x[4]; }
    if (_pimpl->fixcen && (x[0] != 0. || x[1] != 0.)) { 
        calculateZTransform(
            std::complex<double>(x[0],x[1]),_pimpl->bxorder,_pimpl->Taug);
    }
    if (_pimpl->fixgam && (x[2] != 0. || x[3] != 0.)) { 
        calculateGTransform(
            std::complex<double>(x[2],x[3]),_pimpl->bxorder,_pimpl->Saug);
    }
    if (_pimpl->fixmu && x[4] != 0.) { 
        calculateMuTransform(x[4],_pimpl->bxorder,_pimpl->Daug);
    }

    _pimpl->x_short = _pimpl->U * x;

    bool ret;
    try {
        ret = NLSolver::solve(_pimpl->x_short,_pimpl->f_short);
    } catch (...) {
        xdbg<<"Caught exception during NLSolver::solve"<<std::endl;
        ret = false;
    }

    EIGEN_Transpose(x) = EIGEN_Transpose(_pimpl->x_short)*_pimpl->U;
    if (_pimpl->fixcen) { x[0] = _pimpl->fixuc; x[1] = _pimpl->fixvc; }
    if (_pimpl->fixgam) { x[2] = _pimpl->fixg1; x[3] = _pimpl->fixg2; }
    if (_pimpl->fixmu) { x[4] = _pimpl->fixm; }
    EIGEN_Transpose(f) = EIGEN_Transpose(_pimpl->f_short)*_pimpl->U;

    return ret;
}

bool EllipseSolver3::testJ(
    const DVector& x, DVector& f,
    std::ostream* os, double relerr) const 
{
    Assert(x.size() == 5);
    Assert(f.size() == 5);
    _pimpl->xinit = x;
    if (_pimpl->fixcen) { _pimpl->fixuc = x[0]; _pimpl->fixvc = x[1]; }
    if (_pimpl->fixgam) { _pimpl->fixg1 = x[2]; _pimpl->fixg2 = x[3]; }
    if (_pimpl->fixmu) { _pimpl->fixm = x[4]; }

    _pimpl->x_short = _pimpl->U * x;

    return NLSolver::testJ(_pimpl->x_short,_pimpl->f_short,os,relerr);
}

