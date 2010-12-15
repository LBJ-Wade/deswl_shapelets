
#include "EllipseSolver.h"
#include "dbg.h"
#include "PsiHelper.h"

//#define JTEST
//#define ALWAYS_NUMERIC_J

#ifdef JTEST
#include "TestHelper.h"
bool shouldShowTests;
bool shouldThrow;
std::string lastSuccess;
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

    mutable BVec b0;
    mutable BVec bx;
    mutable DMatrix D;
    mutable DMatrix S;
    mutable DMatrix T;
    mutable DMatrix dTdx;
    mutable DMatrix dTdy;
    mutable DMatrix dSdg1;
    mutable DMatrix dSdg2;
    mutable DMatrix dDdmu;
    mutable DVector dbdE;
    mutable DVector Db0;
    mutable DVector GDb0;
    mutable DVector GthDb0;
    mutable DVector SDb0;
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
    const BVec& _b0, int order,
    bool _fixcen, bool _fixgam, bool _fixmu) :
    b0(_b0), bx(order,b0.getSigma()),
    D(bx.size(),bx.size()), S(bx.size(),bx.size()), T(bx.size(),bx.size()),
    dTdx(bx.size(),bx.size()), dTdy(bx.size(),bx.size()),
    dSdg1(bx.size(),bx.size()), dSdg2(bx.size(),bx.size()),
    dDdmu(bx.size(),bx.size()), dbdE(6),
    Db0(bx.size()), GDb0(bx.size()), GthDb0(bx.size()), SDb0(bx.size()),
    fixcen(_fixcen), fixgam(_fixgam), fixmu(_fixmu),
    numeric_j(false), U((fixcen?0:2)+(fixgam?0:2)+(fixmu?0:1),5), 
    xinit(5), xx(5), ff(5), jj(5,5),
    x_short(U.TMV_colsize()), f_short(U.TMV_colsize()),
    Gx(bx.size(),bx.size()), Gy(bx.size(),bx.size()), 
    Gg1(bx.size(),bx.size()), Gg2(bx.size(),bx.size()), 
    Gth(bx.size(),bx.size()), Gmu(bx.size(),bx.size())
{
    xdbg<<"EllipseSolver3: \n";
    xdbg<<"b0 = "<<b0.vec()<<std::endl;
    xdbg<<"b0.order, sigma = "<<b0.getOrder()<<"  "<<b0.getSigma()<<std::endl;
    xdbg<<"order = "<<order<<std::endl;
    xdbg<<"bx.order, sigma = "<<bx.getOrder()<<"  "<<bx.getSigma()<<std::endl;
    xdbg<<"fixcen,gam,mu = "<<fixcen<<"  "<<fixgam<<"  "<<fixmu<<std::endl;
    U.setZero();
    xinit.setZero();

    int k=0;
    if (!fixcen) { U(k++,0) = 1.; U(k++,1) = 1.; }
    if (!fixgam) { U(k++,2) = 1.; U(k++,3) = 1.; }
    if (!fixmu) { U(k,4) = 1.; }

    setupGx(Gx,order,order);
    setupGy(Gy,order,order);
    setupGg1(Gg1,order,order);
    setupGg2(Gg2,order,order);
    setupGth(Gth,order,order);
    setupGmu(Gmu,order,order);

    if (fixcen) { T.setToIdentity(); }
    if (fixgam) { S.setToIdentity(); }
    if (fixmu) { D.setToIdentity(); }
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

    bx = b0;
    //xdbg<<"Start with b0 = "<<b0.vec()<<std::endl;
    //xdbg<<"bx = "<<bx.vec()<<std::endl;
    if (!fixmu) {
        calculateMuTransform(mu,order,D);
        bx.vec() = D * bx.vec();
        //xdbg<<"mu = "<<mu<<" bx => "<<bx.vec()<<std::endl;
    }
    if (!fixgam) {
        calculateGTransform(g,order,S);
        bx.vec() = S * bx.vec();
        //xdbg<<"g = "<<g<<" bx => "<<bx.vec()<<std::endl;
    }
    if (!fixcen) {
        calculateZTransform(zc,order,T);
        bx.vec() = T * bx.vec();
        //xdbg<<"zc = "<<zc<<" bx => "<<bx.vec()<<std::endl;
    }

    f = bx.vec().TMV_subVector(1,6);
    f /= bx(0);
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
    double fact = 1./(1.-norm(g));

    if (!fixgam || !fixmu) {
        Db0 = D.colRange(0,b0.size()) * b0.vec();
    }

    // Leave off one factor of bx(0) for now.
    // We apply it at the end to the whole matrix.
    if (!fixcen) {
        // dT/dx = Gx * T;
        // dbx/dx = Gx T S D b0 = Gx bx
        dbdE = Gx.rowRange(0,6) * bx.vec();
        //dbg<<"dbdE = "<<dbdE<<std::endl;
        J.col(0) = dbdE.subVector(1,6) - dbdE(0) * f;
        //dbg<<"J(0) = "<<J.col(0)<<std::endl;
        dbdE = Gy.rowRange(0,6) * bx.vec();
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
        // dSdg1 = fact * S * (Gg1 + g2 * Gth);
        // dSdg2 = fact * S * (Gg2 - g1 * Gth);
        // For these, S doesn't commute with Gg1, Gg2, or Gth, so it's
        // important that S comes first.  For x,y,mu, it doesn't matter.
        double g1 = real(g);
        double g2 = imag(g);
        Db0 = D.colRange(0,b0.size()) * b0.vec();
        GDb0 = Gg1 * Db0;
        GthDb0 = Gth * Db0;
        GDb0 += g2*GthDb0;
        SDb0 = fact * S * GDb0;
        dbdE = T.rowRange(0,6) * SDb0;
        //dbg<<"dbdE = "<<dbdE<<std::endl;
        J.col(2) = dbdE.subVector(1,6) - dbdE(0) * f;
        //dbg<<"J(2) = "<<J.col(2)<<std::endl;
        GDb0 = Gg2 * Db0;
        GDb0 -= g1*GthDb0;
        SDb0 = fact * S * GDb0;
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
        // dD/dmu = -GmuT D
        GDb0 = -Gmu.transpose() * Db0;
        SDb0 = S * GDb0; 
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
    double mu = x[4];
    int order = bx.getOrder();
    
    // This section was only used for debugging purposes, but
    // I'm leaving it in, since the derivations of Gx[i] are
    // helpful in understanding some of the code following this section.
    DVector f0(5);
    doF(x,f0);
    DMatrix J_num(5,5);
    DVector f1(5);
    DVector f2(5);
    double dE = 1.e-8;
    for(int j=0;j<5;++j) {
        DVector x1 = x;
        x1(j) -= dE;
        doF(x1,f1);
        DVector x2 = x;
        x2(j) += dE;
        doF(x2,f2);
        J_num.col(j) = (f2-f1)/(2.*dE);
        dbg<<"j = "<<j<<std::endl;
        dbg<<"(f2-f0)/dE = "<<(f2-f0)/dE<<std::endl;
        dbg<<"(f0-f1)/dE = "<<(f0-f1)/dE<<std::endl;
        dbg<<"diff = "<<(f2-f0-f0+f1)/dE<<std::endl;
    }

    dbg<<"numerical df: "<<J_num<<std::endl;
    dbg<<"analytic: "<<J.clip(1.e-5)<<std::endl;
    dbg<<"Norm(diff) = "<<Norm(J-J_num)<<std::endl;
    doF(x,f0); // return everything to original values

    //
    // dD/dmu:
    //
    DMatrix D0(bx.size(),bx.size());
    DMatrix D1(bx.size(),bx.size());
    DMatrix D2(bx.size(),bx.size());
    calculateMuTransform(mu,order,D0);
    calculateMuTransform(mu-dE,order,D1);
    calculateMuTransform(mu+dE,order,D2);

    // technically, this needs to be done at order+2, but we 
    // are already bumping up the order pretty high, so I don't 
    // think it's necessary. 
    // But we'll only test the agreement to order-2
    const int n1 = (order-1)*order/2;

    DMatrix dDdmu_num = (D2-D1)/(2.*dE);
    DMatrix d2D_num = (D2+D1-2.*D0)/(dE*dE);
    // I'm not sure why the - and transpose are necessary here.
    // Maybe there is an error in Gmu?  
    dDdmu = -Gmu.transpose() * D;
    dbg<<"Gmu = "<<Gmu.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"D = "<<D.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"dDdmu = "<<dDdmu.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"dDdmu_num = "<<dDdmu_num.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"Norm(dD/dmu - numeric dD/dmu) = "<< 
        Norm((dDdmu-dDdmu_num).subMatrix(0,n1,0,n1))<<std::endl;
    dbg<<"dE*Norm(d2D_num) = "<<
        dE*Norm(d2D_num.subMatrix(0,n1,0,n1))<<std::endl;
    test(Norm((dDdmu-dDdmu_num).subMatrix(0,n1,0,n1)) <
         10.*dE*Norm(d2D_num.subMatrix(0,n1,0,n1)),"dDdmu");

    //
    // dS/dg1
    //
    DMatrix S0(bx.size(),bx.size());
    DMatrix S1(bx.size(),bx.size());
    DMatrix S2(bx.size(),bx.size());
    calculateGTransform(g,order,S0);
    calculateGTransform(g-std::complex<double>(dE,0),order,S1);
    calculateGTransform(g+std::complex<double>(dE,0),order,S2);

    DMatrix dSdg1_num = (S2-S1)/(2.*dE);
    DMatrix d2S_num = (S2+S1-2.*S0)/(dE*dE);
    dSdg1 = fact * S * (Gg1 + imag(g) * Gth);
    dbg<<"Gg1 = "<<Gg1.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"Gth = "<<Gth.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"S = "<<S.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"dSdg1 = "<<dSdg1.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"dSdg1_num = "<<dSdg1_num.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"Norm(dS/dg1 - numeric dS/dg1) = "<< 
        Norm((dSdg1-dSdg1_num).subMatrix(0,n1,0,n1))<<std::endl;
    dbg<<"dE*Norm(d2S_num) = "<<
        dE*Norm(d2S_num.subMatrix(0,n1,0,n1))<<std::endl;
    test(Norm((dSdg1-dSdg1_num).subMatrix(0,n1,0,n1)) <
         10.*dE*Norm(d2S_num.subMatrix(0,n1,0,n1)),"dSdg1");

    //
    // dS/dg2
    //
    calculateGTransform(g-std::complex<double>(0,dE),order,S1);
    calculateGTransform(g+std::complex<double>(0,dE),order,S2);

    DMatrix dSdg2_num = (S2-S1)/(2.*dE);
    d2S_num = (S2+S1-2.*S0)/(dE*dE);
    dSdg2 = fact * S * (Gg2 - real(g) * Gth);
    dbg<<"Gg2 = "<<Gg2.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"Gth = "<<Gth.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"S = "<<S.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"dSdg2 = "<<dSdg2.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"dSdg2_num = "<<dSdg2_num.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"Norm(dS/dg2 - numeric dS/dg2) = "<< 
        Norm((dSdg2-dSdg2_num).subMatrix(0,n1,0,n1))<<std::endl;
    dbg<<"dE*Norm(d2S_num) = "<<
        dE*Norm(d2S_num.subMatrix(0,n1,0,n1))<<std::endl;
    test(Norm((dSdg2-dSdg2_num).subMatrix(0,n1,0,n1)) < 
         10.*dE*Norm(d2S_num.subMatrix(0,n1,0,n1)),"dSdg2");

    //
    // dT/dx
    //
    DMatrix T0(bx.size(),bx.size());
    DMatrix T1(bx.size(),bx.size());
    DMatrix T2(bx.size(),bx.size());
    calculateZTransform(zc,order,T0);
    calculateZTransform(zc-std::complex<double>(dE,0),order,T1);
    calculateZTransform(zc+std::complex<double>(dE,0),order,T2);

    DMatrix dTdx_num = (T2-T1)/(2.*dE);
    DMatrix d2T_num = (T2+T1-2.*T0)/(dE*dE);
    dTdx = Gx * T;
    dbg<<"Gx = "<<Gx.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"T = "<<T.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"dTdx = "<<dTdx.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"dTdx_num = "<<dTdx_num.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"Norm(dT/dx - numeric dT/dx) = "<< 
        Norm((dTdx-dTdx_num).subMatrix(0,n1,0,n1))<<std::endl;
    dbg<<"dE*Norm(d2T_num) = "<<
        dE*Norm(d2T_num.subMatrix(0,n1,0,n1))<<std::endl;
    test(Norm((dTdx-dTdx_num).subMatrix(0,n1,0,n1)) < 
         10.*dE*Norm(d2T_num.subMatrix(0,n1,0,n1)),"dTdx");

    //
    // dT/dy
    //
    calculateZTransform(zc-std::complex<double>(0,dE),order,T1);
    calculateZTransform(zc+std::complex<double>(0,dE),order,T2);

    DMatrix dTdy_num = (T2-T1)/(2.*dE);
    d2T_num = (T2+T1-2.*T0)/(dE*dE);
    dTdy = Gy * T;
    dbg<<"Gy = "<<Gy.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"T = "<<T.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"dTdy = "<<dTdy.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"dTdy_num = "<<dTdy_num.subMatrix(0,10,0,10)<<std::endl;
    dbg<<"Norm(dT/dy - numeric dT/dy) = "<< 
        Norm((dTdy-dTdy_num).subMatrix(0,n1,0,n1))<<std::endl;
    dbg<<"dE*Norm(d2T_num) = "<<
        dE*Norm(d2T_num.subMatrix(0,n1,0,n1))<<std::endl;
    test(Norm((dTdy-dTdy_num).subMatrix(0,n1,0,n1)) < 
         10.*dE*Norm(d2T_num.subMatrix(0,n1,0,n1)),"dTdy");

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

