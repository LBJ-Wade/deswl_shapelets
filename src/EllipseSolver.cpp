
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

class ESImpl 
{

  public :

    friend class EllipseSolver;

    ESImpl(const std::vector<PixelList>& _pix,
	int _order, double _sigma, 
	bool _fixcen, bool _fixgam, bool _fixmu, bool _useflux);

    ESImpl(const std::vector<PixelList>& _pix,
	const std::vector<BVec>& _psf, 
	double _fp, int _order, double _sigma, 
	bool _fixcen, bool _fixgam, bool _fixmu, bool _useflux);

    void F(const tmv::Vector<double>& x, tmv::Vector<double>& f) const;
    void J(const tmv::Vector<double>& x, const tmv::Vector<double>& f,
	tmv::Matrix<double>& df) const;

    void DoF(const tmv::Vector<double>& x, tmv::Vector<double>& f) const;
    void DoJ(const tmv::Vector<double>& x, const tmv::Vector<double>& f,
	tmv::Matrix<double>& df) const;

  private :

    size_t nsize, np2size;
    mutable BVec b;
    const std::vector<PixelList>& pix;
    const std::vector<BVec>* psf;
    double f_psf;
    std::vector<double> sigma_obs;
    tmv::Vector<double> I;
    tmv::Vector<std::complex<double> > Z;
    mutable tmv::Vector<std::complex<double> > Z1;
    tmv::Vector<double> W;
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
    mutable tmv::Matrix<double> dff;
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

EllipseSolver::EllipseSolver(const std::vector<PixelList>& pix,
    int order, double sigma, 
    bool fixcen, bool fixgam, bool fixmu, bool useflux) :
  pimpl(new ESImpl(pix,order,sigma,fixcen,fixgam,fixmu,useflux))
{}

EllipseSolver::EllipseSolver(const std::vector<PixelList>& pix,
    const std::vector<BVec>& psf, double fp,
    int order, double sigma, 
    bool fixcen, bool fixgam, bool fixmu, bool useflux) :
  pimpl(new ESImpl(pix,psf,fp,order,sigma,fixcen,fixgam,fixmu,useflux))
{}

EllipseSolver::~EllipseSolver()
{ delete pimpl; }

void EllipseSolver::F(
    const tmv::Vector<double>& x, tmv::Vector<double>& f) const
{ pimpl->F(x,f); }

void EllipseSolver::J(
    const tmv::Vector<double>& x, const tmv::Vector<double>& f,
    tmv::Matrix<double>& df) const
{ 
  if (pimpl->numeric_j) NLSolver::J(x,f,df);
  else pimpl->J(x,f,df); 
}

static size_t SumSize(const std::vector<PixelList>& v)
{
  size_t sum = 0;
  for(size_t i=0;i<v.size();i++) sum += v[i].size();
  return sum;
}

ESImpl::ESImpl(const std::vector<PixelList>& _pix, 
    int _order, double _sigma, 
    bool _fixcen, bool _fixgam, bool _fixmu, bool _useflux) :
  nsize((_order+1)*(_order+2)/2),
  np2size((_order+3)*(_order+4)/2), b(_order,_sigma),
  pix(_pix), psf(0), f_psf(1.), sigma_obs(pix.size(),_sigma),
  I(SumSize(pix)), Z(I.size()), Z1(I.size()), W(I.size()), 
  psi(I.size(),nsize,0.), A_aux(I.size(),np2size,0.), AC(0,0),
  A(A_aux.Cols(0,nsize)),
  fixcen(_fixcen), fixgam(_fixgam), fixmu(_fixmu), useflux(_useflux),
  numeric_j(false), U((fixcen?0:2)+(fixgam?0:2)+(fixmu?0:1),5,0.), 
  xinit(5,0.), xx(5), ff(6), dff(6,5),
  x_short(U.colsize()), f_short(U.colsize()+(useflux?1:0)),
  Iresid(I.size()), Cb(0), Gb(np2size,5),
  dbdE(nsize,5), dAdEb(I.size(),5), temp(nsize,5), AtI(np2size),
  dCdE(0,0), flux(0), maxpsforder(0),
  _Gx(np2size,nsize), _Gy(np2size,nsize), 
  _Gg1(np2size,nsize), _Gg2(np2size,nsize), 
  _Gth(np2size,nsize), _Gmu(np2size,nsize),
  _dTdE(0,0)
{
  for(size_t k=0,n=0;k<pix.size();k++) {
    for(size_t i=0;i<pix[k].size();i++,n++) {
      I(n) = pix[k][i].I*pix[k][i].wt;
      W(n) = pix[k][i].wt;
      Z(n) = pix[k][i].z;
    }
  }

  A.DivideUsing(tmv::QR);
  A.SaveDiv();
  Assert(A.colsize() >= A.rowsize());
  // Otherwise there are too few pixels to constrain model.

  int k=0;
  if (!fixcen) { U(k++,0) = 1.; U(k++,1) = 1.; }
  if (!fixgam) { U(k++,2) = 1.; U(k++,3) = 1.; }
  if (!fixmu) { U(k,4) = 1.; }

  SetupGx(_Gx,_order+2,_order);
  SetupGy(_Gy,_order+2,_order);
  SetupGg1(_Gg1,_order+2,_order);
  SetupGg2(_Gg2,_order+2,_order);
  SetupGth(_Gth,_order+2,_order);
  SetupGmu(_Gmu,_order+2,_order);
}

static int CalcMaxPSFOrder(const std::vector<BVec>& psf)
{
  int maxpsforder = 0;
  for(size_t k=0;k<psf.size();k++) 
    if (psf[k].GetOrder() > maxpsforder)
      maxpsforder = psf[k].GetOrder();
  return maxpsforder;
}

#define ORDER_1 (std::max(maxpsforder,_order))
#define ORDER_2 (std::max(maxpsforder,_order+2))
#define ORDER_3 (std::max(maxpsforder+2,_order))
#define SIZE_1 (ORDER_1+1)*(ORDER_1+2)/2
#define SIZE_1b (ORDER_1+3)*(ORDER_1+4)/2
#define SIZE_2 (ORDER_2+1)*(ORDER_2+2)/2
#define SIZE_3 (ORDER_3+1)*(ORDER_3+2)/2
#define SIZE_4 (maxpsforder+1)*(maxpsforder+2)/2
#define SIZE_5 (maxpsforder+3)*(maxpsforder+4)/2
ESImpl::ESImpl(const std::vector<PixelList>& _pix, 
    const std::vector<BVec>& _psf,
    double _fp, int _order, double _sigma, 
    bool _fixcen, bool _fixgam, bool _fixmu, bool _useflux) :
  nsize((_order+1)*(_order+2)/2),
  np2size((_order+3)*(_order+4)/2), b(_order,_sigma),
  pix(_pix), psf(&_psf), f_psf(_fp), sigma_obs(pix.size()),
  I(SumSize(pix)), Z(I.size()), Z1(I.size()), W(I.size()), 
  psi(I.size(),nsize,0.), A_aux(I.size(),np2size), AC(I.size(),nsize),
  A(A_aux.Cols(0,nsize)), 
  fixcen(_fixcen), fixgam(_fixgam), fixmu(_fixmu), useflux(_useflux),
  numeric_j(false), U((fixcen?0:2)+(fixgam?0:2)+(fixmu?0:1),5,0.),
  xinit(5,0.), xx(5), ff(6), dff(6,5),
  x_short(U.colsize()), f_short(U.colsize()+(useflux?1:0)),
  Iresid(I.size()), Cb(nsize), Gb(np2size,5),
  dbdE(nsize,5), dAdEb(I.size(),5), temp(nsize,5), AtI(np2size),
  dCdE(nsize,nsize), flux(0), maxpsforder(CalcMaxPSFOrder(_psf)),
  _Gx(np2size,nsize), _Gy(np2size,nsize),
  _Gg1(SIZE_1b,SIZE_1), _Gg2(SIZE_1b,SIZE_1),
  _Gth(SIZE_1b,SIZE_1), _Gmu(SIZE_2,SIZE_3),
  _dTdE(SIZE_4,SIZE_4)
{
  for(size_t k=0,n=0;k<pix.size();k++) {
    for(size_t i=0;i<pix[k].size();i++,n++) {
      I(n) = pix[k][i].I*pix[k][i].wt;
      W(n) = pix[k][i].wt;
      Z(n) = pix[k][i].z;
    }
  }

  int maxpsforder = 0;
  C.reserve(pix.size());
  S.reserve(pix.size());
  D.reserve(pix.size());
  for(size_t k=0;k<pix.size();k++) {
    int psfsize = (*psf)[k].size();
    int psfsize2 = ((*psf)[k].GetOrder()+3)*((*psf)[k].GetOrder()+4)/2;
    C.push_back(tmv::Matrix<double>(nsize,nsize));
    S.push_back(tmv::Matrix<double>(psfsize,psfsize2));
    D.push_back(tmv::Matrix<double>(psfsize2,psfsize));
    if ((*psf)[k].GetOrder() > maxpsforder)
      maxpsforder = (*psf)[k].GetOrder();
    sigma_obs[k] = sqrt(pow(b.GetSigma(),2)+f_psf*pow((*psf)[k].GetSigma(),2));
  }

  AC.DivideUsing(tmv::QR);
  AC.SaveDiv();
  Assert(AC.colsize() >= AC.rowsize());
  // Otherwise there are too few pixels to constrain model.

  int k=0;
  if (!fixcen) { U(k++,0) = 1.; U(k++,1) = 1.; }
  if (!fixgam) { U(k++,2) = 1.; U(k++,3) = 1.; }
  if (!fixmu) { U(k,4) = 1.; }

  SetupGx(_Gx,_order+2,_order);
  SetupGy(_Gy,_order+2,_order);
  SetupGg1(_Gg1,ORDER_1+2,ORDER_1);
  SetupGg2(_Gg2,ORDER_1+2,ORDER_1);
  SetupGth(_Gth,ORDER_1+2,ORDER_1);
  SetupGmu(_Gmu,ORDER_2,ORDER_3);
}
#undef ORDER_1
#undef ORDER_2
#undef ORDER_3klsdfhjlkajsdhflkjashdflkjahsdflkjhasdlkfjhalksjdhflkashdflkashdflkjhasdflkjh
#undef SIZE_1
#undef SIZE_1b
#undef SIZE_2
#undef SIZE_3

void ESImpl::DoF(const tmv::Vector<double>& x, tmv::Vector<double>& f) const
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
  if (NormInf(x-xinit) > 4.) 
  {
    f = b.SubVector(0,6);
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
  if (gsq > 0.99 || (mu < -3. && Norm(x-xinit) > 0.3))
  {
    f = b.SubVector(0,6);
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
  Z1 -= m0*g*Z.Conjugate();
  Z1.AddToAll(m0*(g*std::conj(zc)-zc));
  for(size_t k=0,n=0,nx;k<pix.size();k++,n=nx) {
    nx = n+pix[k].size();
    Z1.SubVector(n,nx) /= sigma_obs[k];
  }
  MakePsi(A,Z1,b.GetOrder(),&W);

  if (psf) {
    for(size_t k=0,n=0,nx;k<pix.size();k++,n=nx) {
      nx = n+pix[k].size();
      BVec newpsf = (*psf)[k];
      const int psfsize = newpsf.size();
      tmv::MatrixView<double> S1 = S[k].Cols(0,psfsize);
      GTransform(g,newpsf.GetOrder(),S1);
      newpsf = S1 * (*psf)[k];
      tmv::MatrixView<double> D1 = D[k].Rows(0,psfsize);
      MuTransform(mu,newpsf.GetOrder(),D1);
      newpsf = D1 * newpsf;
      PSFConvolve(newpsf,b.GetOrder(),b.GetSigma(),C[k]);
      AC.Rows(n,nx) = A.Rows(n,nx) * C[k];
    }
    AC.ReSetDiv();
    if (AC.Singular()) {
      // I used to keep going and just give an error. But these never 
      // eventually work, so just bail on the whole calculation.
      dbg<<"Found singular AC.  Bail.\n";
      throw tmv::Singular();
    }
    b = I/AC;
  } else {
    A.ReSetDiv();
    if (A.Singular()) {
      dbg<<"Found singular A.  Bail.\n";
      throw tmv::Singular();
    }
    b = I/A;
  }

  f = b.SubVector(0,6);
  if (flux == 0.) flux = b(0);
  f /= flux;
  if (useflux) f(0) -= 1.;
}

void ESImpl::F(const tmv::Vector<double>& x, tmv::Vector<double>& f) const
{
  Assert(x.size() == U.colsize());
  if (useflux) Assert(f.size() == U.colsize()+1);
  else Assert(f.size() == U.colsize());

  Assert(xx.size() == U.rowsize());
  xx = x*U;
  if (fixcen) { xx[0] = fixuc;  xx[1] = fixvc; }
  if (fixgam) { xx[2] = fixg1;  xx[3] = fixg2; }
  if (fixmu) { xx[4] = fixm; }

  DoF(xx,ff);

  if (useflux) {
    f(0) = ff(0);
    f.SubVector(1,f.size()) = U*ff.SubVector(1,6);
  } else {
    f = U*ff.SubVector(1,6);
  }
}

void ESImpl::DoJ(const tmv::Vector<double>& x, const tmv::Vector<double>& ,
    tmv::Matrix<double>& df) const
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
  // So dC/dE is obtained from the PSFConvolve routine,
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
  Assert(df.rowsize() == 5);
  Assert(df.colsize() == 6);

  std::complex<double> zc(x[0],x[1]);
  std::complex<double> g(x[2],x[3]);
  double fact = 1./(1.-std::norm(g)); // used below
  double mu = x[4];
  double m0 = exp(-mu)*sqrt(fact);

  if (psf) Iresid = I-AC*b;
  else Iresid = I-A*b;

  AugmentPsi(A_aux,Z1,b.GetOrder());

#ifdef JTEST
  tmv::Matrix<double> A_aux2(Z1.size(),np2size);
  MakePsi(A_aux2.View(),Z1,b.GetOrder()+2,&W);
  Test(Norm(A_aux2-A_aux) < 1.e-8*Norm(A_aux),"A_aux");
#endif

  double g1 = std::real(g);
  double g2 = std::imag(g);

  tmv::ConstMatrixView<double> Gx = _Gx.SubMatrix(0,np2size,0,nsize);
  tmv::ConstMatrixView<double> Gy = _Gy.SubMatrix(0,np2size,0,nsize);
  tmv::ConstMatrixView<double> Gg1 = _Gg1.SubMatrix(0,np2size,0,nsize);
  tmv::ConstMatrixView<double> Gg2 = _Gg2.SubMatrix(0,np2size,0,nsize);
  tmv::ConstMatrixView<double> Gmu = _Gmu.SubMatrix(0,np2size,0,nsize);
  tmv::ConstMatrixView<double> Gth = _Gth.SubMatrix(0,np2size,0,nsize);

#ifdef JTEST
  // This section was only used for debugging purposes, but
  // I'm leaving it in, since the derivations of Gx[i] are
  // helpful in understanding some of the code following this section.
  tmv::Vector<double> f0(6);
  DoF(x,f0);
  tmv::Matrix<double> df_num(6,5);
  tmv::Matrix<double> d2f_num(6,5);
  tmv::Vector<double> f1(6);
  tmv::Vector<double> f2(6);
  double dE = 1.e-8;
  for(int j=0;j<5;j++) {
    tmv::Vector<double> x1 = x;
    x1(j) -= dE;
    DoF(x1,f1);
    tmv::Vector<double> x2 = x;
    x2(j) += dE;
    DoF(x2,f2);
    df_num.col(j) = (f2-f1)/(2.*dE);
    d2f_num.col(j) = (f2+f1-2.*f0)/(dE*dE);
  }
  dbg<<"numerical df: "<<df_num<<std::endl;
  dbg<<"numerical d2f: "<<d2f_num<<std::endl;
  dbg<<"(f2-f0)/dE = "<<(f2-f0)/dE<<std::endl;
  dbg<<"(f0-f1)/dE = "<<(f0-f1)/dE<<std::endl;
  dbg<<"diff = "<<(f2-f0-f0+f1)/dE<<std::endl;
  DoF(x,f0); // return A to original value

  tmv::Matrix<double> dAdEb_num(I.size(),5);
  tmv::Matrix<double> dAdEtI_num(nsize,5);
  for(int j=0;j<5;j++) {
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
    for(size_t k=0,n=0,nx;k<pix.size();k++,n=nx) {
      nx = n+pix[k].size();
      double m = m0/sigma_obs[k];
      z0.SubVector(n,nx) = m*Z.SubVector(n,nx);
      z0.SubVector(n,nx) -= m*g*Z.Conjugate().SubVector(n,nx);
      z0.SubVector(n,nx).AddToAll(m*(g*std::conj(zc)-zc));
      double m_1 = exp(-mu_1)/sigma_obs[k]/sqrt(1.-std::norm(g_1));
      z1.SubVector(n,nx) = m_1*Z.SubVector(n,nx);
      z1.SubVector(n,nx) -= m_1*g_1*Z.Conjugate().SubVector(n,nx);
      z1.SubVector(n,nx).AddToAll(m_1*(g_1*std::conj(zc_1)-zc_1));
      double m_2 = exp(-mu_2)/sigma_obs[k]/sqrt(1.-std::norm(g_2));
      z2.SubVector(n,nx) = m_2*Z.SubVector(n,nx);
      z2.SubVector(n,nx) -= m_2*g_2*Z.Conjugate().SubVector(n,nx);
      z2.SubVector(n,nx).AddToAll(m_2*(g_2*std::conj(zc_2)-zc_2));
    }
    MakePsi(A0.View(),z0,b.GetOrder(),&W);
    MakePsi(A1.View(),z1,b.GetOrder(),&W);
    MakePsi(A2.View(),z2,b.GetOrder(),&W);
    tmv::Matrix<double> dAdE_num = (A2-A1)/(2.*dE);
    tmv::Matrix<double> d2A_num = (A2+A1-2.*A0)/(dE*dE);

    if (psf) {
      for(size_t k=0,n=0,nx;k<pix.size();k++,n=nx) {
	nx = n+pix[k].size();
	A0.Rows(n,nx) *= C[k];

	tmv::Matrix<double> C1(nsize,nsize);
	BVec newpsf1 = (*psf)[k];
	ApplyG(g_1,newpsf1);
	ApplyMu(mu_1,newpsf1);
	PSFConvolve(newpsf1,b.GetOrder(),b.GetSigma(),C1);
	A1.Rows(n,nx) *= C1;

	tmv::Matrix<double> C2(nsize,nsize);
	BVec newpsf2 = (*psf)[k];
	ApplyG(g_2,newpsf2);
	ApplyMu(mu_2,newpsf2);
	PSFConvolve(newpsf2,b.GetOrder(),b.GetSigma(),C2);
	A2.Rows(n,nx) *= C2;
      }
    }

    tmv::Matrix<double> dACdE_num = (A2-A1)/(2.*dE);
    tmv::Matrix<double> d2AC_num = (A2+A1-2.*A0)/(dE*dE);

    tmv::Matrix<double> dACdE(I.size(),nsize);
    tmv::Matrix<double> dAdE(I.size(),nsize);

    for(size_t k=0,n=0,nx;k<pix.size();k++,n=nx) {
      nx = n+pix[k].size();
      double m = m0/sigma_obs[k];

      std::vector<tmv::Matrix<double> > GG(5,
	  tmv::Matrix<double>(np2size,nsize));

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

      dAdE.Rows(n,nx) = A_aux.Rows(n,nx)*GG[j];
      if (n == 0) {
	dbg<<"A_aux = "<<A_aux.Rows(0,5)<<std::endl;
	dbg<<"GG["<<j<<"] = "<<GG[j]<<std::endl;
	dbg<<"dAdE = "<<dAdE.Rows(0,5)<<std::endl;
	dbg<<"j = "<<j<<std::endl;
      }
      if (psf) {
	dACdE.Rows(n,nx) = dAdE.Rows(n,nx)*C[k];

	if (j>1) {
	  tmv::Matrix<double> C1(nsize,nsize);
	  BVec newpsf1 = (*psf)[k];
	  ApplyG(g_1,newpsf1);
	  ApplyMu(mu_1,newpsf1);
	  PSFConvolve(newpsf1,b.GetOrder(),b.GetSigma(),C1);

	  tmv::Matrix<double> C2(nsize,nsize);
	  BVec newpsf2 = (*psf)[k];
	  ApplyG(g_2,newpsf2);
	  ApplyMu(mu_2,newpsf2);
	  PSFConvolve(newpsf2,b.GetOrder(),b.GetSigma(),C2);

	  // For j=0,1 dC/dE = 0, so dA/dE is just AGC
	  // If not, we need to add A(dC/dE)
	  tmv::Matrix<double> dCdE_num = (C2-C1)/(2.*dE);
	  tmv::Matrix<double> d2C_num = (C2+C1-2.*C[k])/(dE*dE);

	  int n1 = ((*psf)[k].GetOrder()+1)*((*psf)[k].GetOrder()+2)/2;
	  int n2 = ((*psf)[k].GetOrder()+3)*((*psf)[k].GetOrder()+4)/2;
	  tmv::Matrix<double> S0(n1,n1);
	  tmv::Matrix<double> D0(n1,n1);
	  GTransform(g,(*psf)[k].GetOrder(),S0);
	  MuTransform(mu,(*psf)[k].GetOrder(),D0);
	  tmv::Matrix<double> T0 = D0*S0;
	  tmv::Matrix<double> S1(n1,n1);
	  tmv::Matrix<double> D1(n1,n1);
	  GTransform(g_1,(*psf)[k].GetOrder(),S1);
	  MuTransform(mu_1,(*psf)[k].GetOrder(),D1);
	  tmv::Matrix<double> T1 = D1*S1;
	  tmv::Matrix<double> D2(n1,n1);
	  tmv::Matrix<double> S2(n1,n1);
	  GTransform(g_2,(*psf)[k].GetOrder(),S2);
	  MuTransform(mu_2,(*psf)[k].GetOrder(),D2);
	  tmv::Matrix<double> T2 = D2*S2;
	  tmv::Matrix<double> dTdE_num = (T2-T1)/(2.*dE);
	  tmv::Matrix<double> d2T_num = (T2+T1-2.*T0)/(dE*dE);

	  tmv::Matrix<double> dTdE(n1,n1);
	  if (j==2) {
	    tmv::Matrix<double> Sx(n2,n2);
	    GTransform(g,(*psf)[k].GetOrder()+2,Sx);
	    tmv::Matrix<double> dSdE_num = (S2-S1)/(2.*dE);
	    tmv::Matrix<double> d2S_num = (S2+S1-2.*S0)/(dE*dE);
	    dbg<<"Gg1 size = "<<_Gg1.colsize()<<"  "<<_Gg1.rowsize()<<std::endl;
	    dbg<<"Gth size = "<<_Gth.colsize()<<"  "<<_Gth.rowsize()<<std::endl;
	    dbg<<"n1 = "<<n1<<", n2 = "<<n2<<std::endl;
	    tmv::Matrix<double> dSdE = fact*Sx.Rows(0,n1)*_Gg1.SubMatrix(0,n2,0,n1)+fact*g2*Sx.Rows(0,n1)*_Gth.SubMatrix(0,n2,0,n1);
	    dbg<<"Norm(dS/dE - numeric dS/dE) = "<<Norm(dSdE-dSdE_num)<<std::endl;
	    Test(Norm(dSdE-dSdE_num) < 10.*dE*Norm(d2S_num),"dSdE");
	    dTdE = fact*D[k].Rows(0,n1)*Sx.Rows(0,n1)*_Gg1.SubMatrix(0,n2,0,n1)+fact*g2*D[k].Rows(0,n1)*Sx.Rows(0,n1)*_Gth.SubMatrix(0,n2,0,n1);
	  } else if (j==3) {
	    tmv::Matrix<double> Sx(n2,n2);
	    GTransform(g,(*psf)[k].GetOrder()+2,Sx);
	    tmv::Matrix<double> dSdE_num = (S2-S1)/(2.*dE);
	    tmv::Matrix<double> d2S_num = (S2+S1-2.*S0)/(dE*dE);
	    tmv::Matrix<double> dSdE = fact*Sx.Rows(0,n1)*_Gg2.SubMatrix(0,n2,0,n1)-fact*g1*Sx.Rows(0,n1)*_Gth.SubMatrix(0,n2,0,n1);
	    dbg<<"Norm(dS/dE - numeric dS/dE) = "<<Norm(dSdE-dSdE_num)<<std::endl;
	    Test(Norm(dSdE-dSdE_num) < 10.*dE*Norm(d2S_num),"dSdE");
	    dTdE = fact*D[k].Rows(0,n1)*Sx.Rows(0,n1)*Gg2.SubMatrix(0,n2,0,n1)-fact*g1*D[k].Rows(0,n1)*Sx.Rows(0,n1)*Gth.SubMatrix(0,n2,0,n1);
	  } else if (j==4) {
	    tmv::Matrix<double> Dx(n2,n2);
	    MuTransform(mu,(*psf)[k].GetOrder()+2,Dx);
	    tmv::Matrix<double> dDdE_num = (D2-D1)/(2.*dE);
	    tmv::Matrix<double> d2D_num = (D2+D1-2.*D0)/(dE*dE);
	    tmv::Matrix<double> dDdE = _Gmu.SubMatrix(0,n1,0,n2)*Dx.Cols(0,n1);
	    dbg<<"Norm(dD/dE - numeric dD/dE) = "<<Norm(dDdE-dDdE_num)<<std::endl;
	    Test(Norm(dDdE-dDdE_num) < 10.*dE*Norm(d2D_num),"dDdE");
	    dTdE = _Gmu.SubMatrix(0,n1,0,n2)*Dx.Cols(0,n1)*S[k].Cols(0,n1);
	  }
	  dbg<<"Norm(dT/dE - numeric dT/dE) = "<<Norm(dTdE-dTdE_num)<<std::endl;
	  Test(Norm(dTdE-dTdE_num) < 10.*dE*Norm(dTdE_num),"dTdE");

	  BVec dbpsfdE((*psf)[k].GetOrder(),(*psf)[k].GetSigma());
	  dbpsfdE = dTdE*(*psf)[k];
	  tmv::Matrix<double> dCdE(nsize,nsize);
	  PSFConvolve(dbpsfdE,b.GetOrder(),b.GetSigma(),dCdE);
	  dbg<<"Norm(dC/dE - numeric dC/dE) = "<<Norm(dCdE-dCdE_num)<<"\n";
	  Test(Norm(dCdE-dCdE_num) < 10.*dE*Norm(d2C_num),"dCdE");
	  dACdE.Rows(n,nx) += A.Rows(n,nx)*dCdE;
	}
      }
    }
    dbg<<"Norm(dA/dE - numeric dA/dE) = "<<Norm(dAdE-dAdE_num)<<std::endl;
    dbg<<"Norm(dA/dE) = "<<Norm(dAdE)<<std::endl;
    dbg<<"Norm(d2A/dE2) = "<<Norm(d2A_num)<<std::endl;
    if (Norm(dAdE-dAdE_num) > 10.*dE*Norm(d2A_num)) {
      dbg<<"dAdE = "<<dAdE.Rows(0,5)<<std::endl;
      dbg<<"dAdE_num = "<<dAdE_num.Rows(0,5)<<std::endl;
      dbg<<"diff = "<<(dAdE-dAdE_num).Rows(0,5)<<std::endl;
    }
    Test(Norm(dAdE-dAdE_num) < 10.*dE*Norm(d2A_num),"dAdE");

    if (psf) {
      dbg<<"Norm(dAC/dE - numeric dAC/dE) = "<<Norm(dACdE-dACdE_num)<<std::endl;
      dbg<<"Norm(dAC/dE) = "<<Norm(dACdE)<<std::endl;
      dbg<<"Norm(d2A/dE2) = "<<Norm(d2A_num)<<std::endl;
      if (Norm(dACdE-dACdE_num) > 10.*dE*Norm(d2A_num)) {
	dbg<<"dACdE = "<<dACdE<<std::endl;
	dbg<<"dACdE_num = "<<dACdE_num<<std::endl;
	dbg<<"diff = "<<dACdE-dACdE_num<<std::endl;
      }
      Test(Norm(dACdE-dACdE_num) < 10.*dE*Norm(d2AC_num),"dACdE");
      dAdE = dACdE;
    }

    dAdEb_num.col(j) = dAdE * b;
    dAdEtI_num.col(j) = dAdE.Transpose() * Iresid;

    // (At A) db/dE = (dA/dE)t (I - A b) - At (dA/dE) b
    tmv::Vector<double> dbdE_calc = dAdE.Transpose() * Iresid;
    dbg<<"dAdEt I = "<<dbdE_calc<<std::endl;
    double kappa;
    if (psf) {
      dbdE_calc -= AC.Transpose() * dAdE * b;
      dbg<<"(ACt) dAdE b = "<<AC.Transpose()*dAdE*b<<std::endl;
      dbg<<"dbdE_calc => "<<dbdE_calc<<std::endl;
      dbg<<"R.diag() = "<<AC.QRD().GetR().diag()<<std::endl;
      dbdE_calc /= AC.QRD().GetR().Transpose();
      dbg<<"(/= Rt) dbdE_calc => "<<dbdE_calc<<std::endl;
      dbdE_calc /= AC.QRD().GetR();
      dbg<<"(/= R) dbdE_calc => "<<dbdE_calc<<std::endl;
      kappa = AC.QRD().GetR().DoCondition();
    } else {
      dbdE_calc -= A.Transpose() * dAdE * b;
      dbdE_calc /= A.QRD().GetR().Transpose();
      dbdE_calc /= A.QRD().GetR();
      kappa = A.QRD().GetR().DoCondition();
    }
    tmv::Vector<double> b0 = I/A0;
    tmv::Vector<double> b1 = I/A1;
    tmv::Vector<double> b2 = I/A2;
    tmv::Vector<double> dbdE_num = (b2-b1)/(2.*dE);
    tmv::Vector<double> d2b_num = (b2+b1-2.*b0)/(dE*dE);
    dbg<<"Norm(db/dE_calc - numeric db/dE) = "<<Norm(dbdE_calc-dbdE_num)<<std::endl;
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
    Test(Norm(dbdE_calc-dbdE_num) < 10.*dE*kappa*kappa*Norm(d2b_num),"dbdE_calc");
    df.col(j) = dbdE_calc.SubVector(0,6);
  }
  dbg<<"Semi-empirical: "<<df<<std::endl;
  dbg<<"Norm(diff) = "<<Norm(df-df_num)<<std::endl;
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

  dbdE.Zero();

  for(size_t k=0,n=0,nx;k<pix.size();k++,n=nx) {
    nx = n + pix[k].size();
    double m = m0/sigma_obs[k];
    if (psf) {
      Cb = C[k] * b;

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

      dAdEb.Rows(n,nx) = A_aux.Rows(n,nx) * Gb;
    } else {
      // GG[0] = -m*((1.-g1)*Gx - g2*Gy);
      Gb.col(0) = -m * (1.-g1) * Gx * b;
      Gb.col(0) += m * g2 * Gy * b;

      // GG[1] = -m*((1.+g1)*Gy - g2*Gx);
      Gb.col(1) = -m * (1.+g1) * Gy * b;
      Gb.col(1) += m * g2 * Gx * b;

      // GG[2] = -fact*(Gg1 - g2*Gth);
      Gb.col(2) = -fact * Gg1 * b;
      Gb.col(2) += fact * g2 * Gth * b;

      // GG[3] = -fact*(Gg2 + g1*Gth);
      Gb.col(3) = -fact * Gg2 * b;
      Gb.col(3) -= fact * g1 * Gth * b;

      // GG[4] = -Gmu;
      Gb.col(4) = -Gmu * b;

      dAdEb.Rows(n,nx) = A_aux.Rows(n,nx) * Gb;
    }

    AtI = A_aux.Rows(n,nx).Transpose()*Iresid.SubVector(n,nx);

    // Here, temp = (dAdE)t * I = GtAtI
    // GtAtI.col(i) = GG[i].Transpose() * A_aux.Transpose() * Iresid
    temp.col(0) = -m*(1.-g1) * Gx.Transpose() * AtI;
    temp.col(0) += m*g2 * Gy.Transpose() * AtI;

    temp.col(1) = -m*(1.+g1) * Gy.Transpose() * AtI;
    temp.col(1) += m*g2 * Gx.Transpose() * AtI;

    temp.col(2) = -fact * Gg1.Transpose() * AtI;
    temp.col(2) += fact*g2 * Gth.Transpose() * AtI;

    temp.col(3) = -fact * Gg2.Transpose() * AtI;
    temp.col(3) -= fact*g1 * Gth.Transpose() * AtI;

    temp.col(4) = -Gmu.Transpose() * AtI;

    if (psf) {
      // with a psf, dAdE = AGC + A(dC/dE)
      // AGCb is already accounted for correctly by doing the Cb thing above.
      // for (dAdE)t, CtGtAtI is effected by the following line:
      size_t psfsize = (*psf)[k].size();

      temp = C[k].Transpose() * temp;

      // The rest of this section does the A(dC/dE) terms, which
      // only exist for E = {g1,g2,mu}.  ie. not E = {x,y}.

      int n2 = ((*psf)[k].GetOrder()+3)*((*psf)[k].GetOrder()+4)/2;
      tmv::MatrixView<double> dTdE = _dTdE.SubMatrix(0,psfsize,0,psfsize);
      tmv::ConstMatrixView<double> D1 = D[k].Rows(0,psfsize);
      tmv::ConstMatrixView<double> S1 = S[k].Cols(0,psfsize);

      AugmentGTransformCols(g,(*psf)[k].GetOrder(),S[k].View());
      AugmentMuTransformRows(mu,(*psf)[k].GetOrder(),D[k].View());

#ifdef JTEST
      tmv::Matrix<double> S2(n2,n2);
      GTransform(g,(*psf)[k].GetOrder()+2,S2.View());
      dbg<<"S[k] = "<<S[k]<<std::endl;
      dbg<<"S2 = "<<S2.Rows(0,psfsize)<<std::endl;
      dbg<<"diff = "<<S[k]-S2.Rows(0,psfsize)<<std::endl;
      dbg<<"Norm(S[k]-S2) = "<<Norm(S[k]-S2.Rows(0,psfsize))<<std::endl;
      dbg<<"Norm(S2) = "<<Norm(S2)<<std::endl;
      Test(Norm(S[k]-S2.Rows(0,psfsize)) < 1.e-8*Norm(S2),"AugmentG");
      tmv::Matrix<double> D2(n2,n2);
      MuTransform(mu,(*psf)[k].GetOrder()+2,D2.View());
      dbg<<"Norm(diff) = "<<Norm(D[k]-D2.Cols(0,psfsize))<<std::endl;
      dbg<<"Norm(D2) = "<<Norm(D2)<<std::endl;
      Test(Norm(D[k]-D2.Cols(0,psfsize)) < 1.e-8*Norm(D2),"AugmentMu");
#endif

      // E = g1
      dTdE = fact*D1*S[k]*_Gg1.SubMatrix(0,n2,0,psfsize);
      dTdE += fact*g2*D1*S[k]*_Gth.SubMatrix(0,n2,0,psfsize);
      BVec dbpsfdE((*psf)[k].GetOrder(),(*psf)[k].GetSigma());
      dbpsfdE = dTdE*(*psf)[k];
      PSFConvolve(dbpsfdE,b.GetOrder(),b.GetSigma(),dCdE);
      dAdEb.col(2,n,nx) += A.Rows(n,nx)*(dCdE*b);
      temp.col(2) += dCdE.Transpose() * AtI.SubVector(0,nsize);

      // E = g2
      dTdE = fact*D1*S[k]*_Gg2.SubMatrix(0,n2,0,psfsize);
      dTdE -= fact*g1*D1*S[k]*_Gth.SubMatrix(0,n2,0,psfsize);
      dbpsfdE = dTdE*(*psf)[k];
      PSFConvolve(dbpsfdE,b.GetOrder(),b.GetSigma(),dCdE);
      dAdEb.col(3,n,nx) += A.Rows(n,nx)*(dCdE*b);
      temp.col(3) += dCdE.Transpose() * AtI.SubVector(0,nsize);

      // E = mu
      dTdE = _Gmu.SubMatrix(0,psfsize,0,n2)*D[k]*S1;
      dbpsfdE = dTdE*(*psf)[k];
      PSFConvolve(dbpsfdE,b.GetOrder(),b.GetSigma(),dCdE);
      dAdEb.col(4,n,nx) += A.Rows(n,nx)*(dCdE*b);
      temp.col(4) += dCdE.Transpose() * AtI.SubVector(0,nsize);
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
    dbdE /= AC.QRD().GetR().Transpose();
    //dbdE -= AC.QRD().GetQ().Transpose() * dAdEb;
    temp = dAdEb / AC.QRD().GetQ();
    dbdE -= temp;
    dbdE /= AC.QRD().GetR();
  } else {
    dbdE /= A.QRD().GetR().Transpose();
    //dbdE -= A.QRD().GetQ().Transpose() * dAdEb;
    temp = dAdEb / A.QRD().GetQ();
    dbdE -= temp;
    dbdE /= A.QRD().GetR();
  }
  df = dbdE.Rows(0,6);
  df /= flux;

#ifdef JTEST
  dbg<<"analytic: "<<df.Clip(1.e-5)<<std::endl;
  dbg<<"Norm(diff) = "<<Norm(df-df_num)<<std::endl;
#endif
}

void ESImpl::J(const tmv::Vector<double>& x, const tmv::Vector<double>& f,
    tmv::Matrix<double>& df) const
{
  Assert(x.size() == U.colsize());
  Assert(df.rowsize() == U.colsize());
  if (useflux) {
    Assert(df.colsize() == U.colsize()+1);
  } else {
    Assert(df.colsize() == U.colsize());
  }

  xx = x*U;
  if (fixcen) { xx[0] = fixuc;  xx[1] = fixvc; }
  if (fixgam) { xx[2] = fixg1;  xx[3] = fixg2; }
  if (fixmu) { xx[4] = fixm; }

  DoJ(xx,f,dff);
  if (useflux) {
    df.row(0) = dff.row(0)*U.Transpose();
    df.Rows(1,df.colsize()) = U*dff.Rows(1,6)*U.Transpose();
  } else {
    df = U*dff.Rows(1,6)*U.Transpose();
  }
}

void EllipseSolver::UseNumericJ() { pimpl->numeric_j = true; }

const BVec& EllipseSolver::GetB() const { return pimpl->b; }

void EllipseSolver::GetBCov(tmv::Matrix<double>& bcov) const 
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
  // This is what AC.InverseATA returns.
  if (pimpl->psf) pimpl->AC.InverseATA(bcov);
  else pimpl->A.InverseATA(bcov);
}

void EllipseSolver::getCovariance(tmv::Matrix<double>& cov) const 
{
  tmv::Matrix<double> cov1(pimpl->x_short.size(),pimpl->x_short.size());
  NLSolver::getCovariance(cov1);
  cov = pimpl->U.Transpose()*cov1*pimpl->U;
  dbg<<"getCovariance:\n";
  dbg<<"cov1 = "<<cov1<<std::endl;
  dbg<<"full cov = "<<cov<<std::endl;
}

void EllipseSolver::getInverseCovariance(tmv::Matrix<double>& invcov) const 
{
  tmv::Matrix<double> invcov1(pimpl->x_short.size(),pimpl->x_short.size());
  NLSolver::getInverseCovariance(invcov1);
  invcov = pimpl->U.Transpose()*invcov1*pimpl->U;
  dbg<<"getCovariance:\n";
  dbg<<"invcov1 = "<<invcov1<<std::endl;
  dbg<<"full invcov = "<<invcov<<std::endl;
}

void EllipseSolver::CallF(const tmv::Vector<double>& x,
    tmv::Vector<double>& f) const
{
  Assert(x.size() == 5);
  Assert(f.size() == 5);
  pimpl->xinit = x;
  if (pimpl->fixcen) { pimpl->fixuc = x[0]; pimpl->fixvc = x[1]; }
  if (pimpl->fixgam) { pimpl->fixg1 = x[2]; pimpl->fixg2 = x[3]; }
  if (pimpl->fixmu) { pimpl->fixm = x[4]; }

  pimpl->x_short = pimpl->U * x;

  pimpl->flux = 0;

  F(pimpl->x_short,pimpl->f_short);

  if (pimpl->useflux) 
    f = pimpl->f_short.SubVector(1,pimpl->f_short.size()) * pimpl->U;
  else f = pimpl->f_short*pimpl->U;
}

bool EllipseSolver::Solve(tmv::Vector<double>& x, tmv::Vector<double>& f) const
{
  Assert(x.size() == 5);
  Assert(f.size() == 5);
  pimpl->xinit = x;
  if (pimpl->fixcen) { pimpl->fixuc = x[0]; pimpl->fixvc = x[1]; }
  if (pimpl->fixgam) { pimpl->fixg1 = x[2]; pimpl->fixg2 = x[3]; }
  if (pimpl->fixmu) { pimpl->fixm = x[4]; }

  pimpl->x_short = pimpl->U * x;

  pimpl->flux = 0;
  bool ret;
  try
  {
    ret = NLSolver::Solve(pimpl->x_short,pimpl->f_short);
  }
  catch (...)
  {
    ret = false;
  }

  x = pimpl->x_short*pimpl->U;
  if (pimpl->fixcen) { x[0] = pimpl->fixuc; x[1] = pimpl->fixvc; }
  if (pimpl->fixgam) { x[2] = pimpl->fixg1; x[3] = pimpl->fixg2; }
  if (pimpl->fixmu) { x[4] = pimpl->fixm; }
  if (pimpl->useflux) 
    f = pimpl->f_short.SubVector(1,pimpl->f_short.size()) * pimpl->U;
  else f = pimpl->f_short*pimpl->U;

  return ret;
}

void EllipseSolver::numericH(
    const tmv::Vector<double>& x, const tmv::Vector<double>& f,
    tmv::SymMatrix<double>& h) const
{
  Assert(x.size() == 5);
  Assert(f.size() == 5);
  Assert(h.size() == 5);
  pimpl->xinit = x;
  if (pimpl->fixcen) { pimpl->fixuc = x[0]; pimpl->fixvc = x[1]; }
  if (pimpl->fixgam) { pimpl->fixg1 = x[2]; pimpl->fixg2 = x[3]; }
  if (pimpl->fixmu) { pimpl->fixm = x[4]; }

  pimpl->x_short = pimpl->U * x;

  tmv::SymMatrix<double> h_short(pimpl->x_short.size(),0.);

  pimpl->flux = 0;
  NLSolver::numericH(pimpl->x_short,pimpl->f_short,h_short);

  h = tmv::SymMatrix<double>(tmv::Matrix<double>(
	pimpl->U.Transpose() * h_short * pimpl->U));
  if (pimpl->fixcen) { h(0,0) = h(1,1) = 1.0; }
  if (pimpl->fixgam) { h(2,2) = h(3,3) = 1.0; }
  if (pimpl->fixmu) { h(4,4) = 1.0; }
}

bool EllipseSolver::TestJ(const tmv::Vector<double>& x, tmv::Vector<double>& f,
    std::ostream* os, double relerr) const 
{
  Assert(x.size() == 5);
  Assert(f.size() == 5);
  pimpl->xinit = x;
  if (pimpl->fixcen) { pimpl->fixuc = x[0]; pimpl->fixvc = x[1]; }
  if (pimpl->fixgam) { pimpl->fixg1 = x[2]; pimpl->fixg2 = x[3]; }
  if (pimpl->fixmu) { pimpl->fixm = x[4]; }

  pimpl->x_short = pimpl->U * x;

  pimpl->flux = 0;
  return NLSolver::TestJ(pimpl->x_short,pimpl->f_short,os,relerr);
}

