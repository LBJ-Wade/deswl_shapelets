
#include "EllipseSolver.h"
#include "TMV.h"
#include "dbg.h"
#include "PsiHelper.h"

static int CalcMaxPSFOrder(const std::vector<BVec>& psf)
{
  int maxpsforder = 0;
  for(size_t k=0;k<psf.size();k++) 
    if (psf[k].GetOrder() > maxpsforder)
      maxpsforder = psf[k].GetOrder();
  return maxpsforder;
}

static size_t SumSize(const std::vector<PixelList>& v)
{
  size_t sum = 0;
  for(size_t i=0;i<v.size();i++) sum += v[i].size();
  return sum;
}

//
// EllipseSolver2
//

class ESImpl2 
{

  friend class EllipseSolver2;

  public :

    ESImpl2(const std::vector<PixelList>& _pix,
	int _order, double _sigma, double _pixscale,
	bool _fixcen, bool _fixgam, bool _fixmu, bool _useflux);

    ESImpl2(const std::vector<PixelList>& _pix,
	const std::vector<BVec>& _psf, double _fp,
	int _order, double _sigma, double _pixscale,
	bool _fixcen, bool _fixgam, bool _fixmu, bool _useflux);

    void F(const tmv::Vector<double>& x, tmv::Vector<double>& f) const;

    void J(const tmv::Vector<double>& x, const tmv::Vector<double>& f,
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
    mutable std::vector<tmv::Matrix<double> > A_aux;
    mutable std::vector<tmv::MatrixView<double> > A;
    mutable std::vector<tmv::Matrix<double> > C;
    mutable std::vector<tmv::Matrix<double> > S;
    mutable std::vector<tmv::Matrix<double> > D;
    bool fixcen, fixgam, fixmu, useflux, numeric_j;
    tmv::Matrix<double> U;
    mutable tmv::Vector<double> xinit;
    mutable tmv::Vector<double> xx;
    mutable tmv::Vector<double> x_short;
    mutable tmv::Vector<double> f_short;
    mutable tmv::Vector<double> b1;
    mutable tmv::Matrix<double> dbdE;
    mutable tmv::Matrix<double> dbdE1;
    mutable tmv::Vector<double> IA;
    mutable tmv::Matrix<double> dCdE;
    mutable double fixuc,fixvc,fixg1,fixg2,fixm,flux;
    double pixscale;
    tmv::DiagMatrix<double> normD;
    int maxpsforder;
    tmv::Matrix<double> _Gx;
    tmv::Matrix<double> _Gy;
    tmv::Matrix<double> _Gg1;
    tmv::Matrix<double> _Gg2;
    tmv::Matrix<double> _Gth;
    tmv::Matrix<double> _Gmu;
    mutable tmv::Matrix<double> _dTdE;
};

EllipseSolver2::EllipseSolver2(const std::vector<PixelList>& pix,
    int order, double sigma, double pixscale,
    bool fixcen, bool fixgam, bool fixmu, bool useflux) :
  pimpl(new ESImpl2(pix,order,sigma,pixscale,
	fixcen,fixgam,fixmu,useflux))
{}

EllipseSolver2::EllipseSolver2(const std::vector<PixelList>& pix,
    const std::vector<BVec>& psf, double fp,
    int order, double sigma, double pixscale,
    bool fixcen, bool fixgam, bool fixmu, bool useflux) :
  pimpl(new ESImpl2(pix,psf,fp,order,sigma,pixscale,
	fixcen,fixgam,fixmu,useflux))
{}

EllipseSolver2::~EllipseSolver2()
{ delete pimpl; }

void EllipseSolver2::F(
    const tmv::Vector<double>& x, tmv::Vector<double>& f) const
{ pimpl->F(x,f); }

void EllipseSolver2::J(
    const tmv::Vector<double>& x, const tmv::Vector<double>& f,
    tmv::Matrix<double>& df) const
{
  if (pimpl->numeric_j) NLSolver::J(x,f,df);
  else pimpl->J(x,f,df); 
}

ESImpl2::ESImpl2(const std::vector<PixelList>& _pix,
    int _order, double _sigma, double _pixscale,
    bool _fixcen, bool _fixgam, bool _fixmu, bool _useflux) :
  nsize((_order+1)*(_order+2)/2),
  np2size((_order+3)*(_order+4)/2), b(_order,_sigma),
  pix(_pix), psf(0), f_psf(1.), sigma_obs(pix.size(),_sigma), 
  I(SumSize(pix)), Z(I.size()), Z1(I.size()),
  fixcen(_fixcen), fixgam(_fixgam), fixmu(_fixmu), useflux(_useflux),
  numeric_j(false), U((fixcen?0:2)+(fixgam?0:2)+(fixmu?0:1),5,0.),
  xinit(5,0.), xx(5), x_short(U.colsize()), f_short(U.colsize()+(useflux?1:0)),
  b1(nsize), dbdE(nsize,5), dbdE1(nsize,5), IA(np2size),
  dCdE(0,0), pixscale(_pixscale), normD(b.size()), maxpsforder(0),
  _Gx(np2size,nsize), _Gy(np2size,nsize),
  _Gg1(np2size,nsize), _Gg2(np2size,nsize),
  _Gth(np2size,nsize), _Gmu(np2size,nsize),
  _dTdE(0,0)
{
  A_aux.reserve(pix.size());
  A.reserve(pix.size());
  for(size_t k=0,n=0;k<pix.size();k++) {
    for(size_t i=0;i<pix[k].size();i++,n++) {
      I(n) = pix[k][i].I;
      Z(n) = pix[k][i].z;
    }
    A_aux.push_back(tmv::Matrix<double>(pix[k].size(),np2size));
    A.push_back(A_aux[k].Cols(0,nsize));
    //dbg<<"Allocated A matrix:"<<std::endl;
    //dbg<<"A_aux[k].ptr = "<<A_aux[k].ptr()<<std::endl;
    //dbg<<"A[k].ptr = "<<A[k].ptr()<<std::endl;
    //dbg<<"A_aux[k].col0 = "<<A_aux[k].col(0)<<std::endl;
    //dbg<<"A[k].col0 = "<<A[k].col(0)<<std::endl;
  }

  int j=0;
  if (!fixcen) { U(j++,0) = 1.; U(j++,1) = 1.; }
  if (!fixgam) { U(j++,2) = 1.; U(j++,3) = 1.; }
  if (!fixmu) { U(j,4) = 1.; }

  SetupGx(_Gx,_order+2,_order);
  SetupGy(_Gy,_order+2,_order);
  SetupGg1(_Gg1,_order+2,_order);
  SetupGg2(_Gg2,_order+2,_order);
  SetupGth(_Gth,_order+2,_order);
  SetupGmu(_Gmu,_order+2,_order);

  double val = pixscale*pixscale/(_sigma*_sigma)/pix.size();
  for(int n=0,k=0;n<=b.GetOrder();++n) {
    for(int p=n,q=0;p>=q;--p,++q,++k) {
      if (p!=q) { normD(k) = val/2.; normD(++k) = val/2.; }
      else normD(k) = val;
    }
  }
}

#define ORDER_1 (std::max(maxpsforder,_order))
#define ORDER_2 (std::max(maxpsforder,_order+2))
#define ORDER_3 (std::max(maxpsforder+2,_order))
#define SIZE_1 (ORDER_1+1)*(ORDER_1+2)/2
#define SIZE_1b (ORDER_1+3)*(ORDER_1+4)/2
#define SIZE_2 (ORDER_2+1)*(ORDER_2+2)/2
#define SIZE_3 (ORDER_3+1)*(ORDER_3+2)/2
#define SIZE_4 (maxpsforder+1)*(maxpsforder+2)/2

ESImpl2::ESImpl2(const std::vector<PixelList>& _pix,
    const std::vector<BVec>& _psf, double _fp, int _order, double _sigma, 
    double _pixscale, bool _fixcen, bool _fixgam, bool _fixmu, bool _useflux) :
  nsize((_order+1)*(_order+2)/2),
  np2size((_order+3)*(_order+4)/2), b(_order,_sigma),
  pix(_pix), psf(&_psf), f_psf(_fp), sigma_obs(pix.size()), 
  I(SumSize(pix)), Z(I.size()), Z1(I.size()),
  fixcen(_fixcen), fixgam(_fixgam), fixmu(_fixmu), useflux(_useflux),
  numeric_j(false), U((fixcen?0:2)+(fixgam?0:2)+(fixmu?0:1),5,0.),
  xinit(5,0.), xx(5), x_short(U.colsize()), f_short(U.colsize()+(useflux?1:0)),
  b1(nsize), dbdE(nsize,5), dbdE1(nsize,5), IA(np2size),
  dCdE(nsize,nsize), pixscale(_pixscale), normD(b.size()),
  maxpsforder(CalcMaxPSFOrder(_psf)),
  _Gx(np2size,nsize), _Gy(np2size,nsize),
  _Gg1(SIZE_1b,SIZE_1), _Gg2(SIZE_1b,SIZE_1),
  _Gth(SIZE_1b,SIZE_1), _Gmu(SIZE_2,SIZE_3),
  _dTdE(SIZE_4,SIZE_4)

{
  int maxpsforder = 0;
  A_aux.reserve(pix.size());
  A.reserve(pix.size());
  C.reserve(pix.size());
  S.reserve(pix.size());
  D.reserve(pix.size());
  for(size_t k=0,n=0,nx;k<pix.size();k++,n=nx) {
    nx = n+pix[k].size();

    for(size_t i=0,n1=n;i<pix[k].size();i++,n1++) {
      I(n1) = pix[k][i].I;
      Z(n1) = pix[k][i].z;
    }

    int psfsize = (*psf)[k].size();
    int psfsize2 = ((*psf)[k].GetOrder()+3)*((*psf)[k].GetOrder()+4)/2;
    A_aux.push_back(tmv::Matrix<double,tmv::ColMajor>(pix[k].size(),np2size));
    A.push_back(A_aux[k].Cols(0,nsize));
    //dbg<<"Allocated A matrix:"<<std::endl;
    //dbg<<"A_aux[k].ptr = "<<A_aux[k].ptr()<<std::endl;
    //dbg<<"A[k].ptr = "<<A[k].ptr()<<std::endl;
    //dbg<<"A_aux[k].col0 = "<<A_aux[k].col(0)<<std::endl;
    //dbg<<"A[k].col0 = "<<A[k].col(0)<<std::endl;
    C.push_back(tmv::Matrix<double>(nsize,nsize));
    S.push_back(tmv::Matrix<double>(psfsize,psfsize2));
    D.push_back(tmv::Matrix<double>(psfsize2,psfsize));
    if ((*psf)[k].GetOrder() > maxpsforder)
      maxpsforder = (*psf)[k].GetOrder();
    sigma_obs[k] = sqrt(pow(_sigma,2)+f_psf*pow((*psf)[k].GetSigma(),2));
    C[k].SaveDiv();
    
    I.SubVector(n,nx) *= pow(_sigma,2) / pow(sigma_obs[k],2);
  }

  int j=0;
  if (!fixcen) { U(j++,0) = 1.; U(j++,1) = 1.; }
  if (!fixgam) { U(j++,2) = 1.; U(j++,3) = 1.; }
  if (!fixmu) { U(j,4) = 1.; }

  SetupGx(_Gx,_order+2,_order);
  SetupGy(_Gy,_order+2,_order);
  SetupGg1(_Gg1,ORDER_1+2,ORDER_1);
  SetupGg2(_Gg2,ORDER_1+2,ORDER_1);
  SetupGth(_Gth,ORDER_1+2,ORDER_1);
  SetupGmu(_Gmu,ORDER_2,ORDER_3);

  double val = pixscale*pixscale/(_sigma*_sigma)/pix.size();
  for(int n=0,k=0;n<=b.GetOrder();++n) {
    for(int p=n,q=0;p>=q;--p,++q,++k) {
      if (p!=q) { normD(k) = val/2.; normD(++k) = val/2.; }
      else normD(k) = val;
    }
  }
}
#undef ORDER_1
#undef ORDER_2
#undef ORDER_3
#undef SIZE_1
#undef SIZE_1b
#undef SIZE_2
#undef SIZE_3

void ESImpl2::F(const tmv::Vector<double>& x, tmv::Vector<double>& f) const
{
  Assert(x.size() == U.colsize());
  if (useflux) Assert(f.size() == U.colsize()+1);
  else Assert(f.size() == U.colsize());
  // x(0),x(1) = u_cen, v_cen
  // x(2),x(3) = gamma_1, gamma_2
  // x(4) = mu
  //
  // ( u )' = exp(-mu)/sqrt(1-gsq) ( 1-g1  -g2  ) ( u-uc )
  // ( v )                         (  -g2  1+g1 ) ( v-vc )
  // 

  xx = x*U;
  if (fixcen) { xx[0] = fixuc;  xx[1] = fixvc; }
  if (fixgam) { xx[2] = fixg1;  xx[3] = fixg2; }
  if (fixmu) { xx[4] = fixm; }

  if (NormInf(xx-xinit) > 4.) 
  {
    if (flux == 0.) flux = b(0);
    if (useflux) {
      f(0) = 2.*(b(0)/flux-1.);
      f.SubVector(1,f.size()) = 2.*U*b.SubVector(1,6)/flux; 
    }
    else f = 2.*U*b.SubVector(1,6)/flux;
    return; 
  }

  std::complex<double> zc(xx[0],xx[1]);
  std::complex<double> g(xx[2],xx[3]);
  double mu = xx[4];
  double gsq = std::norm(g);

  if (gsq > 0.99 || (mu < -2. && Norm(xx-xinit) > 0.3))
  {
    if (flux == 0.) flux = b(0);
    if (flux == 0.) flux = 1.;
    if (useflux) {
      f(0) = 2.*(b(0)/flux-1.);
      f.SubVector(1,f.size()) = 2.*U*b.SubVector(1,6)/flux; 
    }
    else f = 2.*U*b.SubVector(1,6)/flux; 
    return; 
  }


  double m0 = exp(-mu)/sqrt(1.-std::norm(g));
  b.Zero();

  // z' = m*(z-zc) - m*g*conj(z-zc)
  //    = m*z - m*g*conj(z) - m*zc + m*g*conj(zc);
  Z1 = m0*Z;
  Z1 -= m0*g*Z.Conjugate();
  Z1.AddToAll(m0*(g*std::conj(zc)-zc));
  for(size_t k=0,n=0,nx;k<pix.size();k++,n=nx) {
    nx = n+pix[k].size();

    Z1.SubVector(n,nx) /= sigma_obs[k];

    MakePsi(A[k],Z1.SubVector(n,nx),b.GetOrder());
    b1 = I.SubVector(n,nx) * A[k];
    b1 *= normD;

    if (psf) {
      int psfsize = (*psf)[k].size();

      BVec newpsf((*psf)[k].GetOrder(),(*psf)[k].GetSigma());
      tmv::MatrixView<double> S1 = S[k].Cols(0,psfsize);
      GTransform(g,(*psf)[k].GetOrder(),S1);
      newpsf = S1 * (*psf)[k];
      tmv::MatrixView<double> D1 = D[k].Rows(0,psfsize);
      MuTransform(mu,(*psf)[k].GetOrder(),D1);
      newpsf = D1 * newpsf;
      PSFConvolve(newpsf,b.GetOrder(),b.GetSigma(),C[k]);
      C[k].ReSetDiv();
      b1 /= C[k];
    }
    b += b1;
  }

  // First time through, set flux
  if (flux == 0.) {
    flux = b(0);
  }
  if (useflux) {
    f(0) = b(0)/flux-1.;
    f.SubVector(1,U.colsize()+1) = U*b.SubVector(1,6)/flux;
  } else {
    f = U*b.SubVector(1,6)/flux;
  }
  b *= exp(-2.*mu);
}

void ESImpl2::J(const tmv::Vector<double>& x, const tmv::Vector<double>& f,
    tmv::Matrix<double>& df) const
{
  // b = exp(-2 mu) C^-1 normD AT I
  // db/dE = exp(-2 mu) C^-1 normD (dA/dE)T I
  //         - exp(-2 mu) C^-1 (dC/dE) C^-1 normD AT I
  //         [and if E == mu] -2 b
  Assert(x.size() == U.colsize());
  Assert(df.rowsize() == U.colsize());
  if (useflux) {
    Assert(f.size() == U.colsize()+1);
    Assert(df.colsize() == U.colsize()+1);
  } else {
    Assert(f.size() == U.colsize());
    Assert(df.colsize() == U.colsize());
  }

  xx = x*U;
  if (fixcen) { xx[0] = fixuc;  xx[1] = fixvc; }
  if (fixgam) { xx[2] = fixg1;  xx[3] = fixg2; }
  if (fixmu) { xx[4] = fixm; }
  std::complex<double> zc(xx[0],xx[1]);
  std::complex<double> g(xx[2],xx[3]);
  double mu = xx[4];
  double fact = 1./(1.-std::norm(g));

  double m0 = exp(-mu)*sqrt(fact);
  double g1 = std::real(g);
  double g2 = std::imag(g);

  tmv::ConstMatrixView<double> Gx = _Gx.SubMatrix(0,np2size,0,nsize);
  tmv::ConstMatrixView<double> Gy = _Gy.SubMatrix(0,np2size,0,nsize);
  tmv::ConstMatrixView<double> Gg1 = _Gg1.SubMatrix(0,np2size,0,nsize);
  tmv::ConstMatrixView<double> Gg2 = _Gg2.SubMatrix(0,np2size,0,nsize);
  tmv::ConstMatrixView<double> Gmu = _Gmu.SubMatrix(0,np2size,0,nsize);
  tmv::ConstMatrixView<double> Gth = _Gth.SubMatrix(0,np2size,0,nsize);

  dbdE.Zero();
  for(size_t k=0,n=0,nx;k<pix.size();k++,n=nx) {
    nx = n+pix[k].size();
    double m = m0/sigma_obs[k];

    AugmentPsi(A_aux[k],Z1.SubVector(n,nx),b.GetOrder());

    //tmv::Vector<double> IA = exp(-2.*mu) * I.SubVector(n,nx) * A_aux;
    IA = I.SubVector(n,nx) * A_aux[k];

    dbdE1.col(0) = -IA * m*(1.-g1)*Gx;
    dbdE1.col(0) += IA * m*g2*Gy;

    dbdE1.col(1) = -IA * m*(1.+g1)*Gy;
    dbdE1.col(1) += IA * m*g2*Gx;

    dbdE1.col(2) = -IA * fact*Gg1;
    dbdE1.col(2) += IA * fact*g2*Gth;

    dbdE1.col(3) = -IA * fact*Gg2;
    dbdE1.col(3) -= IA * fact*g1*Gth;

    dbdE1.col(4) = -IA * Gmu;

    dbdE1 = normD * dbdE1;
    b1 = normD * IA.SubVector(0,nsize);

    if (psf) {
      // dbdE -= dCdE b
      int n2 = ((*psf)[k].GetOrder()+3)*((*psf)[k].GetOrder()+4)/2;
      size_t psfsize = (*psf)[k].size();

      Assert(C[k].DivIsSet());
      b1 /= C[k];

      tmv::MatrixView<double> dTdE = _dTdE.SubMatrix(0,psfsize,0,psfsize);
      tmv::ConstMatrixView<double> D1 = D[k].Rows(0,psfsize);
      tmv::ConstMatrixView<double> S1 = S[k].Cols(0,psfsize);

      AugmentGTransformCols(g,(*psf)[k].GetOrder(),S[k].View());
      AugmentMuTransformRows(mu,(*psf)[k].GetOrder(),D[k].View());

      // E = g1
      dTdE = fact*D1*S[k]*_Gg1.SubMatrix(0,n2,0,psfsize);
      dTdE += fact*g2*D1*S[k]*_Gth.SubMatrix(0,n2,0,psfsize);
      BVec dbpsfdE((*psf)[k].GetOrder(),(*psf)[k].GetSigma());
      dbpsfdE = dTdE*(*psf)[k];
      PSFConvolve(dbpsfdE,b.GetOrder(),b.GetSigma(),dCdE);
      dbdE1.col(2) -= dCdE*b1;

      // E = g2
      dTdE = fact*D1*S[k]*_Gg2.SubMatrix(0,n2,0,psfsize);
      dTdE -= fact*g1*D1*S[k]*_Gth.SubMatrix(0,n2,0,psfsize);
      dbpsfdE = dTdE*(*psf)[k];
      PSFConvolve(dbpsfdE,b.GetOrder(),b.GetSigma(),dCdE);
      dbdE1.col(3) -= dCdE*b1;

      // E = mu
      dTdE = _Gmu.SubMatrix(0,psfsize,0,n2)*D[k]*S1;
      dbpsfdE = dTdE*(*psf)[k];
      PSFConvolve(dbpsfdE,b.GetOrder(),b.GetSigma(),dCdE);
      dbdE1.col(4) -= dCdE*b1;

      Assert(C[k].DivIsSet());
      dbdE1 /= C[k];
      //dbg<<"C["<<k<<"] = "<<C[k]<<std::endl;
      //tmv::Matrix<double> C2 = C[k];
      //tmv::DiagMatrix<double> S(C2.rowsize());
      //SV_Decompose(C2.View(),S.View(),false);
      //dbg<<"S = "<<S.diag()<<std::endl;
    }
    //dbdE1.col(4) -= 2.*b1;

    dbdE += dbdE1;
  }

  if (useflux) {
    df.row(0) = dbdE.row(0)*U.Transpose();
    df.Rows(1,U.colsize()+1) = U*dbdE.Rows(1,6)*U.Transpose();
  } else {
    df = U*dbdE.Rows(1,6)*U.Transpose();
  }
  //df *= exp(2.*mu)/flux;
  df /= flux;
  //df.col(4) += 2.*f;
}

void EllipseSolver2::CallF(const tmv::Vector<double>& x,
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

bool EllipseSolver2::Solve(tmv::Vector<double>& x, tmv::Vector<double>& f) const
{
  Assert(x.size() == 5);
  Assert(f.size() == 5);
  pimpl->xinit = x;
  if (pimpl->fixcen) { pimpl->fixuc = x[0]; pimpl->fixvc = x[1]; }
  if (pimpl->fixgam) { pimpl->fixg1 = x[2]; pimpl->fixg2 = x[3]; }
  if (pimpl->fixmu) { pimpl->fixm = x[4]; }

  pimpl->x_short = pimpl->U * x;

  pimpl->flux = 0;
  bool ret = NLSolver::Solve(pimpl->x_short,pimpl->f_short);

  x = pimpl->x_short*pimpl->U;
  if (pimpl->fixcen) { x[0] = pimpl->fixuc; x[1] = pimpl->fixvc; }
  if (pimpl->fixgam) { x[2] = pimpl->fixg1; x[3] = pimpl->fixg2; }
  if (pimpl->fixmu) { x[4] = pimpl->fixm; }
  if (pimpl->useflux) 
    f = pimpl->f_short.SubVector(1,pimpl->f_short.size()) * pimpl->U;
  else f = pimpl->f_short*pimpl->U;

  return ret;
}

bool EllipseSolver2::TestJ(const tmv::Vector<double>& x, tmv::Vector<double>& f,
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

void EllipseSolver2::UseNumericJ() { pimpl->numeric_j = true; }

const BVec& EllipseSolver2::GetB() const { return pimpl->b; }

void EllipseSolver2::GetBCov(tmv::Matrix<double>&) const 
{ Assert(false); }
