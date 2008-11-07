#include "Function2D.h"
#include <algorithm>
#include <complex>
using std::complex;
#include "TMV.h"
#include "TMV_Tri.h"
using tmv::Matrix;
using tmv::Vector;
#include <iostream>
using std::ostream;
using std::istream;
using std::endl;
#include <vector>
using std::vector;
#include "dbg.h"

using std::max;
using std::min;

using std::norm;
inline double norm(const double& x) { return x*x; }

template <class T>
Constant2D<T>::Constant2D(istream& fin) : Function2D<T>()
{
  if(!(fin >> (*coeffs)(0,0))) f2d_error("reading constant");
}

template <class T>
void Constant2D<T>::Write(ostream& fout) const
{
  int oldprec = fout.precision(6);
  std::ios_base::fmtflags oldf = fout.setf(std::ios_base::scientific,std::ios_base::floatfield);
  fout << "C " << (*coeffs)(0,0) << endl;
  if (!fout) f2d_error("writing (constant) function");
  fout.precision(oldprec);
  fout.flags(oldf);
}

template <class T>
void Constant2D<T>::operator+=(const Function2D<T>& rhs)
{
  const Constant2D<T>* crhs = dynamic_cast<const Constant2D<T>*>(&rhs);
  Assert(crhs);
  (*coeffs)(0,0) += (*crhs->coeffs)(0,0);
}


template <class T>
void Polynomial2D<T>::SetFunction(int xo, int yo, const Vector<T>& fvect)
{
  if(xorder != xo || yorder != yo) {
    xorder = xo; yorder = yo;
    coeffs.reset(new Matrix<T>(xo+1,yo+1,0.));
  }
  int k=0;
  //for(int m=0;m<=max(xorder,yorder);m++) 
  //  for(int i=min(m,xorder);m-i<=min(m,yorder);i--)
  //    (*coeffs)(i,m-i) = fvect(k++);
  
  int maxorder = max(xo,yo);
  for(int m=0;m<=maxorder;m++) {
    int i0 = min(xo,m);
    int len = min(yo,i0)+1;
    tmv::VectorView<T> mdiag = coeffs->SubVector(i0,m-i0,-1,1,len);
    tmv::ConstVectorView<T> subf = fvect.SubVector(k,k+len);
    mdiag = subf;
    k += len;
  }

  Assert(k==(int)fvect.size());
}

template <class T>
Polynomial2D<T>::Polynomial2D(istream& fin) : Function2D<T>()
// Order of parameters:  (example is for xorder = 2, yorder = 3
// xorder(2) yorder(3) a00 a10 a01 a20 a11 a02 a21 a12 a03
// where f = a00 + a10 x + a01 y + a20 x^2 + a11 xy + a02 y^2
//           + a21 x^2 y + a12 xy^2 + a03 y^3
// Note that total order doesn't go past the max of xorder and yorder.
// Also note that a30 is not listed since xorder is only 2.
// Note that aij are complex numbers so each is listed as real_part imag_part.
{
  int xo,yo;
  if (!(fin >> xo >> yo >> scale)) f2d_error("reading xorder,yorder,scale");
  xorder = xo;
  yorder = yo;
  coeffs.reset(new Matrix<T>(xo+1,yo+1,0.));
  int maxorder = max(xo,yo);
  for(int m=0;m<=maxorder;m++) {
    int i0 = min(xo,m);
    int len = min(yo,i0)+1;
    tmv::VectorView<T> mdiag = coeffs->SubVector(i0,m-i0,-1,1,len);
    for(int i=0;i<len;i++) fin >> mdiag(i);
  }
  if (!fin) f2d_error("reading (polynomial)");
}

template <class T>
void Polynomial2D<T>::Write(ostream& fout) const
{
  int oldprec = fout.precision(6);
  std::ios_base::fmtflags oldf = fout.setf(std::ios_base::scientific,std::ios_base::floatfield);
  int maxorder = max(xorder,yorder);
  if (maxorder == 0) fout << "C " << (*coeffs)(0,0) << endl;
  else {
    fout << "P " << xorder << ' ' << yorder << ' ' << scale << ' ';
    for(int m=0;m<=maxorder;m++) {
      int i0 = min(xorder,m);
      int len = min(yorder,i0)+1;
      tmv::VectorView<T> mdiag = coeffs->SubVector(i0,m-i0,-1,1,len);
      for(int i=0;i<len;i++) fout << mdiag(i) << ' ';
    }
  }
  fout << endl;
  if (!fout) f2d_error("writing (polynomial) function");
  fout.flags(oldf);
  fout.precision(oldprec);
}

template <class T>
void Polynomial2D<T>::AddLinear(T a, T b, T c)
{
  (*coeffs)(0,0) += a;
  (*coeffs)(1,0) += b*scale;
  (*coeffs)(0,1) += c*scale;
}

template <class T>
void Polynomial2D<T>::LinearPreTransform(T a, T b, T c, T d, T e, T f)
  // F(x,y) = Sum_i,j a(i,j) x^i y^j
  // F'(x,y) = F(a+bx+cy,d+ex+fy)
  //         = Sum_i,j a(i,j) (a+bx+cy)^i (d+ex+fy)^j
  //         = Sum_i,j a(i,j) (Sum_kl iCk kCl a^i-k (bx)^k-l (cy)^l) *
  //                Sum_mn jCm mCn d^j-m (ex)^m-n (fy)^n
{
  int maxorder = max(xorder,yorder);
  vector<double> scaletothe(maxorder+1);
  scaletothe[0] = 1.; scaletothe[1] = scale;
  for(int i=2;i<=maxorder;i++) scaletothe[i] = scaletothe[i-1]*scale;
  vector<T> atothe(maxorder+1);
  vector<T> btothe(maxorder+1);
  vector<T> ctothe(maxorder+1);
  vector<T> dtothe(maxorder+1);
  vector<T> etothe(maxorder+1);
  vector<T> ftothe(maxorder+1);
  atothe[0] = 1.; atothe[1] = a;
  btothe[0] = 1.; btothe[1] = b;
  ctothe[0] = 1.; ctothe[1] = c;
  dtothe[0] = 1.; dtothe[1] = d;
  etothe[0] = 1.; etothe[1] = e;
  ftothe[0] = 1.; ftothe[1] = f;
  for(int i=2;i<=maxorder;i++) atothe[i] = atothe[i-1]*a;
  for(int i=2;i<=maxorder;i++) btothe[i] = btothe[i-1]*b;
  for(int i=2;i<=maxorder;i++) ctothe[i] = ctothe[i-1]*c;
  for(int i=2;i<=maxorder;i++) dtothe[i] = dtothe[i-1]*d;
  for(int i=2;i<=maxorder;i++) etothe[i] = etothe[i-1]*e;
  for(int i=2;i<=maxorder;i++) ftothe[i] = ftothe[i-1]*f;
  xorder = maxorder; yorder = maxorder;

  tmv::LowerTriMatrix<T,tmv::UnitDiag> binom(maxorder+1);
  for(int n=1;n<=maxorder;++n) {
    binom(n,0) = T(1);
    for(int m=1;m<n;++m) {
      binom(n,m) = binom(n-1,m-1) + binom(n-1,m);
    }
  }

  std::auto_ptr<Matrix<T> > oldcoeffs = coeffs;
  coeffs.reset(new Matrix<T>(xorder+1,yorder+1,0.));
  for(int i=0;i<=xorder;i++) for(int j=0;j<=yorder&&i+j<=maxorder;j++) {
    for(int k=0;k<=i;k++) for(int l=0;l<=k;l++) {
      for(int m=0;m<=j;m++) for(int n=0;n<=m;n++) {
        (*coeffs)(k-l+m-n,l+n) += 
          binom(i,k)*binom(k,l)*binom(j,m)*binom(m,n)*
          atothe[i-k]*btothe[k-l]*ctothe[l]*
          dtothe[j-m]*etothe[m-n]*ftothe[n]*
          (*oldcoeffs)(i,j)/scaletothe[i+j-k-m];
      }
    }
  }
}

template <class T>
void Polynomial2D<T>::operator+=(const Function2D<T>& rhs)
{
  const Polynomial2D<T>* prhs = dynamic_cast<const Polynomial2D<T>*>(&rhs);
  Assert(prhs);
  Assert(scale == prhs->scale);
  if (xorder == prhs->xorder && yorder == prhs->yorder) {
    *coeffs += *prhs->coeffs;
  } else {
    int newxorder = max(xorder,prhs->xorder);
    int newyorder = max(yorder,prhs->yorder);
    std::auto_ptr<Matrix<T> > newc(new Matrix<T>(newxorder+1,newyorder+1,0.));
    newc->SubMatrix(0,xorder+1,0,yorder+1) = *coeffs;
    newc->SubMatrix(0,prhs->xorder+1,0,prhs->yorder+1) += *prhs->coeffs;
    coeffs = newc;
    xorder = newxorder;
    yorder = newyorder;
  }
}

template <class T>
std::auto_ptr<Function2D<T> > Polynomial2D<T>::DFDX() const
{
  if (xorder == 0) 
    return std::auto_ptr<Function2D<T> >(new Constant2D<T>());
  if (xorder == 1 && yorder == 0) 
    return std::auto_ptr<Function2D<T> >(new Constant2D<T>((*coeffs)(1,0)));

  int newxorder = xorder-1;
  int newyorder = xorder > yorder ? yorder : yorder-1;

  std::auto_ptr<Polynomial2D<T> > temp(
    new Polynomial2D<T>(newxorder,newyorder));

  int maxorder = max(newxorder,newyorder);
  for(int i=newxorder;i>=0;i--) 
    for(int j=min(maxorder-i,newyorder);j>=0;j--) {
      Assert(i+1<=xorder);
      Assert(j<=yorder);
      Assert(i+1+j<=max(xorder,yorder));
      (*temp->coeffs)(i,j) = (*coeffs)(i+1,j)*(i+1.)/scale;
  }
  return std::auto_ptr<Function2D<T> >(temp);
}

template <class T>
std::auto_ptr<Function2D<T> > Polynomial2D<T>::DFDY() const 
{
  if (yorder == 0)
    return std::auto_ptr<Function2D<T> >(new Constant2D<T>());
  if (yorder == 1 && xorder == 0)
    return std::auto_ptr<Function2D<T> >(new Constant2D<T>((*coeffs)(0,1)));

  int newxorder = yorder > xorder ? xorder : xorder-1;
  int newyorder = yorder-1;

  std::auto_ptr<Polynomial2D<T> > temp(
    new Polynomial2D<T>(newxorder,newyorder));

  int maxorder = max(newxorder,newyorder);
  for(int i=newxorder;i>=0;i--) 
    for(int j=min(maxorder-i,newyorder);j>=0;j--) {
      Assert(i<=xorder);
      Assert(j+1<=yorder);
      Assert(i+j+1<=max(xorder,yorder));
      (*temp->coeffs)(i,j) = (*coeffs)(i,j+1)*(j+1.)/scale;
  }
  return std::auto_ptr<Function2D<T> >(temp);
}

template <> std::auto_ptr<Function2D<double> > Function2D<double>::Conj() const
{
  return Copy();
}

template <> std::auto_ptr<Function2D<complex<double> > > Function2D<complex<double> >::Conj() const
{
  std::auto_ptr<Function2D<complex<double> > > temp = Copy();
  temp->coeffs->ConjugateSelf();
  return temp;
}

template <class T>
T Function2D<T>::operator()(double x,double y) const
{
  //xdbg<<"Start op()("<<x<<","<<y<<")\n";
  Vector<double> px = DefinePX(xorder,x);
  //xdbg<<"px = "<<px<<endl;
  Vector<double> py = DefinePY(yorder,y);
  //xdbg<<"py = "<<py<<endl;
  //xdbg<<"coeffs = "<<*coeffs<<endl;
  T result = px * (*coeffs) * py;
  //xdbg<<"result = "<<result<<endl;
  return result;
  //T result = 0.;
  //int maxorder = max(xorder,yorder);
  //for(int i=0;i<=xorder;i++) for(int j=0;j<=min(maxorder-i,yorder);j++) {
    //result += (*coeffs)(i,j)*px(i)*py(j);
  //}
  //return result;
}

template <class T>
std::auto_ptr<Function2D<T> > Function2D<T>::Read(istream& fin) 
{
  char fc,tc;

  fin >> fc >> tc;
  if (tc != 'D' && tc != 'C') fin.putback(tc);
  switch(fc) {
    case 'C' : return std::auto_ptr<Function2D<T> >(new Constant2D<T>(fin));
    case 'P' : return std::auto_ptr<Function2D<T> >(new Polynomial2D<T>(fin));
#ifdef LEGENDRE2D_H
    case 'L' : return std::auto_ptr<Function2D<T> >(new Legendre2D<T>(fin));
#endif
#ifdef CHEBY2D_H
    case 'X' : return std::auto_ptr<Function2D<T> >(new Cheby2D<T>(fin));
#endif
    default: f2d_error("invalid type"); return std::auto_ptr<Function2D<T> >(0);
  }
}

template <class T>
void Function2D<T>::LinearTransform(T a, T b, T c,
	const Function2D<T>& f, const Function2D<T>& g)
{
#ifndef NDEBUG
  if(dynamic_cast<Constant2D<T>*>(this)) {
    Assert(dynamic_cast<const Constant2D<T>*>(&f));
    Assert(dynamic_cast<const Constant2D<T>*>(&g));
  }
  if(dynamic_cast<Polynomial2D<T>*>(this)) {
    Assert(dynamic_cast<const Polynomial2D<T>*>(&f));
    Assert(dynamic_cast<const Polynomial2D<T>*>(&g));
  }
#ifdef LEGENDRE2D_H
  if(dynamic_cast<Legendre2D<T>*>(this)) {
    const Legendre2D<T>* lf = dynamic_cast<const Legendre2D<T>*>(&f);
    const Legendre2D<T>* lg = dynamic_cast<const Legendre2D<T>*>(&g);
    Assert(lf);
    Assert(lg);
    Assert(lf->GetBounds() == lg->GetBounds());
  }
#endif
#ifdef CHEBY2D_H
  if(dynamic_cast<Cheby2D<T>*>(this)) {
    const Cheby2D<T>* cf = dynamic_cast<const Cheby2D<T>*>(&f);
    const Cheby2D<T>* cg = dynamic_cast<const Cheby2D<T>*>(&g);
    Assert(cf);
    Assert(cg);
    Assert(cf->GetBounds() == cg->GetBounds());
  }
#endif
  Assert(f.GetXOrder() == g.GetXOrder());
  Assert(f.GetYOrder() == g.GetYOrder());
#endif

  if (xorder != f.GetXOrder() || yorder != f.GetYOrder()) {
    xorder = f.GetXOrder();
    yorder = f.GetYOrder();
    coeffs.reset(new Matrix<T>(xorder+1,yorder+1,0.));
  } else coeffs->Zero();
  for(int i=0;i<=xorder;i++) for(int j=0;j<=yorder;j++) {
    (*coeffs)(i,j) = a + b*f.GetCoeffs()(i,j) + c*g.GetCoeffs()(i,j);
  }
}

inline int FitSize(const int xorder, const int yorder)
{
  int loworder = min(xorder,yorder);
  int highorder = max(xorder,yorder);
  return (loworder+1)*(loworder+2)/2 + (loworder+1)*(highorder-loworder);
}

template <class T>
void Function2D<T>::DoSimpleFit(size_t xorder, size_t yorder, 
    const vector<Position>& pos, const vector<T>& v, 
    const vector<bool>& use, Vector<T> *f, const vector<double>* siglist,
    int *dof, Vector<T> *diff, Matrix<double>* cov)
// f(x,y) = Sum_pq k_pq px_p py_q
//        = P . K (where each is a vector in pq)
// chisq = Sum_n ((v_n - f(x,y))/s_n)^2
// minchisq => 0 = Sum_n (v_n - f(x,y)) P_pq/s_n^2
// => [Sum_n P_pq P/s_n^2].K = [Sum_n v_n P_pq/s_n^2]
// Using ^ to represent an outer product, this can be written:
//
//    [Sum_n (P/s_n)^(P/s_n)] . K = [Sum_n (v_n/s_n) (P/s_n)]
//
// Or if P' is a matrix with n index for rows, and pq index for columns,
// with each element P'(n,pq) = P_n,pq/s_n
// and v' is a vector with v'(n) = v_n/s_n
//
// Then, this can be written:
//
// P' K = v'
//
// The solution to this equation, then, gives our answer for K.
//
{
  Assert(pos.size() == v.size());
  Assert(use.size() == v.size());
  Assert((siglist==0) || (siglist->size() == v.size()));
  xdbg<<"Start SimpleFit: size = "<<pos.size()<<endl;
  xdbg<<"order = "<<xorder<<','<<yorder<<endl;

  size_t highorder = max(xorder,yorder);
  size_t size = FitSize(xorder,yorder);
  xdbg<<"size = "<<size<<endl;

  Assert(f->size() == size);
  Assert(!diff || diff->size() == v.size());
  Assert(!cov || (cov->colsize() == size && cov->rowsize() == size));

  size_t nuse = 0;
  for(size_t i=0;i<use.size();i++) if (use[i]) nuse++;
  xdbg<<"nuse = "<<nuse<<endl;

  Matrix<double> P(nuse,size,0.);
  Vector<T> V(nuse);

  size_t ii=0;
  for(size_t i=0;i<v.size();i++) if (use[i]) {
    if (siglist) {
      Assert((*siglist)[i] > 0.);
      V(ii) = v[i]/(*siglist)[i];
    } else {
      V(ii) = v[i];
    }

    Vector<double> px = DefinePX(xorder,pos[i].GetX());
    Vector<double> py = DefinePY(yorder,pos[i].GetY());
    size_t pq=0;
    for(size_t pplusq=0;pplusq<=highorder;pplusq++) { 
      for(size_t p=min(pplusq,xorder),q=pplusq-p;
          q<=std::min(pplusq,yorder);p--,q++,pq++) {
        Assert(p<px.size());
        Assert(q<py.size());
        Assert(pq<P.rowsize());
        P(ii,pq) = px[p]*py[q];
      }
    }
    Assert(pq == P.rowsize());
    if (siglist) P.row(ii) /= (*siglist)[i];
    ii++;
  }
  Assert(ii==nuse);

  xdbg<<"Done make V,P\n";
  xdbg<<"V = "<<V<<endl;
  P.DivideUsing(tmv::QR);
  P.SaveDiv();
  *f = V/P;
  xdbg<<"*f = "<<*f<<endl;

  if (diff) {
    size_t k=0;
    for(size_t i=0;i<v.size();i++) {
      if (use[i]) {
	(*diff)(i) = V(k) - P.row(k) * (*f);
	k++;
      }
      else (*diff)(i) = 0.;
    }
  }
  if (dof) {
    *dof = P.colsize() - P.rowsize();
    if (*dof < 0) *dof = 0;
  }
  if (cov) {
    P.InverseATA(*cov);
  }
  xdbg<<"Done simple fit\n";
}
    
template <class T>
void Function2D<T>::SimpleFit(int order, const vector<Position>& pos, 
    const vector<T>& v, const vector<bool>& use, 
    const vector<double>* sig,
    double *chisqout, int *dofout, Matrix<double>* cov) 
{
  Vector<T> fvect(FitSize(order,order));
  if (chisqout) {
    Vector<T> diff(v.size());
    DoSimpleFit(order,order,pos,v,use,&fvect,sig,dofout,&diff,cov);
    *chisqout = NormSq(diff);
  }
  else DoSimpleFit(order,order,pos,v,use,&fvect,sig);
  SetFunction(order,order,fvect);
}

template <class T>
void Function2D<T>::OutlierFit(int order,double nsig,
    const vector<Position>& pos, const vector<T>& v, vector<bool> *use,
    const vector<double>* sig, double *chisqout, int *dofout, 
    Matrix<double> *cov)
{
  xdbg<<"start outlier fit\n";
  bool done=false;
  Vector<T> fvect(FitSize(order,order));
  int dof;
  double chisq=0.;
  double nsigsq = nsig*nsig;
  while (!done) {
    Vector<T> diff(v.size());
    xdbg<<"before dosimple\n";
    DoSimpleFit(order,order,pos,v,*use,&fvect,sig,&dof,&diff,cov);
    xdbg<<"after dosimple\n";
    // Caclulate chisq, keeping the vector diffsq for later when
    // looking for outliers
    chisq = NormSq(diff);
    // If sigmas are given, then chisq should be 1, since diff is
    // already normalized by sig_i.  But if not, then this effectively
    // assumes that all the given errors are off by some uniform factor.

    // Look for outliers, setting done = false if any are found
    done = true;
    if (dof <= 0) break;
    double thresh=nsigsq*chisq/dof;
    xdbg<<"chisq = "<<chisq<<", thresh = "<<thresh<<endl;
    for(size_t i=0;i<v.size();i++) if( (*use)[i]) { 
      xdbg<<"pos ="<<pos[i]<<", v = "<<v[i]<<" diff = "<<diff[i]<<endl;
      if (norm(diff[i]) > thresh) {
        done = false; 
        (*use)[i] = false;
        xdbg<<i<<" ";
      }
    }
    if (!done) xdbg<<" are outliers\n";
  }
  SetFunction(order,order,fvect);
  if (chisqout) *chisqout = chisq;
  if (dofout) *dofout = dof;
  xdbg<<"done outlier fit\n";
}

inline double betai(double a,double b,double x);

inline bool Equivalent(double chisq1,double chisq2, int n1, int n2,
  double equivprob)
{
  if (chisq2 <= chisq1) return true;
  // should only happpen when essentially equal but have rounding errors
  if (chisq1 <= 0.) return (chisq2 <= 0.);

  Assert(n1 < n2);
  if (n1 <= 0) return (n2 <= 0);
  double f = (chisq2-chisq1)/(n2-n1) / (chisq1/n1);
  
  double prob = betai(0.5*n2,0.5*n1,n2/(n2+n1*f));
  // = probability that these chisq would happen by chance if equiv
  // (technically if underlying chisq2 really smaller or = to chisq1)
  // so equiv if prob is large, not equiv if prob is small.
  return (prob > 1.-equivprob);
}

template <class T>
void Function2D<T>::OrderFit(int maxorder, double equivprob,
    const vector<Position>& pos,const vector<T>& v,
    const vector<bool>& use, const vector<double>* sig,
    double *chisqout, int *dofout, Matrix<double>* cov) 
{
  xdbg<<"Start OrderFit\n";
  Vector<T> fvectmax(FitSize(maxorder,maxorder));
  Vector<T> diff(v.size());
  int dofmax;
  DoSimpleFit(maxorder,maxorder,pos,v,use,&fvectmax,sig,&dofmax,&diff,cov);
  double chisqmax = NormSq(diff);
  xdbg<<"chisq,dof(n="<<maxorder<<") = "<<chisqmax<<','<<dofmax<<endl;
  int tryorder;
  double chisq=-1.;
  int dof=-1;
  std::auto_ptr<Vector<T> > fvect(0);
  for(tryorder=0;tryorder<maxorder;tryorder++) {
    fvect.reset(new Vector<T>(FitSize(tryorder,tryorder)));
    DoSimpleFit(tryorder,tryorder,pos,v,use,fvect.get(),sig,&dof,&diff,cov);
    chisq = NormSq(diff);
    xdbg<<"chisq,dof(n="<<tryorder<<") = "<<chisq<<','<<dof<<"....  ";
    if (Equivalent(chisqmax,chisq,dofmax,dof,equivprob)) {
      xdbg<<"equiv\n";
      break;
    }
    xdbg<<"not equiv\n";
  }
  if (tryorder == maxorder) {
    SetFunction(tryorder,tryorder,fvectmax);
    if(chisqout) *chisqout = chisqmax;
    if(dofout) *dofout = dofmax;
  } else {
    Assert(fvect.get());
    SetFunction(tryorder,tryorder,*fvect);
    if(chisqout) *chisqout = chisq;
    if(dofout) *dofout = dof;
  }
}

#include <math.h>
#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

inline double betacf(double a,double b,double x)
{
  double qab = a+b;
  double qap = a+1.0;
  double qam = a-1.0;
  double c = 1.0;
  double d = 1.0 - qab*x/qap;
  if (fabs(d) < FPMIN) d = FPMIN;
  d = 1.0/d;
  double h = d;
  int m;
  for(m=1;m<=MAXIT;m++) {
    int m2 = 2*m;
    double aa = m*(b-m)*x/((qam+m2)*(a+m2));
    d = 1.0+aa*d;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = 1.0+aa/c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d = 1.0+aa*d;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = 1.0+aa/c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1.0/d;
    double del = d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (m>MAXIT) f2d_error("a or b too big in betacf, or MAXIT too small");
  return h;
}

#undef MAXIT
#undef EPS
#undef FPMIN

inline double gammln(double x)
{
  const double cof[6]={76.18009172947146,-86.50532032941677,
    24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
    -0.5395239384953e-5};
  double temp = x+5.5;
  temp -= (x+0.5)*log(temp);
  double ser = 1.000000000190015;
  double y=x;
  for(int j=0;j<6;j++) ser += cof[j]/(y+=1.0);
  return -temp+log(2.5066282746310005*ser/x);
}

inline double betai(double a,double b,double x)
{
  if (x<0.0 || x>1.0) f2d_error("Bad x in betai");
  if (x==0.0) return 0.;
  if (x==1.0) return 1.;
  double bt = exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
  if (x < (a+1.0)/(a+b+2.0)) return bt*betacf(a,b,x)/a;
  else return 1.0-bt*betacf(b,a,1.0-x)/b;
}

template class Function2D<complex<double> >;
template class Function2D<double>;
template class Constant2D<complex<double> >;
template class Constant2D<double>;
template class Polynomial2D<complex<double> >;
template class Polynomial2D<double>;