
#include "BVec.h"
#include "TMV.h"
#include <cmath>
#include "dbg.h"

const double PI = 3.14159265359;
const double sqrtpi = sqrt(PI);

BVec& BVec::operator=(const BVec& rhs)
{ 
  Assert(order >= rhs.order);
  //Assert(sigma == rhs.sigma);
  sigma = rhs.sigma;
  SubVector(0,rhs.size()) = rhs;
  SubVector(rhs.size(),size()).Zero();
  return *this;
}

BVec& BVec::operator=(const tmv::AssignableToVector<double>& rhs)
{ 
  Assert(int(rhs.size()) == (order+1)*(order+2)/2); 
  rhs.AssignToV(this->View());
  return *this;
}

void ZTransform(std::complex<double> z, int order, 
    const tmv::MatrixView<double>& T)
{
  if (z == 0.0) { T.SetToIdentity(); return; }

  Assert(int(T.colsize()) == (order+1)*(order+2)/2);
  Assert(int(T.rowsize()) == (order+1)*(order+2)/2);
  Assert(T.iscm());
  Assert(!T.isconj());

  // T(st,pq) = f(p,s) f*(q,t)
  // f(p,0) = (-z/2)^p / sqrt(p!) exp(-|z|^2/8)
  // f(p,s+1) = (sqrt(p) f(p-1,s) + 1/2 z* f(p,s) )/sqrt(s+1)

  std::complex<double> zo2 = -z/2.;
  tmv::Matrix<std::complex<double>,tmv::ColMajor> f(order+1,order+1);
  std::complex<double>* fcols = f.ptr();
  std::vector<double> isqrt(order+1);

  fcols[0] = exp(-std::norm(z)/8.); // f(0,0)
  isqrt[0] = 0.;
  for(int p=1;p<=order;p++) 
  {
    isqrt[p] = sqrt(double(p));
    fcols[p] = fcols[p-1]*(-zo2)/isqrt[p]; // f(p,0)
  }
  std::complex<double>* fcolsp1 = fcols+f.stepj();
  for(int s=0;s<order;++s) 
  {
    fcolsp1[0] = std::conj(zo2)*fcols[0]/isqrt[s+1]; // f(0,s+1)
    for(int p=1;p<=order;p++) 
    {
      fcolsp1[p] = (isqrt[p]*fcols[p-1] + std::conj(zo2)*fcols[p])/
	isqrt[s+1]; // f(p,s+1)
    }
    fcols = fcolsp1;
    fcolsp1 += f.stepj();
  }

  for(int n=0,pq=0;n<=order;++n) 
  {
    for(int p=n,q=0;p>=q;--p,++q,++pq) 
    {
      double* Tpq = T.col(pq).ptr();
      double* Tpq1 = Tpq + T.stepj();
      for(int nn=0,st=0;nn<=order;++nn) 
      {
	for(int s=nn,t=0;s>=t;--s,++t,++st) 
	{
	  std::complex<double> t1 = f.cref(p,s) * std::conj(f.cref(q,t));
	  if (p==q) 
	  {
	    if (s==t) Tpq[st] = std::real(t1); // T(st,pq)
	    else 
	    {
	      Tpq[st] = std::real(t1);
	      Tpq[++st] = std::imag(t1);
	    }
	  }
	  else if (s==t) 
	  {
	    // b_st = t_stpq b_pq + t_stqp b_qp
	    // In this case with s=t, t_stpq == t_stqp*
	    Tpq[st] = 2.*std::real(t1);
	    Tpq1[st] = -2.*std::imag(t1);
	  }
	  else 
	  {
	    std::complex<double> t2 = f.cref(q,s) * std::conj(f.cref(p,t));
	    Tpq[st] = std::real(t1) + std::real(t2);
	    Tpq1[st] = -std::imag(t1) + std::imag(t2);
	    Tpq[++st]= std::imag(t1) + std::imag(t2);
	    Tpq1[st] = std::real(t1) - std::real(t2);
	  }
	}
      }
      if (p!=q) ++pq;
    }
  }
}

void ApplyZ(std::complex<double> z, BVec& b)
{
  if (z != 0.0) 
  {
    tmv::Matrix<double> T(b.size(),b.size());
    ZTransform(z/b.GetSigma(),b.GetOrder(),T);
    b = T * b;
  }
}

void MuTransform(double mu, int order, 
    const tmv::MatrixView<double>& D)
{
  if (mu == 0.0) { D.SetToIdentity(); return; }

  Assert(int(D.colsize()) == (order+1)*(order+2)/2);
  Assert(int(D.rowsize()) == (order+1)*(order+2)/2);
  // The easiest recursion for the D_mu matrix is not actually
  // the one found in BJ02.  Rather, I use the following:
  // D(00,pq) = e^mu sech(mu) (-tanh(mu))^q delta_pq
  //   where pCq is the binomial coefficient.
  // D(s0,pq) = 1/sqrt(s) sech(mu) sqrt(p) D(s-10,p-1q)
  // D(st,pq) = 1/sqrt(t) [ sech(mu) sqrt(q) D(st-1,pq-1)
  //                        + tanh(mu) sqrt(s) D(s-1t-1,pq) ]
  // Note: Only s-t = p-q terms are nonzero.
  //       This means that this can be significantly sped up
  //       by forming smaller matrices for each m, and using 
  //       permutations to rotabe b so that these elements are
  //       continuous and then rotate back after multiplying.
  //       However, I'll wait to implement this until such speed 
  //       is found to be necessary.

  mu = -mu; 
  double tmu = tanh(mu);
  double smu = 1./cosh(mu);

  std::vector<double> isqrt(order+1);
  for(int i=0;i<=order;++i) isqrt[i] = sqrt(double(i));

  D.Zero();
  double Dqq = exp(mu)*smu;
  // Dqq = exp(mu) sech(mu) (-tanh(mu))^q
  // This variable keeps the latest Dqq value:
  for(int n=0,pq=0;n<=order;++n) 
  {
    for(int p=n,q=0;p>=q;--p,++q,++pq) 
    {
      double* Dpq = D.col(pq).ptr();
      double* Dpq_n = D.col(pq-n).ptr();
      double* Dpq_n_2 = q>0?D.col(pq-n-2).ptr():0;
      double* Dpq1 = p>q?D.col(pq+1).ptr():0;
      if (p==q) 
      {
	if (p > 0) Dqq *= -tmu;
	Dpq[0] = Dqq;  // D(0,pq)
      }
      for(int m=1,st=1;m<=order;++m) 
      {
	for(int s=m,t=0;s>=t;--s,++t,++st) 
	{
	  if (p-q==s-t) 
	  {
	    double temp;
	    if (t == 0) 
	    {
	      //temp = smu*sqrt(double(p))*D(st-m,pq-n)/sqrt(double(s));
	      temp = smu*isqrt[p]*Dpq_n[st-m]/isqrt[s];
	    }
	    else 
	    {
	      temp = tmu*isqrt[s]*Dpq[st-2*m-1];
	      if (q > 0) 
	      {
		temp += smu*isqrt[q]*Dpq_n_2[st-m-2];
	      }
	      temp /= isqrt[t];
	    }
	    Dpq[st] = temp; // D(st,pq)
	    if (s!=t) Dpq1[st+1] = temp; // D(st+1,pq+1)
	  }
	  if (s!=t) ++st;
	}
      }
      if (p!=q) ++pq;
    }
  }
}

void AugmentMuTransformRows(double mu, int order, 
    const tmv::MatrixView<double>& D)
{
  Assert(int(D.colsize()) == (order+3)*(order+4)/2);
  Assert(int(D.rowsize()) == (order+1)*(order+2)/2);

  const int size1 = D.rowsize();
  const int size2 = D.colsize();

  D.Rows(size1,size2).Zero();
  if (mu == 0.0) return;

  // The easiest recursion for the D_mu matrix is not actually
  // the one found in BJ02.  Rather, I use the following:
  // D(00,pq) = e^mu sech(mu) (-tanh(mu))^q delta_pq
  //   where pCq is the binomial coefficient.
  // D(s0,pq) = 1/sqrt(s) sech(mu) sqrt(p) D(s-10,p-1q)
  // D(st,pq) = 1/sqrt(t) [ sech(mu) sqrt(q) D(st-1,pq-1)
  //                        + tanh(mu) sqrt(s) D(s-1t-1,pq) ]
  // Note: Only s-t = p-q terms are nonzero.
  //       This means that this can be significantly sped up
  //       by forming smaller matrices for each m, and using 
  //       permutations to rotabe b so that these elements are
  //       continuous and then rotate back after multiplying.
  //       However, I'll wait to implement this until such speed 
  //       is found to be necessary.

  mu = -mu; 
  double tmu = tanh(mu);
  double smu = 1./cosh(mu);
  std::vector<double> isqrt(order+3);
  for(int i=0;i<=order+2;++i) isqrt[i] = sqrt(double(i));

  for(int n=0,pq=0;n<=order;++n) 
  {
    for(int p=n,q=0;p>=q;--p,++q,++pq) 
    {
      double* Dpq = D.col(pq).ptr();
      double* Dpq_n = D.col(pq-n).ptr();
      double* Dpq_n_2 = q>0?D.col(pq-n-2).ptr():0;
      double* Dpq1 = p>q?D.col(pq+1).ptr():0;
      for(int m=order+1,st=size1;m<=order+2;++m) 
      {
	for(int s=m,t=0;s>=t;--s,++t,++st) 
	{
	  if (p-q==s-t) 
	  {
	    double temp;
	    if (t == 0) 
	    {
	      temp = smu*isqrt[p]*Dpq_n[st-m]/isqrt[s];
	    }
	    else 
	    {
	      temp = tmu*isqrt[s]*Dpq[st-2*m-1];
	      if (q > 0) 
	      {
		temp += smu*isqrt[q]*Dpq_n_2[st-m-2];
	      }
	      temp /= isqrt[t];
	    }
	    Dpq[st] = temp;
	    if (s!=t) Dpq1[st+1] = temp;
	  }
	  if (s!=t) ++st;
	}
      }
      if (p!=q) ++pq;
    }
  }
}

void ApplyMu(double mu, BVec& b)
{
  if (mu != 0.0) 
  {
    tmv::Matrix<double> D(b.size(),b.size());
    MuTransform(mu,b.GetOrder(),D);
    b = D * b;
  }
}

void GTransform(std::complex<double> g, int order,
    const tmv::MatrixView<double>& S)
{
  if (g == 0.0) { S.SetToIdentity(); return; }

  Assert(int(S.colsize()) == (order+1)*(order+2)/2);
  Assert(int(S.rowsize()) == (order+1)*(order+2)/2);
  // S(st,pq) = f(p,s) f(q,t) (eta/|eta|)^(s-t-p+q)
  // f(p,0) = sqrt(p!)/(p/2)! sqrt(sech(|eta|/2)) (-tanh(|eta|/2)/2)^p/2
  // f(p,s+1) = (sqrt(p) sech(|eta|/2) f(p-1,s) +
  //                      sqrt(s) tanh(|eta|/2) f(p,s-1))/sqrt(s+1)
  //
  // tanh(|eta|/2) = |g|
  // sech(|eta|/2) = sqrt(1-|g|^2)
  //
  // Note: Like with the matrix in ApplyMu, this one is also fairly
  // sparse.  Could get a speedup by expoiting that, but currently don't.
  // I'll wait until the speedup is found to be necessary.

  double absg = std::abs(g);
  double normg = std::norm(g);
  std::vector<std::complex<double> > phase(2*order+1);
  phase[0] = 1.;
  std::complex<double> ph = -std::conj(g)/absg;
  // I'm not sure why conj was needed here.  Maybe there is an
  // error in the phase indices below that this corrects.
  for(int i=1;i<=2*order;i++) phase[i] = phase[i-1]*ph;

  tmv::Matrix<double,tmv::ColMajor> f(order+1,order+1,0.);
  double te = absg;
  double se = sqrt(1.-normg);

  double* fcols = f.ptr();
  fcols[0] = sqrt(se); // f(0,0)
  // only terms with p+s even are non-zero.
  for(int p=2;p<=order;p+=2) 
  {
    fcols[p] = fcols[p-2]*(-te/2.)*sqrt(double(p*(p-1)))/double(p/2); // f(p,0)
  }
  double* fcolsm1 = 0;
  double* fcolsp1 = fcols + f.stepj();
  std::vector<double> isqrt(order+1);
  for(int i=0;i<=order;i++) isqrt[i] = sqrt(double(i));

  for(int s=0;s<order;s++) 
  {
    if (s%2==1) 
    {
      fcolsp1[0] = isqrt[s]*te*fcolsm1[0]/isqrt[s+1]; // f(0,s+1)
    }
    for(int p=s%2+1;p<=order;p+=2) 
    {
      double temp = isqrt[p]*se*fcols[p-1]; 
      if (s>0) temp += isqrt[s]*te*fcolsm1[p];
      temp /= isqrt[s+1]; 
      fcolsp1[p] = temp; // f(p,s+1)
    }
    fcolsm1 = fcols;
    fcols = fcolsp1;
    fcolsp1 += f.stepj();
  }

  S.Zero();
  for(int n=0,pq=0;n<=order;++n) 
  {
    for(int p=n,q=0;p>=q;--p,++q,++pq) 
    {
      double* Spq = S.col(pq).ptr();
      double* Spq1 = p>q?S.col(pq+1).ptr():0;
      for(int nn=n%2,st=(nn==0?0:1);nn<=order;nn+=2,st+=nn) 
      {
	for(int s=nn,t=0;s>=t;--s,++t,++st) 
	{

	  double s0 = f.cref(p,s) * f.cref(q,t);

	  int iphase = s-t-p+q;
	  std::complex<double> s1 = s0 * 
	    (iphase >= 0 ? phase[iphase/2] : std::conj(phase[-iphase/2]));

	  if (p==q) 
	  {
	    if (s==t) Spq[st] = std::real(s1);
	    else 
	    {
	      Spq[st] = std::real(s1);
	      Spq[++st] = std::imag(s1);
	    }
	  }
	  else if (s==t) 
	  {
	    // b_st = t_stpq b_pq + t_stqp b_qp
	    // In this case with s=t, t_stpq == t_stqp*
	    Spq[st] = 2.*std::real(s1);
	    Spq1[st] = -2.*std::imag(s1);
	  }
	  else 
	  {
	    s0 = f.cref(q,s) * f.cref(p,t);
	    iphase = s-t-q+p;
	    std::complex<double> s2 = s0 *
	      (iphase >= 0 ? phase[iphase/2] : std::conj(phase[-iphase/2]));
	    Spq[st] = std::real(s1) + std::real(s2);
	    Spq1[st] = -std::imag(s1) + std::imag(s2);
	    Spq[++st] = std::imag(s1) + std::imag(s2);
	    Spq1[st] = std::real(s1) - std::real(s2);
	  }
	}
      }
      if (p!=q) ++pq;
    }
  }
}

void AugmentGTransformCols(std::complex<double> g, int order,
    const tmv::MatrixView<double>& S)
{
  Assert(int(S.colsize()) == (order+1)*(order+2)/2);
  Assert(int(S.rowsize()) == (order+3)*(order+4)/2);

  const int size1 = S.colsize();
  const int size2 = S.rowsize();

  S.Cols(size1,size2).Zero();
  if (g == 0.0) { return; }

  double absg = std::abs(g);
  double normg = std::norm(g);
  std::vector<std::complex<double> > phase(2*order+5);
  phase[0] = 1.;
  std::complex<double> ph = -std::conj(g)/absg;
  for(int i=1;i<=2*order+4;i++) phase[i] = phase[i-1]*ph;

  tmv::Matrix<double> f(order+3,order+1,0.);
  double te = absg;
  double se = sqrt(1.-normg);

  double* fcols = f.ptr();
  fcols[0] = sqrt(se); // f(0,0)
  for(int p=2;p<=order+2;p+=2) 
  {
    fcols[p] = fcols[p-2]*(-te/2.)*sqrt(double(p*(p-1)))/double(p/2); // f(p,0)
  }
  double* fcolsm1 = 0;
  double* fcolsp1 = fcols + f.stepj();
  std::vector<double> isqrt(order+3);
  for(int i=0;i<=order+2;i++) isqrt[i] = sqrt(double(i));

  for(int s=0;s<order;s++) 
  {
    if (s%2==1)
      fcolsp1[0] = isqrt[s]*te*fcolsm1[0]/isqrt[s+1]; // f(0,s+1)
    for(int p=s%2+1;p<=order+2;p+=2) 
    {
      double temp = isqrt[p]*se*fcols[p-1]; 
      if (s>0) temp += isqrt[s]*te*fcolsm1[p];
      temp /= isqrt[s+1]; 
      fcolsp1[p] = temp; // f(p,s+1)
    }
    fcolsm1 = fcols;
    fcols = fcolsp1;
    fcolsp1 += f.stepj();
  }

  for(int n=order+1,pq=size1;n<=order+2;++n) 
  {
    for(int p=n,q=0;p>=q;--p,++q,++pq) 
    {
      double* Spq = S.col(pq).ptr();
      double* Spq1 = p>q?S.col(pq+1).ptr():0;
      for(int nn=n%2,st=(nn==0?0:1);nn<=order;nn+=2,st+=nn) 
      {
	for(int s=nn,t=0;s>=t;--s,++t,++st) 
	{

	  double s0 = f.cref(p,s) * f.cref(q,t);
	  int iphase = s-t-p+q;
	  std::complex<double> s1 = s0 *
	    (iphase >= 0 ? phase[iphase/2] : std::conj(phase[-iphase/2]));

	  if (p==q) 
	  {
	    if (s==t) Spq[st] = std::real(s1);
	    else 
	    {
	      Spq[st] = std::real(s1);
	      Spq[++st] = std::imag(s1);
	    }
	  }
	  else if (s==t) 
	  {
	    Spq[st] = 2.*std::real(s1);
	    Spq1[st] = -2.*std::imag(s1);
	  }
	  else 
	  {
	    s0 = f.cref(q,s) * f.cref(p,t);
	    iphase = s-t-q+p;
	    std::complex<double> s2 = s0 *
	      (iphase >= 0 ? phase[iphase/2] : std::conj(phase[-iphase/2]));
	    Spq[st]= std::real(s1) + std::real(s2);
	    Spq1[st] = -std::imag(s1) + std::imag(s2);
	    Spq[++st] = std::imag(s1) + std::imag(s2);
	    Spq1[st] = std::real(s1) - std::real(s2);
	  }
	}
      }
      if (p!=q) ++pq;
    }
  }
}

void ApplyG(std::complex<double> g, BVec& b)
{
  if (g != 0.0) 
  {
    tmv::Matrix<double> S(b.size(),b.size());
    GTransform(g,b.GetOrder(),S);
    b = S * b;
  }
}

#if 0
void OldPSFConvolve(const BVec& bpsf, int order, double sigma,
    const tmv::MatrixView<double>& C)
{
  Assert(int(bpsf.size()) == (bpsf.GetOrder()+1)*(bpsf.GetOrder()+2)/2);
  Assert(int(C.colsize()) == (order+1)*(order+2)/2);
  Assert(int(C.rowsize()) == (order+1)*(order+2)/2);

  // C(st,pq) = 2sqrt(pi) Sum_uv sqrt(p!u!q!v!/s!t!)/w! 
  //                               G(s,p,u) G(t,q,v) bpsf_uv
  // The sum is only over terms for which  p+u-s == q+v-t,
  // and w = p+u-s = q+v-t >= 0
  // 
  // G(0,p,u) = binom(p+u,u) (-A)^u B^p
  // G(s+1,p,u) = A G(s,p-1,u) + B G(s,p,u-1)
  // where A = sigma_init / sigma_obs, B = sigma_psf / sigma_obs
  // D = A^2 = 1-B^2
  //
  // It is more efficient to combine the sqrt(p!u!/s!w!) into the G(s,p,u).
  // Call the product H(s,p,u).  We need to translate the above formulae:
  // H(0,p,u) = sqrt(p!u!/(p+u)!) (p+u)!/(p!u!) (-A)^u B^p
  //          = sqrt((p+u)!/(p!u!)) (-A)^u B^p
  // H(s+1,p,u) = sqrt(p!u!/(s+1)!(p+u-s-1)!) * [
  //                    A H(s,p-1,u) / sqrt((p-1)!u!/s!(p+u-s-1)!) +
  //                    B H(s,p,u-1) / sqrt(p!(u-1)!/s!(p+u-s-1)!) ]
  //            = A sqrt(p)/sqrt(s+1) H(s,p-1,u) + 
  //              B sqrt(u)/sqrt(s+1) H(s,p,u-1)
  // 
  const double twosqrtpi = 2.*sqrtpi;
  // sigma^2 = exp(mu)
  // D = sigma_i^2 / (sigma_i^2 + sigma_psf^2) 
  //   = 1 / (1 + (sigma_psf/sigma_i)^2)
  double D = 1./(1.+pow(bpsf.GetSigma()/sigma,2));
  double A = sqrt(D);
  double B = sqrt(1-D);

  int ntot = order+bpsf.GetOrder()+1;
  double isqrt[ntot];
  for(int i=0;i<ntot;i++) isqrt[i] = sqrt(double(i));

  std::vector<tmv::Matrix<double,tmv::ColMajor> > H(order+1,
      tmv::Matrix<double,tmv::ColMajor>(order+1,bpsf.GetOrder()+1));
  const int stepj = H[0].stepj();
  double* H0u = H[0].ptr();
  H0u[0] = 1.; // H[0](0,0)
  double* H0um1 = 0;
  for(int u=0;u<=bpsf.GetOrder();u++) 
  {
    if (u>0) H0u[0] = -A * H0um1[0];
    for(int p=1;p<=order;p++) 
      H0u[p] = B*isqrt[p+u]/isqrt[p]*H0u[p-1]; // H[0](p,u)
    H0um1 = H0u;
    H0u += stepj;
  }
  for(int s=0;s<order;s++) 
  {
    double* Hsu = H[s].ptr();
    double* Hsp1u = H[s+1].ptr();
    Hsp1u[0] = 0.;
    for(int p=1;p<=order;p++) 
      Hsp1u[p] = A*isqrt[p]*Hsu[p-1]/isqrt[s+1];
    double* Hsum1 = Hsu;
    Hsu += stepj;
    Hsp1u += stepj;
    for(int u=1;u<=bpsf.GetOrder();u++) 
    {
      Hsp1u[0] = B*isqrt[u]*Hsum1[0]/isqrt[s+1];
      for(int p=1;p<=order;p++) 
	Hsp1u[p] = (A*isqrt[p]*Hsu[p-1] + 
	    B*isqrt[u]*Hsum1[p])/ isqrt[s+1];
      Hsum1 = Hsu;
      Hsu += stepj;
      Hsp1u += stepj;
    }
  }

  C.Zero();
  int pq = 0;
  for(int n=0;n<=order;++n) 
  {
    for(int p=n,q=0;p>=q;(p==q?++pq:pq+=2),--p,++q) 
    {
      double* Cpq = C.col(pq).ptr();
      double* Cpq1 = p>q?C.col(pq+1).ptr():0;
      int st = 0;
      for(int nn=0;nn<=order;++nn) 
      {
	for(int s=nn,t=0;s>=t;(s==t?++st:st+=2),--s,++t) 
	{
	  double Cpqst = 0.;
	  double Cpq1st = 0.;
	  double Cpqst1 = 0.;
	  double Cpq1st1 = 0.;
	  const tmv::Matrix<double,tmv::ColMajor>& Hs = H[s];
	  const tmv::Matrix<double,tmv::ColMajor>& Ht = H[t];
	  int uv = 0;
	  for(int upv=0;upv<=bpsf.GetOrder();++upv) 
	  {
	    for(int u=upv,v=0;u>=v;(u==v?++uv:uv+=2),--u,++v) 
	    {

	      // First do terms with p>=q, u>=v  (always keep s>=t)
	      int w = p+u-s;
	      if (w >= 0 && q+v-t == w) 
	      {
		// s-t = p-q + u-v
		if (u==v) 
		{
		  double temp = Hs.cref(p,u)*Ht.cref(q,v)*bpsf(uv);
		  if (s==t) 
		  {
		    Assert(p==q);
		    Cpqst += temp;
		  }
		  else 
		  {
		    Cpqst += temp;
		    Cpq1st1 += temp;
		  }
		}
		else 
		{
		  Assert(s>t);
		  double tempr = Hs.cref(p,u)*Ht.cref(q,v);
		  double tempi = tempr * bpsf(uv+1);
		  tempr *= bpsf(uv);
		  //std::complex<double> temp = Hs(p,u)*Ht(q,v)*
		    //std::complex<double>(bpsf(uv),bpsf(uv+1));
		  if (p==q) 
		  {
		    Cpqst += tempr;
		    Cpqst1 += tempi;
		  }
		  else 
		  {
		    Assert(p>q);
		    Cpqst += tempr;
		    Cpqst1 += tempi;
		    Cpq1st -= tempi;
		    Cpq1st1 += tempr;
		  }
		}
	      }

	      // Next p<q, u>v.  Implement by swapping p,q
	      // These terms account for the fact that 
	      // b_init_qp = b_init_pq*
	      w = q+u-s;
	      if (w >= 0 && p>q && p+v-t == w) 
	      {
		// s-t = q-p + u-v
		Assert(u>v);
		double tempr = Hs.cref(q,u)*Ht.cref(p,v);
		double tempi = tempr * bpsf(uv+1);
		tempr *= bpsf(uv);
		//std::complex<double> temp = Hs(q,u)*Ht(p,v)*
		  //std::complex<double>(bpsf(uv),bpsf(uv+1));
		if (s==t) 
		{
		  Cpqst += tempr;
		  Cpq1st += tempi;
		}
		else {
		  Cpqst += tempr;
		  Cpqst1 += tempi;
		  Cpq1st += tempi;
		  Cpq1st1 -= tempr;
		}
	      }

	      // Next p>q, u<v.
	      // These terms account for b_psf_vu = b_psf_uv*
	      w = p+v-s;
	      if (w >= 0 && u>v && q+u-t == w) {
		// s-t = p-q + v-u
		Assert(p>q);
		double tempr = Hs.cref(p,v)*Ht.cref(q,u);
		double tempi = -tempr * bpsf(uv+1);
		tempr *= bpsf(uv);
		//std::complex<double> temp = Hs(p,v)*Ht(q,u)*
		  //std::complex<double>(bpsf(uv),-bpsf(uv+1));
		if (s==t) {
		  Cpqst += tempr;
		  Cpq1st -= tempi;
		}
		else {
		  Cpqst += tempr;
		  Cpqst1 += tempi;
		  Cpq1st -= tempi;
		  Cpq1st1 += tempr;
		}
	      }

	    }
	  }
	  Assert(uv == int(bpsf.size()));
	  Cpq[st] = Cpqst;
	  if (s != t) Cpq[st+1] = Cpqst1;
	  if (p != q) {
	    Cpq1[st] = Cpq1st;
	    if (s != t) Cpq1[st+1] = Cpq1st1;
	  }
	}
      }
      Assert(st == int(C.colsize()));
    }
  }
  Assert(pq == int(C.rowsize()));
  C *= twosqrtpi;
}
#endif

void PSFConvolve(const BVec& bpsf, int order, double sigma,
    const tmv::MatrixView<double>& C)
// Just make the matrix which multiplies b to effect the convolution.
// b_obs = C * b_init
// However, we do not actually use it this way.  Rather, we use it to 
// switch the model from:
// I = Sum psi_pq b_obs_pq
// to 
// I = Sum psi_pq C b_init_pq
// We use this to solve for the ML b_init.
{
  Assert(int(bpsf.size()) == (bpsf.GetOrder()+1)*(bpsf.GetOrder()+2)/2);
  Assert(int(C.colsize()) == (order+1)*(order+2)/2);
  Assert(int(C.rowsize()) == (order+1)*(order+2)/2);

  // C(st,pq) = 2sqrt(pi) Sum_uv sqrt(p!u!q!v!/s!t!)/w! 
  //                               G(s,p,u) G(t,q,v) bpsf_uv
  // The sum is only over terms for which  p+u-s == q+v-t,
  // and w = p+u-s = q+v-t >= 0
  // 
  // G(0,p,u) = binom(p+u,u) (-A)^u B^p
  // G(s+1,p,u) = A G(s,p-1,u) + B G(s,p,u-1)
  // where A = sigma_init / sigma_obs, B = sigma_psf / sigma_obs
  // D = A^2 = 1-B^2
  //
  // It is more efficient to combine the sqrt(p!u!/s!w!) into the G(s,p,u).
  // Call the product H(s,p,u).  We need to translate the above formulae:
  // H(0,p,u) = sqrt(p!u!/(p+u)!) (p+u)!/(p!u!) (-A)^u B^p
  //          = sqrt((p+u)!/(p!u!)) (-A)^u B^p
  // H(s+1,p,u) = sqrt(p!u!/(s+1)!(p+u-s-1)!) * [
  //                    A H(s,p-1,u) / sqrt((p-1)!u!/s!(p+u-s-1)!) +
  //                    B H(s,p,u-1) / sqrt(p!(u-1)!/s!(p+u-s-1)!) ]
  //            = A sqrt(p)/sqrt(s+1) H(s,p-1,u) + 
  //              B sqrt(u)/sqrt(s+1) H(s,p,u-1)
  // 
  const double twosqrtpi = 2.*sqrtpi;
  // sigma^2 = exp(mu)
  // D = sigma_i^2 / (sigma_i^2 + sigma_psf^2) 
  //   = 1 / (1 + (sigma_psf/sigma_i)^2)
  double D = 1./(1.+pow(bpsf.GetSigma()/sigma,2));
  double A = sqrt(D);
  double B = sqrt(1-D);

  int ntot = order+bpsf.GetOrder()+1;
  std::vector<double> isqrt(ntot);
  for(int i=0;i<ntot;i++) isqrt[i] = sqrt(double(i));

  std::vector<tmv::Matrix<double,tmv::ColMajor> > H(order+1,
      tmv::Matrix<double,tmv::ColMajor>(order+1,bpsf.GetOrder()+1));
  const int stepj = H[0].stepj();
  double* H0u = H[0].ptr();
  H0u[0] = 1.; // H[0](0,0)
  double* H0um1 = 0;
  for(int u=0;u<=bpsf.GetOrder();u++) {
    if (u>0) H0u[0] = -A * H0um1[0];
    for(int p=1;p<=order;p++) 
      H0u[p] = B*isqrt[p+u]/isqrt[p]*H0u[p-1]; // H[0](p,u)
    H0um1 = H0u;
    H0u += stepj;
  }
  for(int s=0;s<order;s++) {
    double* Hsu = H[s].ptr();
    double* Hsp1u = H[s+1].ptr();
    Hsp1u[0] = 0.;
    for(int p=1;p<=order;p++) 
      Hsp1u[p] = A*isqrt[p]*Hsu[p-1]/isqrt[s+1];
    double* Hsum1 = Hsu;
    Hsu += stepj;
    Hsp1u += stepj;
    for(int u=1;u<=bpsf.GetOrder();u++) {
      Hsp1u[0] = B*isqrt[u]*Hsum1[0]/isqrt[s+1];
      for(int p=1;p<=order;p++) 
	Hsp1u[p] = (A*isqrt[p]*Hsu[p-1] + 
	    B*isqrt[u]*Hsum1[p])/ isqrt[s+1];
      Hsum1 = Hsu;
      Hsu += stepj;
      Hsp1u += stepj;
    }
  }

  C.Zero();
  int pq = 0;
  for(int n=0;n<=order;++n) {
    for(int p=n,q=0;p>=q;(p==q?++pq:pq+=2),--p,++q) {
      int pmq = p-q;
      double* Cpq = C.col(pq).ptr();
      double* Cpq1 = p>q?C.col(pq+1).ptr():0;
      int st = 0;
      for(int nn=0;nn<=order;++nn) {
	for(int s=nn,t=0;s>=t;(s==t?++st:st+=2),--s,++t) {
	  int smt = s-t;
	  double Cpqst = 0.;
	  double Cpq1st = 0.;
	  double Cpqst1 = 0.;
	  double Cpq1st1 = 0.;
	  const tmv::Matrix<double,tmv::ColMajor>& Hs = H[s];
	  const tmv::Matrix<double,tmv::ColMajor>& Ht = H[t];
	  int uv0 = 0;
	  int parity = (n+nn)%2;
	  for(int upv=0;upv<=bpsf.GetOrder();++upv,uv0+=upv) {
	    if (upv % 2 != parity) continue;
	    // There are three values of u-v that are worth considering:
	    // u-v = (s-t) - (p-q) >= 0
	    // u-v = (s-t) + (p-q) > 0
	    // u-v = (p-q) - (s-t) < 0

	    {
	      int umv = smt-pmq;
	      if (umv >= 0 && umv <= upv) {
		// First do terms with p>=q, u>=v  (always keep s>=t)
		// s-t = p-q + u-v
		int u = (upv+umv)/2;
		int v = (upv-umv)/2;

		int w = p+u-s;
		Assert(q+v-t == w);
		//Assert((w >= 0) == (umv >= 0));
		if (w >= 0) {
		  int uv = uv0 + 2*v;
		  if (umv == 0) {
		    double temp = Hs.cref(p,u)*Ht.cref(q,v)*bpsf(uv);
		    if (s==t) {
		      Assert(p==q);
		      Cpqst += temp;
		    }
		    else {
		      Cpqst += temp;
		      Cpq1st1 += temp;
		    }
		  } else {
		    Assert(s>t);
		    double tempr = Hs.cref(p,u)*Ht.cref(q,v);
		    double tempi = tempr * bpsf(uv+1);
		    tempr *= bpsf(uv);
		    //std::complex<double> temp = Hs(p,u)*Ht(q,v)*
		    //std::complex<double>(bpsf(uv),bpsf(uv+1));
		    if (p==q) {
		      Cpqst += tempr;
		      Cpqst1 += tempi;
		    }
		    else {
		      Assert(p>q);
		      Cpqst += tempr;
		      Cpqst1 += tempi;
		      Cpq1st -= tempi;
		      Cpq1st1 += tempr;
		    }
		  }
		}
	      }
	    }

	    {
	      int umv = smt+pmq;
	      if (pmq != 0 && umv <= upv) {
		// Next p<q, u>v.  Implement by swapping p,q
		// These terms account for the fact that 
		// b_init_qp = b_init_pq*
		// s-t = q-p + u-v

		Assert(umv > 0);
		int u = (upv+umv)/2;
		int v = (upv-umv)/2;

		int w = q+u-s;
		Assert(p+v-t == w);
		//Assert((w >= 0) == (umv > 0));
		if (w >= 0) {
		  int uv = uv0 + 2*v;
		  //Assert(w > 0);
		  //Assert((w > 0) == (umv > 0));
		  Assert(u>v);
		  double tempr = Hs.cref(q,u)*Ht.cref(p,v);
		  double tempi = tempr * bpsf(uv+1);
		  tempr *= bpsf(uv);
		  //std::complex<double> temp = Hs(q,u)*Ht(p,v)*
		  //std::complex<double>(bpsf(uv),bpsf(uv+1));
		  if (smt==0) {
		    Cpqst += tempr;
		    Cpq1st += tempi;
		  }
		  else {
		    Cpqst += tempr;
		    Cpqst1 += tempi;
		    Cpq1st += tempi;
		    Cpq1st1 -= tempr;
		  }
		}
	      }
	    }

	    {
	      int umv = pmq-smt;
	      if (umv > 0 && umv <= upv) {
		// Next p>q, u<v.
		// These terms account for b_psf_vu = b_psf_uv*
		int u = (upv+umv)/2;
		int v = (upv-umv)/2;

		int w = p+v-s;
		Assert(q+u-t == w);
		if (w >= 0) {
		  int uv = uv0 + 2*v;
		  // s-t = p-q + v-u
		  Assert(p>q);
		  double tempr = Hs.cref(p,v)*Ht.cref(q,u);
		  double tempi = -tempr * bpsf(uv+1);
		  tempr *= bpsf(uv);
		  //std::complex<double> temp = Hs(p,v)*Ht(q,u)*
		  //std::complex<double>(bpsf(uv),-bpsf(uv+1));
		  if (smt==0) {
		    Cpqst += tempr;
		    Cpq1st -= tempi;
		  }
		  else {
		    Cpqst += tempr;
		    Cpqst1 += tempi;
		    Cpq1st -= tempi;
		    Cpq1st1 += tempr;
		  }
		}
	      }
	    }
	  }
	  Assert(uv0 == int(bpsf.size()));
	  Cpq[st] = Cpqst;
	  if (smt != 0) Cpq[st+1] = Cpqst1;
	  if (pmq != 0) {
	    Cpq1[st] = Cpq1st;
	    if (smt != 0) Cpq1[st+1] = Cpq1st1;
	  }
	}
      }
      Assert(st == int(C.colsize()));
    }
  }
  Assert(pq == int(C.rowsize()));
  C *= twosqrtpi;
}

void ApplyPSF(const BVec& bpsf, BVec& b)
{
  tmv::Matrix<double> C(b.size(),b.size());
  PSFConvolve(bpsf,b.GetOrder(),b.GetSigma(),C);
  b = C * b;
}
