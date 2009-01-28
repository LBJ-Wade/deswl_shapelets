
#include "PsiHelper.h"
#include "TMV.h"
#include "dbg.h"

const double PI = 3.14159265359;
const double sqrtpi = sqrt(PI);

// Here is the order of p,q along the indices of psi:
//
// k 0 1,2 3,4 5 6,7 8,9 10,11 12,13 14 15,16 17,18 19,20 21,22 23,24 25,26 27
// p 0  1   2  1  3   2    4     3    2   5     4     3     6     5     4    3
// q 0  0   0  1  0   1    0     1    2   0     1     2     0     1     2    3
// n 0  1   2  2  3   3    4     4    4   5     5     5     6     6     6    6
// m 0  1   2  0  3   1    4     2    0   5     3     1     6     4     2    0
//


void MakePsi1(const tmv::VectorView<double>& psi, 
    std::complex<double> z, int order)
{
  Assert(int(psi.size()) == (order+1)*(order+2)/2);
  const double invsqrtpi = 1./sqrtpi;

  // coeff(p,q) = sqrt( 1/ (pi p! q!) )
  tmv::DiagMatrix<double> coeff(psi.size());
  std::vector<double> sqrtfact(order+1);
  sqrtfact[0] = 1.;
  if (order > 0) sqrtfact[1] = 1.;
  for(int n=0,k=0;n<=order;++n) {
    if (n > 1) sqrtfact[n] = sqrtfact[n-1] * sqrt(double(n));
    for(int p=n,q=0;p>=q;--p,++q,++k) {
      double temp = invsqrtpi / sqrtfact[p] / sqrtfact[q];
      coeff(k) = temp;
      if (p!=q) coeff(++k) = temp;
    }
  }

  double rsq = std::norm(z);
  // psi(p,q) = coeff(p,q) z^m exp(-r^2/2) [(-)^q q! Lqm(r^2)]
  // We'll multiply everything by exp(-r^2/2) at the end.
  // We do coeff, z^m, and the [] here:
  // Define the [] factor as Kqm = (-)^q q! Lqm(r^2)
  //
  // The recursion relation for the Kqm is: 
  // Kqm = (x-(2q+m-1)) K(q-1)m - (q+m-1)*(q-1) K(q-2)m
  // K0m = 1
  // K1m = x-m-1

  // First m=0:
  psi[0] = 1.;
  if (order >= 2) {
    double Kqm = rsq-1.;
    double Kq1m = 1.;
    psi[5] = Kqm;
    double c1 = 1.; // = 2*q+m-1
    double c2 = 0.; // = (q+m-1)*(q-1)
    int k = 5;
    for(int q=2;q<=order/2;++q) {
      c2 += c1;
      c1 += 2.;
      double K = (rsq-c1) * Kqm - c2 * Kq1m;
      Kq1m = Kqm; Kqm = K;
      k += 4*q+1;
      psi[k] = Kqm;
    }
  }

  // Next m > 0:
  // All m > 0 elements are intrinsically complex.
  // However, we are fitting to a real intensity pattern
  // with complex coefficients of the complex shapelets.
  // e.g. b_20 psi_20
  // Since psi_02 = psi_20* (* = complex conjugate),
  // we know that b_02 must be b_20*
  // b_20 psi_20 + b_20* psi_20*
  // = 2 b_20r psi_20r - 2 b_20i psi_20i
  // So the values we want for the real fitter are
  // really 2 real(psi_pq) and -2 imag(psi_pq)
  std::complex<double> zm = z;
  for(int m=1;m<=order;++m) {
    if (m>1) zm *= z;
    int k = m*(m+1)/2;
    psi[k] = 2.*std::real(zm);
    psi[k+1] = -2.*std::imag(zm);
    if (order >= m+2) {
      double c1 = m+1.;  // = 2*q+m-1
      double c2 = 0.;    // = (q+m-1)*(q-1)
      double Kqm = rsq-c1;
      double Kq1m = 1.;
      k += 2*m+5;
      psi[k] = 2.*Kqm * std::real(zm);
      psi[k+1] = -2.*Kqm * std::imag(zm);
      for(int q=2;2*q+m<=order;++q) {
	c2 += c1;
	c1 += 2.;
	double K = (rsq-c1) * Kqm - c2 * Kq1m;
	Kq1m = Kqm; Kqm = K;
	k += 2*m+4*q+1;
	psi[k] = 2.*Kqm * std::real(zm);
	psi[k+1] = -2.*Kqm * std::imag(zm);
      }
    }
  }

  // Now put in the exp(-r^2/2) factor and all the 1/sqrt(pi p! q!) coeffs:
  psi *= exp(-rsq/2.) * coeff;
}

void MakePsi(const tmv::MatrixView<double>& psi, 
    const tmv::Vector<std::complex<double> >& z, int order,
    const tmv::Vector<double>* coeff)
{
  //dbg<<"Start MakePsi"<<std::endl;
  // For p>=q:
  //
  // psi_pq = (pi p! q!)^-1/2 z^m exp(-r^2/2) K_pq(r^2)
  //
  // where K_pq(r^2) satisfies the recursion relation:
  // K_pq = (r^2 - (N-1)) K_p-1,q-1 - (p-1)(q-1) K_p-2,q-2
  // with K_p0 = 1
  //
  // The recursion for psi_pq can then be derived as:
  //
  // psi_pq = (pq)^-1/2 (r^2-(N-1)) psi_p-1,q-1
  //          - sqrt( (p-1)(q-1)/(pq) ) psi_p-2,q-2
  //
  // psi_p0 = (z/sqrt(p)) psi_p-1,0
  // 
  // psi_00 = 1/sqrt(pi) exp(-r^2/2)

  Assert(int(psi.rowsize()) == (order+1)*(order+2)/2);
  Assert(psi.colsize() == z.size());
  if (coeff) Assert(psi.colsize() == coeff->size());
  Assert(psi.iscm());
  Assert(!psi.isconj());

  const double invsqrtpi = 1./sqrtpi;

  // Setup rsq, z vectors and set psi_00
  tmv::Vector<double> rsq(z.size());
  double* rsqit = rsq.ptr();
  double* psi00it = psi.ptr();
  const std::complex<double>* zit = z.cptr();
  //dbg<<"Before loop1"<<std::endl;
  //dbg<<"psi.col(0) = "<<psi.col(0)<<std::endl;
  for(size_t i=0;i<z.size();++i) {
    //dbg<<"i = "<<i<<std::endl;
    rsqit[i] = std::norm(zit[i]);
    //dbg<<"rsq = "<<rsqit[i]<<std::endl;
    psi00it[i] = invsqrtpi * exp(-(rsqit[i])/2.);
    //dbg<<"psi_i0 = "<<psi00it[i]<<std::endl;
  }
  if (coeff) psi.col(0) *= DiagMatrixViewOf(*coeff);
  //dbg<<"psi.col(0) => "<<psi.col(0)<<std::endl;

  tmv::Vector<double> zr = z.Real();
  tmv::Vector<double> zi = z.Imag();
  if (order >= 1) {
    // Set psi_10
    // All m > 0 elements are intrinsically complex.
    // However, we are fitting to a real intensity pattern
    // with complex coefficients of the complex shapelets.
    // Since psi_pq = psi_qp* (* = complex conjugate),
    // we know that b_pq must be b_qp*
    // b_pq psi_pq + b_pq* psi_pq* = 2 b_pqr psi_pqr - 2 b_pqi psi_pqi
    // So the values we want for the real fitter are
    // really 2 real(psi_pq) and -2 imag(psi_pq)
    // Putting the 2's here carries through to the rest of the 
    // elements via the recursion.
    psi.col(1) = 2. * DiagMatrixViewOf(zr) * psi.col(0);
    psi.col(2) = -2. * DiagMatrixViewOf(zi) * psi.col(0);
  }
  //dbg<<"Before loop2"<<std::endl;
  for(int N=2,k=3;N<=order;++N) {
    //dbg<<"N = "<<N<<", k = "<<k<<std::endl;
    // Set psi_N0
    // The signs of these are not what you naively think due to 
    // the +2, -2 discussed above.  You just have to follow through
    // what the complex psi_N0 is, and what value is stored in the
    // psi_N-1,0 location, and what needs to get stored here.
    psi.col(k) = sqrt(1./N) * DiagMatrixViewOf(zr) * psi.col(k-N);
    psi.col(k) += sqrt(1./N) * DiagMatrixViewOf(zi) * psi.col(k-N+1);
    psi.col(k+1) = -sqrt(1./N) * DiagMatrixViewOf(zi) * psi.col(k-N);
    psi.col(k+1) += sqrt(1./N) * DiagMatrixViewOf(zr) * psi.col(k-N+1);
    k+=2;

    // Set psi_pq with q>0
    // The rsq part of this calculation can be done in batch, which 
    // speeds things up a bit.
    psi.Cols(k,k+N-1) = DiagMatrixViewOf(rsq) * psi.Cols(k-2*N-1,k-N-2);
    psi.Cols(k,k+N-1) -= (N-1.) * psi.Cols(k-2*N-1,k-N-2);
    // The other calculation steps are different for each component:
    for(int m=N-2,p=N-1,q=1;m>=0;--p,++q,m-=2) {
      double pq = p*q;
      if (m==0) {
	psi.col(k) /= sqrt(pq);
	if (q > 1) psi.col(k) -= sqrt(1.-(N-1.)/pq)*psi.col(k+2-4*N);
	++k;
      } else {
	psi.Cols(k,k+2) /= sqrt(pq);
	if (q > 1)
	  psi.Cols(k,k+2) -= sqrt(1.-(N-1.)/pq)*psi.Cols(k+2-4*N,k+4-4*N);
	k+=2;
      }
    }
  }
}

void AugmentPsi(tmv::Matrix<double>& psi,
    const tmv::Vector<std::complex<double> >& z, int order)
{
  Assert(int(psi.rowsize()) == (order+3)*(order+4)/2);
  Assert(psi.colsize() == z.size());
  Assert(order >= 1);
  Assert(psi.iscm());
  Assert(!psi.isconj());

  tmv::Vector<double> rsq(z.size());
  double* rsqit = rsq.ptr();
  const std::complex<double>* zit = z.cptr();
  for(size_t i=0;i<z.size();++i) {
    rsqit[i] = std::norm(zit[i]);
  }

  tmv::Vector<double> zr = z.Real();
  tmv::Vector<double> zi = z.Imag();
  for(int N=order+1,k=N*(N+1)/2;N<=order+2;++N) {
    psi.col(k) = sqrt(1./N) * DiagMatrixViewOf(zr) * psi.col(k-N);
    psi.col(k) += sqrt(1./N) * DiagMatrixViewOf(zi) * psi.col(k-N+1);
    psi.col(k+1) = -sqrt(1./N) * DiagMatrixViewOf(zi) * psi.col(k-N);
    psi.col(k+1) += sqrt(1./N) * DiagMatrixViewOf(zr) * psi.col(k-N+1);
    k+=2;

    psi.Cols(k,k+N-1) = DiagMatrixViewOf(rsq) * psi.Cols(k-2*N-1,k-N-2);
    psi.Cols(k,k+N-1) -= (N-1.) * psi.Cols(k-2*N-1,k-N-2);
    for(int m=N-2,p=N-1,q=1;m>=0;--p,++q,m-=2) {
      double pq = p*q;
      if (m==0) {
	psi.col(k) /= sqrt(pq);
	if (q > 1) psi.col(k) -= sqrt(1.-(N-1.)/pq)*psi.col(k+2-4*N);
	++k;
      } else {
	psi.Cols(k,k+2) /= sqrt(pq);
	if (q > 1)
	  psi.Cols(k,k+2) -= sqrt(1.-(N-1.)/pq)*psi.Cols(k+2-4*N,k+4-4*N);
	k+=2;
      }
    }
  }
}

void SetupGx(tmv::Matrix<double>& Gx, int order1, int order2)
{
  Assert(int(Gx.colsize()) == (order1+1)*(order1+2)/2);
  Assert(int(Gx.rowsize()) == (order2+1)*(order2+2)/2);

  Gx.Zero();
  for(int n=0,k=0;n<=order2;n++) for(int p=n,q=0;p>=q;p--,q++,k++) {
    //std::cout<<"n,k,p,q = "<<n<<','<<k<<','<<p<<','<<q<<std::endl;
    double dp = p;
    double dq = q;
    // d(psi(x,y))/dx = psi(x,y) Gx
    // d/dx = 1/2(ap + aq - apt - aqt)
    // Gx( p-1,q , pq ) = 1/2 sqrt(p)
    // Gx( p,q-1 , pq ) = 1/2 sqrt(q)
    // Gx( p+1,q , pq ) = -1/2 sqrt(p+1)
    // Gx( p,q+1 , pq ) = -1/2 sqrt(q+1)
    if (p==q) {
      if (q>0 && n<=order1+1) Gx(k-n-2,k) = sqrt(dp)/2.;
      if (n+1<=order1) Gx(k+n+1,k) = -sqrt(dp+1.)/2.;
    } else if (p==q+1) {
      if (n<=order1+1) Gx(k-n,k) = sqrt(dp);
      if (q>0 && n<=order1+1) Gx(k-n-2,k) = sqrt(dq)/2.;
      if (n+1<=order1) Gx(k+n+1,k) = -sqrt(dp+1.)/2.;
      if (n+1<=order1) Gx(k+n+3,k) = -sqrt(dq+1.);
      if (q>0 && n<=order1+1) Gx(k-n-1,k+1) = sqrt(dq)/2.;
      if (n+1<=order1) Gx(k+n+2,k+1) = -sqrt(dp+1.)/2.;
    } else {
      if (n<=order1+1) Gx(k-n,k) = sqrt(dp)/2.;
      if (q>0 && n<=order1+1) Gx(k-n-2,k) = sqrt(dq)/2.;
      if (n+1<=order1) Gx(k+n+1,k) = -sqrt(dp+1.)/2.;
      if (n+1<=order1) Gx(k+n+3,k) = -sqrt(dq+1.)/2.;
      if (n<=order1+1) Gx(k-n+1,k+1) = sqrt(dp)/2.;
      if (q>0 && n<=order1+1) Gx(k-n-1,k+1) = sqrt(dq)/2.;
      if (n+1<=order1) Gx(k+n+2,k+1) = -sqrt(dp+1.)/2.;
      if (n+1<=order1) Gx(k+n+4,k+1) = -sqrt(dq+1.)/2.;
    }
    if (p > q) ++k;
  }
  //std::cout<<"after define Gx"<<std::endl;
}

void SetupGy(tmv::Matrix<double>& Gy, int order1, int order2)
{
  Assert(int(Gy.colsize()) == (order1+1)*(order1+2)/2);
  Assert(int(Gy.rowsize()) == (order2+1)*(order2+2)/2);

  Gy.Zero();
  for(int n=0,k=0;n<=order2;n++) for(int p=n,q=0;p>=q;p--,q++,k++) {
    double dp = p;
    double dq = q;
    // d(psi(x,y))/dx = psi(x,y) Gx
    // d/dy = 1/2 i (ap - aq + apt - aqt)
    // Gy( p-1,q , pq ) = 1/2 i sqrt(p)
    // Gy( p,q-1 , pq ) = -1/2 i sqrt(q)
    // Gy( p+1,q , pq ) = 1/2 i sqrt(p+1)
    // Gy( p,q+1 , pq ) = -1/2 i sqrt(q+1)
    if (p==q) {
      if (q>0 && n<=order1+1) Gy(k-n-1,k) = -sqrt(dp)/2.;
      if (n+1<=order1) Gy(k+n+2,k) = sqrt(dp+1.)/2.;
    } else if (p==q+1) {
      if (q>0 && n<=order1+1) Gy(k-n-1,k) = -sqrt(dq)/2.;
      if (n+1<=order1) Gy(k+n+2,k) = sqrt(dp+1.)/2.;
      if (n<=order1+1) Gy(k-n,k+1) = -sqrt(dp);
      if (q>0 && n<=order1+1) Gy(k-n-2,k+1) = sqrt(dq)/2.;
      if (n+1<=order1) Gy(k+n+1,k+1) = -sqrt(dp+1.)/2.;
      if (n+1<=order1) Gy(k+n+3,k+1) = sqrt(dq+1.);
    } else {
      if (n+1<=order1) Gy(k-n+1,k) = sqrt(dp)/2.;
      if (q>0 && n<=order1+1) Gy(k-n-1,k) = -sqrt(dq)/2.;
      if (n+1<=order1) Gy(k+n+2,k) = sqrt(dp+1.)/2.;
      if (n+1<=order1) Gy(k+n+4,k) = -sqrt(dq+1.)/2.;
      if (n<=order1+1) Gy(k-n,k+1) = -sqrt(dp)/2.;
      if (q>0 && n<=order1+1) Gy(k-n-2,k+1) = sqrt(dq)/2.;
      if (n+1<=order1) Gy(k+n+1,k+1) = -sqrt(dp+1.)/2.;
      if (n+1<=order1) Gy(k+n+3,k+1) = sqrt(dq+1.)/2.;
    }

    if (p > q) ++k;
  }
}

void SetupGg1(tmv::Matrix<double>& Gg1, int order1, int order2)
{
  Assert(int(Gg1.colsize()) == (order1+1)*(order1+2)/2);
  Assert(int(Gg1.rowsize()) == (order2+1)*(order2+2)/2);

  Gg1.Zero();
  for(int n=0,k=0;n<=order2;n++) for(int p=n,q=0;p>=q;p--,q++,k++) {
    double dp = p;
    double dq = q;
    // d(psi(x,y))/dg1 = psi(x,y) Gg1
    // Gg1 is G_g1
    // d/dg1 = x d/dx - y d/dy
    //       = 1/2 (ap^2 + aq^2 - apt^2 - aqt^2)
    // Gg1( p-2,q , pq ) = 1/2 sqrt(p(p-1))
    // Gg1( p,q-2 , pq ) = 1/2 sqrt(q(q-1))
    // Gg1( p+2,q , pq ) = -1/2 sqrt((p+1)(p+2))
    // Gg1( p,q+2 , pq ) = -1/2 sqrt((q+1)(q+2))
    if (p==q) {
      if (q>1 && n<=order1+2) Gg1(k-2*n-3,k) = sqrt(dp*(dp-1.))/2.;
      if (n+2<=order1) Gg1(k+2*n+3,k) = -sqrt((dp+1.)*(dp+2.))/2.;
    } else if (p==q+1) {
      if (p>1 && n<=order1+2) Gg1(k-2*n-1,k) = sqrt(dp*(dp-1.))/2.;
      if (q>1 && n<=order1+2) Gg1(k-2*n-3,k) = sqrt(dq*(dq-1.))/2.;
      if (n+2<=order1) Gg1(k+2*n+3,k) = -sqrt((dp+1.)*(dp+2.))/2.;
      if (n+2<=order1) Gg1(k+2*n+5,k) = -sqrt((dq+1.)*(dq+2.))/2.;
      if (p>1 && n<=order1+2) Gg1(k-2*n,k+1) = -sqrt(dp*(dp-1.))/2.;
      if (q>1 && n<=order1+2) Gg1(k-2*n-2,k+1) = sqrt(dq*(dq-1.))/2.;
      if (n+2<=order1) Gg1(k+2*n+4,k+1) = -sqrt((dp+1)*(dp+2.))/2.;
      if (n+2<=order1) Gg1(k+2*n+6,k+1) = sqrt((dq+1)*(dq+2.))/2.;
    } else if (p==q+2) {
      if (n<=order1+2) Gg1(k-2*n+1,k) = sqrt(dp*(dp-1.));
      if (q>1 && n<=order1+2) Gg1(k-2*n-3,k) = sqrt(dq*(dq-1.))/2.;
      if (n+2<=order1) Gg1(k+2*n+3,k) = -sqrt((dp+1.)*(dp+2.))/2.;
      if (n+2<=order1) Gg1(k+2*n+7,k) = -sqrt((dq+1.)*(dq+2.));
      if (q>1 && n<=order1+2) Gg1(k-2*n-2,k+1) = sqrt(dq*(dq-1.))/2.;
      if (n+2<=order1) Gg1(k+2*n+4,k+1) = -sqrt((dp+1.)*(dp+2.))/2.;
    } else {
      if (n<=order1+2) Gg1(k-2*n+1,k) = sqrt(dp*(dp-1.))/2.;
      if (q>1 && n<=order1+2) Gg1(k-2*n-3,k) = sqrt(dq*(dq-1.))/2.;
      if (n+2<=order1) Gg1(k+2*n+3,k) = -sqrt((dp+1.)*(dp+2.))/2.;
      if (n+2<=order1) Gg1(k+2*n+7,k) = -sqrt((dq+1.)*(dq+2.))/2.;
      if (n<=order1+2) Gg1(k-2*n+2,k+1) = sqrt(dp*(dp-1.))/2.;
      if (q>1 && n<=order1+2) Gg1(k-2*n-2,k+1) = sqrt(dq*(dq-1.))/2.;
      if (n+2<=order1) Gg1(k+2*n+4,k+1) = -sqrt((dp+1.)*(dp+2.))/2.;
      if (n+2<=order1) Gg1(k+2*n+8,k+1) = -sqrt((dq+1.)*(dq+2.))/2.;
    }

    if (p > q) ++k;
  }
}

void SetupGg2(tmv::Matrix<double>& Gg2, int order1, int order2)
{
  Assert(int(Gg2.colsize()) == (order1+1)*(order1+2)/2);
  Assert(int(Gg2.rowsize()) == (order2+1)*(order2+2)/2);

  Gg2.Zero();
  for(int n=0,k=0;n<=order2;n++) for(int p=n,q=0;p>=q;p--,q++,k++) {
    double dp = p;
    double dq = q;
    // d(psi(x,y))/dg2 = psi(x,y) Gg2
    // Gg2 is G_g2
    // d/dg2 = y d/dx + x d/dy
    //       = 1/2 i (ap^2 - aq^2 + apt^2 - aqt^2)
    // Gg2( p-2,q , pq ) = 1/2 i sqrt(p(p-1))
    // Gg2( p,q-2 , pq ) = -1/2 i sqrt(q(q-1))
    // Gg2( p+2,q , pq ) = 1/2 i sqrt((p+1)(p+2))
    // Gg2( p,q+2 , pq ) = -1/2 i sqrt((q+1)(q+2))
    if (p==q) {
      if (q>1 && n<=order1+2) Gg2(k-2*n-2,k) = -sqrt(dp*(dp-1.))/2.;
      if (n+2<=order1) Gg2(k+2*n+4,k) = sqrt((dp+1.)*(dp+2.))/2.;
    } else if (p==q+1) {
      if (p>1 && n<=order1+2) Gg2(k-2*n,k) = -sqrt(dp*(dp-1.))/2.;
      if (q>1 && n<=order1+2) Gg2(k-2*n-2,k) = -sqrt(dq*(dq-1.))/2.;
      if (n+2<=order1) Gg2(k+2*n+4,k) = sqrt((dp+1.)*(dp+2.))/2.;
      if (n+2<=order1) Gg2(k+2*n+6,k) = sqrt((dq+1.)*(dq+2.))/2.;
      if (p>1 && n<=order1+2) Gg2(k-2*n-1,k+1) = -sqrt(dp*(dp-1.))/2.;
      if (q>1 && n<=order1+2) Gg2(k-2*n-3,k+1) = sqrt(dq*(dq-1.))/2.;
      if (n+2<=order1) Gg2(k+2*n+3,k+1) = -sqrt((dp+1.)*(dp+2.))/2.;
      if (n+2<=order1) Gg2(k+2*n+5,k+1) = sqrt((dq+1.)*(dq+2.))/2.;
    } else if (p==q+2) {
      if (q>1 && n<=order1+2) Gg2(k-2*n-2,k) = -sqrt(dq*(dq-1.))/2.;
      if (n+2<=order1) Gg2(k+2*n+4,k) = sqrt((dp+1.)*(dp+2.))/2.;
      if (n<=order2+2) Gg2(k-2*n+1,k+1) = -sqrt(dp*(dp-1.));
      if (q>1 && n<=order1+2) Gg2(k-2*n-3,k+1) = sqrt(dq*(dq-1.))/2.;
      if (n+2<=order1) Gg2(k+2*n+3,k+1) = -sqrt((dp+1.)*(dp+2.))/2.;
      if (n+2<=order1) Gg2(k+2*n+7,k+1) = sqrt((dq+1.)*(dq+2.));
    } else {
      Gg2(k-2*n+2,k) = sqrt(dp*(dp-1.))/2.;
      if (q>1 && n<=order1+2) Gg2(k-2*n-2,k) = -sqrt(dq*(dq-1.))/2.;
      if (n+2<=order1) Gg2(k+2*n+4,k) = sqrt((dp+1.)*(dp+2.))/2.;
      if (n+2<=order1) Gg2(k+2*n+8,k) = -sqrt((dq+1.)*(dq+2.))/2.;
      if (n<=order2+2) Gg2(k-2*n+1,k+1) = -sqrt(dp*(dp-1.))/2.;
      if (q>1 && n<=order1+2) Gg2(k-2*n-3,k+1) = sqrt(dq*(dq-1.))/2.;
      if (n+2<=order1) Gg2(k+2*n+3,k+1) = -sqrt((dp+1.)*(dp+2.))/2.;
      if (n+2<=order1) Gg2(k+2*n+7,k+1) = sqrt((dq+1.)*(dq+2.))/2.;
    }

    if (p > q) ++k;
  }
}

void SetupGmu(tmv::Matrix<double>& Gmu, int order1, int order2)
{
  Assert(int(Gmu.colsize()) == (order1+1)*(order1+2)/2);
  Assert(int(Gmu.rowsize()) == (order2+1)*(order2+2)/2);

  Gmu.Zero();
  for(int n=0,k=0;n<=order2;n++) for(int p=n,q=0;p>=q;p--,q++,k++) {
    double dp = p;
    double dq = q;
    // d(psi(x,y))/dmu = psi(x,y) Gmu
    // d/dmu = x d/dx + y d/dy
    //         ap aq - apt aqt - 1
    // Gmu( p-1,q-1 , pq ) = sqrt(pq)
    // Gmu( p+1,q+1 , pq ) = -sqrt((p+1)(q+1))
    // Gmu( pq , pq ) = -1
    if (q>0 && n<=order1+2) Gmu(k-2*n-1,k) = sqrt(dp*dq);
    if (n+2<=order1) Gmu(k+2*n+5,k) = -sqrt((dp+1.)*(dq+1.));
    if (n<=order1) Gmu(k,k) = -1.;
    if (p > q) {
      if (q>0 && n<=order1+2) Gmu(k-2*n,k+1) = sqrt(dp*dq);
      if (n+2<=order1) Gmu(k+2*n+6,k+1) = -sqrt((dp+1.)*(dq+1.));
      if (n<=order1) Gmu(k+1,k+1) = -1.;
    }

    if (p > q) ++k;
  }
}

void SetupGth(tmv::Matrix<double>& Gth, int order1, int order2)
{
  Assert(int(Gth.colsize()) == (order1+1)*(order1+2)/2);
  Assert(int(Gth.rowsize()) == (order2+1)*(order2+2)/2);

  Gth.Zero();
  for(int n=0,k=0;n<=order2;n++) for(int p=n,q=0;p>=q;p--,q++,k++) {
    double dp = p;
    double dq = q;
    // d(psi(x,y))/dtheta = psi(x,y) Gth
    // d/dtheta = x d/dy - y d/dx
    //          = m i
    // Gth( pq , pq ) = m i
    if (p>q && n<=order1) {
      Gth(k+1,k) = (dp-dq);
      Gth(k,k+1) = -(dp-dq);
    }

    if (p > q) ++k;
  }
}

