#include <cmath>
#include <iostream>
#include <string>
#include <ios>
#include <assert.h>

#include "Matrix.h"
#include "BinomFact.h"
#include "Laguerre.h"
#include "Solve.h"

using mv::DMatrix;
using mv::CMatrix;
using mv::SqDMatrix;
using mv::SqCMatrix;
using mv::DVector;
using mv::CVector;
using mv::Vector;
using mv::Matrix;
using mv::SqMatrix;

namespace laguerre {

// Build an LVector from a DVector for the real degrees of freedom.
// DVector must have same dimension as needed for LVector of chosen order
LVector::LVector(const DVector& rhs, int order_): 
  order(new int(order_)), 
  pcount(new int(1)),
  v(new DVector(rhs)) {
  if (v->size()!=PQIndex::size(*order)) {
      delete v;
      delete pcount;
      delete order;
      throw LaguerreError("Input to LVector(DVector) is wrong size for order");
  }
}

LVector::LVector(const CVector& cv) {
  int ord = floor(0.5*(sqrt(8*cv.size()+1)-3));
  if (PQIndex::size(ord) != cv.size()) 
    throw LaguerreError("CVector size incompatible with LVector");
  // ??? should also check for Hermiticity
  //**/cerr << "order is: " << ord << endl;
  order = new int(ord);
  pcount = new int(1);
  v = new DVector(cv.size(),0.);
  for (PQIndex pq(0,0); !pq.pastOrder(ord); ++pq) {
    if (pq.needsConjugation()) continue;
    this->set(pq, cv[pq.cIndex()]);
    //**/cerr << "pq: " << pq << "\t value: " << cv[pq.cIndex()] 
    //**/     << "\t value: " << (*this)[pq] <<endl;
  }
}

void
LVector::rotate(double theta) {
  DComplex z(cos(theta), -sin(theta));
  DComplex imz(1., 0.);
  for (int m=1; m<=*order; m++) {
    imz *= z;
    for (PQIndex pq(m,0); !pq.pastOrder(*order); pq.incN()) {
      DComplex newb = getElement(v,pq) * imz;
      set(pq, newb);
    }
  }
}

// Build an LCovar from a SqDMatrix for the real degrees of freedom.
// SqDMatrix must have same dimension as needed for LCovar of chosen order
// **Note we are not checking for symmetry of input matrix.
LCovar::LCovar(const SqDMatrix& rhs, int order_): 
  order(new int(order_)), 
  pcount(new int(1)),
  m(new SqDMatrix(rhs)) {
  if (m->size()!=PQIndex::size(*order)) {
      delete m;
      delete pcount;
      delete order;
      throw LaguerreError("Input to LMatrix(SqDMatrix) is wrong size for order");
  }
}


// routines to retrieve and save complex elements of LCovar:
DComplex 
LCovar::operator()(PQIndex pq1, PQIndex pq2) const {
  assert(pq1.pqValid() && !pq1.pastOrder(*order));
  assert(pq2.pqValid() && !pq2.pastOrder(*order));
  int r1index=pq1.rIndex();
  int r2index=pq2.rIndex();
  int i1index=(pq1.isReal()? r1index: r1index+1);
  int i2index=(pq2.isReal()? r2index: r2index+1);

  return DComplex( (*m)(r1index,r2index) 
		   - pq1.iSign()*pq2.iSign()*(*m)(i1index,i2index),
		   pq2.iSign()*(*m)(r1index,i2index)
		   + pq1.iSign()*(*m)(i1index,r2index) );
}

void 
LCovar::set(PQIndex pq1, PQIndex pq2,
	    DComplex covpq1pq2,
	    DComplex covpq1qp2) {
  assert(pq1.pqValid() && !pq1.pastOrder(*order));
  assert(pq2.pqValid() && !pq2.pastOrder(*order));

  if (pq2.needsConjugation()) {
    pq2 = pq2.swapPQ();
    DComplex tmp=covpq1qp2;
    covpq1qp2 = covpq1pq2;
    covpq1pq2 = tmp;
  }
  if (pq1.needsConjugation()) {
    pq1 = pq1.swapPQ();
    DComplex tmp=conj(covpq1qp2);
    covpq1qp2 = conj(covpq1pq2);
    covpq1pq2 = tmp;
  }

  int rIndex1 = pq1.rIndex();
  int rIndex2 = pq2.rIndex();
  int iIndex1 = rIndex1 + 1;
  int iIndex2 = rIndex2 + 1;

  DComplex z=covpq1pq2 + covpq1qp2;
  (*m)(rIndex1, rIndex2) = 0.5*z.real();
  (*m)(rIndex2, rIndex1) = 0.5*z.real();
  if (pq1.isReal()) {
    if (z.imag()!=0.) {
      std::ostringstream oss;
      oss << "Imaginary part " << z.imag()
	  << " in assignment to real LCovar element";
      throw LaguerreError(oss.str());
    }
  } else {
    (*m)(iIndex1, rIndex2) = 0.5*z.imag();
    (*m)(rIndex2, iIndex1) = 0.5*z.imag();
  }

  z=covpq1pq2 - covpq1qp2;
  if (pq2.isReal()) {
    if (z.imag()!=0.) {
      std::ostringstream oss;
      oss << "Imaginary part " << z.imag()
	  << " in assignment to real LCovar element";
      throw LaguerreError(oss.str());
    }
  } else {
    (*m)(rIndex1, iIndex2) = 0.5*z.imag();
    (*m)(iIndex2, rIndex1) = 0.5*z.imag();
  }

  if (pq1.isReal() || pq2.isReal()) {
    if (z.real()!=0.) {
      std::ostringstream oss;
      oss << "Imaginary part " << z.imag()
	  << " in assignment to real LCovar element";
      throw LaguerreError(oss.str());
    }
  } else {
    (*m)(iIndex1, iIndex2) = -0.5*z.real();
    (*m)(iIndex2, iIndex1) = -0.5*z.real();
  }
}

void
LTransform::identity() {
  *m = 0.;
  for (int i=0; i<MIN(m->getM(), m->getN()); i++)
    (*m)(i,i)=1.;
}

//
LTransform::LTransform(const DMatrix& rhs, int orderOut_, int orderIn_):
  orderIn(new int(orderIn_)),
  orderOut(new int(orderOut_)),
  pcount(new int(1)),
  m(new DMatrix(rhs))  {
  if (m->getN()!=PQIndex::size(*orderIn)
      || m->getM()!=PQIndex::size(*orderOut)) {
    /**cerr << "size desired: " << PQIndex::size(*orderIn)
	     << " x " << PQIndex::size(*orderOut)
	     << " size given: " << rhs.getN()
	     << " x " << rhs.getM()
	     << endl; //**/
      delete m;
      delete pcount;
      delete orderIn;
      delete orderOut;
      throw LaguerreError("Input to LTransform(DMatrix) is wrong size for orders");
  }
}

// routines to retrieve and save complex elements of LTransform:
DComplex 
LTransform::operator()(PQIndex pq1, PQIndex pq2) const {
  assert(pq1.pqValid() && !pq1.pastOrder(*orderOut));
  assert(pq2.pqValid() && !pq2.pastOrder(*orderIn));
  int r1index=pq1.rIndex();
  int r2index=pq2.rIndex();
  int i1index=(pq1.isReal()? r1index: r1index+1);
  int i2index=(pq2.isReal()? r2index: r2index+1);

  DComplex z( (*m)(r1index,r2index) 
	      + pq1.iSign()*pq2.iSign()*(*m)(i1index,i2index),
	      pq1.iSign()*(*m)(i1index,r2index)
	      - pq2.iSign()*(*m)(r1index,i2index) );

  return (pq2.isReal()? z : 0.5*z);
}

void 
LTransform::set(PQIndex pq1, PQIndex pq2,
		DComplex Cpq1pq2,
		DComplex Cqp1pq2) {
  assert(pq1.pqValid() && !pq1.pastOrder(*orderOut));
  assert(pq2.pqValid() && !pq2.pastOrder(*orderIn));

  const double RoundoffTolerance=1.e-15;
  DComplex Cpq1qp2;

  if (pq2.needsConjugation()) {
    pq2 = pq2.swapPQ();
    DComplex tmp=conj(Cqp1pq2);
    Cqp1pq2 = conj(Cpq1pq2);
    Cpq1pq2 = tmp;
  }
  if (pq1.needsConjugation()) {
    pq1 = pq1.swapPQ();
    DComplex tmp=Cqp1pq2;
    Cqp1pq2 = Cpq1pq2;
    Cpq1pq2 = tmp;
  }

  int rIndex1 = pq1.rIndex();
  int rIndex2 = pq2.rIndex();
  int iIndex1 = rIndex1+1;
  int iIndex2 = rIndex2+1;

  if (pq1.isReal()) {
    if (Cpq1pq2!=Cqp1pq2) {
      std::ostringstream oss;
      oss << "Invalid LTransform elements for p1=q1, " << Cpq1pq2
	  << " != " << Cqp1pq2;
      throw LaguerreError(oss.str());
    }
    (*m)(rIndex1,rIndex2) = Cpq1pq2.real() * (pq2.isReal()? 1. : 2.);
    if (pq2.isReal()) {
      if (abs(Cpq1pq2.imag()) > RoundoffTolerance) {
	std::ostringstream oss;
	oss << "Nonzero imaginary LTransform elements for p1=q1, p2=q2: " 
	    << Cpq1pq2;
	throw LaguerreError(oss.str());
      }
    } else {
      (*m)(rIndex1,iIndex2) = -2.*Cpq1pq2.imag();
    }
    return;
  } else if (pq2.isReal()) {
    // Here we know p1!=q1:
    if (norm(Cpq1pq2-conj(Cqp1pq2))>RoundoffTolerance) {
      std::ostringstream oss;
      oss << "Inputs to LTransform.set are not conjugate for p2=q2: "
	  << Cpq1pq2 << " vs " << Cqp1pq2 ;
      throw LaguerreError(oss.str());
    }
    (*m)(rIndex1, rIndex2) = Cpq1pq2.real();
    (*m)(iIndex1, rIndex2) = Cpq1pq2.imag();
  } else {
    // Neither pq is real:
    DComplex z=Cpq1pq2 + Cqp1pq2;
    (*m)(rIndex1, rIndex2) = z.real();
    (*m)(rIndex1, iIndex2) = -z.imag();
    z=Cpq1pq2 - Cqp1pq2;
    (*m)(iIndex1, rIndex2) = z.imag();
    (*m)(iIndex1, iIndex2) = z.real();
  }
}

LVector
LTransform::operator*(const LVector rhs) const {
  if (*orderIn != rhs.getOrder()) {
    std::ostringstream oss;
    oss << "Order mismatch between LTransform [" << *orderIn
	<< "] and LVector [" << rhs.getOrder()
	<< "]";
    throw LaguerreError(oss.str());
  }
  DVector out = (*m) * rhs.rVector();
  return LVector(out, *orderOut);
}

LCovar
LTransform::operator*(const LCovar rhs) const {
  if (*orderIn != rhs.getOrder()) {
    std::ostringstream oss;
    oss << "Order mismatch between LTransform [" << *orderIn
	<< "] and LCovar [" << rhs.getOrder()
	<< "]";
    throw LaguerreError(oss.str());
  }
  DMatrix tmp = (*m) * rhs.rMatrix();
  SqDMatrix out(*orderOut,0.);
  for (int i=0; i<*orderOut; i++)
    for (int j=0; j<=i; j++) {
      double sum=0.;
      for (int k=0; k<*orderIn; k++) sum+=tmp(i,k)*(*m)(j,k);
      out(i,j)=sum;
      out(j,i)=sum;
    }
  return LCovar(out, *orderOut);
}

LTransform
LTransform::operator*(const LTransform rhs) const {
  if (*orderIn != rhs.getOrderOut()) {
    std::ostringstream oss;
    oss << "Order mismatch between LTransform [" << *orderIn
	<< "] and LTransform [" << rhs.getOrderOut()
	<< "]";
    throw LaguerreError(oss.str());
  }
  DMatrix out = (*m) * (*rhs.m);
  return LTransform(out, *orderOut, *(rhs.orderIn));
}

LTransform&
LTransform::operator*=(const LTransform rhs) {
  if (*orderIn != rhs.getOrderOut()) {
    std::ostringstream oss;
    oss << "Order mismatch between LTransform [" << *orderIn
	<< "] and LTransform [" << rhs.getOrderOut()
	<< "]";
    throw LaguerreError(oss.str());
  }
  (*m) = (*m) * (*rhs.m);
  *orderIn = *(rhs.orderOut);
  return *this;
}

void
LTransform::preRotation(double theta) {

  int maxM = getOrderIn();

  CVector e_imtheta(maxM+1);  // vector of e^{-im theta}
  const DComplex ONE(1.,0.);
  const DComplex ZERO(0.,0.);
  e_imtheta[0] = ONE;
  DComplex e_itheta(cos(theta), -sin(theta));
  for (int i=1; i<=maxM; i++)
    e_imtheta[i] = e_imtheta[i-1] * e_itheta;

  // Advance through p's
  for (PQIndex pqprime(0,0); 
       !pqprime.pastOrder(maxM);
       pqprime.nextDistinct()) {
    if (pqprime.m()==0) continue;
    DComplex ph= e_imtheta[pqprime.m()];
    
    for (PQIndex pq(0,0); 
	 !pq.pastOrder(getOrderOut());
	 pq.nextDistinct()) {
      DComplex pqpq=(*this)(pq,pqprime);
      DComplex qppq=(*this)(pq.swapPQ(),pqprime);
      set(pq, pqprime, pqpq*ph, qppq*ph);
    }
  }
}

//----------------------------------------------------------------
//----------------------------------------------------------------
// Calculate Laguerre polynomials and wavefunctions:

//----------------------------------------------------------------
// Fill the object's b vector with the values of Laguerre polynomials.
void LVector::fillLaguerre(const double x) {
  double t0, t1, t2, tqmx, tqm;
  int N = *order;
  PQIndex pq(0,0);
  // Do m=0 first:
  t0 = 1.;
  set(pq, t0); //q=0
  pq.incN();
  if (!pq.pastOrder(N)) {
    tqm = 0.;
    t1 = tqmx = 1.-x;	//q=1
    set(pq, t1);
    
    for (pq.incN(); !pq.pastOrder(N); pq.incN()) {
      tqm += 1.;
      tqmx += 2.;
      t2 = (tqmx*t1 - tqm*t0)/pq.getQ();	//q=2,3,...
      set(pq,t2);
      t0 = t1;
      t1 = t2; 
    }
  }
  //proceed with higher m's
  for (int m=1; m<=N; m++) {
    pq.setPQ(m,0);
    t0 = 1.;	//q=0
    set(pq, t0);
    pq.incN();
    if (!pq.pastOrder(N)) {
      tqm = m;
      t1 = tqmx = m+1.-x;	//q=1
      set(pq, t1);
      
      for (pq.incN(); !pq.pastOrder(N); pq.incN()) {
	tqm += 1.;
	tqmx += 2.;
	t2 = (tqmx*t1 - tqm*t0)/pq.getQ();	//q=2,3,...
	set(pq, t2);
	t0 = t1;
	t1 = t2;
      }
    }
  }
}


void 
LVector::fillPsi(double xraw, double yraw,
		  const Ellipse& e) {
  Position<double> xy = e.inv(Position<double>(xraw,yraw));
  fillPsi(xy.x, xy.y, exp(e.getMu()));
}

// **** Note that the following routine, which is in the middle of many pixel
// **** loops, has been rewritten with direct access into the real vector.
// **** So it is dependent upon the storage order NOT changing.
// **** The "old" version below shows what it would be using calls to PQIndex
// **** to generate indices in a more OO way.
void 
LVector::fillPsi(double xunit, double yunit, 
		  double sigma) {
  // fill with psi_pq(z), where psi now defined to have 1/sigma^2 in
  // front.
  double t0, t1, t2, tqmx, tqm;
  static DVector factorials;
  DComplex z(xunit,yunit);

  // First set up the factorials[] vector, which is filled with 
  // (-1)^q sqrt(q!/p!) / 2 pi
  int N = *order;
  if (factorials.size() < PQIndex::size(N)) {
    factorials.resize(PQIndex::size(N));
    PQIndex pq;
    double accum;
    for (int p=0; p<=N; p++) {
      pq.setPQ(p,0);
      accum = 1./(2.*PI*sqrtfact(p));
      int iR=pq.rIndex();
      if (p!=0) factorials[iR++] = accum;
      factorials[iR] = accum;
      pq.incQ();
      for (int q=1; !pq.pastOrder(N) && p>=q; pq.incQ(), q++) {
	accum *= -sqrtn(q);
	int iR=pq.rIndex();
	if (p!=q) factorials[iR++] = accum;
	factorials[iR] = accum;
      }
    }
  }

  PQIndex pq(0,0);
  DComplex zm=1.;	//z^m
  double x = norm(z);
  double prefactor = exp(-0.5*x) / (sigma*sigma);
  // Do m=0 first:
  t0 = prefactor;	//q=0
  int iR = pq.rIndex();
  double* vptr = &((*v)[iR]);
  *vptr = t0;

  int iN=2;
  if (iN<=N) {
    tqm = 0.;
    tqmx = 1.-x;	//q=1
    t1 = prefactor*tqmx;
    vptr += 2*iN+1;
    *vptr = t1;

    int iQ=2;
    for (iN+=2; iN<=N; iN+=2, iQ++) {
      tqm += 1.;
      tqmx += 2.;
      t2 = (tqmx*t1 - tqm*t0)/iQ;	//q=2,3,...
      vptr += 2*iN+1;
      *vptr = t2;
      t0 = t1;
      t1 = t2; 
    }
  }
  //proceed with higher m's
  for (int m=1; m<=N; m++) {
    zm *= z;
    double zr=zm.real();
    double zi=zm.imag();
    pq.setPQ(m,0);
    t0 = prefactor;	//q=0
    int iR = pq.rIndex();
    vptr = &((*v)[iR]);
    *(vptr++) = t0*zr;
    *vptr = t0*zi;
    int iN = pq.N()+2;
    if (iN<=N) {
      tqm = m;
      tqmx = m+1.-x;	//q=1
      t1 = prefactor * tqmx;
      vptr += 2*iN;
      *(vptr++) = t1*zr;
      *vptr = t1*zi;
      
      int iQ = 2;
      for (iN+=2; iN<=N; iN+=2, iQ++) {
	tqm += 1.;
	tqmx += 2.;
	t2 = (tqmx*t1 - tqm*t0)/iQ;	//q=2,3,...
	vptr += 2*iN;
	*(vptr++) = t2*zr;
	*vptr = t2*zi;
	t0 = t1;
	t1 = t2;
      }
    }
  }
  // Now multiply in the factorial array
  vptr = &((*v)[0]);
  double* fptr = &(factorials[0]);
  for (int i=0; i<v->size(); i++)
    *(vptr++) *= *(fptr++);
}

/* OLD version works but is slower:
void 
LVector::fillPsi(double xunit, double yunit, 
		  double sigma) {
  // fill with psi_pq(z), where psi now defined to have 1/sigma^2 in
  // front.
  double t0, t1, t2, tqmx, tqm;
  static DVector factorials;
  DComplex z(xunit,yunit);

  // First set up the factorials[] vector, which is filled with 
  // (-1)^q sqrt(q!/p!) / 2 pi
  int N = *order;
  if (factorials.size() < PQIndex::size(N)) {
    factorials.resize(PQIndex::size(N));
    PQIndex pq;
    double accum;
    for (int p=0; p<=N; p++) {
      pq.setPQ(p,0);
      accum = 1./(2.*PI*sqrtfact(p));
      factorials[pq.cIndex()] = accum;
      pq.incQ();
      for (int q=1; !pq.pastOrder(N) && p>=q; pq.incQ(), q++) {
	accum *= -sqrtn(q);
	factorials[pq.cIndex()] = accum;
      }
    }
  }

  PQIndex pq(0,0);
  DComplex zm=1.;	//z^m
  double x = norm(z);
  double prefactor = exp(-0.5*x) / (sigma*sigma);
  // Do m=0 first:
  t0 = prefactor;	//q=0
  set(pq, t0* factorials[pq.cIndex()]);	//q=0
  pq.incN();
  if (!pq.pastOrder(N)) {
    tqm = 0.;
    tqmx = 1.-x;	//q=1
    t1 = prefactor*tqmx;
    set(pq, t1 * factorials[pq.cIndex()]);
    
    for (pq.incN(); !pq.pastOrder(N); pq.incN()) {
      tqm += 1.;
      tqmx += 2.;
      t2 = (tqmx*t1 - tqm*t0)/pq.getQ();	//q=2,3,...
      set(pq, t2 * factorials[pq.cIndex()]);
      t0 = t1;
      t1 = t2; 
    }
  }
  //proceed with higher m's
  for (int m=1; m<=N; m++) {
    zm *= z;
    pq.setPQ(m,0);
    t0 = prefactor;	//q=0
    set(pq, t0*factorials[pq.cIndex()]*zm);
    pq.incN();
    if (!pq.pastOrder(N)) {
      tqm = m;
      tqmx = m+1.-x;	//q=1
      t1 = prefactor * tqmx;
      set(pq, t1*factorials[pq.cIndex()]*zm);
      
      for (pq.incN(); !pq.pastOrder(N); pq.incN()) {
	tqm += 1.;
	tqmx += 2.;
	t2 = (tqmx*t1 - tqm*t0)/pq.getQ();	//q=2,3,...
	set(pq, t2*factorials[pq.cIndex()]*zm);
	t0 = t1;
	t1 = t2;
      }
    }
  }
}
*/
DVector
LVector::dotVector(0);
int
LVector::dotVectorOrder=-1;

void
LVector::checkDotVector(int order) {
  if (order<=dotVectorOrder) return;
  dotVector.resize(PQIndex::size(order));
  dotVectorOrder = order;
  for (PQIndex pq(0,0); !pq.pastOrder(order); pq.nextDistinct()) {
    int ri = pq.rIndex();
    if (pq.isReal()) dotVector[ri] = 1.;
    else {
      dotVector[ri]=2.;
      dotVector[ri+1]=-2.;
    }
  }
}

double 
LVector::dot(const LVector& rhs) const {
  double sum=0.;
  checkDotVector(*order);
  for (int j=0; j<v->size(); j++)
    sum += (*v)[j]* (*(rhs.v))[j] * dotVector[j];
  return sum;
}

DVector
LVector::realRHS() const {
  checkDotVector(*order);
  DVector X(rVector());	//new copy of the psi array
  for (int j=0; j<v->size(); j++)
    // ??? speed this up!
    X[j] *= dotVector[j];
  return X;
}

DVector
LVector::realPsi(int maxN, double x, double y,
			   const Ellipse& e) {
  int j;
  LVector lv(maxN);
  lv.fillPsi(x, y, e);
  return lv.realRHS();
}

void 
LVector::realPsi(int N, const DVector& xx, const DVector& yy, 
			   const Ellipse& e, 
			   DMatrix* psi) {
  // This is a dumb implementation:
  assert(xx.size()==yy.size());
  psi->resize(PQIndex::size(N), xx.size());
  LVector lv(N);
  for (int j=0; j<xx.size(); j++) {
    lv.fillPsi(xx[j], yy[j], e);
    DVector xx(lv.realRHS());
    for (int i=0; i<xx.size(); i++)
      (*psi)(i,j) = xx[i];
  }
}

void
LVector::fillKPsi(double xunit, double yunit) {
  this->fillPsi(xunit, yunit, 1.);  // Fourier[Psi_pq] is unitless
  // rotate kvalues of Psi with (-i)^(p+q)
  for (PQIndex pq(0,0); !pq.pastOrder(*order); ++pq) {
    if (pq.needsConjugation()) continue;
    switch (pq.N() % 4) {
    case 0: this->set(pq, (*this)[pq]);  break;
    case 1: this->set(pq, (*this)[pq]*DComplex(0,-1));  break;
    case 2: this->set(pq, (*this)[pq]*DComplex(-1,0));  break;
    case 3: this->set(pq, (*this)[pq]*DComplex(0,1));  break;
    }
  }
}

CVector*
LVector::cVector() const {
  CVector *cv = new CVector(PQIndex::size(*order), 0.);
  for (PQIndex pq(0,0); !pq.pastOrder(*order); ++pq) {
    (*cv)[pq.cIndex()] = (*this)[pq];
  }
  return cv;
}

CVector*
LVector::cVectorConjugate() const {
  CVector *cv = new CVector(PQIndex::size(*order), 0.);
  DComplex c;
  for (PQIndex pq(0,0); !pq.pastOrder(*order); ++pq) {
    (*cv)[pq.cIndex()] = conj((*this)[pq]);
  }
  return cv;
}

CVector*
LVector::kSpaceCVector() const {
  CVector *cv = new CVector(PQIndex::size(*order), 0.);
  for (PQIndex pq(0,0); !pq.pastOrder(*order); ++pq) {
    if (!pq.needsConjugation())
      (*cv)[pq.cIndex()] = (*this)[pq];
    else
      if (pq.N()%2 == 0)
	(*cv)[pq.cIndex()] = (*this)[pq];
      else
	(*cv)[pq.cIndex()] = -(*this)[pq];
  }
  return cv;
}

CVector*
LVector::kSpaceCVectorConjugate() const {
  CVector *cv = new CVector(PQIndex::size(*order), 0.);
  DComplex c;
  for (PQIndex pq(0,0); !pq.pastOrder(*order); ++pq) {
    if (!pq.needsConjugation())
      (*cv)[pq.cIndex()] = conj((*this)[pq]);
    else 
      if (pq.N()%2 == 0)
	(*cv)[pq.cIndex()] = conj((*this)[pq]);
      else
	(*cv)[pq.cIndex()] = -conj((*this)[pq]);
  }
  return cv;
}



//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Flux determinations
double 
LVector::flux(int maxP) const {
  if (maxP<0) maxP = getOrder()/2;
  if (maxP > getOrder()/2) maxP=getOrder()/2;
  double retval=0.;
  for (int p=0; p<=maxP; p++)
    retval += getElement(v,PQIndex(p,p)).real();
  return retval;
}
double 
LVector::apertureFlux(double R_, int maxP) const {
  static DVector fp;
  static double R=-1.;
  static double psize=-1;

  assert(R_>=0.);

  if (maxP<0) maxP= getOrder()/2;
  if (maxP > getOrder()/2) maxP=getOrder()/2;

  if (R_ != R || maxP>psize) {
    if (maxP != psize) {
      fp.resize(maxP+1);
      psize = maxP;
    }
    R = R_;
    DVector Lp(maxP+1);
    DVector Qp(maxP+1);
    double x = R*R;
    double efact = exp(-0.5*x);
    Lp[0] = Qp[0]=1.;
    if (maxP>0) {
      Lp[1] = 1. - x;
      Qp[1] = -1. - x;
    } 
    for (int p=1; p<maxP; p++) {
      Lp[p+1] = ((2*p+1-x)*Lp[p]-p*Lp[p-1])/(p+1);
      Qp[p+1] = (-x*Lp[p]-Qp[p]+p*Qp[p-1])/(p+1);
    }
    for (int p=0; p<=maxP; p++)
      fp[p] = 1. - efact*Qp[p]*(p%2==0 ? 1. : -1.);
  }

  double flux = 0.;
  for (int p=0; p<=maxP; p++)
    flux += getElement(v,PQIndex(p,p)).real() * fp[p];
  return flux;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// I/O Routines

ostream& operator<<(ostream& os, const LVector& lv) {
  lv.write(os); return os;
}

void LVector::write(ostream& os, int maxorder) const
{
  int oldprec = os.precision(8);
  std::ios::fmtflags oldf = os.setf(std::ios::scientific,std::ios::floatfield);
  if (maxorder < 0 || maxorder > *order)
    maxorder = *order;
  os << *order << endl;
  for (int n=0; n<=maxorder; n++) 
    for(PQIndex pq(n,0); !pq.needsConjugation(); pq.decm()) {
      os << " " << std::setw(2) << pq.getP() 
	 << " " << std::setw(2) << pq.getQ() ;
      if (pq.isReal()) 
	os << " " << std::setw(15) << (*this)[pq].real() << endl;
      else
	os << " " << std::setw(15) << (*this)[pq].real() 
	   << " " << std::setw(15) << (*this)[pq].imag() << endl;
    }
  os.precision(oldprec);
  os.flags(oldf);
}

istream& operator>>(istream& is, LVector& lv) {
  lv.read(is); return is;
}

void LVector::read(istream& is)
{
  // get order
  int order;
  is >> order;
  resize(order);
  // discard p,q info, read into rVector
  int p, q; 
  double re, im;
  for (int n=0; n<=order; n++) 
    for(PQIndex pq(n,0); !pq.needsConjugation(); pq.decm()) {
      is >> p >> q;
      if (pq.isReal()) {
	is >> re;  im = 0.;
	setElement(v, pq, DComplex(re,im));
      } else {
	is >> re >> im;
	setElement(v, pq, DComplex(re,im));
      }
    }
}

void PQIndex::write(ostream& os) const
{
  os << std::setw(2) << getP() 
     << "," << std::setw(2) << getQ() ;
}

// Transformation generators - these are static quantities:
const DMatrix& 
LVector::Generator(int iparam, int order) {
  static DMatrix gmu;
  static DMatrix gx;
  static DMatrix gy;
  static DMatrix ge1;
  static DMatrix ge2;

  if (iparam==iMu) {
    if (gmu.getM()<PQIndex::size(order)) {
      LTransform lt(order, order+2);

      for (PQIndex pq(0,0); !pq.pastOrder(order); pq.nextDistinct()) {
	int p=pq.getP();
	int q=pq.getQ();
	DComplex zz(-1.,0.);
	if (pq.isReal())
	  lt.set(pq,pq,zz, zz);
	else
	  lt.set(pq,pq,zz, 0.);
	PQIndex pqprime(p+1, q+1);
	zz = DComplex(sqrtn(p+1)*sqrtn(q+1), 0.);
	if (pq.isReal())
	  lt.set(pq,pqprime,zz, zz);
	else
	  lt.set(pq,pqprime,zz, 0.);
	
	if (q>0) {
	  pqprime.setPQ(p-1,q-1);
	  zz = DComplex(-sqrtn(p)*sqrtn(q), 0.);
	  if (pq.isReal())
	    lt.set(pq,pqprime,zz, zz);
	  else
	    lt.set(pq,pqprime,zz, 0.);
	}
      }
      gmu = lt.rMatrix();
      LVector::checkDotVector(order+2);
      for (int i=0; i<gmu.getM(); i++)
	for (int j=0; j<gmu.getN(); j++)
	  gmu(i,j) *= LVector::dotVector[i]/LVector::dotVector[j];
    }
    return gmu;
  }
  if (iparam==iX) {
    if (gx.getM()<PQIndex::size(order)) {
      LTransform lt(order, order+2);

      for (PQIndex pq(0,0); !pq.pastOrder(order); pq.nextDistinct()) {
	int p=pq.getP();
	int q=pq.getQ();
	PQIndex pqprime(p+1, q);
	DComplex zz(0.5*sqrtn(p+1),0.);
	if (pq.isReal()) {
	  lt.set(pq,pqprime,zz, zz);
	  if (p>0) {
	    zz = DComplex(-0.5*sqrtn(p), 0.);
	    pqprime.setPQ(p-1,q);
	    lt.set(pq,pqprime,zz, zz);
	  }
	} else {
	  lt.set(pq,pqprime,zz, 0.);
	  pqprime.setPQ(p, q+1);
	  zz = DComplex(0.5*sqrtn(q+1),0.);
	  if (pq.m()==1) {
	    lt.set(pq,pqprime, zz, zz);
	  } else {
	    lt.set(pq,pqprime, zz, 0.);
	  }
	  pqprime.setPQ(p-1,q);
	  zz = DComplex(-0.5*sqrtn(p), 0.);
	  if (pq.m()==1) {
	    lt.set(pq,pqprime, zz, zz);
	  } else {
	    lt.set(pq,pqprime, zz, 0.);
	  }
	  if (q>0) {
	    pqprime.setPQ(p,q-1);
	    zz = DComplex(-0.5*sqrtn(q), 0.);
	    lt.set(pq,pqprime, zz, 0.);
	  }
	}
      }
      gx = lt.rMatrix();
      LVector::checkDotVector(order+2);
      for (int i=0; i<gx.getM(); i++)
	for (int j=0; j<gx.getN(); j++)
	  gx(i,j) *= LVector::dotVector[i]/LVector::dotVector[j];
    }
    return gx;
  }

  if (iparam==iY) {
    if (gy.getM()<PQIndex::size(order)) {
      LTransform lt(order, order+2);

      for (PQIndex pq(0,0); !pq.pastOrder(order); pq.nextDistinct()) {
	int p=pq.getP();
	int q=pq.getQ();
	PQIndex pqprime(p+1, q);
	DComplex zz(0.,-0.5*sqrtn(p+1));
	if (pq.isReal()) {
	  lt.set(pq,pqprime,zz, zz);
	  if (p>0) {
	    zz = DComplex(0.,-0.5*sqrtn(p));
	    pqprime.setPQ(p-1,q);
	    lt.set(pq,pqprime,zz, zz);
	  }
	} else {
	  lt.set(pq,pqprime,zz, 0.);
	  pqprime.setPQ(p, q+1);
	  zz = DComplex(0.,0.5*sqrtn(q+1));
	  if (pq.m()==1) {
	    lt.set(pq,pqprime, zz, conj(zz));
	  } else {
	    lt.set(pq,pqprime, zz, 0.);
	  }
	  pqprime.setPQ(p-1,q);
	  zz = DComplex(0.,-0.5*sqrtn(p));
	  if (pq.m()==1) {
	    lt.set(pq,pqprime, zz, conj(zz));
	  } else {
	    lt.set(pq,pqprime, zz, 0.);
	  }
	  if (q>0) {
	    pqprime.setPQ(p,q-1);
	    zz = DComplex(0.,0.5*sqrtn(q));
	    lt.set(pq,pqprime, zz, 0.);
	  }
	}
      }
      gy = lt.rMatrix();
      LVector::checkDotVector(order+2);
      for (int i=0; i<gy.getM(); i++)
	for (int j=0; j<gy.getN(); j++)
	  gy(i,j) *= LVector::dotVector[i]/LVector::dotVector[j];
    }
    return gy;
  }

  if (iparam==iE1) {
    if (ge1.getM()<PQIndex::size(order)) {
      LTransform lt(order, order+2);

      for (PQIndex pq(0,0); !pq.pastOrder(order); pq.nextDistinct()) {
	int p=pq.getP();
	int q=pq.getQ();
	PQIndex pqprime(p+2, q);
	DComplex zz(0.25*sqrtn(p+1)*sqrtn(p+2),0.);
	if (pq.isReal()) {
	  lt.set(pq,pqprime,zz, zz);
	  if (p>1) {
	    zz = DComplex(-0.25*sqrtn(p)*sqrtn(p-1),0.);
	    pqprime.setPQ(p-2,q);
	    lt.set(pq,pqprime,zz, zz);
	  }
	} else {
	  lt.set(pq,pqprime,zz, 0.);
	  pqprime.setPQ(p, q+2);
	  zz = DComplex(0.25*sqrtn(q+1)*sqrtn(q+2),0.);
	  if (pq.m()==2) {
	    lt.set(pq,pqprime, zz, zz);
	  } else {
	    lt.set(pq,pqprime, zz, 0.);
	  }
	  if (p>1) {
	    pqprime.setPQ(p-2,q);
	    zz = DComplex(-0.25*sqrtn(p)*sqrtn(p-1),0.);
	    if (pq.m()==2) {
	      lt.set(pq,pqprime, zz, zz);
	    } else {
	      lt.set(pq,pqprime, zz, 0.);
	    }
	    if (q>1) {
	      pqprime.setPQ(p,q-2);
	      zz = DComplex(-0.25*sqrtn(q)*sqrtn(q-1),0.);
	      lt.set(pq,pqprime, zz, 0.);
	    }
	  }
	}
      }
      ge1 = lt.rMatrix();
      LVector::checkDotVector(order+2);
      for (int i=0; i<ge1.getM(); i++)
	for (int j=0; j<ge1.getN(); j++)
	  ge1(i,j) *= LVector::dotVector[i]/LVector::dotVector[j];
    }
    return ge1;
  }

  if (iparam==iE2) {
    if (ge2.getM()<PQIndex::size(order)) {
      LTransform lt(order, order+2);

      for (PQIndex pq(0,0); !pq.pastOrder(order); pq.nextDistinct()) {
	int p=pq.getP();
	int q=pq.getQ();
	PQIndex pqprime(p+2, q);
	DComplex zz(0., -0.25*sqrtn(p+1)*sqrtn(p+2));
	if (pq.isReal()) {
	  lt.set(pq,pqprime,zz, zz);
	  if (p>1) {
	    zz = DComplex(0.,-0.25*sqrtn(p)*sqrtn(p-1));
	    pqprime.setPQ(p-2,q);
	    lt.set(pq,pqprime,zz, zz);
	  }
	} else {
	  lt.set(pq,pqprime,zz, 0.);
	  pqprime.setPQ(p, q+2);
	  zz = DComplex(0.,0.25*sqrtn(q+1)*sqrtn(q+2));
	  if (pq.m()==2) {
	    lt.set(pq,pqprime, zz, conj(zz));
	  } else {
	    lt.set(pq,pqprime, zz, 0.);
	  }
	  if (p>1) {
	    pqprime.setPQ(p-2,q);
	    zz = DComplex(0.,-0.25*sqrtn(p)*sqrtn(p-1));
	    if (pq.m()==2) {
	      lt.set(pq,pqprime, zz, conj(zz));
	    } else {
	      lt.set(pq,pqprime, zz, 0.);
	    }
	    if (q>1) {
	      pqprime.setPQ(p,q-2);
	      zz = DComplex(0.,0.25*sqrtn(q)*sqrtn(q-1));
	      lt.set(pq,pqprime, zz, 0.);
	    }
	  }
	}
      }
      ge2 = lt.rMatrix();
      checkDotVector(order+2);
      for (int i=0; i<ge2.getM(); i++)
	for (int j=0; j<ge2.getN(); j++)
	  ge2(i,j) *= dotVector[i]/dotVector[j];
    }
    return ge2;
  }
  else
    throw (LaguerreError("Unknown parameter for LVector::Generator()"));
}

// Function to solve for radius enclosing a specified flux.
// Return negative radius if no root is apparent.
class FRSolve {
public:
  FRSolve(const LVector& lv_, double thresh_, int maxP_): 
    lv(lv_), maxP(maxP_), thresh(thresh_) {
    assert(lv.getOrder() >= 2*maxP);
  }
  double operator()(double u) const {return lv.apertureFlux(u,maxP)-thresh;}
private:
  const LVector& lv;
  double thresh;
  int maxP;
};

double
fluxRadius(const LVector& lv, double threshold, int maxP) {
  if (maxP<0) maxP= lv.getOrder()/2;
  if (maxP > lv.getOrder()/2) maxP=lv.getOrder()/2;
  FRSolve func(lv, threshold, maxP);

  // First we step through manually at intervals roughly the smallest that 
  // a function of this order can oscillate, in order to bracket the root
  // closest to the origin.

  const double TOLERANCE=0.001;	//radius accuracy required
  const double maxR = 5.;
  double ustep=0.5/sqrt(static_cast<double> (maxP+1));
  double u1 = 0.0001;
  double f1 = func(u1);
  double u2;
  while (u1<maxR) {
    u2 = u1 + ustep;
    double f2 = func(u2);
    if ( f1*f2<=0.) break;
    u1 = u2;
    f1 = f2;
  }
  if (u1>=maxR) {
    u2 = 2*maxR;
    double f2 = func(u2);
    if ( f1*f2>0.) {
      // At this point there appears to be no root.
      return -1.;
    }
  }

  // Now a bisection solution for the root
  solve::Solve<FRSolve> s(func, u1, u2);
  s.setXTolerance(TOLERANCE);
  return s.root();
}

}
