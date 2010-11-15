// $Id: LTransforms.cpp,v 1.10 2007-03-14 15:52:38 reiko3 Exp $
#include "Array3d.h"
#include "Laguerre.h"
#include "BinomFact.h"
#include <math.h>
#include <assert.h>

namespace laguerre {

//-------------------------------------------------------------------
// TRANSFORMATION MATRICES: Shear
//-------------------------------------------------------------------

// Make the F matrix for shear.  
void ShearFMatrix(CMatrix* &Fptr,
		  const Shear eta, 
		  int orderOut, int orderIn) {

  double secheta = 1./cosh(0.5*eta.getEta());
  
  double twobeta = eta.getBeta()*2.;
  DComplex phase(cos(twobeta),-sin(twobeta));
  DComplex conjphase = conj(phase);

  int maxOrder = MAX(orderOut, orderIn);
  int sumOrder = orderOut+orderIn;

  if (!Fptr) Fptr = new CMatrix(orderOut+1, orderIn+1);
  else Fptr->resize(orderOut+1, orderIn+1);
  *Fptr = DComplex(0.,0.);

  static valarray<double> sqrtp;	//commonly used vector of sqrt(p)
  if (sqrtp.size() < sumOrder+1) {
    sqrtp.resize(sumOrder+1);
    for (int i=0; i<=sumOrder; i++) sqrtp[i]=sqrtn(i);
  }
  // also static work vectors needed
  static valarray<complex<double> > *gp, *gpm1, *gpm2, *gtmp;
  if (gp==0) gp = new valarray<complex<double> >(sumOrder+1);
  if (gpm1==0) gpm1 = new valarray<complex<double> >(sumOrder+1);
  if (gpm2==0) gpm2 = new valarray<complex<double> >(sumOrder+1);
  gp->resize(sumOrder+1);
  gpm1->resize(sumOrder+1);
  gpm2->resize(sumOrder+1);
  (*gp) = DComplex(0.,0.);
  (*gpm1) = DComplex(0.,0.);
  (*gpm2) = DComplex(0.,0.);

  // An array for powers of tanh
  static valarray<double> tanh_powers;
  if (tanh_powers.size() < maxOrder/2+1) 
    tanh_powers.resize(maxOrder/2+1);
  {
    double accum=1.;
    double tanheta = tanh(0.5*eta.getEta());
    for (int i=0; i<tanh_powers.size(); i++) {
      tanh_powers[i] = accum;
      accum *= tanheta;
    }
  }

  // Begin the regression process.  Element p of the vectors 
  // gp, gpm1, gpm2 will hold the quantities
  // g(p',p), g(p'-1,p), and g(p'-2,p)
  // defined in notes of 1/2/04.  These need to be computed up to
  // p+p' <= sumOrder
  // but only for p'<=p.

  // The matrix element F(p', p) will be defined as
  // F(p', p) = g(p',p) * tanh(eta/2)^(p-p')/2
  // With this definition, the transformation matrix is
  // C(p'q' pq) = F(p',p) * conj(F(q',q))

  // Create (0,p) entries - only even ones are non-zero
  DComplex accum = sqrt(secheta);

  (*gp)[0] = accum;
  for (int p2=1; 2*p2<=sumOrder; p2++) {
    accum *= -sqrtp[2*p2-1]*conjphase / sqrtp[2*p2];
    (*gp)[2*p2] = accum;
  }
  // Fill in the F matrix elements
  for (int p2=0; 2*p2<=orderIn; p2++) {
    (*Fptr)(0,2*p2) = (*gp)[2*p2] * tanh_powers[p2];
  }
  for (int pp2=2; 2*pp2<=orderOut; pp2+=2) {
    (*Fptr)(2*pp2,0) = conj((*gp)[2*pp2]) * tanh_powers[pp2];
  }
  for (int pp2=1; 2*pp2<=orderOut; pp2+=2) {
    (*Fptr)(2*pp2,0) = -conj((*gp)[2*pp2]) * tanh_powers[pp2];
  }

  // Use recursion formula to generate higher p'
  DComplex pm1fac, pm2fac;
  DComplex ztmp;
  for (int pprime=1; pprime<=maxOrder; pprime++) {
    // cycle the work vectors:
    gtmp = gpm2; gpm2=gpm1; gpm1=gp; gp=gtmp;
    (*gp)=DComplex(0.,0.);

    pm1fac = -phase*secheta / sqrtp[pprime];
    pm2fac = phase* sqrtp[pprime-1] / sqrtp[pprime];

    for (int p=pprime; p<=sumOrder-pprime; p+=2) {
      (*gp)[p] = pm1fac*sqrtp[p+1]*(*gpm1)[p+1] + pm2fac * (*gpm2)[p];
    }

    // Fill in F matrix elements, first (p',p)
    int dp2, p;
    if (pprime <= orderOut)
      for (dp2=0, p=pprime; p<=orderIn; dp2++, p+=2) 
	(*Fptr)(pprime,p) = (*gp)[p] * tanh_powers[dp2];
    // Then (p,p'); note sign flip when (p-prime)=2 (mod 4):
    if (pprime <= orderIn) {
      for (dp2=2, p=pprime+4; p<=orderOut; dp2+=2, p+=4) 
	//	(*Fptr)(p,pprime) = conj((*gp)[p]) * tanh_powers[dp2];
	(*Fptr)(p,pprime) = conj((*gp)[p]) * tanh_powers[dp2];
      for (dp2=1, p=pprime+2; p<=orderOut; dp2+=2, p+=4) 
	(*Fptr)(p,pprime) = -conj((*gp)[p]) * tanh_powers[dp2];
    }
  }
}

// Make full matrix for a shear
LTransform MakeLTransform(Shear eta, 
			  int orderOut, int orderIn,
			  bool shiftCoords) {

  assert(orderOut>=0);
  assert(orderIn>=0);

  static CMatrix* F=0;
  static int oOut=-1;
  static int oIn=-1;
  static Shear oldEta;

  if (shiftCoords) eta = -eta;  //coord shift needs inverse 

  // Build the F matrix for this order & shear
  if ( (!F) 
       || oOut<orderOut
       || oIn<orderIn
       || !(oldEta==eta) ) {
    ShearFMatrix(F,eta, orderOut, orderIn);
    oOut = orderOut;
    oIn  = orderIn;
    oldEta=eta;
  }

  // Create the output matrix
  LTransform result(orderOut, orderIn);

  // And fill it:
  for (PQIndex pqi; pqi.N()<=orderIn; pqi.nextDistinct()) {
    for (PQIndex pqo; pqo.N()<=orderOut; pqo.nextDistinct()) {
      // Zero matrix elements when delta m is odd:
      if ( (pqi.m()+pqo.m())%2 ) continue;
      result.set(pqo, pqi,
		 (*F)(pqo.getP(),pqi.getP()) * 
		 conj((*F)(pqo.getQ(),pqi.getQ())),
		 (*F)(pqo.getQ(),pqi.getP()) * 
		 conj((*F)(pqo.getP(),pqi.getQ())));
    }
  }

  return result;
}

//-------------------------------------------------------------------
// TRANSFORMATION MATRICES: Translation
//-------------------------------------------------------------------

void
TranslationFMatrix(CMatrix* &Fptr,
		   const Position<double> x0, 
		   int orderOut, int orderIn) {

  int maxOrder = MAX(orderOut, orderIn);
  if (!Fptr) Fptr = new CMatrix(orderOut+1, orderIn+1);
  else Fptr->resize(orderOut+1, orderIn+1);
  *Fptr = DComplex(0.,0.);

  DComplex z(x0.x,x0.y);
  DComplex z2=0.5*z;
  DComplex zbar2 = 0.5*conj(z);

  static valarray<double> sqrtp;	//commonly used vector of sqrt(p)
  if (sqrtp.size() < maxOrder+1) {
    sqrtp.resize(maxOrder+1);
    for (int i=0; i<=maxOrder; i++) sqrtp[i]=sqrtn(i); 
  }

  // Build first row of f
  (*Fptr)(0,0) = exp(-z*conj(z)/8.);
  for (int p=1; p<=orderIn; p++) 
    (*Fptr)(0,p) = -(*Fptr)(0,p-1)*z2/sqrtp[p];

  // Recursion for pprime>0 rows
  for (int pprime=1; pprime<=orderOut; pprime++) {
    (*Fptr)(pprime,0) = zbar2*(*Fptr)(pprime-1,0)/sqrtp[pprime];
    for (int p=1; p<=orderIn; p++)
      (*Fptr)(pprime,p) = ( (*Fptr)(pprime-1,p-1)*sqrtp[p] 
			    + zbar2* (*Fptr)(pprime-1,p) )
	/ sqrtp[pprime];
  }
}

LTransform 
MakeLTransform(Position<double> x0, 
	       int orderOut, int orderIn,
		   bool shiftCoords) 
{
  assert(orderOut>=0);
  assert(orderIn>=0);

  static CMatrix* F=0;
  static int oOut=-1;
  static int oIn=-1;
  static Position<double> oldX0;

  if (shiftCoords) {
    //invert xfrm for coordinate shift
    x0.x = -x0.x; 
    x0.y = -x0.y;
  }
  // Build the F matrix for this order & shear
  if ( (!F) 
       || oOut<orderOut
       || oIn<orderIn
       || oldX0.x!=x0.x
       || oldX0.y!=x0.y ) {
    TranslationFMatrix(F, x0, orderOut, orderIn);
    oOut = orderOut;
    oIn  = orderIn;
    oldX0= x0;
  }

  // Create the output matrix
  LTransform result(orderOut, orderIn);

  // And fill it:
  for (PQIndex pqi; pqi.N()<=orderIn; pqi.nextDistinct()) {
    for (PQIndex pqo; pqo.N()<=orderOut; pqo.nextDistinct()) {
      result.set(pqo, pqi,
		 (*F)(pqo.getP(), pqi.getP()) 
		   * conj((*F)(pqo.getQ(), pqi.getQ())),
		 (*F)(pqo.getQ(), pqi.getP()) 
		   * conj((*F)(pqo.getP(), pqi.getQ())));
    }
  }

  return result;
}


//-------------------------------------------------------------------
// TRANSFORMATION MATRICES: Dilation
//-------------------------------------------------------------------

LTransform
MakeLTransform(double mu,
	       int orderOut, int orderIn,
	       bool coordShift) {

  assert(orderOut>=0);
  assert(orderIn>=0);

  int sumOrder = orderOut + orderIn;
  int maxQp = sumOrder/2;

  LTransform D(orderOut, orderIn);

  if (coordShift) mu=-mu;	//invert for coord transformation
  double tanhmu = tanh(mu);
  double cmu = cosh(mu);
  double smu = sinh(mu);

  int p,q,pprime,qprime;

  // matrices to hold m, m-1 submatrices
  static DVector *rowp0=0, *rowpm10=0, *rowpq=0, *rowpqm1=0, *rowtemp=0;

  // The temporary vectors for recursion are
  // rowp0[q']    = C(p'q', p   0), p'= p+q'
  // rowpm10[q']  = C(p'q', p-1 0), p'= p+q'-1
  // rowpq[q']    = C(p'q', p   q), p'= p+q'-q
  // rowpqm1[q']  = C(p'q', p q-1), p'= p+q'+1

  if (rowp0==0) rowp0 = new DVector(maxQp+1);
  if (rowpm10==0) rowpm10 = new DVector(maxQp+1);
  if (rowpq==0) rowpq = new DVector(maxQp+1);
  if (rowpqm1==0) rowpqm1 = new DVector(maxQp+1);
  rowp0->Resize(maxQp+1);
  rowpm10->Resize(maxQp+1);
  rowpq->Resize(maxQp+1);
  rowpqm1->Resize(maxQp+1);

  static valarray<double> sqrtp;	//commonly used vector of sqrt(p)
  if (sqrtp.size() < sumOrder+1) {
    sqrtp.resize(sumOrder+1);
    for (int i=0; i<=sumOrder; i++) sqrtp[i]=sqrtn(i); 
  }

  // Start with m=0, p=0;
  double factor = tanh(mu);
  // This is the flux-conserving prefactor, if psi's are defined to have
  // unit flux in any basis.
  double accum=exp(-mu)/cmu;	
  for (qprime = 0; qprime<=maxQp; qprime++) {
    (*rowp0)[qprime] = accum;
    if (2*qprime<=orderOut) 
      D.set(PQIndex(qprime,qprime), PQIndex(0,0),
	    DComplex(accum), DComplex(accum));
    accum *= factor;
  }

  // Advance through p's
  for (p=1; p<=orderIn; p++) {
    // cycle buffers:
    rowtemp = rowpm10; rowpm10 = rowp0; rowp0 = rowtemp;
    // generate D[pprime qprime][p0] vector.  Recursion is
    // sqrt(p) D[p'q'][pq] = cosh(mu) sqrt(p') D[p'-1 q'][p-1 q]
    //		-sinh(mu) sqrt(q'+1) D[p' q'+1][p-1 q]
    double factor1=cmu/sqrtp[p];
    double factor2=-smu/sqrtp[p];
    for (qprime=0, pprime=p; 
	 qprime+pprime<=sumOrder-p; 
	 qprime++, pprime++) {
      (*rowp0)[qprime] = factor1 * sqrtp[pprime] * (*rowpm10)[qprime]
	+ factor2* sqrtp[qprime+1] * (*rowpm10)[qprime+1];
      // Fill spots in the matrix:
      if (pprime+qprime<=orderOut) {
	D.set(PQIndex(pprime, qprime), PQIndex(p,0),
	      DComplex( (*rowp0)[qprime]),
	      DComplex(0.,0.));
      }
    }

    // Now use recursion to climb the q ladder (for q<=p).
    (*rowpq) = (*rowp0);
    for (q=1; q<=p && q<=orderIn-p; q++) {
      PQIndex pqi(p,q);
      // cycle buffers
      rowtemp = rowpqm1; rowpqm1 = rowpq; rowpq = rowtemp;
      // generate [pprime qprime][pq] vector.  Recursion is
      // sqrt(q) D[p'q'][pq] = cosh(mu) sqrt(q') D[p' q'-1][p q-1]
      //		-sinh(mu) sqrt(p'+1) D[p'+1 q'][p q-1].
      double factor1=cmu/sqrtp[q];
      double factor2=-smu/sqrtp[q];
      PQIndex pqo(pqi.m(), 0);
      qprime=pqo.getQ();
      (*rowpq)[qprime] = factor2* sqrtp[pqo.getP()+1] * (*rowpqm1)[qprime];
      if (pqo.N()<=orderOut) {
	  if (pqo.isReal())
	    D.set(pqo, pqi,
		  DComplex((*rowpq)[qprime]),
		  DComplex((*rowpq)[qprime]));
	  else
	    D.set(pqo, pqi,
		  DComplex((*rowpq)[qprime]),
		  DComplex(0.,0.) );
      }


      for (pqo.incN(); pqo.N()<=sumOrder-pqi.N(); pqo.incN()) {
	qprime=pqo.getQ();
	(*rowpq)[qprime] = factor1 * sqrtp[qprime] * (*rowpqm1)[qprime-1]
	  + factor2* sqrtp[pqo.getP()+1] * (*rowpqm1)[qprime];
	// Fill spots in the matrix:
	if (pqo.N()<=orderOut) {
	  if (pqo.isReal())
	    D.set(pqo, pqi,
		  DComplex((*rowpq)[qprime]),
		  DComplex((*rowpq)[qprime]));
	  else
	    D.set(pqo, pqi,
		  DComplex((*rowpq)[qprime]),
		  DComplex(0.,0.) );
	}
      }
    }
  }
  
  return D;
}

//-------------------------------------------------------------------
// TRANSFORMATION MATRICES: General Ellipse (Translate+Dilate+Shear)
//-------------------------------------------------------------------

LTransform 
MakeLTransform(const Ellipse& e, 
	       int orderOut, int orderIn,
	       bool shiftCoords,
	       int orderIntermediate) {
  // Recall that the Ellipse transform is defined to be
  // translate, then scale, then shear.
  // Create the matrices in reverse order (left to right)

  // Default is to have intermediate results be larger of in/out orders
  if (orderIntermediate<0)
    orderIntermediate=MAX(orderOut,orderIn);

  bool identityShear = (e.getS().getE()==0.);
  bool identityDilate = e.getMu()==0.;
  bool identityTranslate = (e.getX0().x==0. && e.getX0().y==0.);

  if (identityShear && identityDilate && identityTranslate) {
    // Return an identity transformation
    LTransform mOut(orderOut, orderIn);
    for (PQIndex pq; 
	 pq.N()<=MIN(orderOut,orderIn);
	 pq.nextDistinct())
      mOut.set(pq,pq,
	       DComplex(1.,0.),
	       DComplex(1.,0.));
    return mOut;
  }

  // Now start accumulating the three matrices as needed
  LTransform mOut(0,0);

  if (shiftCoords) {
    // compound in order SDT to describe same flux with new coords:
    if (!identityShear) {
      int thisIn = ( (identityDilate && identityTranslate) ? 
		     orderIn : orderIntermediate);
      mOut = MakeLTransform(e.getS(), orderOut, thisIn,shiftCoords);
    }

    if (!identityDilate) {
      int thisIn = ( identityTranslate ? 
		     orderIn : orderIntermediate);
      if (identityShear) {
	// Dilation is the first operation:
	mOut = MakeLTransform(e.getMu(), orderOut, thisIn, shiftCoords);
      } else {
	// Compound dilation with shear:
	mOut = mOut * MakeLTransform(e.getMu(), orderIntermediate, thisIn,
				     shiftCoords);
      }
    }

    if (!identityTranslate) {
      if (identityShear && identityDilate) {
	// Translation is the first operation:
	mOut = MakeLTransform(e.getX0(), orderOut, orderIn,
			      shiftCoords);
      } else {
	// Compound shear with dilation+shear
	mOut = mOut * MakeLTransform(e.getX0(), orderIntermediate, orderIn,
				     shiftCoords);
      }
    }

  } else {
    // compound in order TDS to move flux on fixed coords:

    if (!identityTranslate) {
      int thisIn = ( (identityDilate && identityShear) ? 
		     orderIn : orderIntermediate);
      mOut = MakeLTransform(e.getX0(), orderOut, thisIn, shiftCoords);
    }

    if (!identityDilate) {
      int thisIn = ( identityShear ? 
		     orderIn : orderIntermediate);
      if (identityTranslate) {
	// Dilation is the first operation:
	mOut = MakeLTransform(e.getMu(), orderOut, thisIn, shiftCoords);
      } else {
	// Compound dilation with shear:
	mOut = mOut * MakeLTransform(e.getMu(), orderIntermediate, thisIn,
				     shiftCoords);
      }
    }

    if (!identityShear) {
      if (identityTranslate && identityDilate) {
	// Translation is the first operation:
	mOut = MakeLTransform(e.getS(), orderOut, orderIn, shiftCoords);
      } else {
	// Compound shear with dilation+translation
	mOut = mOut * MakeLTransform(e.getS(), orderIntermediate, orderIn,
				     shiftCoords);
      }
    }
  }
  return mOut;
}

//-------------------------------------------------------------------
// TRANSFORMATION MATRICES: Rotation
//-------------------------------------------------------------------

// Note that executing a rotation with the matrix is a bit of a waste
// since operation should be O(N) instead of O(N^2)
LTransform
RotationLTransform(double theta,
		   int orderOut, int orderIn,
		   bool coordShift) {

  assert(orderOut>=0);
  assert(orderIn>=0);

  int maxM = MIN(orderIn,orderOut);

  if (coordShift) theta=-theta;	//invert for coord transformation

  CVector e_imtheta(maxM+1);  // vector of e^{-im theta}
  const DComplex ONE(1.,0.);
  const DComplex ZERO(0.,0.);
  e_imtheta[0] = ONE;
  DComplex e_itheta(cos(theta), -sin(theta));
  for (int i=1; i<=maxM; i++)
    e_imtheta[i] = e_imtheta[i-1] * e_itheta;

  LTransform D(orderOut, orderIn);

  // Advance through p's
  for (PQIndex pq(0,0); 
       !pq.pastOrder(maxM);
       pq.nextDistinct()) {
    if (pq.m()==0) 
      D.set(pq, pq, ONE, ONE);
    else
      D.set(pq, pq, e_imtheta[pq.m()], ZERO);
  }
  return D;
}

//-------------------------------------------------------------------
// CONVOLUTION MATRIX:
//-------------------------------------------------------------------

// Make the "F matrix" for convolutions, which is defined to be
// sqrt( pi! p*! / po! Delta!) * G(po, pi, p*) from Appendix B of BJ02.

// If pointer is null, new matrix is created, else resize & fill.
// D is defined as sig_i^2 / sig_o^2.

// The 3d Array will be calculated only to the desired orders in the PSF,
// intrinsic, and observed (Out) b vectors, but a cubic array will be
// allocated for convenience.
void ConvolveFMatrix(mv::Cube<double>* &fptr,
		     const double D,
		     const int orderOut,
		     const int orderIn,
		     const int orderStar)
{

  double sigi = sqrt(D);
  double sigstar = sqrt(1-D);

  int maxOrder = MAX(orderIn, orderStar);
  maxOrder = MAX(maxOrder, orderOut);

  // Get the matrix
  if (fptr==0) fptr = new mv::Cube<double>(maxOrder+1); 
  else fptr->resize(maxOrder+1);
  *fptr = 0.;

  static valarray<double> sqrtp;	//commonly used vector of sqrt(p)
  {
    int topSqrt = orderStar + orderIn;
    if (sqrtp.size() < topSqrt+1) {
      sqrtp.resize(topSqrt+1);
      for (int i=0; i<=topSqrt; i++) sqrtp[i]=sqrtn(i); 
    }
  }

  int po, pi, pstar;

  // Set up the po=0 elements
  pi = 0;
  for (double ti=1. ; pi<=orderIn; pi++, ti*=sigstar) {
    (*fptr)(0,pi,0) = ti;
    pstar = 1;
    for (double tstar=ti; 
	 pstar<=orderStar; 
	 pstar++) {
      tstar *= -sigi*sqrtp[pstar+pi]/sqrtp[pstar];
      (*fptr)(0,pi,pstar) = tstar;
    }
  }

  // Now po>0:
  //double next;
  for (po=1; po<=MIN(orderOut,orderIn+orderStar); po++) {
    // Do the pi=0 separately
    for (pstar=po; pstar<=orderStar; pstar++) {
      (*fptr)(po,0,pstar) = sqrtp[pstar]*sigstar*(*fptr)(po-1,0,pstar-1)
	/ sqrtp[po];
    }
    // Now the pi>0:
    for (pi=1; pi<=orderIn; pi++) {
      if (po-pi<1) {
	(*fptr)(po,pi,0) = sqrtp[pi]*sigi*(*fptr)(po-1,pi-1,0)
	/ sqrtp[po];
	pstar = 1;
      } else {
	pstar = po-pi;
      }
      for ( ; pstar<=orderStar; pstar++)
	(*fptr)(po,pi,pstar) = 
	  ( sqrtp[pi]*sigi*(*fptr)(po-1,pi-1,pstar) +
	    sqrtp[pstar]*sigstar*(*fptr)(po-1,pi,pstar-1) )
	  / sqrtp[po];
    }
  }
}
    
// Create the matrix representation of convolution
// psf is the input PSF, and 
// D = "deconvolution factor" = sig_i^2 / sig_o^2

LTransform
MakeLTransform(const LVector psf,
	       const double D,
	       const int orderOut,
	       const int orderIn,
	       const int orderStar)
{
  // 3d Array of p-only coefficients:
  static mv::Cube<double>* F=0;
  static int oOut = -1;
  static int oIn = -1;
  static int oStar = -1;
  static double oldD = -1.;

  assert(psf.getOrder() >= orderStar);

  // (Re-)build the F matrix if needed
  if (F==0
      || orderOut > oOut
      || orderIn > oIn
      || orderStar != oStar
      || D != oldD) {
    ConvolveFMatrix(F, D, orderOut, orderIn, orderStar);
    oOut = orderOut;
    oIn = orderIn;
    oStar = orderStar;
    oldD = D;
  }

  //Make the output transform
  LTransform result(orderOut, orderIn);

  //int po, qo, pi, qi, pstar, qstar, Delta;
  DComplex b;

  for (PQIndex pqi; pqi.N()<=orderIn; pqi.nextDistinct()) {
    for (PQIndex pqo; pqo.N()<=MIN(orderIn,orderOut); pqo.nextDistinct()) {
      DComplex Cpqpq(0.,0.);
      DComplex Cqppq(0.,0.);

      for (PQIndex pqStar(pqo.getP() - pqi.getP(), pqo.getQ() - pqi.getQ());
	   pqStar.N() <= orderStar;
	   pqStar.incN()) {
	if (!pqStar.pqValid()) continue;
	/*cerr << "Computing element for pqi, star, o: " << pqi
		 << " " << pqStar << " " << pqo
		 << " terms " << (*F)(pqo.getP(), pqi.getP(), pqStar.getP())
		 << " , " <<  (*F)(pqo.getQ(), pqi.getQ(), pqStar.getQ())
		 << endl;*/

	Cpqpq += psf[pqStar] *
	  (*F)(pqo.getP(), pqi.getP(), pqStar.getP()) *
	  (*F)(pqo.getQ(), pqi.getQ(), pqStar.getQ());
      }

      // Now calculate for swap of p0 & q0:
      if (pqo.isReal()) Cqppq = Cpqpq;
      else {
	for (PQIndex pqStar(pqo.getQ() - pqi.getP(), pqo.getP() - pqi.getQ());
	     pqStar.N() <= orderStar;
	     pqStar.incN()) {
	  if (!pqStar.pqValid()) continue;
	  /*cerr << "Computing element for pqi, star, o: " << pqi
		   << " " << pqStar << " *" << pqo << endl;*/
	  Cqppq += psf[pqStar] *
	    (*F)(pqo.getQ(), pqi.getP(), pqStar.getP()) *
	    (*F)(pqo.getP(), pqi.getQ(), pqStar.getQ());
	}
      }
      result.set(pqo, pqi, Cpqpq, Cqppq);
    }
  }
  return result;
}

} // namespace Laguerre

