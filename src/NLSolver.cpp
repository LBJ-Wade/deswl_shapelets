// The algorithms contained in this file are taken from the paper
// "Methods for Nonlinear Least-Squares Problems", by Madsen, Nielsen,
// and Tingleff (2004).  
// A copy of this paper should be included with the code in the file
// madsen04.pdf.  Please refer to this paper for more details about
// how the algorithms work.


#include "NLSolver.h"
#include <iostream>
#include <limits>
#include <algorithm>

template <class T> inline T SQR(T x) { return x*x; }
const double sqrteps = sqrt(std::numeric_limits<double>::epsilon());

#define dbg if(nlout) (*nlout)
#define xdbg if(verbose && nlout) (*nlout)

void NLSolver::J(
    const tmv::Vector<double>& x, const tmv::Vector<double>& f, 
    tmv::Matrix<double>& df) const
{
  // Do a finite difference calculation for J.
  // This function is virtual, so if there is a better way to 
  // calculate J, then you should override this version.

  tmv::Vector<double> x2 = x;
  tmv::Vector<double> f2(f.size());
  for(size_t j=0;j<x.size();j++) {
    const double dx = sqrteps * (Norm(x) + 1.);
    x2(j) += dx;
    this->F(x2,f2);
    df.col(j) = (f2-f)/dx;
    x2(j) = x(j);
  }
}

bool NLSolver::TestJ(
    const tmv::Vector<double>& x, tmv::Vector<double>& f,
    std::ostream* os, double relerr) const
{
  this->F(x,f);
  _pJ.reset(new tmv::Matrix<double>(f.size(),x.size()));
  tmv::Matrix<double>& J = *_pJ;
  this->J(x,f,J);
  tmv::Matrix<double> Jn(f.size(),x.size());
  NLSolver::J(x,f,Jn);
  double err = MaxAbsElement(J-Jn) / Norm(Jn);
  if (!relerr) relerr = 10.*sqrteps;
  if (os) {
    *os << "TestJ:\n";
    if (verbose) {
      *os << "x = "<<x<<std::endl;
      *os << "f = "<<f<<std::endl;
      *os << "Direct J = "<<J<<std::endl;
      *os << "Numeric J = "<<Jn<<std::endl;
    }
    *os << "MaxAbsElement(J-J_num) / Norm(J) = "<<err<<std::endl;
    *os << "cf. relerr = "<<relerr<<std::endl;
    if (err >= relerr) {
      tmv::Matrix<double> diff = J-Jn;
      *os << "J-J_num = "<<diff;
      double maxel = diff.MaxAbsElement();
      *os << "Max element = "<<maxel<<std::endl;
      for(size_t i=0;i<diff.colsize();i++) {
	for(size_t j=0;j<diff.rowsize();j++) {
	  if (std::abs(diff(i,j)) > 0.9*maxel) {
	    *os<<"J("<<i<<','<<j<<") = "<<J(i,j)<<"  ";
	    *os<<"J_num("<<i<<','<<j<<") = "<<Jn(i,j)<<"  ";
	    *os<<"diff = "<<J(i,j)-Jn(i,j)<<std::endl;
	  }
	}
      }
    }
  }
  return err < relerr;
}

// H(i,j) = d^2 Q / dx_i dx_j
// where Q = 1/2 Sum_k |f_k|^2
// H = JT J + Sum_k f_k d^2f_k/(dx_i dx_j)
void NLSolver::numericH(
    const tmv::Vector<double>& x,
    const tmv::Vector<double>& f, 
    tmv::SymMatrix<double>& h) const
{
  // Do a finite difference calculation for H.

  const double sqrteps = sqrt(std::numeric_limits<double>::epsilon());
  const double dx = sqrt(sqrteps) * (Norm(x) + 1.);
  double q0 = 0.5 * NormSq(f);

  tmv::Vector<double> x2 = x;
  tmv::Vector<double> f2(f.size());
  for(size_t i=0;i<x.size();i++) {
    x2(i) = x(i) + dx;
    this->F(x2,f2);
    double q2a = 0.5*NormSq(f2);
    x2(i) = x(i) - dx;
    this->F(x2,f2);
    double q2b = 0.5*NormSq(f2);
    h(i,i) = (q2a + q2b - 2.*q0) / (dx*dx);
    x2(i) = x(i);

    for(size_t j=i+1;j<x.size();j++) {
      x2(i) = x(i) + dx;
      x2(j) = x(j) + dx;
      this->F(x2,f2);
      double q2a = 0.5*NormSq(f2);

      x2(i) = x(i) + dx;
      x2(j) = x(j) - dx;
      this->F(x2,f2);
      double q2b = 0.5*NormSq(f2);

      x2(i) = x(i) - dx;
      x2(j) = x(j) + dx;
      this->F(x2,f2);
      double q2c = 0.5*NormSq(f2);

      x2(i) = x(i) - dx;
      x2(j) = x(j) - dx;
      this->F(x2,f2);
      double q2d = 0.5*NormSq(f2);

      h(i,j) = (q2a - q2b - q2c + q2d) / (4.*dx*dx);
      x2(i) = x(i);
      x2(j) = x(j);
    }
  }
}


#define CHECKF(norminff) \
  do { \
    double checkf_temp = (norminff); \
    if (checkf_temp < ftol) { \
      dbg<<"Found ||f|| ~= 0\n"; \
      dbg<<"||f||_inf = "<<checkf_temp<<" < "<<ftol<<std::endl; \
      return true; \
    } \
  } while (false)

#define CHECKG(norminfg) \
  do { \
    double checkg_temp = (norminfg); \
    if (checkg_temp < gtol) { \
      dbg<<"Found local minimum of ||f||\n"; \
      dbg<<"||g||_inf = "<<checkg_temp<<" < "<<gtol<<std::endl; \
      return true; \
    } \
  } while (false)

#define SHOWFAILFG \
  do { \
    dbg<<"||f||_inf = "<<NormInf(f)<<" !< "<<ftol<<std::endl; \
    dbg<<"||g||_inf = "<<NormInf(g)<<" !< "<<gtol<<std::endl; \
  } while (false)

#define CHECKSTEP(normh) \
  do { \
    double checkstep_temp1 = (normh); \
    double checkstep_temp2 = min_step*(Norm(x)+min_step); \
    if (checkstep_temp1 < checkstep_temp2) { \
      dbg<<"Step size became too small\n"; \
      dbg<<"||h|| = "<<checkstep_temp1<<" < "<<checkstep_temp2<<std::endl; \
      SHOWFAILFG; \
      return false; \
    } \
  } while (false)

bool NLSolver::Solve_Newton(
    tmv::Vector<double>& x, tmv::Vector<double>& f) const
// This is a simple descent method which uses either the 
// Newton direction or the steepest descent direction.
{
  const double gamma1 = 0.1;
  const double gamma2 = 0.5;
  dbg<<"Start Solve_Newton\n";

  _pJ.reset(new tmv::Matrix<double>(f.size(),x.size()));
  tmv::Matrix<double>& J = *_pJ;
  tmv::Vector<double> g(x.size());
  tmv::Vector<double> h(x.size());
  tmv::Vector<double> xnew(x.size());
  tmv::Vector<double> fnew(f.size());
  tmv::Vector<double> gnew(x.size());

  xdbg<<"x = "<<x<<std::endl;
  this->F(x,f);
  xdbg<<"f = "<<f<<std::endl;
  CHECKF(NormInf(f));
  double Q = 0.5*NormSq(f);
  xdbg<<"Q = "<<Q<<std::endl;
  this->J(x,f,J);
  if (usesvd) J.DivideUsing(tmv::SV);
  J.SaveDiv();
  xdbg<<"J = "<<J<<std::endl;
  g = J.Transpose() * f;
  xdbg<<"g = "<<g<<std::endl;
  CHECKG(NormInf(g));
  double alpha = Q/NormSq(g);
  bool newton = true;

  dbg<<"iter   |f|inf   Q   |g|inf   alpha\n";
  for(int k=0;k<max_iter;k++) {
    newton = true;
    dbg<<k<<"   "<<NormInf(f)<<"   "<<Q<<"   "<<NormInf(g)<<"   "<<alpha<<std::endl;

    h = -f/J;
    xdbg<<"h = "<<h<<std::endl;
    double normh = Norm(h);
    CHECKSTEP(normh);

    // phi(alpha) = Q(x + alpha h)
    // phi'(alpha) = fT J h = gT h
    // where g is measured at xnew, not x
    // m = phi'(0)
    double m = h*g;
    double normg = Norm(g);
    double m2 = -normg*normg;

    if ((k%5 == 0 && m >= 0.) || (k%5 != 0 && m/normh >= -0.01*normg)) {
      // i.e. either m >= 0 or |m/normh| < 0.01 * |m2/normh2|
      newton = false;
      xdbg<<"Newton is not a good descent direction - use steepest descent\n";
      h = -g;
      CHECKSTEP(normg);
      m = m2;
    } else {
      xdbg<<"m = "<<m<<", Steepest m = "<<m2<<std::endl;
      xdbg<<"m/Norm(h) = "<<m/normh<<", Steepest m/Norm(h) = "<<-normg<<std::endl;
    }

    if (newton && alpha > 0.1) alpha = 1.0;
    for(int k2=0;k2<=max_iter;k2++) {
      if (k2 == max_iter) { 
	dbg<<"Maximum iterations exceeded in subloop of Newton method\n";
	dbg<<"This can happen when there is a singularity (or close to it)\n";
	dbg<<"along the gradient direction:\n";
	if (usesvd) {
	  dbg<<"J Singular values = \n"<<J.SVD().GetS().diag()<<std::endl;
	  dbg<<"V = \n"<<J.SVD().GetV()<<std::endl;
	}
	SHOWFAILFG; 
	return false;
      }
      xnew = x + alpha*h;
      if (alpha < min_step) {
	dbg<<"alpha became too small ("<<alpha<<" < "<<min_step<<")\n";
	SHOWFAILFG; 
	return false;
      }
      xdbg<<"xnew = "<<xnew<<std::endl;
      this->F(xnew,fnew);
      xdbg<<"fnew = "<<fnew<<std::endl;
      double Qnew = 0.5*NormSq(fnew);
      xdbg<<"Qnew = "<<Qnew<<std::endl;

      // Check that phi has decreased significantly
      // Require phi(alpha) <= phi(0) + gamma1 phi'(0) alpha
      if (Qnew > Q + gamma1 * m * alpha) {
	alpha /= 2;
	newton = false;
	xdbg<<"Qnew not small enough: alpha => "<<alpha<<std::endl;
	continue;
      }
      this->J(xnew,fnew,J);
      J.UnSetDiv();
      xdbg<<"Jnew = "<<J<<std::endl;
      gnew = J.Transpose() * fnew;
      xdbg<<"gnew = "<<gnew<<std::endl;

      // Check that alpha is not too small
      // Require phi'(alpha) >= gamma2 phi'(0)
      double mnew = h*gnew;
      if (mnew < gamma2 * m) {
	alpha *= 3.;
	newton = false;
	xdbg<<"New slope not shallow enough: alpha => "<<alpha<<std::endl;
	xdbg<<"(m = "<<m<<", mnew = "<<mnew<<")\n";
	continue;
      }
      xdbg<<"Good choice\n";
      x = xnew; f = fnew; Q = Qnew; g = gnew;
      break;
    }
    CHECKF(NormInf(f));
    CHECKG(NormInf(g));
  }
  dbg<<"Maximum iterations exceeded in Newton method\n";
  SHOWFAILFG; 
  return false;
}

bool NLSolver::Solve_LM(
    tmv::Vector<double>& x, tmv::Vector<double>& f) const
// This is the Levenberg-Marquardt method
{
  dbg<<"Start Solve_LM\n";

  _pJ.reset(new tmv::Matrix<double>(f.size(),x.size()));
  tmv::Matrix<double>& J = *_pJ;
  tmv::Vector<double> h(x.size());
  tmv::Vector<double> xnew(x.size());
  tmv::Vector<double> fnew(f.size());
  tmv::Vector<double> gnew(x.size());

  xdbg<<"x = "<<x<<std::endl;
  this->F(x,f);
  xdbg<<"f = "<<f<<std::endl;
  CHECKF(NormInf(f));
  double Q = 0.5*NormSq(f);
  xdbg<<"Q = "<<Q<<std::endl;
  this->J(x,f,J);
  if (usesvd) J.DivideUsing(tmv::SV);
  J.SaveDiv();
  xdbg<<"J = "<<J<<std::endl;
  tmv::Vector<double> g = J.Transpose() * f;
  xdbg<<"g = "<<g<<std::endl;
  CHECKG(NormInf(g));

  tmv::SymMatrix<double> A = J.Transpose() * J;
  xdbg<<"JT J = "<<A<<std::endl;
  if (usesvd) A.DivideUsing(tmv::SV);
  else if (startwithch) A.DivideUsing(tmv::CH);
  else A.DivideUsing(tmv::LU);
  double mu = tau * NormInf(A.diag());
  xdbg<<"initial mu = "<<tau<<" * "<<NormInf(A.diag())<<" = "<<mu<<std::endl;
  A += mu;
  A.SaveDiv();
  double nu = 2.;

  dbg<<"iter   |f|inf   Q   |g|inf   mu\n";
  for(int k=0;k<max_iter;k++) {
    dbg<<k<<"   "<<NormInf(f)<<"   "<<Q<<"   "<<NormInf(g)<<"   "<<mu<<std::endl;
    xdbg<<"k = "<<k<<std::endl;
    xdbg<<"mu = "<<mu<<std::endl;
    xdbg<<"A = "<<A<<std::endl;
#ifndef NOTHROW
    try {
#endif
      //dbg<<"Before h=-g/A"<<std::endl;
      h = -g/A;
      //dbg<<"After h=-g/A"<<std::endl;
#ifndef NOTHROW
    } 
    catch (tmv::NonPosDef) {
      xdbg<<"NonPosDef caught - switching division to LU method.\n";
      // Once the Cholesky decomp fails, just use LU from that point on.
      A.DivideUsing(tmv::LU);
      h = -g/A;
    }
#endif
    xdbg<<"h = "<<h<<std::endl;
    CHECKSTEP(Norm(h));

    xnew = x + h;
    xdbg<<"xnew = "<<xnew<<std::endl;
    this->F(xnew,fnew);
    xdbg<<"fnew = "<<fnew<<std::endl;
    double Qnew = 0.5*NormSq(fnew);
    xdbg<<"Qnew = "<<Qnew<<std::endl;

    if (Qnew < Q) {
      xdbg<<"improved\n";
      x = xnew; f = fnew; 
      CHECKF(NormInf(f));

      this->J(x,f,J);
      J.UnSetDiv();
      A = J.Transpose() * J;
      gnew = J.Transpose() * f;
      xdbg<<"gnew = "<<gnew<<std::endl;
      CHECKG(NormInf(gnew));

      // Use g as a temporary for (g - mu*h)
      g -= mu*h;
      double rho = (Q-Qnew) / (-0.5*h*g);
      xdbg<<"rho = "<<Q-Qnew<<" / "<<(-0.5*h*g)<<" = "<<rho<<std::endl;
      mu *= std::max(1./3.,1.-std::pow(2.*rho-1.,3)); nu = 2.;
      xdbg<<"mu *= "<<std::max(1./3.,1.-std::pow(2.*rho-1.,3))<<" = "<<mu<<std::endl;
      A += mu;
      A.UnSetDiv();
      Q = Qnew; g = gnew;
    } else {
      xdbg<<"not improved\n";
      A += mu*(nu-1.); mu *= nu; nu *= 2.;
      A.UnSetDiv();
      xdbg<<"mu *= (nu = "<<nu<<") = "<<mu<<std::endl;
    }
  }
  dbg<<"Maximum iterations exceeded in LM method\n";
  SHOWFAILFG; 
  return false;
}

bool NLSolver::Solve_Dogleg(
    tmv::Vector<double>& x, tmv::Vector<double>& f) const
// This is the Dogleg method
{
  dbg<<"Start Solve_Dogleg\n";
  _pJ.reset(new tmv::Matrix<double>(f.size(),x.size()));
  tmv::Matrix<double>& J = *_pJ;
  tmv::Vector<double> h(x.size());
  tmv::Vector<double> temp(x.size());
  tmv::Vector<double> xnew(x.size());
  tmv::Vector<double> fnew(f.size());

  xdbg<<"x = "<<x<<std::endl;
  this->F(x,f);
  xdbg<<"f = "<<f<<std::endl;
  CHECKF(NormInf(f));
  double Q = 0.5*NormSq(f);
  xdbg<<"Q = "<<Q<<std::endl;
  this->J(x,f,J);
  if (usesvd) J.DivideUsing(tmv::SV);
  J.SaveDiv();
  xdbg<<"J = "<<J<<std::endl;
  xdbg<<"J.svd = "<<J.SVD().GetS().diag()<<std::endl;

  tmv::Vector<double> g = J.Transpose() * f;
  xdbg<<"g = "<<g<<std::endl;
  CHECKG(NormInf(g));

  double delta = delta0;
  int maxnsing = std::min(f.size(),x.size());
  int nsing = maxnsing;

  dbg<<"iter   |f|inf   Q   |g|inf   delta\n";
  for(int k=0;k<max_iter;k++) {
    dbg<<k<<"   "<<NormInf(f)<<"   "<<Q<<"   "<<NormInf(g)<<"   "<<delta<<std::endl;
    if (usesvd && nsing == maxnsing && J.Singular() && nsing > 1) {
      dbg<<"Singular J, so try lowering number of singular values.\n";
      nsing = J.SVD().GetKMax();
      dbg<<"J Singular values = \n"<<J.SVD().GetS().diag()<<std::endl;
      dbg<<"nsing -> "<<nsing<<std::endl;
    }
    h = -f/J;
    xdbg<<"h = "<<h<<std::endl;

    double normsqg = NormSq(g);
    double normh = Norm(h);
    double normh1 = normh;
    double rhodenom;

    if (normh <= delta) {
      xdbg<<"|h| < delta\n";
      rhodenom = Q;
      xdbg<<"rhodenom = "<<rhodenom<<std::endl;
    } else {
      double alpha = normsqg / NormSq(J*g);
      double normg = sqrt(normsqg);
      if (normg >= delta / alpha) {
	xdbg<<"|g| > delta/alpha\n";
	h = -(delta / normg) * g;
	xdbg<<"h = "<<h<<std::endl;
	rhodenom = delta*(2.*alpha*normg-delta)/(2.*alpha);
	xdbg<<"rhodenom = "<<rhodenom<<std::endl;
      }
      else {
	xdbg<<"dogleg\n";
	temp = h + alpha*g;
	double a = NormSq(temp);
	double b = -alpha * g * temp;
	double c = alpha*alpha*NormSq(g)-delta*delta;
	// beta is the solution of 0 = a beta^2 + 2b beta + c
	xdbg<<"a,b,c = "<<a<<" "<<b<<" "<<c<<std::endl;
	double beta = (b <= 0) ?
	  (-b + sqrt(b*b - a*c)) / a :
	  -c / (b + sqrt(b*b - a*c));
	xdbg<<"alpha = "<<alpha<<std::endl;
	xdbg<<"beta = "<<beta<<std::endl;
	h = -alpha*g + beta*temp;
	xdbg<<"h = "<<h<<std::endl;
	xdbg<<"Norm(h) = "<<Norm(h)<<"  delta = "<<delta<<std::endl;
	rhodenom = 0.5*alpha*SQR((1.-beta)*normg) + beta*(2.-beta)*Q;
	xdbg<<"rhodenom = "<<rhodenom<<std::endl;
      }
      normh = Norm(h);
    }

    CHECKSTEP(normh);

    xnew = x + h;
    xdbg<<"xnew = "<<xnew<<std::endl;
    this->F(xnew,fnew);
    xdbg<<"fnew = "<<fnew<<std::endl;
    double Qnew = 0.5*NormSq(fnew);
    xdbg<<"Qnew = "<<Qnew<<std::endl;

    bool deltaok = false;
    if (Qnew < Q) {
      double rho = (Q-Qnew) / rhodenom;
      xdbg<<"rho = "<<Q-Qnew<<" / "<<rhodenom<<" = "<<rho<<std::endl;
      x = xnew; f = fnew; Q = Qnew;
      CHECKF(NormInf(f));
      this->J(x,f,J);
      J.UnSetDiv();
      g = J.Transpose() * f;
      xdbg<<"g = "<<g<<std::endl;
      CHECKG(NormInf(g));
      if (rho > 0.75) {
	delta = std::max(delta,3.*normh);
	deltaok = true;
      }
    }
    if (deltaok) {
      nsing = maxnsing;
    } else {
      double normsqh = normh1*normh1;
      if (usesvd && delta < normh1 && normsqg < 0.01 * normsqh && nsing > 1) {
	dbg<<"normsqg == "<<normsqg/normsqh<<
	  " * normsqh, so try lowering number of singular values.\n";
        --nsing;
	dbg<<"nsing -> "<<nsing<<std::endl;
	dbg<<"J Singular values = \n"<<J.SVD().GetS().diag()<<std::endl;
	J.SVD().Top(nsing);
      } else {
	delta /= 2.;
	double min_delta = min_step * (Norm(x)+min_step);
	if (delta < min_delta) {
	  dbg<<"delta became too small ("<<delta<<" < "<<min_delta<<")\n";
	  SHOWFAILFG; 
	  return false;
	}
      }
    }
  }
  dbg<<"Maximum iterations exceeded in Dogleg method\n";
  SHOWFAILFG; 
  return false;
}

bool NLSolver::Solve_Hybrid(
    tmv::Vector<double>& x, tmv::Vector<double>& f) const
// This is the Hybrid method which starts with the L-M method,
// but switches to a quasi-newton method if ||f|| isn't approaching 0.
{
  dbg<<"Start Solve_Hybrid\n";
  _pJ.reset(new tmv::Matrix<double>(f.size(),x.size()));
  tmv::Matrix<double>& J = *_pJ;
  tmv::Vector<double> h(x.size());
  tmv::Vector<double> xnew(x.size());
  tmv::Vector<double> fnew(f.size());
  tmv::Vector<double> gnew(x.size());
  tmv::Matrix<double> Jnew(f.size(),x.size());
  tmv::Vector<double> y(x.size());
  tmv::Vector<double> v(x.size());

  xdbg<<"x = "<<x<<std::endl;
  this->F(x,f);
  xdbg<<"f = "<<f<<std::endl;
  double norminff = NormInf(f);
  CHECKF(norminff);
  double Q = 0.5*NormSq(f);
  xdbg<<"Q = "<<Q<<std::endl;
  this->J(x,f,J);
  if (usesvd) J.DivideUsing(tmv::SV);
  J.SaveDiv();
  xdbg<<"J = "<<J<<std::endl;

  tmv::SymMatrix<double> A = J.Transpose()*J;
  if (usesvd) A.DivideUsing(tmv::SV);
  else if (startwithch) A.DivideUsing(tmv::CH);
  else A.DivideUsing(tmv::LU);
  A.SaveDiv();
  tmv::SymMatrix<double> H(x.size());
  if (usesvd) H.DivideUsing(tmv::SV);
  else if (startwithch) H.DivideUsing(tmv::CH);
  else H.DivideUsing(tmv::LU);
  H.SaveDiv();
  bool directH = hasdirecth;
  if (directH) {
#ifndef NOTHROW
    try {
#endif
      this->H(x,f,J,H);
#ifndef NOTHROW
    }
    catch(int) {
      dbg<<"No direct H calculation - calculate on the fly\n";
      directH = false;
      H.SetToIdentity();
    }
#endif
  } else {
    H.SetToIdentity();
  }

  tmv::Vector<double> g = J.Transpose() * f;
  xdbg<<"g = "<<g<<std::endl;
  double norminfg = NormInf(g);
  CHECKG(norminfg);

  double mu = tau * NormInf(A.diag());
  A += mu;
  double nu = 2.;
  double delta = delta0;
  bool quasinewton = false;
  int count = 0;

  dbg<<"iter   |f|inf   Q   |g|inf   mu   delta  LM/QN\n";
  for(int k=0;k<max_iter;k++) {
    dbg<<k<<"   "<<norminff<<"   "<<Q<<"   "<<norminfg<<"   "<<mu<<"   "<<delta<<"   "<<(quasinewton?"QN":"LM")<<std::endl;
    xdbg<<"k = "<<k<<std::endl;
    xdbg<<"mu = "<<mu<<std::endl;
    xdbg<<"delta = "<<delta<<std::endl;
    xdbg<<"A = "<<A<<std::endl;
    xdbg<<"H = "<<H<<std::endl;
    xdbg<<"method = "<<(quasinewton ? "quasinewton\n" : "LM\n");
    bool better = false;
    bool switchmethod = false;

    if (quasinewton) {
#ifndef NOTHROW
      try { 
#endif
	//dbg<<"Before h=-g/H"<<std::endl;
	h = -g/H; 
	//dbg<<"After h=-g/H"<<std::endl;
#ifndef NOTHROW
      }
      catch (tmv::NonPosDef) {
	xdbg<<"NonPosDef caught - switching division to LU method for H\n";
	H.DivideUsing(tmv::LU); 
	h = -g/H;
      }
#endif
    } else {
#ifndef NOTHROW
      try { 
#endif
	//dbg<<"Before h=-g/A"<<std::endl;
	h = -g/A; 
	//dbg<<"After h=-g/A"<<std::endl;
#ifndef NOTHROW
      }
      catch (tmv::NonPosDef) {
	xdbg<<"NonPosDef caught - switching division to LU method for A\n";
	A.DivideUsing(tmv::LU);
	h = -g/A; 
      }
#endif
    }

    xdbg<<"h = "<<h<<std::endl;
    double normh = Norm(h);
    CHECKSTEP(normh);
    if (quasinewton && normh > delta) h *= delta/normh;

    xnew = x + h;
    xdbg<<"xnew = "<<xnew<<std::endl;
    this->F(xnew,fnew);
    xdbg<<"fnew = "<<fnew<<std::endl;
    double Qnew = 0.5*NormSq(fnew);
    xdbg<<"Qnew = "<<Qnew<<std::endl;

    if (!directH || quasinewton || Qnew < Q) {
      this->J(xnew,fnew,Jnew);
      xdbg<<"Jnew = "<<Jnew<<std::endl;
    }
    if (quasinewton || Qnew < Q) {
      gnew = Jnew.Transpose() * fnew;
      xdbg<<"gnew = "<<gnew<<std::endl;
    }
    double norminfgnew = NormInf(gnew);

    if (quasinewton) {
      xdbg<<"quasinewton\n";
      better = (Qnew < Q) || (Qnew <= (1.+sqrteps)*Q && norminfgnew < norminfg);
      xdbg<<"better = "<<better<<std::endl;
      switchmethod = (norminfgnew >= norminfg);
      xdbg<<"switchmethod = "<<switchmethod<<std::endl;
      if (Qnew < Q) {
	double rho = (Q-Qnew) / (-h*g-0.5*NormSq(J*h));
	if (rho > 0.75) {
	  delta = std::max(delta,3.*normh);
	}
	else if (rho < 0.25) {
	  delta /= 2.;
	  double min_delta = min_step * (Norm(x)+min_step);
	  if (delta < min_delta) {
	    dbg<<"delta became too small ("<<delta<<" < "<<min_delta<<")\n";
	    SHOWFAILFG; 
	    return false;
	  }
	}
      } else {
	delta /= 2.;
	double min_delta = min_step * (Norm(x)+min_step);
	if (delta < min_delta) {
	  dbg<<"delta became too small ("<<delta<<" < "<<min_delta<<")\n";
	  SHOWFAILFG; 
	  return false;
	}
      }
    } else {
      xdbg<<"LM\n";
      if (Qnew <= Q) {
	better = true;
	// we don't need the g vector anymore, so use this space
	// to calculate g-mu*h
	//double rho = (Q-Qnew) / (0.5*h*(mu*h-g));
	g -= mu*h;
	double rho = (Q-Qnew) / (-0.5*h*g);
	mu *= std::max(1./3.,1.-std::pow(2.*rho-1.,3)); nu = 2.;
	xdbg<<"check1: "<<norminfgnew<<" <? "<<0.02*Qnew<<std::endl;
	xdbg<<"check2: "<<Q-Qnew<<" <? "<<0.02*Qnew<<std::endl;
	if (std::min(norminfgnew,Q-Qnew) < 0.02 * Qnew) {
	  count++;
	  if (count == 3) switchmethod = true;
	} else {
	  count = 0;
	}
	if (count != 3) {
	  A = Jnew.Transpose() * Jnew;
	  A += mu;
	}
      } else {
	A += mu*(nu-1.); mu *= nu; nu *= 2.;
	count = 0;
	// MJ: try this?
	switchmethod = (nu >= 32.);
      }
      A.UnSetDiv();
      xdbg<<"better = "<<better<<std::endl;
      xdbg<<"switchmethod = "<<switchmethod<<std::endl;
      xdbg<<"count = "<<count<<std::endl;
    }

    if (!directH) {
      y = Jnew.Transpose()*(Jnew*h) + (Jnew-J).Transpose()*fnew;
      double hy = h*y;
      xdbg<<"hy = "<<hy<<std::endl;
      if (hy > 0.) {
	v = H*h;
        xdbg<<"v = "<<v<<std::endl;
        xdbg<<"y = "<<y<<std::endl;
	xdbg<<"hv = "<<h*v<<std::endl;
	H -= (1./(h*v)) * (v^v);
	xdbg<<"H -> "<<H<<std::endl;
	H += (1./hy) * (y^y);
	H.UnSetDiv();
	xdbg<<"H -> "<<H<<std::endl;
      }
    }
    
    if (better) {
      xdbg<<"better"<<std::endl;
      x = xnew; f = fnew; Q = Qnew; J = Jnew; g = gnew; 
      J.UnSetDiv();
      norminff = NormInf(f); norminfg = norminfgnew;
      if (directH && quasinewton && !switchmethod) this->H(x,f,J,H);
      CHECKF(norminff);
      CHECKG(norminfg);
    }
    if (switchmethod) {
      if (quasinewton) {
	xdbg<<"switch to LM\n";
	A = J.Transpose() * J;
	//mu = tau * NormInf(A.diag());
	A += mu;
	A.UnSetDiv();
	quasinewton = false;
	count = 0;
      } else {
	xdbg<<"switch to quasinewton\n";
	delta = std::max(1.5*min_step*(Norm(x)+min_step),0.2*normh);
	if (directH) {
	  this->H(x,f,J,H);
	  H.UnSetDiv();
	}
	quasinewton = true;
      }
    }
  }
  dbg<<"Maximum iterations exceeded in Hybrid method\n";
  SHOWFAILFG; 
  return false;
}

bool NLSolver::Solve_SecantLM(
    tmv::Vector<double>& x, tmv::Vector<double>& f) const
// This is the Secant version of the Levenberg-Marquardt method
{
  dbg<<"Start Solve_SecantLM\n";
  _pJ.reset(new tmv::Matrix<double>(f.size(),x.size()));
  tmv::Matrix<double>& J = *_pJ;
  tmv::Vector<double> h(x.size());
  tmv::Vector<double> xnew(x.size());
  tmv::Vector<double> fnew(f.size());
  tmv::Vector<double> gnew(x.size());

  xdbg<<"x = "<<x<<std::endl;
  this->F(x,f);
  xdbg<<"f = "<<f<<std::endl;
  CHECKF(NormInf(f));
  double Q = 0.5*NormSq(f);
  xdbg<<"Q = "<<Q<<std::endl;
  this->J(x,f,J);
  if (usesvd) J.DivideUsing(tmv::SV);
  xdbg<<"J = "<<J<<std::endl;
  tmv::SymMatrix<double> A = J.Transpose() * J;
  if (usesvd) A.DivideUsing(tmv::SV);
  else if (startwithch) A.DivideUsing(tmv::CH);
  else A.DivideUsing(tmv::LU);
  tmv::Vector<double> g = J.Transpose() * f;
  xdbg<<"g = "<<g<<std::endl;
  CHECKG(NormInf(g));

  double mu = tau * NormInf(A.diag());
  A += mu;
  double nu = 2.;

  dbg<<"iter   |f|inf   Q   |g|inf   mu\n";
  for(int k=0,j=0;k<max_iter;k++) {
    dbg<<k<<"   "<<NormInf(f)<<"   "<<Q<<"   "<<NormInf(g)<<"   "<<mu<<std::endl;
    xdbg<<"k = "<<k<<std::endl;
    xdbg<<"mu = "<<mu<<std::endl;
    xdbg<<"J = "<<J<<std::endl;
#ifndef NOTHROW
    try {
#endif
      //dbg<<"Before h=-g/A"<<std::endl;
      h = -g/A;
      //dbg<<"After h=-g/A"<<std::endl;
#ifndef NOTHROW
    }
    catch (tmv::NonPosDef) {
      xdbg<<"NonPosDef caught - switching division to LU method.\n";
      A.DivideUsing(tmv::LU);
      h = -g/A;
    }
#endif
    xdbg<<"h = "<<h<<std::endl;
    double normh = Norm(h);
    CHECKSTEP(normh);

    xdbg<<"j = "<<j<<std::endl;
    if (h(j) < 0.8 * normh) {
      xnew = x; 
      //double eta = sqrteps * (std::abs(x(j)) + sqrteps);
      double eta = min_step * (Norm(x) + 1.);
      xnew(j) += eta;
      this->F(xnew,fnew);
      J.col(j) = (fnew-f)/eta;
      xdbg<<"J -> "<<J<<std::endl;
    }
    j = (j+1)%J.ncols();

    xnew = x + h;
    xdbg<<"xnew = "<<xnew<<std::endl;
    this->F(xnew,fnew);
    xdbg<<"fnew = "<<fnew<<std::endl;
    double Qnew = 0.5*NormSq(fnew);
    xdbg<<"Qnew = "<<Qnew<<std::endl;
    J += (1./NormSq(h)) * ((fnew - f - J*h) ^ h);
    xdbg<<"J -> "<<J<<std::endl;

    if (Qnew < Q) {
      x = xnew; f = fnew; 
      CHECKF(NormInf(f));

      A = J.Transpose() * J;
      gnew = J.Transpose() * f;
      CHECKG(NormInf(g));

      g -= mu*h;
      double rho = (Q-Qnew) / (-0.5*h*g);
      xdbg<<"rho = "<<Q-Qnew<<" / "<<(-0.5*h*g)<<" = "<<rho<<std::endl;
      mu *= std::max(1./3.,1.-std::pow(2.*rho-1.,3)); nu = 2.;
      xdbg<<"mu = "<<mu<<std::endl;
      A += mu;
      Q = Qnew; g = gnew;
    } else {
      A += mu*(nu-1.); mu *= nu; nu *= 2.;
    }
  }
  dbg<<"Maximum iterations exceeded in Secant LM method\n";
  SHOWFAILFG; 
  return false;
}

bool NLSolver::Solve_SecantDogleg(
    tmv::Vector<double>& x, tmv::Vector<double>& f) const
// This is the Secant version of the Dogleg method
{
  dbg<<"Start Solve_SecantDogleg\n";
  _pJ.reset(new tmv::Matrix<double>(f.size(),x.size()));
  tmv::Matrix<double>& J = *_pJ;
  tmv::Vector<double> h(x.size());
  tmv::Vector<double> temp(x.size());
  tmv::Vector<double> xnew(x.size());
  tmv::Vector<double> fnew(f.size());
  tmv::Vector<double> y(f.size());
  tmv::Vector<double> djodjy(f.size());

  xdbg<<"x = "<<x<<std::endl;
  this->F(x,f);
  xdbg<<"f = "<<f<<std::endl;
  CHECKF(NormInf(f));
  double Q = 0.5*NormSq(f);
  xdbg<<"Q = "<<Q<<std::endl;
  this->J(x,f,J);
  if (usesvd) J.DivideUsing(tmv::SV);
  //xdbg<<"J = "<<J<<std::endl;
  tmv::Matrix<double> D = J.Inverse();
  //xdbg<<"D = "<<D<<std::endl;

  tmv::Vector<double> g = J.Transpose() * f;
  xdbg<<"g = "<<g<<std::endl;
  CHECKG(NormInf(g));
  double delta = delta0;

  dbg<<"iter   |f|inf   Q   |g|inf   delta\n";
  for(int k=0,j=0;k<max_iter;k++) {
    dbg<<k<<"   "<<NormInf(f)<<"   "<<Q<<"   "<<NormInf(g)<<"   "<<delta<<std::endl;
    //xdbg<<"k = "<<k<<std::endl;
    //xdbg<<"delta = "<<delta<<std::endl;
    h = -D*f;
    xdbg<<"h = "<<h<<std::endl;

    double normsqg = NormSq(g);
    double alpha = normsqg / NormSq(J*g);
    double normh = Norm(h);
    double rhodenom;

    if (normh <= delta) {
      xdbg<<"|h| < delta \n";
      rhodenom = Q;
    } else {
      double normg = sqrt(normsqg);
      if (normg >= delta / alpha) {
	xdbg<<"|g| > delta/alpha \n";
	h = -(delta / normg) * g;
	xdbg<<"h = "<<h<<std::endl;
	rhodenom = delta*(2.*alpha*normg-delta)/(2.*alpha);
      }
      else {
	xdbg<<"dogleg\n";
	temp = h + alpha*g;
	double a = NormSq(temp);
	double b = -alpha * g * temp;
	double c = alpha*alpha*NormSq(g)-delta*delta;
	// beta is the solution of 0 = a beta^2 + 2b beta + c
	//xdbg<<"a,b,c = "<<a<<" "<<b<<" "<<c<<std::endl;
	double beta = (b <= 0) ?
	  (-b + sqrt(b*b - a*c)) / a :
	  -c / (b + sqrt(b*b - a*c));
	xdbg<<"alpha = "<<alpha<<std::endl;
	xdbg<<"beta = "<<beta<<std::endl;
	h = -alpha*g + beta*temp;
	xdbg<<"h = "<<h<<std::endl;
	rhodenom = 0.5*alpha*SQR((1.-beta)*normg) + beta*(2.-beta)*Q;
      }
      normh = Norm(h);
    }

    CHECKSTEP(normh);

    //xdbg<<"j = "<<j<<std::endl;
    bool resetd = false;
    if (h(j) < 0.8 * normh) {
      xnew = x; 
      //double eta = sqrteps * (std::abs(x(j)) + sqrteps);
      double eta = min_step * (Norm(x) + 1.);
      xnew(j) += eta;
      this->F(xnew,fnew);
      y = fnew-f;
      J.col(j) = y/eta;
      double djy = D.row(j)*y;
      if (djy < sqrteps*eta) resetd = true;
      else {
	djodjy = D.row(j)/djy;
        D -= ((D*y) ^ djodjy);
	D.row(j) += eta*djodjy;
      }
      //xdbg<<"J -> "<<J<<std::endl;
      //xdbg<<"D -> "<<D<<std::endl;
      //xdbg<<"D*J -> "<<D*J<<std::endl;
    }
    j = (j+1)%J.ncols();

    xnew = x + h;
    this->F(xnew,fnew);
    double Qnew = 0.5*NormSq(fnew);

    y = fnew - f;
    J += (1./NormSq(h)) * ((fnew - f - J*h) ^ h);
    double hDy = h*D*y;
    if (resetd || hDy < sqrteps*Norm(h)) {
      D = J.Inverse();
    } else {
      D += 1./(hDy) * ((h-D*y) ^ (h*D));
    }
    //xdbg<<"J -> "<<J<<std::endl;
    //xdbg<<"D -> "<<D<<std::endl;
    //xdbg<<"D*J -> "<<D*J<<std::endl;

    if (Qnew < Q) {
      double rho = (Q-Qnew) / rhodenom;
      xdbg<<"rho = "<<Q-Qnew<<" / "<<rhodenom<<" = "<<rho<<std::endl;
      x = xnew; f = fnew; Q = Qnew;
      CHECKF(NormInf(f));
      g = J.Transpose() * f;
      xdbg<<"g = "<<g<<std::endl;
      CHECKG(NormInf(g));
      if (rho > 0.75) {
	delta = std::max(delta,3.*normh);
      }
      else if (rho < 0.25) {
	delta /= 2.;
	double min_delta = min_step * (Norm(x)+min_step);
	if (delta < min_delta) {
	  dbg<<"delta became too small ("<<delta<<" < "<<min_delta<<")\n";
	  SHOWFAILFG; 
	  return false;
	}
      }
    } else {
      delta /= 2.;
      double min_delta = min_step * (Norm(x)+min_step);
      if (delta < min_delta) {
	dbg<<"delta became too small ("<<delta<<" < "<<min_delta<<")\n";
	SHOWFAILFG; 
	return false;
      }
    }
  }
  dbg<<"Maximum iterations exceeded in Secant Dogleg method\n";
  SHOWFAILFG; 
  return false;
}

bool NLSolver::Solve(
    tmv::Vector<double>& x, tmv::Vector<double>& f) const
  // On input, x is the initial guess
  // On output, if return is true, then
  // x is the solution for which either Norm(f) ~= 0
  // or f is a local minimum.
{
#ifndef NOTHROW
  try {
#endif
    switch (method) {
      case Hybrid : return Solve_Hybrid(x,f);
      case Dogleg : return Solve_Dogleg(x,f);
      case LM : return Solve_LM(x,f);
      case Newton : return Solve_Newton(x,f);
      case SecantLM : return Solve_SecantLM(x,f);
      case SecantDogleg : return Solve_SecantDogleg(x,f);
      default : dbg<<"Unknown method\n"; return false;
    }
#ifndef NOTHROW
  } 
#if 0
  catch (int) {}
#else
  catch (tmv::Singular& e) {
    dbg<<"Singular matrix encountered in NLSolver::Solve\n";
    dbg<<e<<std::endl;
  }
  catch (tmv::Error& e) {
    dbg<<"TMV error encountered in NLSolver::Solve\n";
    dbg<<e<<std::endl;
  }
  catch (...) {
    dbg<<"Error encountered in NLSolver::Solve\n";
  }
#endif
  return false;
#endif
}

void NLSolver::getCovariance(tmv::Matrix<double>& cov) const
{
  if (!_pJ.get()) {
    throw std::runtime_error(
	"J not set before calling getCovariance");
  }
  tmv::Matrix<double>& J = *_pJ;
  // This might have changed between solve and getCovariance:
  // And we need to set the threshold to sqrt(eps) rather than eps
  if (usesvd) {
    J.DivideUsing(tmv::SV);
    J.SVD().Thresh(sqrteps); 
  }
  J.InverseATA(cov);
}

void NLSolver::getInverseCovariance(tmv::Matrix<double>& invcov) const
{
  if (!_pJ.get()) {
    throw std::runtime_error(
      "J not set before calling getInverseCovariance");
  }
  tmv::Matrix<double>& J = *_pJ;
  invcov = J.Transpose() * J;
}
