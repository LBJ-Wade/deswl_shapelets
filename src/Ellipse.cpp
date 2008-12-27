
#include "Ellipse.h"
#include "EllipseSolver.h"
#include "TMV.h"
#include "TMV_Tri.h"
#include <cmath>
#include "dbg.h"
#include <fstream>
#include "TimeVars.h"
#include "PsiHelper.h"

#define N_FLUX_ATTEMPTS 0
#define MAXITER 3

static std::complex<double> GammaAdd(const std::complex<double> g1,
    const std::complex<double> g2)
{
  double absg1 = std::abs(g1);
  if (absg1 == 0.) return g2;
  double absd1 = tanh(atanh(absg1)*2.);
  std::complex<double> d1 = absd1 / absg1 * g1;
  double absg2 = std::abs(g2);
  if (absg2 == 0.) return g1;
  double absd2 = tanh(atanh(absg2)*2.);
  std::complex<double> d2 = absd2 / absg2 * g2;

  double d2sq = absd2*absd2;
  double x = (1.-std::sqrt(1.-d2sq))/d2sq;
  std::complex<double> d3 = d1+d2*(1.+x*(std::real(d1)*d2-std::real(d2)*d1));
  d3 /= 1. + std::real(d1*std::conj(d2));
  double absd3 = std::abs(d3);
  double absg3 = tanh(atanh(absd3)/2.);
  std::complex<double> g3 = absg3/absd3*d3;
  xdbg<<"GammaAdd: "<<g1<<" + "<<g2<<" = "<<g3<<std::endl;
  return g3;
}

bool Ellipse::DoMeasure(const std::vector<std::vector<Pixel> >& pix,
    const std::vector<BVec>* psf,
    int order, double sigma, bool use_integ, tmv::Matrix<double>* cov,
    BVec* bret, tmv::Matrix<double>* bcov)
{
  timeval tp;
  double t1=0.,t2=0.;

  // The non-linear solver is pretty sensitive to having a good
  // initial estimate.  So start with a simple estimate.
  if (order > 3) {
    //if (!DoMeasure(pix,psf,order-2,sigma,use_integ)) return false;
    DoMeasure(pix,psf,order-2,sigma,use_integ);
  }

  xdbg<<"Start DoMeasure: order = "<<order<<", psf = "<<bool(psf)<<std::endl;
  for(size_t i=0;i<pix.size();i++) xdbg<<"npix["<<i<<"] = "<<pix[i].size()<<std::endl;

  std::auto_ptr<BaseEllipseSolver> solver;

  tmv::Vector<double> x(5,0.);
  x(0) = std::real(cen);  x(1) = std::imag(cen); 
  x(2) = std::real(gamma); x(3) = std::imag(gamma);
  x(4) = mu;
  tmv::Vector<double> xinit = x;
  tmv::Vector<double> f(5);

  if (use_integ && order <= 3) {
    if (dotimings) {
      gettimeofday(&tp,0);
      t1 = tp.tv_sec + tp.tv_usec/1.e6;
    }
    // First we make the approximations the the galaxy is well-sampled
    // and has uniform variance.
    // Also, we start by just fitting the centroid.
    if (!fixcen && (!fixgam || !fixmu)) {
#ifdef _OPENMP
//#pragma omp critical
#endif
      {
	if (psf)
	  // pixscale doesn't really matter unless we want an accurate B in
	  // the end, so just use 1 here.
	  solver.reset(new EllipseSolver2(pix,*psf,f_psf,order,sigma,1.,
		false,true,true));
	else
	  solver.reset(new EllipseSolver2(pix,order,sigma,1.,
		false,true,true));
      }
      xdbg<<"xinit (int fixgam) = "<<x<<std::endl;
      solver->method = NLSolver::Dogleg;
#ifdef NOTHROW
      solver->startwithch = false;
#endif
      solver->ftol = 1.e-3;
      if (XDEBUG) solver->nlout = dbgout;
      solver->delta0 = 0.01;
      solver->min_step = 1.e-12;
      solver->max_iter = 60;
      xdbg<<"Integrating solver, centroid only:\n";
      //if (XDEBUG) if (!solver->TestJ(x,f,dbgout,1.e-5)) exit(1);
      if (!(solver->Solve(x,f)))
      {
	dbg<<"failed integrating solver, centroid - x = "<<x<<std::endl;
	dbg<<"f = "<<f<<"  Norm(f) = "<<Norm(f)<<std::endl;
	dbg<<"b = "<<solver->GetB()<<std::endl;
	return false;
      }
      xdbg<<"Done: x (Integ, centroid only) = "<<x<<std::endl;
      xdbg<<"f = "<<f<<std::endl;
      xdbg<<"b = "<<solver->GetB()<<std::endl;
    }

    // Next allow the shear and/or mu to be fit as well:
#ifdef N_FLUX_ATTEMPTS
#ifdef _OPENMP
//#pragma omp critical
#endif
    {
      if (psf)
	// pixscale doesn't really matter unless we want an accurate B in
	// the end, so just use 1 here.
	solver.reset(
	    new EllipseSolver2(pix,*psf,f_psf,order,sigma,1.,
	      fixcen,fixgam,fixmu,true));
      else
	solver.reset(new EllipseSolver2(pix,order,sigma,1.,
	      fixcen,fixgam,fixmu,true));
    }
    xdbg<<"xinit (integ) = "<<x<<std::endl;
    solver->method = NLSolver::Hybrid;
#ifdef NOTHROW
    solver->startwithch = false;
    solver->hasdirecth = false;
#endif
    solver->tau = 1.;
    solver->ftol = 3.e-3;
    solver->gtol = 1.e-5;
    if (XDEBUG) solver->nlout = dbgout;
    solver->min_step = 1.e-12;
    solver->max_iter = 10;
    xdbg<<"Integrating solver:\n";
    for(int iter = 0; iter < N_FLUX_ATTEMPTS; iter++) {
      xdbg<<"Attempt #"<<iter<<std::endl;
      //if (XDEBUG) if (!solver->TestJ(x,f,dbgout,1.e-5)) exit(1);
      solver->Solve(x,f);
      xdbg<<"x => "<<x<<std::endl;
      xdbg<<"f => "<<f<<"  Norm(f) = "<<Norm(f)<<std::endl;
      xdbg<<"b => "<<solver->GetB()<<std::endl;
      if (Norm(f) < 1.e-2) break;
    } 
#endif
    // Repeat, but allow flux to change.
#ifdef _OPENMP
//#pragma omp critical
#endif
    {
      if (psf)
	solver.reset(
	    new EllipseSolver2(pix,*psf,f_psf,order,sigma,1.,
	      fixcen,fixgam,fixmu));
      else
	solver.reset(new EllipseSolver2(pix,order,sigma,1.,
	      fixcen,fixgam,fixmu));
    }
    xdbg<<"xinit (integ) = "<<x<<std::endl;
    solver->method = NLSolver::Dogleg;
#ifdef NOTHROW
      solver->startwithch = false;
#endif
    solver->ftol = 1.e-3;
    if (XDEBUG) solver->nlout = dbgout;
    solver->delta0 = 0.01;
    solver->min_step = 1.e-12;
    solver->max_iter = 10;
    xdbg<<"Final Integrating solver:\n";
    //if (XDEBUG) if (!solver->TestJ(x,f,dbgout,1.e-5)) exit(1);
    //solver->verbose = true;
    solver->Solve(x,f);
    xdbg<<"Done: x (Integrating) = "<<x<<std::endl;
    xdbg<<"f (integ) = "<<f<<std::endl;
    xdbg<<"b = "<<solver->GetB()<<std::endl;
    if (dotimings) {
      gettimeofday(&tp,0);
      t2 = tp.tv_sec + tp.tv_usec/1.e6;
      t_integ += t2-t1;
    }
  }

  // Finally, we should be close enough for the exact solver
  // to work.
#if 0
  // But just to be safe, fit each of centroid, gamma, mu separately first:
  if (dotimings) {
    gettimeofday(&tp,0);
    t1 = tp.tv_sec + tp.tv_usec/1.e6;
  }
  if (!fixcen) {
#ifdef _OPENMP
//#pragma omp critical
#endif
    {
      if (psf)
	solver.reset(new EllipseSolver(pix,*psf,f_psf,order,sigma,
	      false,true,true));
      else
	solver.reset(new EllipseSolver(pix,order,sigma,
	      false,true,true));
    }
    xdbg<<"xinit = "<<x<<std::endl;
    solver->method = NLSolver::Dogleg;
#ifdef NOTHROW
      solver->startwithch = false;
#endif
    solver->ftol = 1.e-3;
    if (XDEBUG) solver->nlout = dbgout;
    solver->delta0 = 0.01;
    solver->min_step = 1.e-12;
    solver->max_iter = 50;
    xdbg<<"ML solver centroid:\n";
    if (!(solver->Solve(x,f))) {
      dbg<<"ML solver, centroid, failed - x = "<<x<<std::endl;
      dbg<<"f = "<<f<<std::endl;
      dbg<<"b = "<<solver->GetB()<<std::endl;
      if (XDEBUG) solver->TestJ(x,f,dbgout,1.e-5);
      if (dotimings) {
	gettimeofday(&tp,0);
	t2 = tp.tv_sec + tp.tv_usec/1.e6;
	t_centroid += t2-t1;
      }
      return false;
    }
    xdbg<<"Done: x (ML centroid only) = "<<x<<std::endl;
    xdbg<<"f = "<<f<<std::endl;
    xdbg<<"b = "<<solver->GetB()<<std::endl;
  }
  if (dotimings) {
    gettimeofday(&tp,0);
    t2 = tp.tv_sec + tp.tv_usec/1.e6;
    t_centroid += t2-t1;
  }

  if (dotimings) {
    gettimeofday(&tp,0);
    t1 = tp.tv_sec + tp.tv_usec/1.e6;
  }
  if (!fixgam) {
    // start with an initial guess from the latest b vector:
    if (solver.get()) {
      BVec b1 = solver->GetB();
      std::complex<double> g1(b1[3]/b1[0]*sqrt(2.),-b1[4]/b1[0]*sqrt(2.));
      dbg<<"g1 = "<<g1<<std::endl;
      if (std::abs(g1) < 1.) {
	std::complex<double> g0(x[2],x[3]);
	std::complex<double> g = GammaAdd(g0,g1);
	dbg<<"Initial gamma = "<<g<<std::endl;
	x[2] = std::real(g);
	x[3] = std::imag(g);
      }
    }

#ifdef _OPENMP
//#pragma omp critical
#endif
    {
      if (psf)
	solver.reset(new EllipseSolver(pix,*psf,f_psf,order,sigma,
	      true,false,true));
      else
	solver.reset(new EllipseSolver(pix,order,sigma,
	      true,false,true));
    }
    xdbg<<"xinit = "<<x<<std::endl;
    solver->method = NLSolver::Dogleg;
#ifdef NOTHROW
      solver->startwithch = false;
#endif
    solver->ftol = 1.e-3;
    if (XDEBUG) solver->nlout = dbgout;
    solver->delta0 = 0.01;
    solver->min_step = 1.e-12;
    solver->max_iter = 50;
    xdbg<<"ML solver gamma:\n";
    if (!(solver->Solve(x,f))) {
      dbg<<"ML solver, gamma, failed - x = "<<x<<std::endl;
      dbg<<"f = "<<f<<std::endl;
      dbg<<"b = "<<solver->GetB()<<std::endl;
      //if (XDEBUG) solver->TestJ(x,f,dbgout,1.e-5);
      if (dotimings) {
	gettimeofday(&tp,0);
	t2 = tp.tv_sec + tp.tv_usec/1.e6;
	t_gamma += t2-t1;
      }
      return false;
    }
    xdbg<<"Done: x (ML gamma only) = "<<x<<std::endl;
    xdbg<<"f = "<<f<<std::endl;
    xdbg<<"b = "<<solver->GetB()<<std::endl;
  }
  if (dotimings) {
    gettimeofday(&tp,0);
    t2 = tp.tv_sec + tp.tv_usec/1.e6;
    t_centroid += t2-t1;
  }

  if (dotimings) {
    gettimeofday(&tp,0);
    t1 = tp.tv_sec + tp.tv_usec/1.e6;
  }
  if (!fixmu) {
#ifdef _OPENMP
//#pragma omp critical
#endif
    {
      if (psf)
	solver.reset(new EllipseSolver(pix,*psf,f_psf,order,sigma,
	      true,true,false));
      else
	solver.reset(new EllipseSolver(pix,order,sigma,
	      true,true,false));
    }
    xdbg<<"xinit = "<<x<<std::endl;
    solver->method = NLSolver::Dogleg;
#ifdef NOTHROW
      solver->startwithch = false;
#endif
    solver->ftol = 1.e-3;
    if (XDEBUG) solver->nlout = dbgout;
    solver->delta0 = 0.01;
    solver->min_step = 1.e-12;
    solver->max_iter = 50;
    xdbg<<"ML solver mu:\n";
    if (!(solver->Solve(x,f))) {
      dbg<<"ML solver, mu, failed - x = "<<x<<std::endl;
      dbg<<"f = "<<f<<std::endl;
      dbg<<"b = "<<solver->GetB()<<std::endl;
      //if (XDEBUG) solver->TestJ(x,f,dbgout,1.e-5);
      if (dotimings) {
	gettimeofday(&tp,0);
	t2 = tp.tv_sec + tp.tv_usec/1.e6;
	t_gamma += t2-t1;
      }
      return false;
    }
    xdbg<<"Done: x (ML mu only) = "<<x<<std::endl;
    xdbg<<"f = "<<f<<std::endl;
    xdbg<<"b = "<<solver->GetB()<<std::endl;
  }
  if (dotimings) {
    gettimeofday(&tp,0);
    t2 = tp.tv_sec + tp.tv_usec/1.e6;
    t_centroid += t2-t1;
  }
#endif

  // Now fit everything, but first time, try to maintain the flux level
#ifdef N_FLUX_ATTEMPTS
  if (dotimings) {
    gettimeofday(&tp,0);
    t1 = tp.tv_sec + tp.tv_usec/1.e6;
  }
#ifdef _OPENMP
//#pragma omp critical
#endif
  {
    if (psf)
      solver.reset(new EllipseSolver(pix,*psf,f_psf,order,sigma,
	    fixcen,fixgam,fixmu,true));
    else
      solver.reset(new EllipseSolver(pix,order,sigma,
	    fixcen,fixgam,fixmu,true));
  }
  xdbg<<"xinit = "<<x<<std::endl;
  solver->method = NLSolver::Hybrid;
#ifdef NOTHROW
  solver->startwithch = false;
  solver->hasdirecth = false;
#endif
  solver->tau = 1.;
  solver->ftol = 3.e-3;
  solver->gtol = 1.e-5;
  if (XDEBUG) solver->nlout = dbgout;
  solver->delta0 = 0.01;
  solver->max_iter = 10;
  solver->min_step = 1.e-12;
  xdbg<<"ML solver use flux:\n";
  for(int iter = 0; iter < N_FLUX_ATTEMPTS; iter++) {
    xdbg<<"Attempt #"<<iter<<std::endl;
    //if (XDEBUG) if (!solver->TestJ(x,f,dbgout,1.e-5)) exit(1);
    solver->Solve(x,f);
    xdbg<<"x => "<<x<<std::endl;
    xdbg<<"f => "<<f<<"  Norm(f) = "<<Norm(f)<<std::endl;
    xdbg<<"b = "<<solver->GetB()<<std::endl;
    if (Norm(f) < 1.e-2) break;
  }
  xdbg<<"Done: x (ML fixed flux) = "<<x<<std::endl;
  xdbg<<"f = "<<f<<std::endl;
  xdbg<<"b = "<<solver->GetB()<<std::endl;
  if (dotimings) {
    gettimeofday(&tp,0);
    t2 = tp.tv_sec + tp.tv_usec/1.e6;
    t_fixflux += t2-t1;
  }
#endif

  if (dotimings) {
    gettimeofday(&tp,0);
    t1 = tp.tv_sec + tp.tv_usec/1.e6;
  }
  // Finally allow the flux change as needed.
#ifdef _OPENMP
//#pragma omp critical
#endif
  {
    if (psf)
      solver.reset(new EllipseSolver(pix,*psf,f_psf,order,sigma,
	    fixcen,fixgam,fixmu));
    else
      solver.reset(new EllipseSolver(pix,order,sigma,
	    fixcen,fixgam,fixmu));
  }
  xdbg<<"xinit = "<<x<<std::endl;
  solver->method = NLSolver::Dogleg;
#ifdef NOTHROW
  solver->startwithch = false;
#endif
  solver->ftol = 1.e-8;
  solver->gtol = 1.e-25;
  if (XDEBUG) solver->nlout = dbgout;
  solver->delta0 = 0.01;
  solver->min_step = 1.e-15;
  solver->max_iter = 50;
  xdbg<<"Final solver:\n";
  for(int iter = 1;iter<=MAXITER;iter++) {
    //if (XDEBUG) if (!solver->TestJ(x,f,dbgout,1.e-6)) exit(1);
    if (solver->Solve(x,f,cov)) break;
    else if (iter == MAXITER) {
      if (dotimings) {
	gettimeofday(&tp,0);
	t2 = tp.tv_sec + tp.tv_usec/1.e6;
	t_final += t2-t1;
      }
      return false;
    }
    dbg<<"Final solver pass failed - x = "<<x<<std::endl;
    dbg<<"b = "<<solver->GetB()<<std::endl;
    //if (XDEBUG) if (!solver->TestJ(x,f,dbgout,1.e-5)) exit(1);
#if 1 
    // Try bumping it to get out of a local well and continue:
    BVec b1 = solver->GetB();
    if (!fixcen) {
      dbg<<"Current z = "<<std::complex<double>(x[0],x[1])<<std::endl;
      std::complex<double> z1(2.*b1[1],-2.*b1[2]);
      z1 /= b1[0] - b1[5];
      if (std::abs(z1) < 1.) {
	dbg<<"Add "<<z1<<std::endl;
	x[0] += std::real(z1); x[1] += std::imag(z1); 
	dbg<<"New z = "<<std::complex<double>(x[0],x[1])<<std::endl;
      }
    }
#endif
#if 1
    if (!fixgam) {
      std::complex<double> g0(x[2],x[3]);
      dbg<<"Current gamma = "<<g0<<std::endl;
      std::complex<double> g1(std::sqrt(2.)*b1[3],-std::sqrt(2.)*b1[4]);
      g1 /= b1[0] - (b1.GetOrder() >= 4 ? b1[14] : 0.);
      if (std::abs(g1) < 1.) {
	dbg<<"Add "<<g1<<std::endl;
	g0 = GammaAdd(g1,g0);
	x[2] = std::real(g0); x[3] = std::imag(g0); 
	dbg<<"New gamma = "<<std::complex<double>(x[2],x[3])<<std::endl;
      }
    }
#endif
#if 0
    if (!fixmu) {
      dbg<<"Current mu = "<<x[4]<<std::endl;
      double m1 = -b1[5];
      m1 /= b1[0] - (b1.GetOrder() >= 4 ? 2.*b1[14] : 0.);
      if (std::abs(m1) < 1.) {
	x[4] += m1; 
	dbg<<"New mu = "<<x[4]<<std::endl;
      }
    }
    dbg<<"New x = "<<x<<std::endl;
#endif

#ifdef _OPENMP
//#pragma omp critical
#endif
    {
      if (psf)
	solver.reset(new EllipseSolver(pix,*psf,f_psf,order,sigma,
	      fixcen,fixgam,true));
      else
	solver.reset(new EllipseSolver(pix,order,sigma,
	      fixcen,fixgam,true));
    }

    solver->method = NLSolver::Hybrid;
#ifdef NOTHROW
    solver->startwithch = false;
#endif
    solver->ftol = 1.e-8;
    solver->gtol = 1.e-10;
    if (XDEBUG) solver->nlout = dbgout;
    solver->delta0 = 0.01;
    solver->min_step = 1.e-15;
    solver->max_iter = 50;
  }
  xdbg<<"Done: Final solver pass successful: x = "<<x<<std::endl;
  xdbg<<"f = "<<f<<std::endl;
  xdbg<<"b = "<<solver->GetB()<<std::endl;
  if (cov) xdbg<<"cov = "<<*cov<<std::endl;
  if (XDEBUG && psf) {
    for(size_t k=0;k<pix.size();k++) {
      double sig_psf = (*psf)[k].GetSigma();
      double D = sigma*sigma / (sig_psf*sig_psf);
      xdbg<<"D = "<<D<<std::endl;
    }
  }
  if (dotimings) {
    gettimeofday(&tp,0);
    t2 = tp.tv_sec + tp.tv_usec/1.e6;
    t_final += t2-t1;
  }

  SetCen(std::complex<double>(x(0),x(1)));
  SetGamma(std::complex<double>(x(2),x(3))); 
  SetMu(x(4));
  if (bret) {
    *bret = solver->GetB();
    xdbg<<"bret = "<<*bret<<std::endl;
    if (bcov) solver->GetBCov(*bcov);
  }
  return true;
}

bool Ellipse::Measure(const std::vector<std::vector<Pixel> >& pix,
    const std::vector<BVec>& psf,
    int order, double sigma, bool use_integ, tmv::Matrix<double>* cov,
    BVec* bret, tmv::Matrix<double>* bcov)
{ 
  didstatus3output = false;
  return DoMeasure(pix,&psf,order,sigma,use_integ,cov,bret,bcov); 
}

bool Ellipse::Measure(const std::vector<std::vector<Pixel> >& pix,
    int order, double sigma, bool use_integ, tmv::Matrix<double>* cov,
    BVec* bret, tmv::Matrix<double>* bcov)
{
  didstatus3output = false;
  return DoMeasure(pix,0,order,sigma,use_integ,cov,bret,bcov);
}

void Ellipse::DoMeasureShapelet(const std::vector<std::vector<Pixel> >& pix,
    const std::vector<BVec>* psf, BVec& b, tmv::Matrix<double>* bcov) const
{
  xdbg<<"Start MeasureShapelet: order = "<<b.GetOrder()<<std::endl;
  xdbg<<"sigma = "<<b.GetSigma()<<std::endl;
  xdbg<<"el = "<<*this<<std::endl;
  
  // ( u )' = exp(-mu)/sqrt(1-gsq) ( 1-g1  -g2  ) ( u-uc )
  // ( v )                         ( -g2   1+g1 ) ( v-vc )
  // 
  // z' = u' + I v' 
  //    = exp(-mu)/sqrt(1-gsq)
  //             [ (1-g1)(u-uc)-g2(v-vc)-Ig2(u-uc)+I(1+g1)(v-vc) ]
  //    = exp(-mu)/sqrt(1-gsq) [ z-zc - g1(z-zc)* -Ig2(z-zc)* ]
  //    = exp(-mu)/sqrt(1-gsq) ( z-zc - g (z-zc)* )

  int order = b.GetOrder();
  double sigma = b.GetSigma();
  double gsq = std::norm(gamma);
  Assert(gsq < 1.);

  double m = exp(-mu)/sqrt(1.-gsq);

  size_t ntot = 0;
  for(size_t i=0;i<pix.size();i++) ntot += pix[i].size();
  dbg<<"ntot = "<<ntot<<" in "<<pix.size()<<" images\n";

  tmv::Vector<double> I(ntot);
  tmv::Vector<double> W(ntot);
  tmv::Vector<std::complex<double> > Z(ntot);

  for(size_t k=0,n=0;k<pix.size();k++) {
    double sigma_obs = psf ?
      sqrt(pow(sigma,2)+f_psf*pow((*psf)[k].GetSigma(),2)) :
      sigma;
    for(size_t i=0;i<pix[k].size();i++,n++) {
      I(n) = pix[k][i].I*pix[k][i].wt;
      W(n) = pix[k][i].wt;
      std::complex<double> z1 = pix[k][i].z;
      std::complex<double> z2 = m*(z1-cen) - m*gamma*conj(z1-cen);
      Z(n) = z2 / sigma_obs;
    }
  }

  size_t nsize = (order+1)*(order+2)/2;
  tmv::Matrix<double> A(ntot,nsize);
  MakePsi(A.View(),Z,order,&W);
  xxdbg<<"after makepsi\n";
  xxdbg<<"A.rows(0,10) = "<<A.Rows(0,10);

  if (psf) {
    for(size_t k=0,n=0,nx;k<pix.size();k++,n=nx) {
      size_t psfsize = (*psf)[k].size();
      BVec newpsf = (*psf)[k];
      tmv::Matrix<double> S(psfsize,psfsize);
      GTransform(gamma,newpsf.GetOrder(),S);
      newpsf = S * (*psf)[k];
      tmv::Matrix<double> D(psfsize,psfsize);
      MuTransform(mu,newpsf.GetOrder(),D);
      newpsf = D * newpsf;
      tmv::Matrix<double> C(nsize,nsize);
      PSFConvolve(newpsf,b.GetOrder(),b.GetSigma(),C);
      nx = n+pix[k].size();
      A.Rows(n,nx) *= C;
    }
    xxdbg<<"after psf correction\n";
  }
  b = I/A;
  if (bcov) {
    tmv::UpperTriMatrix<double> Rinv = A.QRD().GetR();
    Rinv.InvertSelf();
    *bcov = Rinv * Rinv.Transpose();
  }
}

void Ellipse::MeasureShapelet(const std::vector<std::vector<Pixel> >& pix,
    const std::vector<BVec>& psf, BVec& b) const
{ DoMeasureShapelet(pix,&psf,b); }

void Ellipse::MeasureShapelet(const std::vector<std::vector<Pixel> >& pix,
    BVec& b) const
{ DoMeasureShapelet(pix,0,b); }

