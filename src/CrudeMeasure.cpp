
#include "Ellipse.h"
#include "dbg.h"
#include "NLSolver.h"
#include <vector>
#include "TMV.h"

class CrudeSolver : public NLSolver 
{
  public :
    CrudeSolver(const PixelList& pix, double _sigma, double _I1,
	tmv::Vector<double>& _xinit);

    void F(const tmv::Vector<double>& x, tmv::Vector<double>& f) const;
    void J(const tmv::Vector<double>& x, const tmv::Vector<double>& f,
	tmv::Matrix<double>& df) const;

  private : 
    double sigma;
    tmv::Vector<double> I;
    tmv::Vector<std::complex<double> > Z;
    tmv::Vector<double> W;
    mutable tmv::Vector<std::complex<double> > Z1;
    mutable tmv::Vector<double> E;
    mutable tmv::Vector<double> Rsq;
    mutable tmv::Vector<double> f1;
    double I1;
    tmv::Vector<double>& xinit;
};

CrudeSolver::CrudeSolver(const PixelList& pix,
    double _sigma, double _I1, tmv::Vector<double>& _xinit) :
  sigma(_sigma), I(pix.size()), Z(pix.size()), W(pix.size()), 
  Z1(pix.size()), E(pix.size()), Rsq(pix.size()), f1(pix.size()), I1(_I1),
  xinit(_xinit)
{
  for(size_t i=0;i<pix.size();i++) {
    Z(i) = pix[i].z;
    I(i) = pix[i].I;
    double mask = exp(-std::norm(pix[i].z/sigma)/2.);
    W(i) = pix[i].wt*mask;
  }
}

void CrudeSolver::F(const tmv::Vector<double>& params,
    tmv::Vector<double>& f) const
{
  xdbg<<"Start F\n";
  xxdbg<<"params = "<<params<<std::endl;
  xxdbg<<"sigma = "<<sigma<<std::endl;
  xxdbg<<"I1 = "<<I1<<std::endl;
  Assert(params.size() == 4);
  Assert(f.size() == Z.size());

  if (NormInf((params-xinit).SubVector(0,3)) > 2.) 
  { f = 2.e10*f1; f.AddToAll(1.); return; }
  if (params(3) < 0.) { f = 2.e10*f1; f.AddToAll(1.); return; }

  std::complex<double> zc(params[0],params[1]);
  double mu = params[2];
  double I0 = I1 * params[3];
  xxdbg<<"I0 = "<<I0<<std::endl;

  double m0 = exp(-mu)/sigma;
  // z' = m*(z-zc)
  //    = m*z - m*zc
  // x' = m x - m xc
  // y' = m y - m yc
  Z1 = m0*Z;
  Z1.AddToAll(-m0*zc);

  xxdbg<<"In CrudeMeasure::F\n";
  xxdbg<<"sigma = "<<sigma<<std::endl;
  xxdbg<<"mu = "<<mu<<std::endl;
  xdbg<<"m0 = "<<m0<<std::endl;
  for(size_t i=0;i<Z.size();i++) {
    double rsq = std::norm(Z1[i]);
    Rsq[i] = rsq;
    E[i] = exp(-rsq/2.);
    if (rsq < 2.) {
      xxdbg<<Z[i]-zc<<"  "<<Z1[i]<<"  "<<I[i]<<"  "<<I0*E[i]<<std::endl;
    }
  }
  xdbg<<"After Set elements\n";
  f1 = I - I0*E;
  xdbg<<"After f1 = I-I0*E\n";
  ElementProd(1.,W,f1.View());
  xdbg<<"After ElementProd\n";
  f = f1;

  xdbg<<"Done F\n";
  xdbg<<"norm(f) = "<<Norm(f)<<std::endl;
}

void CrudeSolver::J(const tmv::Vector<double>& params, const tmv::Vector<double>& f,
    tmv::Matrix<double>& df) const
{
  xdbg<<"Start J\n";
  xxdbg<<"params = "<<params<<std::endl;
  Assert(params.size() == 4);
  Assert(f.size() == Z.size());
  Assert(df.colsize() == Z.size());
  Assert(df.rowsize() == 4);

  //double xc = params[0];
  //double yc = params[1];
  double mu = params[2];
  double I0 = I1 * params[3];
  double m0 = exp(-mu)/sigma;

  // fi = Wi * (Ii - I0 Ei)
  // dfi/dI0 = -Wi Ei
  // dfi/dmu = Wi I0 Ei (x1 dx1/dxc + y1 dy1/dxc) 
  // likewise for the other 4
  //
  // dx1/dxc = -m 		dy1/dxc = 0
  // dx1/dyc = 0 		dy1/dyc = -m
  // dx1/dmu = -x1		dy1/dmu = -y1
  //

  df.col(0) = -m0 * Z1.Real();
  df.col(1) = -m0 * Z1.Imag();
  df.col(2) = -Rsq;
  df.Cols(0,3) = I0 * DiagMatrixViewOf(E) * df.Cols(0,3);
  df.col(3) = -I1 * E;
  df = DiagMatrixViewOf(W) * df;

  xdbg<<"Done J\n";
  xxdbg<<"J = "<<df<<std::endl;
}

void Ellipse::CrudeMeasure(const PixelList& pix, double sigma)
{
  // We use as our initial estimate the exact value for a 
  // well-sampled, uniform-variance, undistorted Gaussian intensity pattern
  // with no PSF convolution.
  //
  // That is, we assume the model:
  //
  // I(x,y) = I0 exp( -(|z-zc|^2/ (2 sigma'^2) )
  // where zz = (z-zc) 
  // and sigma' = exp(mu) sigma
  // 
  xdbg<<"Current centroid = "<<cen<<std::endl;
  xdbg<<"Current mu = "<<mu<<std::endl;
  xdbg<<"sigma = "<<sigma<<std::endl;

#if 1
  // With a weight of exp(-|z|^2/(2 sigma^2)),
  // the weighted moments of this function are:
  // Iz/I = zc / (1+exp(2mu))
  // Irr/I = ( |zc|^2 + 2 exp(2mu) sigma^2 ) / (1+exp(2mu))

  std::complex<double> Iz = 0.;
  double I = 0.;
  double sig2 = sigma * exp(mu);
  for(size_t i=0;i<pix.size();i++) {
    double wt = exp(-std::norm((pix[i].z-cen)/sig2)/2.);
    Iz += wt * pix[i].I * (pix[i].z-cen);
    I += wt * pix[i].I;
    if (std::abs(pix[i].z-cen) < 2.) 
      xdbg<<pix[i].z<<"  "<<pix[i].I<<std::endl;
  }
  xdbg<<"Iz = "<<Iz<<", I = "<<I<<std::endl;

  // If I <= 0 then this isn't going to work.  Just return and hope
  // the regular measure method might do better.
  if (!(I > 0.)) return;

  std::complex<double> zc = Iz / I;
  // If zc is more than 2, something is probably wrong, so abort now.
  if (!(std::abs(zc) < 2.)) return;
  
  xdbg<<"Initial offset to centroid = "<<zc<<std::endl;
  zc += cen;
  xdbg<<"zc = "<<zc<<std::endl;

  double Irr = 0.;
  double W = 0.;
  I = 0.;
  Iz = 0.;
  for(size_t i=0;i<pix.size();i++) {
    double wt = exp(-std::norm((pix[i].z-zc)/sig2)/2.);
    Iz += wt * pix[i].I * (pix[i].z-zc);
    Irr += wt * pix[i].I * std::norm(pix[i].z-zc);
    I += wt * pix[i].I;
    W += wt;
#if 0
    if (std::abs(wt * pix[i].I * std::norm(pix[i].z-zc)) > 1.) {
      dbg<<"pix.I = "<<pix[i].I<<std::endl;
      dbg<<"pix.z = "<<pix[i].z<<std::endl;
      dbg<<"pix.wt = "<<pix[i].wt<<std::endl;
      dbg<<"wt = "<<wt<<std::endl;
    }
#endif
  }
  xdbg<<"Iz = "<<Iz<<", Irr = "<<Irr<<", I = "<<I<<", W = "<<W<<std::endl;

  std::complex<double> zc1 = Iz/I;
  double S = Irr/I - norm(zc1);
  // S is now 2 exp(2mu) sigma^2 / (1 + exp(2mu))
  xdbg<<"S = "<<S<<std::endl;
  double exp2mu = S / 2. / (sig2*sig2);
  xdbg<<"exp2mu/(1+exp2mu) = "<<exp2mu<<std::endl;
  // It's actually exp(2mu) / (1+exp(2mu)) at this point
  if (exp2mu < 0.2)
    exp2mu = 0.25;
  else if (exp2mu < 0.8)
    exp2mu = 1./(1./exp2mu-1.); // Now it is really exp(2mu)
  else 
    // The above formula is unstable, and probably inappropriate, since
    // we probably have a failure of our model approximation.
    // So just multiply it by 5 -- the correct factor for exp2mu = 0.8
    exp2mu *= 5.;
  xdbg<<"exp2mu = "<<exp2mu<<std::endl;
  if (exp2mu <= 0.) exp2mu = 1.;

  double m = log(exp2mu)/2.;
  xdbg<<"mu = "<<m<<std::endl;

  zc += zc1 * (1.+exp2mu);
  m += mu;

  xdbg<<"Approx cen = "<<zc<<std::endl;
  xdbg<<"Approx mu = "<<m<<std::endl;
  
  // I/W = I0 exp(2mu) / (1+exp(2mu)) * exp(-|zc1|^2/2sigma^2*(1+exp(2mu)))
  double I0 = (I/W)*(1.+exp2mu)/exp2mu /
    exp(-norm(zc1)*(1.+exp2mu)/(2.*sig2*sig2));
#else
  std::complex<double> zc = cen;
  double m = mu;

  double model = 0.;
  double obs = 0.;
  double minx = 1.e100, maxx = -1.e100, miny = 1.e100, maxy = -1.e100;
  for(size_t i=0;i<pix.size();i++) {
    if (real(pix[i].z) < minx) minx = real(pix[i].z);
    if (real(pix[i].z) > maxx) maxx = real(pix[i].z);
    if (imag(pix[i].z) < miny) miny = imag(pix[i].z);
    if (imag(pix[i].z) > maxy) maxy = imag(pix[i].z);
    double wt = exp(-std::norm(pix[i].z/sigma)/2.);
    model += mask;
    obs += pix[i].I * mask;
  }
  double pixarea = (maxx-minx)*(maxy-miny)/pix.size();
  xdbg<<"pixarea = "<<pixarea<<std::endl; 
  xdbg<<"obs = "<<obs<<std::endl;
  xdbg<<"model = "<<model<<std::endl;
  double I0 = 2.*obs / model;
#endif

  xdbg<<"Initial I0 estimate = "<<I0<<std::endl;
  if (I0 < 1.e-6) {
    xdbg<<"Warning: small or negative I0: "<<I0<<" -- Use 1.0\n";
    I0 = 1.;
  }

  //if (std::abs(zc) > 1.0) zc /= std::abs(zc);
  //if (std::abs(m) > 1.0) m /= std::abs(m);
  tmv::Vector<double> x(4);
  x[0] = std::real(zc); x[1] = std::imag(zc); 
  x[2] = m; 
  x[3] = 1.;
  tmv::Vector<double> f(pix.size());
  CrudeSolver s(pix,sigma,I0,x);

  s.method = NLSolver::Hybrid;
  s.ftol = 1.e-4;
  s.gtol = 1.e-8;
  s.min_step = 1.e-15;
  s.tau = 1.0;
  s.delta0 = 0.05;
#ifdef __PGI
  s.startwithch = false;
#endif
  if (XDEBUG) s.nlout = dbgout;
  //s.verbose = true;
  xdbg<<"Before CrudeSolver: x = "<<x<<std::endl;
  //s.TestJ(x,f,dbgout);
  //if (!s.TestJ(x,f,dbgout,1.e-5)) exit(1);
  s.Solve(x,f);
  //s.TestJ(x,f,dbgout);
  //if (!s.TestJ(x,f,dbgout,1.e-5)) exit(1);
  xdbg<<"After CrudeSolver: x = "<<x<<std::endl;
  s.ftol = 1.e-4 * (1.e-4 + std::abs(x[3]));
  s.Solve(x,f);
  xdbg<<"After 2nd CrudeSolver: x = "<<x<<std::endl;

  std::complex<double> newcen(x[0],x[1]);
  double newmu = x[2];

  if (std::abs(newcen-cen) > 2.) {
    dbg<<"Warning: large centroid shift in CrudeMeasure\n";
    dbg<<"Old centroid = "<<cen<<", new centroid = "<<newcen<<std::endl;
    newcen = cen + 2.*(newcen - cen)/std::abs(newcen-cen);
    dbg<<"Scaling back to "<<newcen<<std::endl;
  }

  if (std::abs(newmu-mu) > 2.) {
    dbg<<"Warning: large scale change in CrudeMeasure\n";
    dbg<<"Old mu = "<<mu<<", new mu = "<<newmu<<std::endl;
    newmu = mu + 2.*(newmu - mu)/std::abs(newmu-mu);
    dbg<<"Scaling back to "<<newmu<<std::endl;
  }
  
  xdbg<<"Crude cen = "<<newcen<<std::endl;
  xdbg<<"Crude mu = "<<newmu<<std::endl;

  if (!fixcen) cen = newcen;
  if (!fixmu) mu = newmu;
}

void Ellipse::CrudeMeasure(
    const std::vector<PixelList>& pix, double sigma)
{
  int n = 0;
  for(size_t i=0;i<pix.size();i++) n += pix[i].size();
  PixelList allpix(n);
  for(size_t i=0,k=0;i<pix.size();i++)
    for(size_t j=0;j<pix[i].size();j++,k++)
      allpix[k] = pix[i][j];
  CrudeMeasure(allpix,sigma);
}

void Ellipse::PeakCentroid(const PixelList& pix, double maxr)
{
  double peakI = 0.;
  std::complex<double> peakz = 0.;
  for(size_t i=0;i<pix.size();i++) if (std::abs(pix[i].z) < maxr) {
    if (pix[i].I > peakI) { 
      peakz = pix[i].z;
      peakI = pix[i].I;
    }
  }
  cen = peakz;
}

