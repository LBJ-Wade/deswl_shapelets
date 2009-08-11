
#include "NLSolver.h"
#include <cstdlib>
#include "TMV.h"
#include <iostream>

#define DOLINEAR
#define DONONLINEAR
#define DOPOWELL
#define DOROSENBROCK

#define ROSENBROCKDIRECTH

#define SETMETHOD(solver) \
  solver.method = NLSolver::Dogleg; \
  solver.nlout = &std::cout; solver.verbose = false;

#define USEBEST

const size_t N=2;
const size_t M=5;

class Linear : public NLSolver 
{
  public :

    Linear() : b(M,0.), A(M,N,0.)
    {
      tmv::Vector<double> x(N,1.);
      for(size_t i=0;i<std::min(M,N);i++) {
	A(i,i) = 0.5*i - 6.;
	if (i>0 && 2*i<M) A(2*i,i) = 9.+i;
	if (i>0 && 2*i<M) A(2*i-1,i) = 1.-2*i;
	if (i>0 && 2*i<N) A(i,2*i) = -3.-i;
	if (i>0 && 2*i<N) A(i,2*i-1) = -12.+3*i;
      }
      A.col(0).AddToAll(2.);
      A.row(0).AddToAll(3.);
      b = A*x;
      //b.AddToAll(100.);
    }

    void F(const tmv::Vector<double>& x, tmv::Vector<double>& f) const
    { f = A*x-b; }

    void J(const tmv::Vector<double>& , const tmv::Vector<double>& ,
	tmv::Matrix<double>& j) const
    { j = A; }
  
  protected:

    tmv::Vector<double> b;
    tmv::Matrix<double> A;
};

const size_t MM=50;
class NonLinear : public NLSolver 
{
  public : 

    NonLinear(double sigma) : b(MM), t(MM)
    {
      //double xa[] = {-4.,-5.,4.,-4.};
      tmv::Vector<double> x(4);
      x = tmv::ListInit, -4., -5., 4., -4.;

      for(size_t i=0;i<MM;i++) {
	t(i) = (i+0.5)*(0.5/MM);
	b(i) = func(t(i),x) + sigma*(double(rand())/RAND_MAX -0.5);
      }
      //std::cout<<"t = "<<t<<std::endl;
      //std::cout<<"b = "<<b<<std::endl;
    }

    double func(double t, const tmv::Vector<double>& x) const
    {
      return x(2)*exp(x(0)*t) + x(3)*exp(x(1)*t);
    }

    void F(const tmv::Vector<double>& x, tmv::Vector<double>& f) const
    { for(size_t i=0;i<MM;i++) f(i) = func(t(i),x)-b(i); }

    void J(const tmv::Vector<double>& x, const tmv::Vector<double>& ,
	tmv::Matrix<double>& j) const
    {
      for(size_t i=0;i<MM;i++) {
	j(i,2) = exp(x(0)*t(i));
	j(i,3) = exp(x(1)*t(i));
	j(i,0) = x(2)*t(i)*j(i,2);
	j(i,1) = x(3)*t(i)*j(i,3);
      }
    }
  
  protected:

    tmv::Vector<double> b;
    tmv::Vector<double> t;
};

class NonLinear2 : public NLSolver 
{
  public : 

    NonLinear2(double sigma) : b(MM), t(MM), A(MM,2), y(2)
    {
      double xx[2] = {-4.,-5.};
      double yy[2] = {4.,-4.};

      for(size_t i=0;i<MM;i++) {
	t(i) = (i+0.5)*(0.5/MM);
	b(i) = yy[0]*exp(xx[0]*t(i)) + yy[1]*exp(xx[1]*t(i)) +
	  sigma*(2.*double(rand())/RAND_MAX - 1.);
      }
      A.SaveDiv();
    }

    double func(double t, const tmv::Vector<double>& x, 
	const tmv::Vector<double>& y) const
    {
      return y(0)*exp(x(0)*t) + y(1)*exp(x(1)*t);
    }

    void F(const tmv::Vector<double>& x, tmv::Vector<double>& f) const
    {
      for(size_t i=0;i<MM;i++) {
	if (x(0)*t(i) > 20.) A(i,0) = 5.e8;
	else if (x(0)*t(i) < -20.) A(i,0) = 2.e-9;
	else A(i,0) = exp(x(0)*t(i));
	if (x(1)*t(i) > 20.) A(i,1) = 5.e8;
	else if (x(1)*t(i) < -20.) A(i,1) = 2.e-9;
	else A(i,1) = exp(x(1)*t(i));
      }
      A.ReSetDiv();
      y = b/A;
      f = A*y-b;
    }

    void J(const tmv::Vector<double>& , const tmv::Vector<double>& f,
	tmv::Matrix<double>& j) const
    {
      // At A y = At b
      // df = A dy + dA y
      // (A+dA)t (A+dA) (y+dy) = (A+dA)t b
      // dAt A y + At dA y + At A dy = dAt b
      // At A dy = -dAt (Ay-b) - At dA y
      // At A = Rt R where A = QR is the QR decomposition of A
      // dy = R^-1 Rt^-1 (-dAt f - At dA y)
      // df = A dy + dA y
      //    = A R^-1 Rt^-1 (-dAt f - At dA y) + dA y
      tmv::ConstUpperTriMatrixView<double> R = A.QRD().GetR();
      tmv::Vector<double> dA = DiagMatrixViewOf(t)*A.col(0);
      tmv::Vector<double> dy = -A.Transpose()*dA*y(0);
      dy(0) -= dA*f;
      dy /= R.Transpose();
      dy /= R;
      j.col(0) = A*dy + dA*y(0);
      dA = DiagMatrixViewOf(t)*A.col(1);
      dy = -A.Transpose()*dA*y(1);
      dy(1) -= dA*f;
      dy /= R.Transpose();
      dy /= R;
      j.col(1) = A*dy + dA*y(1);
    }

    const tmv::Vector<double>& GetY() const { return y; }
  
  protected:

    tmv::Vector<double> b;
    tmv::Vector<double> t;
    mutable tmv::Matrix<double> A;
    mutable tmv::Vector<double> y;
};

class Powell : public NLSolver 
{
  public :

    Powell() {}

    void F(const tmv::Vector<double>& x, tmv::Vector<double>& f) const
    {
      f(0) = x(0);
      f(1) = 10.*x(0)/(0.1+x(0)) + 2*x(1)*x(1);
    }
    void J(const tmv::Vector<double>& x, const tmv::Vector<double>& ,
	tmv::Matrix<double>& j) const
    {
      j(0,0) = 1.;
      j(0,1) = 0.;
      j(1,0) = 1./(0.1+x(0))/(0.1+x(0));
      j(1,1) = 4*x(1);
    }
};

class Rosen : public NLSolver 
{
  public :

    Rosen(double _lam) : lambda(_lam) {}

    void F(const tmv::Vector<double>& x, tmv::Vector<double>& f) const
    {
      f(0) = 10.*(x(1) - x(0)*x(0));
      f(1) = 1.-x(0);
      f(2) = lambda;
    }
    void J(const tmv::Vector<double>& x, const tmv::Vector<double>& ,
	tmv::Matrix<double>& j) const
    {
      j(0,0) = -20.*x(0);
      j(1,0) = -1.;
      j(2,0) = 0.;
      j(0,1) = 10.;
      j(1,1) = 0.;
      j(2,1) = 0.;
    }
    
#ifdef ROSENBROCKDIRECTH
    void H(const tmv::Vector<double>& , const tmv::Vector<double>& f,
	const tmv::Matrix<double>& j, tmv::SymMatrix<double>& h) const
    {
      // H = JT J + Sum_k f(k) d^2f/(dxi dxj)
      h = j.Transpose() * j;
      h(0,0) += -20.*f(0);
    }
#endif
    


  private :

    double lambda;
};

int main() try 
{
#ifdef DOLINEAR
  {
    std::cout<<"Test Linear:\n";
    Linear lin;
    SETMETHOD(lin);
#ifdef USEBEST
    lin.method = NLSolver::Newton;
#endif
    tmv::Vector<double> x(N,0.);
    tmv::Vector<double> f(M);

    //std::cout<<"xinit = "<<x<<std::endl;
    //if (!(lin.TestJ(x,f,&std::cout))) return 1;
    if (lin.Solve(x,f))
      std::cout<<"Success: \n";
    else 
      std::cout<<"Failure: \n";
    std::cout<<"xfinal = "<<x<<std::endl;
    //std::cout<<"f = "<<f<<std::endl;
    std::cout<<"Norm(f) = "<<Norm(f)<<std::endl;
  }
#endif

#ifdef DONONLINEAR
  {
    std::cout<<"Test NonLinear:\n";
    NonLinear nlin(1.e-2);
    SETMETHOD(nlin);
#ifdef USEBEST
    nlin.method = NLSolver::Newton;
#endif
    nlin.min_step = 1.e-12;
    //double x0a[] = {-1., -2., 1., -1.};
    tmv::Vector<double> x(4);
    x = tmv::ListInit, -1., -2., 1., -1.;
    tmv::Vector<double> f(MM);

    //std::cout<<"xinit = "<<x<<std::endl;
    //if (!(nlin.TestJ(x,f,&std::cout))) return 1;
    if (nlin.Solve(x,f))
      std::cout<<"Success: \n";
    else 
      std::cout<<"Failure: \n";
    std::cout<<"xfinal = "<<x<<std::endl;
    //std::cout<<"f = "<<f<<std::endl;
    std::cout<<"Norm(f) = "<<Norm(f)<<std::endl;
  }
  {
    std::cout<<"Test NonLinear2:\n";
    NonLinear2 nlin(1.e-2);
    SETMETHOD(nlin);
#ifdef USEBEST
    nlin.method = NLSolver::LM;
#endif
    //double x0a[] = {-1., -2.};
    tmv::Vector<double> x(2);
    x = tmv::ListInit, -1., -2.;
    tmv::Vector<double> f(MM);

    //std::cout<<"xinit = "<<x<<std::endl;
    //if (!(nlin.TestJ(x,f,&std::cout))) return 1;
    if (nlin.Solve(x,f))
      std::cout<<"Success: \n";
    else 
      std::cout<<"Failure: \n";
    tmv::Vector<double> xtot(4);
    xtot.SubVector(0,2) = x;
    xtot.SubVector(2,4) = nlin.GetY();
    std::cout<<"xfinal = "<<xtot<<std::endl;
    //std::cout<<"f = "<<f<<std::endl;
    std::cout<<"Norm(f) = "<<Norm(f)<<std::endl;
  }
#endif

#ifdef DOPOWELL
  {
    std::cout<<"Test Powell\n";
    Powell pow;
    pow.ftol = 1.e-20;
    pow.gtol = 1.e-15;
    pow.min_step = 1.e-15;
    pow.tau = 1.;
    SETMETHOD(pow);
#ifdef USEBEST
    pow.method = NLSolver::Dogleg;
#endif
    tmv::Vector<double> x(2); x(0) = 3.; x(1) = 1.;
    tmv::Vector<double> f(2);

    //std::cout<<"xinit = "<<x<<std::endl;
    if (pow.Solve(x,f))
      std::cout<<"Success: \n";
    else 
      std::cout<<"Failure: \n";
    std::cout<<"xfinal = "<<x<<std::endl;
    //std::cout<<"f = "<<f<<std::endl;
    std::cout<<"Norm(f) = "<<Norm(f)<<std::endl;
  }
#endif

#ifdef DOROSENBROCK
  {
    std::cout<<"Test Rosenbrock\n";
    const int ntest = 4;
    double lambda[ntest] = {0., 1.e-3, 1., 1.e4};
    for(int i=0;i<ntest;i++) {
      Rosen ros(lambda[i]);
      ros.ftol = 1.e-20;
      ros.gtol = 1.e-10;
      ros.min_step = 1.e-14;
      ros.tau = 1.e-3;
      SETMETHOD(ros);
#ifdef USEBEST
      if (i == 0) ros.method = NLSolver::Dogleg;
      else ros.method = NLSolver::Hybrid;
#endif
      tmv::Vector<double> x(2); x(0) = -1.2; x(1) = 1.;
      tmv::Vector<double> x0(2); x0(0) = 1.; x0(1) = 1.;
      tmv::Vector<double> f(3);

      //std::cout<<"xinit = "<<x<<std::endl;
      if (ros.Solve(x,f))
	std::cout<<"Success: \n";
      else 
	std::cout<<"Failure: \n";
      std::cout<<"xfinal = "<<x<<std::endl;
      //std::cout<<"Norm(xf-x0) = "<<Norm(x0-x)<<std::endl;
      //std::cout<<"f = "<<f<<std::endl;
      std::cout<<"Norm(f) = "<<Norm(f)<<std::endl;
    }
  }
#endif
  return 0;
} 
catch (tmv::Error& e)
{
  std::cout<<"Caught "<<e<<std::endl;
  return 1;
}
catch (std::exception& e)
{
  std::cout<<"Caught "<<e.what()<<std::endl;
  return 1;
}
catch (...)
{
  std::cout<<"Caught unknown error\n";
  return 1;
}
