#include "TMV.h"
#include "TMV_Sym.h"

// The class NLSolver is designed as a base class.  Therefore, you should
// define a class with your particular function as a derived class of
// NLSolver:
//
// class MyNLFunction : public NLSolver { ... }
//
// You need to define the virtual functions F and J.
// J = dF/dx is a MxN matrix where M is the number of elements
// in the F vector and N is the number of elements in the x vector.
// This allows you to put any kind of global information that
// might be needed for evaluating the function into the class.
//
// Then:
//
// MyNLSolver nls;
// tmv::Vector<double> x(n);
// tmv::Vector<double> f(m);
// x = [initial guess]
// bool success = nls.Solve(x,f);
// [ if success, x is now a solution to the function. ]
// [ f is returned as F(x), regardless of success. ]
//
// Note: for some functions, the calculation of J is 
// partially degenerate with the calculations needed for F.
// Therefore, we do two things to help out.
// First, the function call for J includes the value for
// f that has previously been calculated for that x.
// Second, we guarantee that every call to J will be for the 
// same x argument as the previous call to F.  This allows
// the class to potentially store useful auxiliary information
// which might make the J calculation easier.
//
// If m > n, then the solution sought is the 
// least squares solution which minimizes 
// Q = 1/2 Sum |F_i|^2
//
// If you want the covariance matrix of the variables at the 
// solution location, you can do:
// tmv::Matrix<double> invcov;
// nls.Solve(x,f,&invcov);
// tmv::Matrix<double> cov = invcov.Inverse();
//
// There are 5 methods implemented here, all of which are based on
// "Methods for Non-Linear Least Squares Problems", 2nd edition,
// April, 2004, K. Madsen, H.B. Nielsen, O. Tingleff,
// Informatics and Mathematical Modelling, Technical University of Denmark.
//
// The Newton method is a simple method which can be very fast for cases
// where the solution is at f=0, but it is also more liable to fail 
// than the below methods.
//
// The Hybrid method is usually the best choice when the number of 
// equations (m) > the number of variables (n).
//
// The Dogleg method is usually best when n == m and you are expecting an
// exact solution.
//
// The LM method is often good in both cases, but sometimes has more trouble 
// converging than the above choices.  For some situations, though, it may
// be faster than these.
//
// The SecantLM method may be best when m > n and a direct calculation of 
// J is not possible.
//
// The SecantDogleg method may be best when m == n and a direct calculation of 
// J is not possible.
//
// To set which method to use for a particular solver objects type:
// solver.method = NLSolver::Newton;
// solver.method = NLSolver::Hybrid;
// solver.method = NLSolver::Dogleg;
// solver.method = NLSolver::LM;
// solver.method = NLSolver::SecantLM;
// solver.method = NLSolver::SecantDogleg;
//
//
// There are two ways that the algorithms can seccessfully converge:
//
// success if NormInf(f) < ftol
// success if NormInf(grad(1/2 NormSq(f)) < gtol
//
// There are also a number of ways they can fail.  
// Two important ones which you should set appropriately are:
//
// fail if Norm2(delta x) < min_step * (Norm2(x) + epsilon2)
// fail if number of iterations > max_iter
//
// The defaults are:
//
// ftol = 1.e-8
// gtol = 1.e-8
// min_step = 1.e-8
// max_iter = 200
//
// There are other failure modes for the various algorithms which use
// the above parameters for their tests, so no need to set anything else.
// However, there are a couple more parameters which affect how 
// various algorithms are initialized:
// 
// tau is a parameter used in LM, Hybrid, and SecantLM.  
// If the initial guess is thought to be very good, use tau = 1.e-6. 
// If it's somewhat reasonable, use tau = 1.e-3. (default)
// If it's bad, use tau = 1.
//
// delta0 is a parameter used in Dogleg and SecantDogleg.
// It gives the initial scale size allowed for steps in x.
// The algorithm will increase or decrease this as appropriate.
// The default value is 1.
//
// You can have progress text sent to an ostream by defining
// nlout = &os;  (e.g. nlout = &std::cout;)
// This will print basic information about each iteration.
// If you want much more information about what's going on, 
// you can set verbose = true;

class NLSolver 
{
  public :

    enum Method { Newton, Hybrid, Dogleg, LM, SecantLM, SecantDogleg };

    NLSolver() : 
      method(Newton),
      ftol(1.e-8), gtol(1.e-8), min_step(1.e-8), max_iter(200),
      tau(1.e-3), delta0(1.), nlout(0), 
      verbose(false), hasdirecth(false), startwithch(true) {}

    virtual ~NLSolver() {}

    // This is the basic function that needs to be overridden.
    virtual void F(const tmv::Vector<double>& x,
	tmv::Vector<double>& f) const =0;

    // J(i,j) = df_i/dx_j
    // If you don't overload the J function, then a finite
    // difference calculation will be performed.
    virtual void J(const tmv::Vector<double>& x,
	const tmv::Vector<double>& f, 
	tmv::Matrix<double>& j) const;

    // Try to solve for F(x) = 0.
    // Returns whether it succeeded.
    // On exit, x is the best solution found.
    // Also, f is returned as the value of F(x) for the best x,
    // regardless of whether the fit succeeded or not.
    virtual bool Solve(tmv::Vector<double>& x, tmv::Vector<double>& f,
	tmv::Matrix<double>* invcov = 0) const;

    // H(i,j) = d^2 Q / dx_i dx_j
    // where Q = 1/2 Sum_k |f_k|^2
    // H = JT J + Sum_k f_k d^2f_k/(dx_i dx_j)
    //
    // This is only used for the Hybrid method, and if it is not
    // overloaded, then an approximation is calculated on the fly.
    // It's not really important to overload this, but in case 
    // the calculation is very easy, a direct calculation would
    // be faster, so we allow for that possibility.
    virtual inline void H(const tmv::Vector<double>& x,
	const tmv::Vector<double>& f, 
	const tmv::Matrix<double>& j, 
	tmv::SymMatrix<double>& h) const;
    // Defined below...

    // This tests whether the direct calculation of J matches
    // the numerical approximation of J calculated from finite differences.
    // It is useful as a test that your analytic formula was coded correctly.
    //
    // If you pass it an ostream (e.g. os = &cout), then debugging
    // info will be sent to that stream.
    //
    // The second optional parameter, relerr, is for functions which 
    // are merely good approximations of the correct J, rather than
    // exact.  Normally, the routine tests whether the J function
    // calculates the same jacobian as a numerical approximation to 
    // within the numerical precision possible for doubles.
    // If J is only approximate, you can set relerr as the relative
    // error in J to be tested for.
    // If relerr == 0 then sqrt(numeric_limits<double>::epsilon())
    // = 1.56e-7 is used.  This is the default.
    //
    // Parameters are x, f, os, relerr.
    virtual bool TestJ(const tmv::Vector<double>& , tmv::Vector<double>& ,
	std::ostream* os=0, double relerr=0.) const;

    // Note: I just left these public, so they can be modified directly.
    // I could have put in Set and Get methods, but I didn't bother.
    Method method;
    double ftol;
    double gtol;
    double min_step;
    int max_iter;
    double tau;
    double delta0;
    std::ostream* nlout;
    bool verbose;
    bool hasdirecth;
    bool startwithch;

  private :

    bool Solve_Newton(tmv::Vector<double>& x, tmv::Vector<double>& f,
	tmv::Matrix<double>* invcov) const;
    bool Solve_LM(tmv::Vector<double>& x, tmv::Vector<double>& f,
	tmv::Matrix<double>* invcov) const;
    bool Solve_Dogleg(tmv::Vector<double>& x, tmv::Vector<double>& f,
	tmv::Matrix<double>* invcov) const;
    bool Solve_Hybrid(tmv::Vector<double>& x, tmv::Vector<double>& f,
	tmv::Matrix<double>* invcov) const;
    bool Solve_SecantLM(tmv::Vector<double>& x, tmv::Vector<double>& f,
	tmv::Matrix<double>* invcov) const;
    bool Solve_SecantDogleg(tmv::Vector<double>& x, tmv::Vector<double>& f,
	tmv::Matrix<double>* invcov) const;

};

inline void NLSolver::H(const tmv::Vector<double>& ,
    const tmv::Vector<double>& , const tmv::Matrix<double>& , 
    tmv::SymMatrix<double>& ) const
{ 
#ifdef NOTHROW
  std::cerr<<"H is undefined\n";
  exit(1);
#else
  throw 0; 
#endif
}


