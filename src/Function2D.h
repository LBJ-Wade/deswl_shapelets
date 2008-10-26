//---------------------------------------------------------------------------
#ifndef Function2D_H
#define Function2D_H
//---------------------------------------------------------------------------

#include <iostream>
#include <functional>
#include "Bounds.h"
#include "TMV.h"
#include "dbg.h"

inline void f2d_error(const char *s)
{ std::cerr<<s<<std::endl; exit(1); }

struct Range_error {
  Position p;
  Bounds b;

  Range_error(const Position& _p, const Bounds& _b) : p(_p),b(_b) {}
};

template <class T> class Function2D 
  //: public std::binary_function<double,double,T> 
{
  // This class heirarchy really isn't organized very well.
  // It's sort of legacy code from when I was first learning C++.
  // But I haven't taken the time to clean it up, since I would also
  // have to change a lot of other code that uses it.
  //
  // For example, Function2D should really be an abstract base class.
  // The order, coeffs, etc. should be all in the derived classes.
  // Likewise, some of the functions here aren't really appropriate
  // for all possible 2d functions that I might want to write
  // (e.g. the fitting routines, addlinear, etc.)
  //

  public:

    Function2D() : xorder(0),yorder(0),coeffs(new tmv::Matrix<T>(1,1,0.)) {}
    Function2D(int xo, int yo, const tmv::Matrix<T>& c) :
      xorder(xo),yorder(yo),coeffs(new tmv::Matrix<T>(c)) {}
    Function2D(const Function2D<T>& rhs) :
      xorder(rhs.xorder), yorder(rhs.yorder),
      coeffs(new tmv::Matrix<T>(*rhs.coeffs)) {}
    virtual ~Function2D() { }

    virtual void Write(std::ostream& fout) const =0;

    static std::auto_ptr<Function2D<T> > Read(std::istream& fin);

    virtual std::auto_ptr<Function2D<T> > Copy() const =0;

    virtual void AddLinear(T a, T b, T c) = 0;
    // Adds a + bx + cy to function
    virtual void LinearPreTransform(T a, T b, T c, T d, T e, T f) = 0;
    // Converts function to f(a+bx+cy,d+ex+fy)
    virtual std::auto_ptr<Function2D<T> > DFDX() const = 0;
    // returns new function x derivative of f
    virtual std::auto_ptr<Function2D<T> > DFDY() const = 0;
    // returns new function y derivative of f

    virtual std::auto_ptr<Function2D<T> > Conj() const;
    virtual void operator*=(double mult)
    { *coeffs *= mult; }
    virtual T operator()(double x,double y) const;
    T operator()(const Position& p) const 
    { return operator()(p.GetX(),p.GetY()); }
    virtual void SetTo(T value) 
    {
      if (xorder || yorder) {
	xorder = 0; yorder = 0; 
	coeffs.reset(new tmv::Matrix<T>(1,1,value));
      } else { // xorder,yorder already 0
	(*coeffs)(0,0) = value;
      }
    }

    bool NonZero() const 
    { return (*coeffs)(0,0) != 0.0 || xorder!=0 || yorder!=0; }
    int GetXOrder() const { return xorder; }
    int GetYOrder() const { return yorder; }
    const tmv::Matrix<T>& GetCoeffs() const { return *coeffs; }

    virtual void SimpleFit(int order, const std::vector<Position>& pos, 
	const std::vector<T>& v, const std::vector<bool>& use,
	const std::vector<double>* siglist=0,
	double* chisqout = 0, int* dofout=0, tmv::Matrix<double>* cov=0);
    // Sets function to fit of f(pos_i) = v_i using i only if use[i] = true

    virtual void OutlierFit(int order,double nsig,
	const std::vector<Position>& pos, const std::vector<T>& v,
	std::vector<bool>* use,
	const std::vector<double>* siglist=0, 
	double* chisqout = 0, int* dofout = 0, tmv::Matrix<double>* cov=0);
    // Sets function to fit of f(pos_i) = v_i using i if fit is within 
    // nsig sigma.  *use is returned as list of good i's.

    virtual void OrderFit(int maxorder, double equivprob,
	const std::vector<Position>& pos, const std::vector<T>& v,
	const std::vector<bool>& use, const std::vector<double>* siglist=0, 
	double *chisqout = 0, int *dofout = 0, tmv::Matrix<double>* cov=0);
    // Sets function to fit of f(pos_i) = v_i reducing the order as far
    // as possible keeping quality of fit the same as for maxorder
    // equivprob = rejection percentile.  eg. 0.9 means a low order fit is
    // rejected if it is statistically rejected at the 90th percentile
    // Higher values of equivprob result in lower order fits.

    virtual void LinearTransform(T a, T b, T c, 
	const Function2D<T>& f, const Function2D<T>& g);
    // Sets function to h(x,y) = a + b*f(x,y) + c*g(x,y)

    virtual void LinearTransform(T a, T b, const Function2D<T>& f)
    { LinearTransform(a,b,0.,f,f); }
    // Sets function to g(x,y) = a + b*f(x,y)

    virtual void operator+=(const Function2D<T>& rhs) = 0;

    virtual void SetFunction(int xorder, int yorder, 
	const tmv::Vector<T>& fvect) = 0;

  protected:

    virtual tmv::Vector<double> DefinePX(int order, double x) const = 0;
    virtual tmv::Vector<double> DefinePY(int order, double y) const = 0;

    int xorder,yorder;
    std::auto_ptr<tmv::Matrix<T> > coeffs;

  private:

    void DoSimpleFit(size_t xorder, size_t yorder, 
	const std::vector<Position>& pos, const std::vector<T>& v,
	const std::vector<bool>& use, tmv::Vector<T> *f, 
	const std::vector<double>* siglist=0, int *dof=0,
	tmv::Vector<T> *diff=0, tmv::Matrix<double>* cov=0);
};

template<class T> inline std::ostream& operator<<(
    std::ostream& fout, const Function2D<T>& f)
{ f.Write(fout); return fout; }

template <class T> inline std::istream& operator>>(
    std::istream& fin, std::auto_ptr<Function2D<T> >& f)
{ f = Function2D<T>::Read(fin); return fin; }

template <class T> class Constant2D : public Function2D<T> 
{

  public:

    Constant2D() : Function2D<T>() {}
    Constant2D(const Constant2D<T>& rhs) : Function2D<T>(rhs) {}
    Constant2D(T value) : Function2D<T>(0,0,tmv::Matrix<T>(1,1,value)) {}
    Constant2D(std::istream& fin);
    virtual ~Constant2D() {}

    virtual T operator()(double ,double ) const 
    { return (*this->coeffs)(0,0); }
    virtual void Write(std::ostream& fout) const;

    virtual std::auto_ptr<Function2D<T> > DFDX() const 
    { return std::auto_ptr<Function2D<T> >(new Constant2D<T>()); }
    virtual std::auto_ptr<Function2D<T> > DFDY() const 
    { return std::auto_ptr<Function2D<T> >(new Constant2D<T>()); }

    std::auto_ptr<Function2D<T> > Copy() const 
    { return std::auto_ptr<Function2D<T> >(new Constant2D<T>(*this)); }

    virtual void AddLinear(T a, T b, T c) 
    { 
      Assert(b == T(0));
      Assert(c == T(0));
      (*this->coeffs)(0,0) += a; 
    }
    virtual void LinearPreTransform(T , T , T , T , T , T )
    { Assert(false); return; }
    virtual void operator+=(const Function2D<T>& rhs);

    virtual void SetFunction(int xorder, int yorder, 
	const tmv::Vector<T>& fvect)
    { 
      Assert(xorder == 0);
      Assert(yorder == 0);
      Assert(fvect.size()==1);
      (*this->coeffs)(0,0) = fvect(0); 
    }

  private:

    virtual tmv::Vector<double> DefinePX(int order, double ) const
    { 
      Assert(order == 0);
      return tmv::Vector<double>(1,1.); 
    }
    virtual tmv::Vector<double> DefinePY(int order, double y) const
    { return DefinePX(order,y); }

    using Function2D<T>::coeffs;
    using Function2D<T>::xorder;
    using Function2D<T>::yorder;
};

template <class T> class Polynomial2D : public Function2D<T> 
{

  public:

    Polynomial2D(double _scale=1.) : Function2D<T>(),scale(_scale) {}
    Polynomial2D(const Polynomial2D<T>& rhs) :
      Function2D<T>(rhs),scale(rhs.scale) {}
    Polynomial2D(const tmv::Matrix<T>& a,double _scale=1.) : 
      Function2D<T>(a.colsize()-1,a.rowsize()-1,a),scale(_scale) {}
    Polynomial2D(std::istream& fin);
    Polynomial2D(int xo,int yo,double _scale=1.) :
      Function2D<T>(xo,yo,tmv::Matrix<T>(xo+1,yo+1,0.)),scale(_scale) {}

    virtual ~Polynomial2D() {}

    virtual void Write(std::ostream& fout) const;

    virtual std::auto_ptr<Function2D<T> > DFDX() const;
    virtual std::auto_ptr<Function2D<T> > DFDY() const;

    std::auto_ptr<Function2D<T> > Copy() const 
    { return std::auto_ptr<Function2D<T> >(new Polynomial2D<T>(*this)); }

    virtual void AddLinear(T a, T b, T c);
    virtual void LinearPreTransform(T a, T b, T c, T d, T e, T f);
    virtual void operator+=(const Function2D<T>& rhs);

    virtual void SetFunction(int _xorder, int _yorder, const tmv::Vector<T>& fvect);

  private:

    double scale;

    virtual tmv::Vector<double> DefinePX(int order, double x) const
    {
      tmv::Vector<double> temp(order+1,1.);
      for(int i=1;i<=order;i++) temp(i) = temp(i-1)*x/scale;
      return temp;
    }
    virtual tmv::Vector<double> DefinePY(int order, double y) const
    { return DefinePX(order,y); }

    using Function2D<T>::coeffs;
    using Function2D<T>::xorder;
    using Function2D<T>::yorder;
};

#endif
