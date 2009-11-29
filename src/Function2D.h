#ifndef Function2D_H
#define Function2D_H

#include <iostream>
#include <functional>
#include "Bounds.h"
#include "TMV.h"
#include "dbg.h"

struct RangeException : public std::runtime_error
{
public :
    RangeException(const Position& p, const Bounds& b);
    ~RangeException() throw() {} 

    const Bounds& getBounds() const { return _b; }
    const Position& getPosition() const { return _p; }

private : 
    Position _p;
    Bounds _b;

};

template <typename T> 
class Function2D 
{
public:

    Function2D() :
        _xOrder(0), _yOrder(0), _coeffs(new tmv::Matrix<T>(1,1,0.)) {}

    Function2D(int xo, int yo, const tmv::Matrix<T>& c) :
        _xOrder(xo),_yOrder(yo),_coeffs(new tmv::Matrix<T>(c)) {}

    Function2D(const Function2D<T>& rhs) :
        _xOrder(rhs._xOrder), _yOrder(rhs._yOrder),
        _coeffs(new tmv::Matrix<T>(*rhs._coeffs)) {}

    virtual ~Function2D() {}

    virtual void write(std::ostream& fout) const =0;

    static std::auto_ptr<Function2D<T> > read(std::istream& fin);

    virtual std::auto_ptr<Function2D<T> > copy() const =0;

    // Adds a + bx + cy to function
    virtual void addLinear(T a, T b, T c) = 0;

    // Converts function to f(a+bx+cy,d+ex+fy)
    virtual void linearPreTransform(T a, T b, T c, T d, T e, T f) = 0;

    // returns new function x derivative of f
    virtual std::auto_ptr<Function2D<T> > dFdX() const = 0;

    // returns new function y derivative of f
    virtual std::auto_ptr<Function2D<T> > dFdY() const = 0;

    virtual std::auto_ptr<Function2D<T> > conj() const;

    virtual void operator*=(double scale)
    { *_coeffs *= scale; }

    virtual T operator()(double x,double y) const;

    T operator()(const Position& p) const 
    { return operator()(p.getX(),p.getY()); }

    virtual void setTo(T value) 
    {
        if (_xOrder || _yOrder) {
            _xOrder = 0; _yOrder = 0; 
            _coeffs.reset(new tmv::Matrix<T>(1,1,value));
        } else { // _xOrder,_yOrder already 0
            (*_coeffs)(0,0) = value;
        }
    }

    bool isNonZero() const 
    { return (*_coeffs)(0,0) != 0.0 || _xOrder!=0 || _yOrder!=0; }

    int getXOrder() const { return _xOrder; }

    int getYOrder() const { return _yOrder; }

    const tmv::Matrix<T>& getCoeffs() const { return *_coeffs; }

    // Sets function to fit of f(pos_i) = v_i using i only if 
    // shouldUse[i] = true
    virtual void simpleFit(
        int order, const std::vector<Position>& pos, 
        const std::vector<T>& v, const std::vector<bool>& shouldUse,
        const std::vector<double>* sigList=0,
        double* chisqOut = 0, int* dofOut=0, tmv::Matrix<double>* cov=0);

    // Sets function to fit of f(pos_i) = v_i using i if fit is within 
    // nsig sigma.  *shouldUse is returned as list of good i's.
    virtual void outlierFit(
        int order,double nsig,
        const std::vector<Position>& pos, const std::vector<T>& v,
        std::vector<bool>* shouldUse,
        const std::vector<double>* sigList=0, 
        double* chisqout = 0, int* dofout = 0, tmv::Matrix<double>* cov=0);

    // Sets function to fit of f(pos_i) = v_i reducing the order as far
    // as possible keeping quality of fit the same as for maxOrder
    // equivProb = rejection percentile.  eg. 0.9 means a low order fit is
    // rejected if it is statistically rejected at the 90th percentile
    // Higher values of equivProb result in lower order fits.
    virtual void orderFit(
        int maxOrder, double equivProb,
        const std::vector<Position>& pos, const std::vector<T>& v,
        const std::vector<bool>& shouldUse, const std::vector<double>* sigList=0, 
        double *chisqout = 0, int *dofout = 0, tmv::Matrix<double>* cov=0);

    // Sets function to h(x,y) = a + b*f(x,y) + c*g(x,y)
    virtual void linearTransform(
        T a, T b, T c, 
        const Function2D<T>& f, const Function2D<T>& g);

    // Sets function to g(x,y) = a + b*f(x,y)
    virtual void linearTransform(T a, T b, const Function2D<T>& f)
    { linearTransform(a,b,0.,f,f); }

    virtual void operator+=(const Function2D<T>& rhs) = 0;

    virtual void setFunction(
        int _xOrder, int _yOrder, const tmv::Vector<T>& fVect) = 0;

protected:

    virtual tmv::Vector<double> definePX(int order, double x) const = 0;
    virtual tmv::Vector<double> definePY(int order, double y) const = 0;

    int _xOrder,_yOrder;
    std::auto_ptr<tmv::Matrix<T> > _coeffs;

private:

    void doSimpleFit(
        int xOrder, int yOrder, 
        const std::vector<Position>& pos, const std::vector<T>& v,
        const std::vector<bool>& shouldUse, tmv::Vector<T> *f, 
        const std::vector<double>* sigList=0, int *dof=0,
        tmv::Vector<T> *diff=0, tmv::Matrix<double>* cov=0);
};

template<typename T> 
inline std::ostream& operator<<(std::ostream& fout, const Function2D<T>& f)
{ f.write(fout); return fout; }

template <typename T> 
inline std::istream& operator>>(
    std::istream& fin, std::auto_ptr<Function2D<T> >& f)
{ f = Function2D<T>::read(fin); return fin; }

template <typename T> 
class Constant2D : public Function2D<T> 
{

public:

    Constant2D() : Function2D<T>() {}

    Constant2D(const Constant2D<T>& rhs) : Function2D<T>(rhs) {}

    Constant2D(T value) : Function2D<T>(0,0,tmv::Matrix<T>(1,1,value)) {}

    Constant2D(std::istream& fin);

    virtual ~Constant2D() {}

    virtual T operator()(double ,double ) const 
    { return (*this->_coeffs)(0,0); }

    virtual void write(std::ostream& fout) const;

    virtual std::auto_ptr<Function2D<T> > dFdX() const 
    { return std::auto_ptr<Function2D<T> >(new Constant2D<T>()); }

    virtual std::auto_ptr<Function2D<T> > dFdY() const 
    { return std::auto_ptr<Function2D<T> >(new Constant2D<T>()); }

    std::auto_ptr<Function2D<T> > copy() const 
    { return std::auto_ptr<Function2D<T> >(new Constant2D<T>(*this)); }

    virtual void addLinear(T a, T b, T c) 
    { 
        Assert(b == T(0));
        Assert(c == T(0));
        (*this->_coeffs)(0,0) += a; 
    }

    virtual void linearPreTransform(T , T , T , T , T , T )
    { Assert(false); return; }

    virtual void operator+=(const Function2D<T>& rhs);

    virtual void setFunction(
        int xOrder, int yOrder, const tmv::Vector<T>& fVect)
    { 
        Assert(_xOrder == 0);
        Assert(_yOrder == 0);
        Assert(fVect.size()==1);
        (*this->_coeffs)(0,0) = fVect(0); 
    }

private:

    virtual tmv::Vector<double> definePX(int order, double ) const
    { 
        Assert(order == 0);
        return tmv::Vector<double>(1,1.); 
    }
    virtual tmv::Vector<double> definePY(int order, double y) const
    { return definePX(order,y); }

    using Function2D<T>::_coeffs;
    using Function2D<T>::_xOrder;
    using Function2D<T>::_yOrder;
};

template <class T> class Polynomial2D : public Function2D<T> 
{

public:

    Polynomial2D(double scale=1.) : Function2D<T>(),_scale(scale) {}

    Polynomial2D(const Polynomial2D<T>& rhs) :
        Function2D<T>(rhs),_scale(rhs._scale) {}

    Polynomial2D(const tmv::Matrix<T>& a,double scale=1.) : 
        Function2D<T>(a.colsize()-1,a.rowsize()-1,a),_scale(scale) {}

    Polynomial2D(std::istream& fin);

    Polynomial2D(int xo,int yo,double scale=1.) :
        Function2D<T>(xo,yo,tmv::Matrix<T>(xo+1,yo+1,0.)),_scale(scale) {}

    virtual ~Polynomial2D() {}

    virtual void write(std::ostream& fout) const;

    virtual std::auto_ptr<Function2D<T> > dFdX() const;

    virtual std::auto_ptr<Function2D<T> > dFdY() const;

    std::auto_ptr<Function2D<T> > copy() const 
    { return std::auto_ptr<Function2D<T> >(new Polynomial2D<T>(*this)); }

    virtual void addLinear(T a, T b, T c);

    virtual void linearPreTransform(T a, T b, T c, T d, T e, T f);

    virtual void operator+=(const Function2D<T>& rhs);

    void makeProductOf(const Polynomial2D<T>& f, const Polynomial2D<T>& g);

    virtual void setFunction(
        int xOrder, int yOrder, const tmv::Vector<T>& fVect);

private:

    double _scale;

    virtual tmv::Vector<double> definePX(int order, double x) const
    {
        tmv::Vector<double> temp(order+1,1.);
        for(int i=1;i<=order;++i) temp(i) = temp(i-1)*x/_scale;
        return temp;
    }
    virtual tmv::Vector<double> definePY(int order, double y) const
    { return definePX(order,y); }

    using Function2D<T>::_coeffs;
    using Function2D<T>::_xOrder;
    using Function2D<T>::_yOrder;
};

extern template class Function2D<std::complex<double> >;
extern template class Function2D<double>;
extern template class Constant2D<std::complex<double> >;
extern template class Constant2D<double>;
extern template class Polynomial2D<std::complex<double> >;
extern template class Polynomial2D<double>;

#endif
