#ifndef LEGENDRE2D_H
#define LEGENDRE2D_H

#include <iostream>
#include "Function2D.h"

template <typename T>
class Legendre2D : public Function2D<T> 
{

public:

    Legendre2D() {}

    Legendre2D(double xMin,double xMax,double yMin,double yMax) :
        Function2D<T>(), _bounds(xMin,xMax,yMin,yMax) {}

    Legendre2D(const Bounds& b) : Function2D<T>(), _bounds(b) {}

    Legendre2D(const Legendre2D<T>& rhs) :
        Function2D<T>(rhs), _bounds(rhs._bounds) {}

    Legendre2D(const Bounds& b, const tmv::Matrix<T>& a) : 
        Function2D<T>(a.colsize()-1,a.rowsize()-1,a), _bounds(b) {}

    Legendre2D(int xo, int yo, const Bounds& b) :
        Function2D<T>(xo,yo,tmv::Matrix<T>(xo+1,yo+1,0.)), _bounds(b) {}

    Legendre2D(std::istream& fin);

    virtual ~Legendre2D() {}

    virtual void write(std::ostream& fout) const;

    virtual std::auto_ptr<Function2D<T> > dFdX() const;

    virtual std::auto_ptr<Function2D<T> > dFdY() const;

    virtual std::auto_ptr<Function2D<T> > copy() const
    { return std::auto_ptr<Function2D<T> >(new Legendre2D<T>(*this)); }

    virtual void addLinear(T a, T b, T c);

    virtual void linearPreTransform(T a, T b, T c, T d, T e, T f);

    virtual void operator+=(const Function2D<T>& rhs);

    double getXMin() const {return _bounds.getXMin();}

    double getXMax() const {return _bounds.getXMax();}

    double getYMin() const {return _bounds.getYMin();}

    double getYMax() const {return _bounds.getYMax();}

    const Bounds& getBounds() const {return _bounds;}

    virtual void setFunction(
        int xorder, int yorder, const tmv::Vector<T>& fvect);

private:

    Bounds _bounds;
    using Function2D<T>::_xOrder;
    using Function2D<T>::_yOrder;
    using Function2D<T>::_coeffs;

    tmv::Vector<double> definePXY(
        int order,double xy,double min,double max) const;

    virtual tmv::Vector<double> definePX(int order, double x) const;

    virtual tmv::Vector<double> definePY(int order, double y) const;

};

extern template class Legendre2D<std::complex<double> >;
extern template class Legendre2D<double>;

#endif
