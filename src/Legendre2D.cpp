#include <algorithm>

#include "Legendre2D.h"
#include "Function2D.cpp"

static std::string makeRangeExceptionMessage(
    const Position& p, const Bounds& b)
{
    std::ostringstream s;
    s << "Range error: Position "<<p<<" is not in Bounds "<<b;
    return s.str();
}

RangeException::RangeException(const Position& p, const Bounds& b) : 
    std::runtime_error(makeRangeExceptionMessage(p,b)), _p(p), _b(b) 
{}

template <typename T>
void Legendre2D<T>::setFunction(
    int xOrder, int yOrder, const tmv::Vector<T>& fVect)
{
    if (_xOrder != xOrder || _yOrder != yOrder) {
        _xOrder = xOrder; _yOrder = yOrder;
        _coeffs.reset(new tmv::Matrix<T>(_xOrder+1,_yOrder+1,0.));
    }
    int k=0;
    for(int m=0; m <= std::max(_xOrder,_yOrder); ++m) {
        for(int i=std::min(m,_xOrder);m-i<=std::min(m,_yOrder);--i) { 
            (*_coeffs)(i,m-i) = fVect(k++);
        }
    }
    Assert(k==(int)fVect.size());
}

template <typename T>
Legendre2D<T>::Legendre2D(std::istream& fin) : Function2D<T>()
{
    // Order of parameters is same as for Polynomial2D.  Difference
    // is that instead of x,x^2,x^3,etc., we use P1(x),P2(x),P3(x),etc.
    fin >> _xOrder >> _yOrder >> _bounds;
    if (!fin) throw std::runtime_error("reading order, bounds");
    _coeffs.reset(new tmv::Matrix<T>(_xOrder+1,_yOrder+1,0.));
    int maxOrder = std::max(_xOrder,_yOrder);
    for(int m=0; m<=maxOrder; ++m) {
        for(int i=std::min(m,_xOrder); m-i<=std::min(m,_yOrder); --i) {
            fin >> (*_coeffs)(i,m-i);
        }
    }
    if (!fin) throw std::runtime_error("reading (legendre) function");
}

template <typename T>
void Legendre2D<T>::write(std::ostream& fout) const
{
    int oldprec = fout.precision(6);
    std::ios::fmtflags oldf = 
        fout.setf(std::ios::scientific,std::ios::floatfield);
    int maxOrder = std::max(_xOrder,_yOrder);
    if (maxOrder == 0) {
        fout << "C " << (*_coeffs)(0,0) << std::endl;
    } else {
        fout << "L " << _xOrder << ' ' << _yOrder << ' ' << _bounds << ' ';
        for(int m=0; m<=maxOrder; ++m) {
            for(int i=std::min(m,_xOrder);m-i<=std::min(m,_yOrder);--i) {
                fout << (*_coeffs)(i,m-i) << ' ';
            }
        }
        fout << std::endl;
    }
    if (!fout) throw std::runtime_error("writing (legendre) function");
    fout.precision(oldprec);
    fout.flags(oldf);
}

template <typename T>
void Legendre2D<T>::addLinear(T a, T b, T c)
{
    double xAve = (getXMin() + getXMax())/2.;
    double yAve = (getYMin() + getYMax())/2.;
    double xRange = getXMax() - getXMin();
    double yRange = getYMax() - getYMin();

    (*_coeffs)(0,0) += a + b*xAve + c*yAve;
    (*_coeffs)(1,0) += b*xRange/2.;
    (*_coeffs)(0,1) += c*yRange/2.;
}

template <typename T>
void Legendre2D<T>::linearPreTransform(T , T , T , T , T , T )
{
    // Not implemented yet.
    Assert(false);
}

template <typename T>
void Legendre2D<T>::operator+=(const Function2D<T>& rhs)
{
    const Legendre2D<T>* lrhs = dynamic_cast<const Legendre2D<T>*>(&rhs);
    Assert(lrhs);
    Assert(getBounds() == lrhs->getBounds());
    if (_xOrder == lrhs->_xOrder && _yOrder == lrhs->_yOrder) {
        *_coeffs += *lrhs->_coeffs;
    } else {
        int newXOrder = std::max(_xOrder,lrhs->_xOrder);
        int newYOrder = std::max(_yOrder,lrhs->_yOrder);
        std::auto_ptr<tmv::Matrix<T> > newc(
            new tmv::Matrix<T>(newXOrder+1,newYOrder+1,0.));
        newc->SubMatrix(0,_xOrder+1,0,_yOrder+1) = *_coeffs;
        newc->SubMatrix(0,lrhs->_xOrder+1,0,lrhs->_yOrder+1) += *lrhs->_coeffs;
        _coeffs = newc;
        _xOrder = newXOrder;
        _yOrder = newYOrder;
    }
}

// dP_2n(x)/dx = Sum_k=0..n-1 (4k+3) P_2k+1(x)
// dP_2n+1(x)/dx = Sum_k=0..n (4k+1) P_2k(x)
template <typename T>
std::auto_ptr<Function2D<T> > Legendre2D<T>::dFdX() const 
{
    if (_xOrder == 0) {
        return std::auto_ptr<Function2D<T> >(new Constant2D<T>());
    }
    if (_xOrder == 1 && _yOrder == 0) {
        return std::auto_ptr<Function2D<T> >(
            new Constant2D<T>((*_coeffs)(1,0)*2./(getXMax()-getXMin())));
    }

    int newXOrder = _xOrder-1;
    int newYOrder = _xOrder > _yOrder ? _yOrder : _yOrder-1;

    std::auto_ptr<Legendre2D<T> > temp(
        new Legendre2D<T>(newXOrder,newYOrder,_bounds));
    // initialized to 0's

    int maxOrder = std::max(_xOrder,_yOrder);
    for(int i=_xOrder;i>=1;--i) {
        for(int j=std::min(maxOrder-i,_yOrder);j>=0;--j) {
            if (i%2 == 0) {
                for(int k=0;k<=i/2-1;++k) 
                    (*temp->_coeffs)(2*k+1,j) += (4.*k+3.)*(*_coeffs)(i,j);
            } else  {
                for(int k=0;k<=(i-1)/2;++k) 
                    (*temp->_coeffs)(2*k,j) += (4.*k+1.)*(*_coeffs)(i,j);
            }
        }
    }
    maxOrder = std::max(newXOrder,newYOrder);
    for(int i=newXOrder;i>=0;--i) {
        for(int j=std::min(maxOrder-i,newYOrder);j>=0;--j) {
            (*temp->_coeffs)(i,j) *= 2./(getXMax()-getXMin());
        }
    }
    return std::auto_ptr<Function2D<T> >(temp);
}

// dP_2n(x)/dx = Sum_k=0..n-1 (4k+3) P_2k+1(x)
// dP_2n+1(x)/dx = Sum_k=0..n-1 (4k+1) P_2k(x)
template <typename T>
std::auto_ptr<Function2D<T> > Legendre2D<T>::dFdY() const 
{
    if (_yOrder == 0) {
        return std::auto_ptr<Function2D<T> >(new Constant2D<T>());
    }
    if (_yOrder == 1 && _xOrder == 0) {
        return std::auto_ptr<Function2D<T> >(new Constant2D<T>(
                (*_coeffs)(0,1)*2./(getYMax()-getYMin())));
    }

    int newXOrder = _yOrder > _xOrder ? _xOrder : _xOrder-1;
    int newYOrder = _yOrder-1;

    std::auto_ptr<Legendre2D<T> > temp(
        new Legendre2D<T>(newXOrder,newYOrder,_bounds));
    // initialized to 0's

    int maxOrder = std::max(_xOrder,_yOrder);
    for(int j=_yOrder;j>=1;--j) {
        for(int i=std::min(maxOrder-j,_xOrder);i>=0;--i) {
            if (j%2 == 0) {
                for(int k=0;k<=(j-2)/2;++k) 
                    (*temp->_coeffs)(i,2*k+1) += (4.*k+3.)*(*_coeffs)(i,j);
            } else {
                for(int k=0;k<=(j-1)/2;++k) 
                    (*temp->_coeffs)(i,2*k) += (4.*k+1.)*(*_coeffs)(i,j);
            }
        }
    }
    maxOrder = std::max(newXOrder,newYOrder);
    for(int j=newYOrder;j>=0;--j) {
        for(int i=std::min(maxOrder-j,newXOrder);i>=0;--i) {
            (*temp->_coeffs)(i,j) *= 2./(getYMax()-getYMin());
        }
    }
    return std::auto_ptr<Function2D<T> >(temp);
}

template <typename T>
tmv::Vector<double> Legendre2D<T>::definePXY(
    int order, double xy, double min, double max) const
{
    tmv::Vector<double> temp(order+1);
    double xyNew = (2.*xy-min-max)/(max-min);
    temp[0] = 1.; if(order>0) temp[1] = xyNew;
    for(int i=2;i<=order;++i)
        temp[i] = ((2.*i-1.)*xyNew*temp[i-1] - (i-1.)*temp[i-2])/i;
    return temp;
}

template <typename T>
tmv::Vector<double> Legendre2D<T>::definePX(int order, double x) const
{
    if (x < getXMin() || x > getXMax()) {
        dbg<<"Error: x = "<<x<<
            ", min,max = "<<getXMin()<<','<<getXMax()<<std::endl;
        throw RangeException(Position(x,getYMin()),_bounds);
    }
    return definePXY(order,x,getXMin(),getXMax());
}

template <typename T>
tmv::Vector<double> Legendre2D<T>::definePY(int order, double y) const
{
    if (y < getYMin() || y > getYMax()) {
        dbg<<"Error: y = "<<y<<
            ", min,max = "<<getYMin()<<','<<getYMax()<<std::endl;
        throw RangeException(Position(getXMin(),y),_bounds);
    }
    return definePXY(order,y,getYMin(),getYMax());
}


template class Legendre2D<std::complex<double> >;
template class Legendre2D<double>;
