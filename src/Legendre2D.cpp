#include <algorithm>
#include <stdexcept>

#include "Legendre2D.h"

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

void Legendre2D::setFunction(
    int xorder, int yorder, const DVector& fVect)
{
    if (_xorder != xorder || _yorder != yorder) {
        _xorder = xorder; _yorder = yorder;
        _coeffs.reset(new DMatrix(_xorder+1,_yorder+1));
        _coeffs->setZero();
    }
    int k=0;
    for(int m=0; m <= std::max(_xorder,_yorder); ++m) {
        for(int i=std::min(m,_xorder);m-i<=std::min(m,_yorder);--i) { 
            (*_coeffs)(i,m-i) = fVect(k++);
        }
    }
    Assert(k==(int)fVect.size());
}

Legendre2D::Legendre2D(std::istream& fin) : Function2D()
{
    // Order of parameters is same as for Polynomial2D.  Difference
    // is that instead of x,x^2,x^3,etc., we use P1(x),P2(x),P3(x),etc.
    fin >> _xorder >> _yorder >> _bounds;
    if (!fin) throw std::runtime_error("reading order, bounds");
    _coeffs.reset(new DMatrix(_xorder+1,_yorder+1));
    _coeffs->setZero();
    int max_order = std::max(_xorder,_yorder);
    for(int m=0; m<=max_order; ++m) {
        for(int i=std::min(m,_xorder); m-i<=std::min(m,_yorder); --i) {
            fin >> (*_coeffs)(i,m-i);
        }
    }
    if (!fin) throw std::runtime_error("reading (legendre) function");
}

void Legendre2D::write(std::ostream& fout) const
{
    int oldprec = fout.precision(6);
    std::ios::fmtflags oldf = 
        fout.setf(std::ios::scientific,std::ios::floatfield);
    int max_order = std::max(_xorder,_yorder);
    if (max_order == 0) {
        fout << "C " << (*_coeffs)(0,0) << std::endl;
    } else {
        fout << "L " << _xorder << ' ' << _yorder << ' ' << _bounds << ' ';
        for(int m=0; m<=max_order; ++m) {
            for(int i=std::min(m,_xorder);m-i<=std::min(m,_yorder);--i) {
                fout << (*_coeffs)(i,m-i) << ' ';
            }
        }
        fout << std::endl;
    }
    if (!fout) throw std::runtime_error("writing (legendre) function");
    fout.precision(oldprec);
    fout.flags(oldf);
}

void Legendre2D::addLinear(double a, double b, double c)
{
    double xAve = (getXMin() + getXMax())/2.;
    double yAve = (getYMin() + getYMax())/2.;
    double xRange = getXMax() - getXMin();
    double yRange = getYMax() - getYMin();

    (*_coeffs)(0,0) += a + b*xAve + c*yAve;
    (*_coeffs)(1,0) += b*xRange/2.;
    (*_coeffs)(0,1) += c*yRange/2.;
}

void Legendre2D::linearPreTransform(
    double , double , double , double , double , double )
{
    // Not implemented yet.
    Assert(false);
}

void Legendre2D::operator+=(const Function2D& rhs)
{
    const Legendre2D* lrhs = dynamic_cast<const Legendre2D*>(&rhs);
    Assert(lrhs);
    Assert(getBounds() == lrhs->getBounds());
    if (_xorder == lrhs->_xorder && _yorder == lrhs->_yorder) {
        *_coeffs += *lrhs->_coeffs;
    } else {
        int new_xorder = std::max(_xorder,lrhs->_xorder);
        int new_yorder = std::max(_yorder,lrhs->_yorder);
        std::auto_ptr<DMatrix > newc(
            new DMatrix(new_xorder+1,new_yorder+1));
        newc->setZero();
        newc->TMV_subMatrix(0,_xorder+1,0,_yorder+1) = *_coeffs;
        newc->TMV_subMatrix(0,lrhs->_xorder+1,0,lrhs->_yorder+1) += *lrhs->_coeffs;
        _coeffs = newc;
        _xorder = new_xorder;
        _yorder = new_yorder;
    }
}

// dP_2n(x)/dx = Sum_k=0..n-1 (4k+3) P_2k+1(x)
// dP_2n+1(x)/dx = Sum_k=0..n (4k+1) P_2k(x)
std::auto_ptr<Function2D > Legendre2D::dFdX() const 
{
    if (_xorder == 0) {
        return std::auto_ptr<Function2D >(new Constant2D());
    }
    if (_xorder == 1 && _yorder == 0) {
        return std::auto_ptr<Function2D >(
            new Constant2D((*_coeffs)(1,0)*2./(getXMax()-getXMin())));
    }

    int new_xorder = _xorder-1;
    int new_yorder = _xorder > _yorder ? _yorder : _yorder-1;

    std::auto_ptr<Legendre2D > temp(
        new Legendre2D(new_xorder,new_yorder,_bounds));
    // initialized to 0's

    int max_order = std::max(_xorder,_yorder);
    for(int i=_xorder;i>=1;--i) {
        for(int j=std::min(max_order-i,_yorder);j>=0;--j) {
            if (i%2 == 0) {
                for(int k=0;k<=i/2-1;++k) 
                    (*temp->_coeffs)(2*k+1,j) += (4.*k+3.)*(*_coeffs)(i,j);
            } else  {
                for(int k=0;k<=(i-1)/2;++k) 
                    (*temp->_coeffs)(2*k,j) += (4.*k+1.)*(*_coeffs)(i,j);
            }
        }
    }
    max_order = std::max(new_xorder,new_yorder);
    for(int i=new_xorder;i>=0;--i) {
        for(int j=std::min(max_order-i,new_yorder);j>=0;--j) {
            (*temp->_coeffs)(i,j) *= 2./(getXMax()-getXMin());
        }
    }
    return std::auto_ptr<Function2D >(temp);
}

// dP_2n(x)/dx = Sum_k=0..n-1 (4k+3) P_2k+1(x)
// dP_2n+1(x)/dx = Sum_k=0..n-1 (4k+1) P_2k(x)
std::auto_ptr<Function2D > Legendre2D::dFdY() const 
{
    if (_yorder == 0) {
        return std::auto_ptr<Function2D >(new Constant2D());
    }
    if (_yorder == 1 && _xorder == 0) {
        return std::auto_ptr<Function2D >(new Constant2D(
                (*_coeffs)(0,1)*2./(getYMax()-getYMin())));
    }

    int new_xorder = _yorder > _xorder ? _xorder : _xorder-1;
    int new_yorder = _yorder-1;

    std::auto_ptr<Legendre2D > temp(
        new Legendre2D(new_xorder,new_yorder,_bounds));
    // initialized to 0's

    int max_order = std::max(_xorder,_yorder);
    for(int j=_yorder;j>=1;--j) {
        for(int i=std::min(max_order-j,_xorder);i>=0;--i) {
            if (j%2 == 0) {
                for(int k=0;k<=(j-2)/2;++k) 
                    (*temp->_coeffs)(i,2*k+1) += (4.*k+3.)*(*_coeffs)(i,j);
            } else {
                for(int k=0;k<=(j-1)/2;++k) 
                    (*temp->_coeffs)(i,2*k) += (4.*k+1.)*(*_coeffs)(i,j);
            }
        }
    }
    max_order = std::max(new_xorder,new_yorder);
    for(int j=new_yorder;j>=0;--j) {
        for(int i=std::min(max_order-j,new_xorder);i>=0;--i) {
            (*temp->_coeffs)(i,j) *= 2./(getYMax()-getYMin());
        }
    }
    return std::auto_ptr<Function2D >(temp);
}

DVector Legendre2D::definePXY(
    int order, double xy, double min, double max) const
{
    DVector temp(order+1);
    double xyNew = (2.*xy-min-max)/(max-min);
    temp[0] = 1.; if(order>0) temp[1] = xyNew;
    for(int i=2;i<=order;++i)
        temp[i] = ((2.*i-1.)*xyNew*temp[i-1] - (i-1.)*temp[i-2])/i;
    return temp;
}

DVector Legendre2D::definePX(int order, double x) const
{
    if (x < getXMin() || x > getXMax()) {
        dbg<<"Error: x = "<<x<<
            ", min,max = "<<getXMin()<<','<<getXMax()<<std::endl;
        throw RangeException(Position(x,getYMin()),_bounds);
    }
    return definePXY(order,x,getXMin(),getXMax());
}

DVector Legendre2D::definePY(int order, double y) const
{
    if (y < getYMin() || y > getYMax()) {
        dbg<<"Error: y = "<<y<<
            ", min,max = "<<getYMin()<<','<<getYMax()<<std::endl;
        throw RangeException(Position(getXMin(),y),_bounds);
    }
    return definePXY(order,y,getYMin(),getYMax());
}


