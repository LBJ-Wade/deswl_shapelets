#include <vector>
#include "Bounds.h"

void Position::read(std::istream& fin)
{ double x,y; fin >> x >> y; _z = std::complex<double>(x,y); }

void Position::write(std::ostream& fout) const
{ fout << getX() << " " << getY() <<" "; }

void Bounds::operator+=(const Position& pos)
    // Expand the bounds to include the given position.
{
    if (_defined) {
        if (pos.getX() < _xmin) _xmin = pos.getX();
        else if (pos.getX() > _xmax) _xmax = pos.getX();
        if (pos.getY() < _ymin) _ymin = pos.getY();
        else if (pos.getY() > _ymax) _ymax = pos.getY();
    } else {
        _xmin = _xmax = pos.getX();
        _ymin = _ymax = pos.getY();
        _defined = true;
    }
}

void Bounds::operator+=(const Bounds& rhs)
    // Expand the bounds to include the given bounds
{
    if (!rhs.isDefined()) return;
    if (_defined) {
        if (rhs.getXMin() < _xmin) _xmin = rhs.getXMin();
        if (rhs.getXMax() > _xmax) _xmax = rhs.getXMax();
        if (rhs.getYMin() < _ymin) _ymin = rhs.getYMin();
        if (rhs.getYMax() > _ymax) _ymax = rhs.getYMax();
    } else {
        *this = rhs;
        _defined = true;
    }
}

bool Bounds::operator==(const Bounds& rhs) const
{
    if (!_defined) return (!rhs._defined);
    else return (_xmin == rhs._xmin && _xmax == rhs._xmax &&
                 _ymin == rhs._ymin && _ymax == rhs._ymax);
}


Position Bounds::getCenter() const
{
    return Position((_xmin + _xmax)/2.,(_ymin + _ymax)/2.);
}

Bounds Bounds::operator&(const Bounds& rhs) const
{
    Bounds temp(_xmin<rhs._xmin ? rhs._xmin : _xmin,
                _xmax>rhs._xmax ? rhs._xmax : _xmax,
                _ymin<rhs._ymin ? rhs._ymin : _ymin,
                _ymax>rhs._ymax ? rhs._ymax : _ymax);
    if (temp._xmin>temp._xmax || temp._ymin>temp._ymax) return Bounds();
    else return temp;
}

void Bounds::snap(double d)
{
    if (_defined) {
        _xmax = d*ceil(_xmax/d);
        _xmin = d*floor(_xmin/d);
        _ymax = d*ceil(_ymax/d);
        _ymin = d*floor(_ymin/d);
    }
}

void Bounds::addXBorder(double d)
{
    if (_defined) { _xmax += d; _xmin -= d; }
}

void Bounds::addYBorder(double d)
{
    if (_defined) { _ymax += d; _ymin -= d; }
}

void Bounds::addBorder(double d)
{
    if (_defined) { _xmax += d; _xmin -= d; _ymax += d; _ymin -= d; }
}

bool Bounds::includes(const Position& pos) const
{
    return (_defined && pos.getX()<=_xmax && pos.getX()>=_xmin &&
            pos.getY()<=_ymax && pos.getY()>=_ymin);
}

bool Bounds::includes(double x, double y) const
{ return includes(Position(x,y)); }

bool Bounds::includes(const Bounds& b2) const
{
    if (!b2.isDefined()) return true;
    else return (
        _defined &&
        b2._xmin >= _xmin && b2._xmax <= _xmax &&
        b2._ymin >= _ymin && b2._ymax <= _ymax);
}

bool Bounds::intersects(const Bounds& b2) const
{
    return (
        _defined && b2._defined &&
        !(b2._xmin >= _xmax) &&
        !(b2._xmax <= _xmin) &&
        !(b2._ymin >= _ymax) &&
        !(b2._ymax <= _ymin) );
}

double Bounds::getArea() const
{ return _defined ? (_xmax-_xmin)*(_ymax-_ymin) : 0.; }

std::vector<Bounds> Bounds::quarter() const
{ return divide(2,2); }

std::vector<Bounds> Bounds::divide(int nx, int ny) const
{
    if (!_defined) return std::vector<Bounds>(nx*ny);
    std::vector<Bounds> temp;
    temp.reserve(nx*ny);
    std::vector<double> x(nx+1);
    std::vector<double> y(ny+1);
    x[0] = _xmin;  x[nx] = _xmax;
    y[0] = _ymin;  y[ny] = _ymax;
    double xstep = (_xmax-_xmin)/nx;
    double ystep = (_ymax-_ymin)/ny;
    for(int i=1;i<nx;++i) x[i] = x[0]+i*xstep;
    for(int j=1;j<ny;++j) y[j] = y[0]+j*ystep;
    for(int i=0;i<nx;++i) for(int j=0;j<ny;++j)
        temp.push_back(Bounds(x[i],x[i+1],y[j],y[j+1]));
    return temp;
}

void Bounds::write(std::ostream& fout) const
{ fout << _xmin << ' ' << _xmax << ' ' << _ymin << ' ' << _ymax << ' '; }

void Bounds::read(std::istream& fin)
{ fin >> _xmin >> _xmax >> _ymin >> _ymax; _defined = true; }

