#ifndef BoundsH
#define BoundsH

#include <vector>
#include <complex>
#include "dbg.h"

class Position
{

public:

    Position() : _z(0.,0.) {}

    Position(const Position& rhs) : _z(rhs._z) {}

    Position(const std::complex<double>& rhs) : _z(rhs) {}

    ~Position() {}

    Position(double x, double y) : _z(x,y) {}

    Position& operator=(const Position& rhs) 
    { _z = rhs._z; return *this; }

    Position& operator=(const std::complex<double>& rhs) 
    { _z=rhs; return *this; }

    operator std::complex<double>() { return _z; }

    std::complex<double> operator-(const Position& p2) const 
    { return _z - p2._z; }

    //std::complex<double> operator-(const std::complex<double>& z2) const 
    //{ return _z - z2; }

    Position& operator*=(double x) 
    { _z *= x; return *this; }

    Position& operator/=(double x) 
    { _z /= x; return *this; }

    Position& operator+=(const std::complex<double>& z2) 
    { _z += z2; return *this; }
    Position operator+(const std::complex<double>& z2) const
    { return Position(*this) += z2; }
    Position& operator-=(const std::complex<double>& z2) 
    { _z -= z2; return *this; }
    Position operator-(const std::complex<double>& z2) const
    { return Position(*this) -= z2; }

    double getX() const { return(_z.real()); }

    double getY() const { return(_z.imag()); }

    void read(std::istream& fin);

    void write(std::ostream& fout) const;

private :

    std::complex<double> _z;

}; // Position

inline std::complex<double> operator-(
    const std::complex<double>& z1, const Position& p2)
{ return (p2-z1); }

inline std::complex<double> operator/(
    const Position& p1, double x)
{ Position p2 = p1;  p2 /= x; return p2; }

inline std::ostream& operator<<(std::ostream& os, const Position& pos)
{ pos.write(os); return os; }

inline std::istream& operator>>(std::istream& os, Position& pos)
{ pos.read(os); return os; }

class Bounds 
{
    // Basically just a rectangle.  This is used to keep track of the bounds of
    // catalogs and fields.  You can set values, but generally you just keep
    // including positions of each galaxy or the bounds of each catalog
    // respectively using the += operators

public:

    Bounds(double x1, double x2, double y1, double y2) :
        _defined(true),_xmin(x1),_xmax(x2),_ymin(y1),_ymax(y2) {}

    Bounds(const Position& pos) :
        _defined(true), _xmin(pos.getX()), _xmax(pos.getX()),
        _ymin(pos.getY()), _ymax(pos.getY()) {}

    Bounds(): _defined(false),_xmin(0.),_xmax(0.),_ymin(0.),_ymax(0.) {}

    ~Bounds() {}

    void setXMin(double x) { _xmin = x; _defined = true; }

    void setXMax(double x) { _xmax = x; _defined = true; }

    void setYMin(double y) { _ymin = y; _defined = true; }

    void setYMax(double y) { _ymax = y; _defined = true; }

    double getXMin() const { return _xmin; }

    double getXMax() const { return _xmax; }

    double getYMin() const { return _ymin; }

    double getYMax() const { return _ymax; }

    bool isDefined() const { return _defined; }

    Position getCenter() const;

    void operator+=(const Position& pos);

    void operator+=(const Bounds& rec);

    bool operator==(const Bounds& rhs) const;

    void snap(double d);

    void addBorder(double d);

    void addXBorder(double d);

    void addYBorder(double d);

    Bounds operator&(const Bounds& rhs) const; // Finds intersection

    bool includes(const Position& pos) const;

    bool includes(double x, double y) const;

    bool includes(const Bounds& b2) const;

    bool intersects(const Bounds& b2) const;

    double getArea() const;

    std::vector<Bounds> quarter() const;

    std::vector<Bounds> divide(int nx, int ny) const;

    void write(std::ostream& fout) const;

    void read(std::istream& fin);

    Position get00() const { return Position(_xmin,_ymin); }

    Position get01() const { return Position(_xmin,_ymax); }

    Position get10() const { return Position(_xmax,_ymin); }

    Position get11() const { return Position(_xmax,_ymax); }

    bool isWide() const { return (_xmax-_xmin > _ymax-_ymin); }

    bool isTall() const { return (_xmax-_xmin < _ymax-_ymin); }

private:

    bool _defined;
    double _xmin,_xmax,_ymin,_ymax;

};

inline std::ostream& operator<<(std::ostream& fout, const Bounds& b)
{ b.write(fout); return fout;}

inline std::istream& operator>>(std::istream& fin, Bounds& b)
{ b.read(fin); return fin;}


#endif
