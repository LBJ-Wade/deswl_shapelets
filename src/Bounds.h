//---------------------------------------------------------------------------
#ifndef BoundsH
#define BoundsH

#include <vector>
#include <complex>

//---------------------------------------------------------------------------

class Position
{

  public:

    Position() : z(0.,0.) {}
    Position(const Position& rhs) : z(rhs.z) {}
    Position(const std::complex<double>& rhs) : z(rhs) {}
    ~Position() {}
    Position(double xin,double yin) : z(xin,yin) {}
    Position& operator=(const Position& rhs) 
    { z = rhs.z; return *this; }
    Position& operator=(const std::complex<double>& rhs) 
    { z=rhs; return *this; }

    operator std::complex<double>() { return z; }
    std::complex<double> operator-(const Position& p2) const 
    { return z - p2.z; }
    std::complex<double> operator-(const std::complex<double>& z2) const 
    { return z - z2; }
    Position& operator*=(double x) 
    { z *= x; return *this; }
    Position& operator/=(double x) 
    { z /= x; return *this; }

    double GetX() const {return(z.real());}
    double GetY() const {return(z.imag());}

    void Read(std::istream& fin) 
    { double x,y; fin >> x >> y; z = std::complex<double>(x,y); }
    void Write(std::ostream& fout) const
    { fout << GetX() << " " << GetY() <<" "; }

  private :

    std::complex<double> z;

}; // Position

inline std::complex<double> operator-(
    const std::complex<double>& z1, const Position& p2)
{ return (p2-z1); }

inline std::ostream& operator<<(std::ostream& os, const Position& pos)
{ pos.Write(os); return os; }

inline std::istream& operator>>(std::istream& os, Position& pos)
{ pos.Read(os); return os; }

class Bounds 
{
  // Basically just a rectangle.  This is used to keep track of the bounds of
  // catalogs and fields.  You can set values, but generally you just keep
  // including positions of each galaxy or the bounds of each catalog
  // respectively using the += operators

  public:

    Bounds(double x1, double x2, double y1, double y2) :
      defined(true),xmin(x1),xmax(x2),ymin(y1),ymax(y2) {}
    Bounds(const Position& pos) :
      defined(true), xmin(pos.GetX()), xmax(pos.GetX()),
      ymin(pos.GetY()), ymax(pos.GetY()) {}
    Bounds(): defined(false),xmin(0.),xmax(0.),ymin(0.),ymax(0.) {}
    ~Bounds() {}
    void SetXMin(double x) { xmin = x; defined = true; }
    void SetXMax(double x) { xmax = x; defined = true; }
    void SetYMin(double y) { ymin = y; defined = true; }
    void SetYMax(double y) { ymax = y; defined = true; }
    double GetXMin() const { return xmin; }
    double GetXMax() const { return xmax; }
    double GetYMin() const { return ymin; }
    double GetYMax() const { return ymax; }
    bool IsDefined() const { return defined; }
    Position Center() const;
    void operator+=(const Position& pos);
    void operator+=(const Bounds& rec);

    bool operator==(const Bounds& rhs) const
    { 
      if (!defined) return (!rhs.defined);
      else return (xmin == rhs.xmin && xmax == rhs.xmax &&
	  ymin == rhs.ymin && ymax == rhs.ymax); 
    }

    void Snap(double d);
    void AddBorder(double d);
    void AddXBorder(double d);
    void AddYBorder(double d);
    Bounds operator&(const Bounds& rhs) const; // Finds intersection
    bool Includes(const Position& pos) const
    { 
      return (defined && pos.GetX()<=xmax && pos.GetX()>=xmin &&
	  pos.GetY()<=ymax && pos.GetY()>=ymin); 
    }
    bool Includes(double x,double y) const
    { return Includes(Position(x,y)); }
    bool Includes(const Bounds& b2) const
    {
      if (!b2.IsDefined()) return true;
      else return (
	  defined && 
	  b2.xmin >= xmin && b2.xmax <= xmax &&
	  b2.ymin >= ymin && b2.ymax <= ymax); 
    }
    bool Intersects(const Bounds& b2) const
    {
      return (
	  defined && b2.defined &&
	  !(b2.xmin >= xmax) &&
	  !(b2.xmax <= xmin) &&
	  !(b2.ymin >= ymax) &&
	  !(b2.ymax <= ymin) );
    }
    double Area() const
    { return defined ? (xmax-xmin)*(ymax-ymin) : 0.; }
    std::vector<Bounds> Quarter() const 
    { return Divide(2,2); }
    std::vector<Bounds> Divide(size_t nx, size_t ny) const;
    void Write(std::ostream& fout) const
    { fout << xmin << ' ' << xmax << ' ' << ymin << ' ' << ymax << ' '; }
    void Read(std::istream& fin)
    { fin >> xmin >> xmax >> ymin >> ymax; defined = true; }
    Position Get00() const {return Position(xmin,ymin);}
    Position Get01() const {return Position(xmin,ymax);}
    Position Get10() const {return Position(xmax,ymin);}
    Position Get11() const {return Position(xmax,ymax);}

    bool IsWide() const { return (xmax-xmin > ymax-ymin); }
    bool IsTall() const { return (xmax-xmin < ymax-ymin); }

  private:
    bool defined;
    double xmin,xmax,ymin,ymax;

};

inline std::ostream& operator<<(std::ostream& fout, const Bounds& b)
{ b.Write(fout); return fout;}

inline std::istream& operator>>(std::istream& fin, Bounds& b)
{ b.Read(fin); return fin;}


#endif
