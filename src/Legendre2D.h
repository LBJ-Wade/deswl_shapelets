//---------------------------------------------------------------------------
#ifndef LEGENDRE2D_H
#define LEGENDRE2D_H
//---------------------------------------------------------------------------

#include <iostream>
#include "Function2D.h"

#define XMINDEF 0.
#define XMAXDEF 2048.
#define YMINDEF 0.
#define YMAXDEF 2048.

template <class T>
class Legendre2D : public Function2D<T> {

  public:

    Legendre2D() : Function2D<T>(), bounds(XMINDEF,XMAXDEF,YMINDEF,YMAXDEF) {}
    Legendre2D(double xmin,double xmax,double ymin,double ymax) :
      Function2D<T>(), bounds(xmin,xmax,ymin,ymax) {}
    Legendre2D(const Bounds& b) : Function2D<T>(),bounds(b) {}
    Legendre2D(const Legendre2D<T>& rhs) :
      Function2D<T>(rhs),bounds(rhs.bounds) {}
    Legendre2D(const Bounds& b, const tmv::Matrix<T>& a) : 
      Function2D<T>(a.colsize()-1,a.rowsize()-1,a),bounds(b) 
    { }
    Legendre2D(int xo, int yo, const Bounds& b) :
      Function2D<T>(xo,yo,tmv::Matrix<T>(xo+1,yo+1,0.)), bounds(b) {}
    Legendre2D(std::istream& fin);
    virtual ~Legendre2D() {}

    virtual void Write(std::ostream& fout) const;

    virtual std::auto_ptr<Function2D<T> > DFDX() const;
    virtual std::auto_ptr<Function2D<T> > DFDY() const;

    virtual std::auto_ptr<Function2D<T> > Copy() const
    { return std::auto_ptr<Function2D<T> >(new Legendre2D<T>(*this)); }

    virtual void AddLinear(T a, T b, T c);
    virtual void LinearPreTransform(T a, T b, T c, T d, T e, T f);
    virtual void operator+=(const Function2D<T>& rhs);

    double GetXMin() const {return bounds.GetXMin();}
    double GetXMax() const {return bounds.GetXMax();}
    double GetYMin() const {return bounds.GetYMin();}
    double GetYMax() const {return bounds.GetYMax();}
    const Bounds& GetBounds() const {return bounds;}

    virtual void SetFunction(int _xorder, int _yorder, const tmv::Vector<T>& fvect);
  private:

    Bounds bounds;
    using Function2D<T>::xorder;
    using Function2D<T>::yorder;
    using Function2D<T>::coeffs;

    tmv::Vector<double> DefinePXY(int order,double xy,double min,double max) const
    {
      //xdbg<<"Start DefPXY: order = "<<order<<endl;
      tmv::Vector<double> temp(order+1);
      //xdbg<<"temp = "<<temp<<endl;
      double newxy = (2.*xy-min-max)/(max-min);
      temp[0] = 1.; if(order>0) temp[1] = newxy;
      for(int i=2;i<=order;i++)
	temp[i] = ((2.*i-1.)*newxy*temp[i-1] - (i-1.)*temp[i-2])/i;
      //xdbg<<"temp = "<<temp<<endl;
      return temp;
    }
    virtual tmv::Vector<double> DefinePX(int order, double x) const
    { 
      if (x < GetXMin() || x > GetXMax()) {
	dbg<<"Error: x = "<<x<<", min,max = "<<GetXMin()<<','<<GetXMax()<<std::endl;
#ifndef NOTHROW
	throw Range_error(Position(x,GetYMin()),bounds);
#endif
      }
      return DefinePXY(order,x,GetXMin(),GetXMax()); 
    }
    virtual tmv::Vector<double> DefinePY(int order, double y) const
    { 
      if (y < GetYMin() || y > GetYMax()) {
	dbg<<"Error: y = "<<y<<", min,max = "<<GetYMin()<<','<<GetYMax()<<std::endl;
#ifndef NOTHROW
	throw Range_error(Position(GetXMin(),y),bounds);
#endif
      }
      return DefinePXY(order,y,GetYMin(),GetYMax()); 
    }

};

#undef XMINDEF
#undef XMAXDEF
#undef YMINDEF
#undef YMAXDEF

#endif
