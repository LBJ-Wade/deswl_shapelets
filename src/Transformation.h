#ifndef Transformation_H
#define Transformation_H

#include "Function2D.h"
#include "Bounds.h"
#include "TMV.h"
#include "TMV_Small.h"

class Transformation {

  public :

    // Start with identity transformation.
    // u(x,y) = x, v(x,y) = y
    Transformation();

    // I/O
    void Read(std::istream& is);
    void Write(std::ostream& os) const;

    // Read WCS information - this is based on code I wrote back in 
    // 2001, so I think the WCS standard has changed since them.
    // But it might provide a guide for how to update this to the 
    // new WCS standard.
    void ReadWCS(std::string fitsfile);

    // Calculate u,v = u(x,y),v(x,y)
    void Transform(Position pxy, Position& puv) const;

    // Update second moments ixx -> iuu, etc. according to the 
    // jacobian of the transformation
    void Distort(Position p, double& ixx, double& ixy, double& iyy) const;

    // Get the jacobian at a particular x,y
    // J = ( dudx  dudy )
    //     ( dvdx  dvdy )
    void GetDistortion(Position p, tmv::Matrix<double>& J) const;
    void GetDistortion(Position p, tmv::SmallMatrix<double,2,2>& J) const;
    void GetDistortion(Position p, 
	double& dudx, double& dudy, double& dvdx, double& dvdy) const;

    // Calculate x,y such that u,v = u(x,y),v(x,y)
    // This uses a non-linear solver to solver for x,y.
    // The value of puv on input is used as an initial guess.
    // The return value indicates whether a solution was found.
    bool InverseTransform(Position pxy, Position& puv) const;

    // Make this Transformation an approximate inverse of t2 using 
    // Legendre polynomials up to the given order in x,y.
    // The resulting transformation will be defined over a square
    // region in (u,v) that corresponds to the region given by
    // bounds in (x,y).
    void MakeInverseOf(const Transformation& t2, 
	const Bounds& bounds, int order);

  private :

    std::auto_ptr<Function2D<double> > up;
    std::auto_ptr<Function2D<double> > vp;
    std::auto_ptr<Function2D<double> > dudxp;
    std::auto_ptr<Function2D<double> > dudyp;
    std::auto_ptr<Function2D<double> > dvdxp;
    std::auto_ptr<Function2D<double> > dvdyp;
};

inline std::ostream& operator<<(std::ostream& os, const Transformation& t)
{ t.Write(os); return os; }

inline std::istream& operator>>(std::istream& is, Transformation& t)
{ t.Read(is); return is; }

#endif
