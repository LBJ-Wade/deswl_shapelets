#ifndef Transformation_H
#define Transformation_H

#include "Function2D.h"
#include "Bounds.h"
#include "TMV.h"
#include "TMV_Small.h"
#include "ConfigFile.h"

class Transformation {

  public :

    // Start with identity transformation.
    // u(x,y) = x, v(x,y) = y
    Transformation();

    // Read parameter file and load specified transformation function
    Transformation(const ConfigFile& params);

    // I/O
    void ReadFunc2D(std::istream& is);
    void ReadWCS(std::string fitsfile, int hdu);
    void SetToScale(double pixel_scale);
    void SetToJacobian(
	double dudx, double dudy, double dvdx, double dvdy);
    void WriteFunc2D(std::ostream& os) const;

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
    // It returns the valid bounds for the new inverse transformation.
    Bounds MakeInverseOf(const Transformation& t2, 
	const Bounds& bounds, int order);



  private :

    std::auto_ptr<Function2D<double> > up;
    std::auto_ptr<Function2D<double> > vp;
    std::auto_ptr<Function2D<double> > dudxp;
    std::auto_ptr<Function2D<double> > dudyp;
    std::auto_ptr<Function2D<double> > dvdxp;
    std::auto_ptr<Function2D<double> > dvdyp;


};

#endif
