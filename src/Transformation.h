#ifndef Transformation_H
#define Transformation_H

#include <stdexcept>
#include "Function2D.h"
#include "Bounds.h"
#include "MyMatrix.h"
#include "ConfigFile.h"

class Transformation 
{

public :

    // Start with identity transformation.
    // u(x,y) = x, v(x,y) = y
    Transformation();

    // Read parameter file and load specified transformation function
    Transformation(const ConfigFile& params);

    // I/O
    void readFunc2D(std::istream& is);
    void readWCS(std::string fitsfile, int hdu);
    void setToScale(double pixel_scale);
    void setToJacobian(
        double dudx, double dudy, double dvdx, double dvdy);
    void writeFunc2D(std::ostream& os) const;
    bool isRaDec() const { return _isRaDec; }

    // Calculate u,v = u(x,y),v(x,y)
    void transform(Position pxy, Position& puv) const;

    Position operator()(Position pxy) const
    { Position puv; transform(pxy,puv); return puv; }

    // Update second moments ixx -> iuu, etc. according to the 
    // jacobian of the transformation
    void distort(Position p, double& ixx, double& ixy, double& iyy) const;

    // Get the jacobian at a particular x,y
    // J = ( dudx  dudy )
    //     ( dvdx  dvdy )
    void getDistortion(Position p, DMatrix& J) const;
    void getDistortion(Position p, DSmallMatrix22& J) const;
    void getDistortion(
        Position p, 
        double& dudx, double& dudy, double& dvdx, double& dvdy) const;

    // Calculate x,y such that u,v = u(x,y),v(x,y)
    // This uses a non-linear solver to solver for x,y.
    // The value of puv on input is used as an initial guess.
    // The return value indicates whether a solution was found.
    bool inverseTransform(Position pxy, Position& puv) const;

    // Make this Transformation an approximate inverse of t2 using 
    // Legendre polynomials u to the given order in x,y.
    // The resulting transformation will be defined over a square
    // region in (u,v) that corresponds to the region given by
    // bounds in (x,y).
    // It returns the valid bounds for the new inverse transformation.
    Bounds makeInverseOf(const Transformation& t2, 
                         const Bounds& bounds, int order);

private :

    bool _isRaDec;

    std::auto_ptr<Function2D> _u;
    std::auto_ptr<Function2D> _v;
    std::auto_ptr<Function2D> _dudx;
    std::auto_ptr<Function2D> _dudy;
    std::auto_ptr<Function2D> _dvdx;
    std::auto_ptr<Function2D> _dvdy;


};

#endif
