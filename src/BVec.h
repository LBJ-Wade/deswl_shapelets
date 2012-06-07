#ifndef BVEC_H
#define BVEC_H

#include <complex>
#include <vector>
#include "dbg.h"
#include "MyMatrix.h"
#include "Bounds.h"

// Just do declarations for Image, Transformation:
template <typename T>
class Image;

class Transformation;

class BVec;

class AssignableToBVec
{
public:
    virtual void assignTo(BVec& bvec) const = 0;
    virtual int getOrder() const = 0;
    virtual double getSigma() const = 0;
    virtual ~AssignableToBVec() {}
};

class BVec : public AssignableToBVec
{

public :

    BVec(int order, double sigma) :
        _order(order), _sigma(sigma),
        _b((_order+1)*(_order+2)/2) 
    { _b.setZero(); }

    BVec(int order, double sigma, double* bvec) :
        _order(order), _sigma(sigma),
#ifdef USE_TMV
        _b( (_order+1)*(_order+2)/2 )
#else
        _b(DVector::Map(bvec,(_order+1)*(_order+2)/2))
#endif
    {
#ifdef USE_TMV
        std::copy(bvec,bvec+_b.size(),_b.begin());
#endif
    }

    BVec(int order, double sigma, const DVector& bvec) :
        _order(order), _sigma(sigma), _b(bvec) 
    { Assert(int(bvec.size()) == (order+1)*(order+2)/2); }

    BVec(const BVec& rhs) :
        _order(rhs._order), _sigma(rhs._sigma), _b(rhs._b)
    {}

    BVec(const AssignableToBVec& rhs) :
        _order(rhs.getOrder()), _sigma(rhs.getSigma()), 
        _b((_order+1)*(_order+2)/2) 
    { rhs.assignTo(*this); }

    ~BVec() {}

    BVec& operator=(const AssignableToBVec& rhs);
    BVec& operator=(const BVec& rhs);

    void assignTo(BVec& bvec) const;

    DVector& vec() { return _b; }
    const DVector& vec() const { return _b; }
    double& operator()(int i) { return _b(i); }
    double operator()(int i) const { return _b(i); }

    int getOrder() const { return _order; }
    double getSigma() const { return _sigma; }
    const DVector& getValues() const { return _b; }
    int size() const { return _b.size(); }

    void setSigma(double sigma) { _sigma = sigma; }

    // For this one v is allowed to be smaller than size, 
    // and the rest of the values are filled in as 0's.
    void setValues(const DVector& v);

    // For this one, the sizes must match.
    // b.setValues(v) is equivalent to b.vec() = v.
    template <typename V>
    void setValues(const V& v) 
    { 
        Assert(v.size() == size());
        _b = v; 
    }

    void normalize() { _b /= _b(0); }

    void conjugateSelf();

    void makeImage(
        Image<double>& im,
        const Position cen, double sky, const Transformation& trans,
        double x_offset, double y_offset) const;

private :

    int _order;
    double _sigma;
    DVector _b;

};

inline std::ostream& operator<<(std::ostream& os, const BVec& b)
{ os << b.getOrder()<<"  "<<b.getSigma()<<"  "<<b.vec(); return os; }

// NB: All of the following calculate and augment functions assume that
// the input matrix has already been zeroed before calling the function.
// This is because the sparsity of the matrices maintains its structure
// for different values of mu, g, theta, and z.  So it is faster to 
// just overwrite the locations that need to be written and skip the zeroes
// each time you call these functions.
void CalculateZTransform(
    std::complex<double> z, int order1, int order2, DMatrix& T);
inline void CalculateZTransform(std::complex<double> z, int order, DMatrix& T)
{ CalculateZTransform(z,order,order,T); }
void AugmentZTransformCols(
    std::complex<double> z, int order1, int order2, DMatrix& T);
inline void AugmentZTransformCols(
    std::complex<double> z, int order, DMatrix& T)
{ AugmentZTransformCols(z,order,order,T); }
void ApplyZ(std::complex<double> z, BVec& b);

void CalculateMuTransform(double mu, int order1, int order2, DMatrix& D);
inline void CalculateMuTransform(double mu, int order, DMatrix& D)
{ CalculateMuTransform(mu,order,order,D); }
void AugmentMuTransformRows(double mu, int order1, int order2, DMatrix& D);
void AugmentMuTransformCols(double mu, int order1, int order2, DMatrix& D);
inline void AugmentMuTransformRows(double mu, int order, DMatrix& D)
{ AugmentMuTransformRows(mu,order,order,D); }
inline void AugmentMuTransformCols(double mu, int order, DMatrix& D)
{ AugmentMuTransformCols(mu,order,order,D); }
void ApplyMu(double mu, BVec& b);

void CalculateThetaTransform(
    double theta, int order1, int order2, DBandMatrix& R);
inline void CalculateThetaTransform(double theta, int order, DBandMatrix& R)
{ CalculateThetaTransform(theta,order,order,R); }
void ApplyTheta(double theta, BVec& b);

void CalculateGTransform(
    std::complex<double> g, int order1, int order2, DMatrix& S);
inline void CalculateGTransform(std::complex<double> g, int order, DMatrix& S)
{ CalculateGTransform(g,order,order,S); }
void AugmentGTransformCols(
    std::complex<double> g, int order1, int order2, DMatrix& S);
inline void AugmentGTransformCols(
    std::complex<double> g, int order, DMatrix& S)
{ AugmentGTransformCols(g,order,order,S); }
void ApplyG(std::complex<double> g, BVec& b);

void CalculatePsfConvolve(
    const BVec& bpsf, int order1, int order2, double sigma, DMatrix& C);
inline void CalculatePsfConvolve(
    const BVec& bpsf, int order, double sigma, DMatrix& C)
{ CalculatePsfConvolve(bpsf,order,order,sigma,C); }
void ApplyPsf(const BVec& bpsf, BVec& b);

#endif
