#ifndef BVEC_H
#define BVEC_H

#ifdef USE_EIGEN
#include "Eigen/Core"
typedef Eigen::VectorXd DVector;
typedef Eigen::Block<DVector,Dynamic,1> DVectorView;
typedef Eigen::MatrixXd DMatrix;
typedef Eigen::Block<DMatrix> DMatrixView;
#else
#include "TMV.h"
typedef tmv::Vector<double> DVector;
typedef tmv::VectorView<double> DVectorView;
typedef tmv::Matrix<double> DMatrix;
typedef tmv::MatrixView<double> DMatrixView;
#endif
#include <complex>
#include <vector>
#include "dbg.h"

class BVec;

class AssignableToBVec
{
public:
    virtual void AssignTo(BVec& bvec) const = 0;
    virtual int getOrder() const = 0;
    virtual double getSigma() const = 0;
};

class BVec : public AssignableToBVec
{

public :

    BVec(int order, double sigma) :
        _order(order), _sigma(sigma),
        _b((_order+1)*(_order+2)/2,0.) 
    {}

    BVec(int order, double sigma, double* bvec) :
        _order(order), _sigma(sigma),
        _b((_order+1)*(_order+2)/2,bvec) 
    {} 

    BVec(const BVec& rhs) :
        _order(rhs._order), _sigma(rhs._sigma), _b(rhs._b)
    {}

    ~BVec() {}

    BVec& operator=(const AssignableToBVec& rhs);
    BVec& operator=(const BVec& rhs);

    void AssignTo(BVec& bvec) const;

    DVector& vec() { return _b; }
    const DVector& vec() const { return _b; }
    double& operator()(int i) { return _b(i); }
    double operator()(int i) const { return _b(i); }

    int getOrder() const { return _order; }
    double getSigma() const { return _sigma; }
    const DVector& getValues() const { return _b; }
    size_t size() const { return _b.size(); }

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

private :

    int _order;
    double _sigma;
    DVector _b;

};

void calculateZTransform(
    std::complex<double> z, int order, const DMatrixView& T);
inline void calculateZTransform(
    std::complex<double> z, int order, DMatrix& T)
{ calculateZTransform(z,order,T.View()); }
void applyZ(std::complex<double> z, BVec& b);

void calculateMuTransform(double mu, int order, const DMatrixView& D);
void augmentMuTransformRows(double mu, int order, const DMatrixView& D);
inline void calculateMuTransform(double mu, int order, DMatrix& D)
{ calculateMuTransform(mu,order,D.View()); }
void applyMu(double mu, BVec& b);

void calculateGTransform(
    std::complex<double> g, int order, const DMatrixView& S);
void augmentGTransformCols(
    std::complex<double> g, int order, const DMatrixView& S);
inline void calculateGTransform(std::complex<double> g, int order, DMatrix& S)
{ calculateGTransform(g,order,S.View()); }
void applyG(std::complex<double> g, BVec& b);

void calculatePsfConvolve(
    const BVec& bpsf, int order, double sigma, const DMatrixView& C);
inline void calculatePsfConvolve(
    const BVec& bpsf, int order, double sigma, DMatrix& C)
{ calculatePsfConvolve(bpsf,order,sigma,C.View()); }
void applyPsf(const BVec& bpsf, BVec& b);

#endif
