#ifndef BVEC_H
#define BVEC_H

#include "TMV.h"
#include <complex>
#include <vector>

class BVec : public tmv::Vector<double> {

  public :

    BVec(int _order, double _sigma) :
      tmv::Vector<double>((_order+1)*(_order+2)/2,0.), 
      order(_order), sigma(_sigma)  {}

    BVec(int _order, double _sigma, double* bvec) :
      tmv::Vector<double>((_order+1)*(_order+2)/2,bvec), 
      order(_order), sigma(_sigma)  {}

    BVec(const BVec& b) :
      tmv::Vector<double>(b), order(b.order), sigma(b.sigma)  {}

    ~BVec() {}

    BVec& operator=(const BVec& rhs);
    BVec& operator=(const tmv::AssignableToVector<double>& rhs);

    int GetOrder() const { return order; }
    double GetSigma() const { return sigma; }
    void SetSigma(double _sigma) { sigma = _sigma; }
    void Normalize() { *this /= (*this)[0]; }

  private :

    int order;
    double sigma;
};

void ZTransform(std::complex<double> z, int order, 
    const tmv::MatrixView<double>& T);
inline void ZTransform(std::complex<double> z, int order,
    tmv::Matrix<double>& T)
{ ZTransform(z,order,T.View()); }
void ApplyZ(std::complex<double> z, BVec& b);

void MuTransform(double mu, int order, const tmv::MatrixView<double>& D);
void AugmentMuTransformRows(double mu, int order,
    const tmv::MatrixView<double>& D);
inline void MuTransform(double mu, int order, tmv::Matrix<double>& D)
{ MuTransform(mu,order,D.View()); }
void ApplyMu(double mu, BVec& b);

void GTransform(std::complex<double> g, int order,
    const tmv::MatrixView<double>& S);
void AugmentGTransformCols(std::complex<double> g, int order,
    const tmv::MatrixView<double>& S);
inline void GTransform(std::complex<double> g, int order,
    tmv::Matrix<double>& S)
{ GTransform(g,order,S.View()); }
void ApplyG(std::complex<double> g, BVec& b);

void PSFConvolve(const BVec& bpsf, int order, double sigma,
    const tmv::MatrixView<double>& C);
inline void PSFConvolve(const BVec& bpsf, int order, double sigma,
    tmv::Matrix<double>& C)
{ PSFConvolve(bpsf,order,sigma,C.View()); }
void ApplyPSF(const BVec& bpsf, BVec& b);

#endif
