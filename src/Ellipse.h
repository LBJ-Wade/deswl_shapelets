#ifndef ELLIPSE_H
#define ELLIPSE_H

#include <complex>
#include <vector>
#include <ostream>
#include "Pixel.h"
#include "BVec.h"
#include "MyMatrix.h"

class Ellipse 
{

public :

    Ellipse() :
        _cen(0.), _gamma(0.), _mu(0.),
        _isFixedCen(false), _isFixedGamma(false), _isFixedMu(false), 
        _shouldDoTimings(false) {}

    Ellipse(double x, double y, double g1, double g2, double mu) :
        _cen(x,y), _gamma(g1,g2), _mu(mu), 
        _isFixedCen(false), _isFixedGamma(false), _isFixedMu(false), 
        _shouldDoTimings(false) {}

    Ellipse(double vals[]) :
        _cen(vals[0],vals[1]), _gamma(vals[2],vals[3]), _mu(vals[4]),
        _isFixedCen(false), _isFixedGamma(false), _isFixedMu(false),
        _shouldDoTimings(false) {}

    bool measure(const std::vector<PixelList>& pix, 
                 const std::vector<BVec>& psf,
                 int order, double sigma, bool shouldUseInteg, long& flag, 
                 DMatrix* cov=0, 
                 BVec* bret=0, DMatrix* bcov=0);
    bool measure(const std::vector<PixelList>& pix, 
                 int order, double sigma, bool shouldUseInteg, long& flag, 
                 DMatrix* cov=0, 
                 BVec* bret=0, DMatrix* bcov=0);

    bool findRoundFrame(const BVec& b, int order, long& flag, DMatrix* cov=0);

    void crudeMeasure(const PixelList& pix, double sigma);
    void crudeMeasure(
        const std::vector<PixelList>& pix, double sigma);

    void peakCentroid(const PixelList& pix, double maxr);

    // order parameter is allowed to be less than the order of bret,
    // in which case the rest of the vector is set to zeros.
    // order2 is the order for the intermediate steps in the calculation.
    void measureShapelet(
        const std::vector<PixelList>& pix, 
        const std::vector<BVec>& psf, BVec& bret, int order, int order2,
        DMatrix* bcov=0) const;
    void measureShapelet(
        const std::vector<PixelList>& pix, BVec& bret, int order, int order2,
        DMatrix* bcov=0) const;

    void altMeasureShapelet(
        const std::vector<PixelList>& pix, 
        const std::vector<BVec>& psf, BVec& bret, int order, int order2,
        double pixScale, DMatrix* bcov=0) const;
    void altMeasureShapelet(
        const std::vector<PixelList>& pix, BVec& bret, int order, int order2,
        double pixScale, DMatrix* bcov=0) const;

    void write(std::ostream& os) const
    { os << _cen<<" "<<_gamma<<" "<<_mu; }

    std::complex<double> getCen() const 
    { return _cen; }
    std::complex<double> getGamma() const 
    { return _gamma; }
    double getMu() const { return _mu; }

    void setCen(const std::complex<double>& cen)
    { _cen = cen; }
    void setGamma(const std::complex<double>& gamma)
    { _gamma = gamma; }
    void setMu(const double mu) { _mu = mu; }

    void shiftBy(const std::complex<double>& cen,
                 const std::complex<double>& gamma,
                 const double mu);

    void fixCen() { _isFixedCen = true; }
    void fixGam() { _isFixedGamma = true; }
    void fixMu() { _isFixedMu = true; }

    void unfixCen() { _isFixedCen = false; }
    void unfixGam() { _isFixedGamma = false; }
    void unfixMu() { _isFixedMu = false; }

    bool isFixedCen() const { return _isFixedCen; }
    bool isFixedGamma() const { return _isFixedGamma; }
    bool isFixedMu() const { return _isFixedMu; }

    void doTimings() { _shouldDoTimings = true; }

private :

    bool doMeasure(
        const std::vector<PixelList>& pix, 
        const std::vector<BVec>* psf, int order, double sigma,
        bool shouldUseInteg, long& flag, DMatrix* cov=0, 
        BVec* bret=0, DMatrix* bcov=0);

    void doMeasureShapelet(
        const std::vector<PixelList>& pix, 
        const std::vector<BVec>* psf, BVec& bret, int order, int order2,
        DMatrix* bcov=0) const;

    void doAltMeasureShapelet(
        const std::vector<PixelList>& pix, 
        const std::vector<BVec>* psf, BVec& bret, int order, int order2,
        double pixScale, DMatrix* bcov=0) const;

    std::complex<double> _cen;
    std::complex<double> _gamma;
    double _mu; 

    bool _isFixedCen,_isFixedGamma,_isFixedMu;

    bool _shouldDoTimings;

};

inline std::ostream& operator<<(std::ostream& os, const Ellipse& s)
{ s.write(os); return os; }

#endif
