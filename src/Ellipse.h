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

    Ellipse(std::complex<double> cen, std::complex<double> gamma,
            std::complex<double> mu) :
        _cen(cen), _gamma(gamma), _mu(mu), 
        _isFixedCen(false), _isFixedGamma(false), _isFixedMu(false), 
        _shouldDoTimings(false) {}

    Ellipse(double vals[]) :
        _cen(vals[0],vals[1]), _gamma(vals[2],vals[3]), _mu(vals[4]),
        _isFixedCen(false), _isFixedGamma(false), _isFixedMu(false),
        _shouldDoTimings(false) {}

    Ellipse(const Ellipse& e2) :
        _cen(e2.getCen()), _gamma(e2.getGamma()),
        _mu(e2.getMu(),e2.getTheta()),
        _isFixedCen(false), _isFixedGamma(false), _isFixedMu(false),
        _shouldDoTimings(false) {}

    Ellipse& operator=(const Ellipse& e2)
    { 
        _cen = e2.getCen();
        _gamma = e2.getGamma();
        _mu = std::complex<double>(e2.getMu(),e2.getTheta());
        return *this;
    }

    bool measure(
        const std::vector<PixelList>& pix, 
        const std::vector<BVec>& psf,
        int order, int order2, double sigma, long& flag, double thresh,
        DMatrix* cov=0, BVec* bret=0, DMatrix* bcov=0,
        std::vector<PixelList>* pixels_model=0);
    bool measure(
        const std::vector<PixelList>& pix, 
        int order, int order2, double sigma, long& flag, double thresh,
        DMatrix* cov=0, BVec* bret=0, DMatrix* bcov=0,
        std::vector<PixelList>* pixels_model=0);

    void crudeMeasure(const PixelList& pix, double sigma);
    void crudeMeasure(
        const std::vector<PixelList>& pix, double sigma);

    void peakCentroid(const PixelList& pix, double maxr);

    // order parameter is allowed to be less than the order of bret,
    // in which case the rest of the vector is set to zeros.
    // order2 is the order for the intermediate steps in the calculation.
    bool measureShapelet(
        const std::vector<PixelList>& pix, 
        const std::vector<BVec>& psf, BVec& bret, int& order, int order2,
        DMatrix* bcov=0) const;
    bool measureShapelet(
        const std::vector<PixelList>& pix, BVec& bret, int& order, int order2,
        DMatrix* bcov=0) const;

    bool altMeasureShapelet(
        const std::vector<PixelList>& pix, 
        const std::vector<BVec>& psf, BVec& bret, int order, int order2,
        double pixScale, DMatrix* bcov=0) const;
    bool altMeasureShapelet(
        const std::vector<PixelList>& pix, BVec& bret, int order, int order2,
        double pixScale, DMatrix* bcov=0) const;

    void write(std::ostream& os) const
    { os << _cen<<" "<<_gamma<<" "<<_mu; }

    std::complex<double> getCen() const 
    { return _cen; }
    std::complex<double> getGamma() const 
    { return _gamma; }
    double getMu() const { return real(_mu); }
    double getTheta() const { return imag(_mu); }

    void setCen(const std::complex<double>& cen) { _cen = cen; }
    void setGamma(const std::complex<double>& gamma) { _gamma = gamma; }
    void setMu(const std::complex<double>& mu) { _mu = mu; }

    void preShiftBy(const std::complex<double>& cen,
                    const std::complex<double>& gamma,
                    const std::complex<double>& mu);
    void postShiftBy(const std::complex<double>& cen,
                     const std::complex<double>& gamma,
                     const std::complex<double>& mu);

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
        const std::vector<BVec>* psf, int order, int order2, double sigma,
        long& flag, double thresh, DMatrix* cov=0, 
        BVec* bret=0, DMatrix* bcov=0, std::vector<PixelList>* pixels_model=0);

    bool doMeasureShapelet(
        const std::vector<PixelList>& pix, 
        const std::vector<BVec>* psf, BVec& bret, int& order, int order2,
        DMatrix* bcov=0) const;

    bool doAltMeasureShapelet(
        const std::vector<PixelList>& pix, 
        const std::vector<BVec>* psf, BVec& bret, int order, int order2,
        double pixScale, DMatrix* bcov=0) const;

    std::complex<double> _cen;
    std::complex<double> _gamma;
    std::complex<double> _mu; 

    bool _isFixedCen,_isFixedGamma,_isFixedMu;

    bool _shouldDoTimings;

};

inline std::ostream& operator<<(std::ostream& os, const Ellipse& s)
{ s.write(os); return os; }

#endif
