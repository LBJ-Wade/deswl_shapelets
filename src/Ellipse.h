#ifndef ELLIPSE_H
#define ELLIPSE_H

#include <complex>
#include <vector>
#include <ostream>
#include "Pixel.h"
#include "BVec.h"
#include "TimeVars.h"

class Ellipse {

  public :

    Ellipse() :
      cen(0.), gamma(0.), mu(0.),
      fixcen(false), fixgam(false), fixmu(false), f_psf(1.0),
      _shouldDoTimings(false) {}

    Ellipse(double _x, double _y, double _g1, double _g2, double _mu) :
      cen(_x,_y), gamma(_g1,_g2), mu(_mu), 
      fixcen(false), fixgam(false), fixmu(false), f_psf(1.0),
      _shouldDoTimings(false) {}

    Ellipse(double vals[]) :
      cen(vals[0],vals[1]), gamma(vals[2],vals[3]), mu(vals[4]),
      fixcen(false), fixgam(false), fixmu(false), f_psf(1.0),
      _shouldDoTimings(false) {}

    bool Measure(const std::vector<PixelList>& pix, 
	const std::vector<BVec>& psf,
	int order, double sigma, bool use_iteg, long& flag, 
	tmv::Matrix<double>* cov=0, 
	BVec* bret=0, tmv::Matrix<double>* bcov=0);
    bool Measure(const std::vector<PixelList>& pix, 
	int order, double sigma, bool use_integ, long& flag, 
	tmv::Matrix<double>* cov=0, 
	BVec* bret=0, tmv::Matrix<double>* bcov=0);

    void CrudeMeasure(const PixelList& pix, double sigma);
    void CrudeMeasure(
	const std::vector<PixelList>& pix, double sigma);

    void PeakCentroid(const PixelList& pix, double maxr);

    void MeasureShapelet(const std::vector<PixelList>& pix, 
	const std::vector<BVec>& psf, BVec& bret,
	tmv::Matrix<double>* bcov=0) const;
    void MeasureShapelet(const std::vector<PixelList>& pix, 
	BVec& bret, tmv::Matrix<double>* bcov=0) const;

    void Write(std::ostream& os) const
    { os << cen<<" "<<gamma<<" "<<mu; }

    std::complex<double> GetCen() const 
    { return cen; }
    std::complex<double> GetGamma() const 
    { return gamma; }
    double GetMu() const { return mu; }

    void SetCen(const std::complex<double>& _cen)
    { cen = _cen; }
    void SetGamma(const std::complex<double>& _gamma)
    { gamma = _gamma; }
    void SetMu(const double _mu) { mu = _mu; }

    void FixCen() { fixcen = true; }
    void FixGam() { fixgam = true; }
    void FixMu() { fixmu = true; }

    void UnFixCen() { fixcen = false; }
    void UnFixGam() { fixgam = false; }
    void UnFixMu() { fixmu = false; }

    bool isFixCen() const { return fixcen; }
    bool isFixGam() const { return fixgam; }
    bool isFixMu() const { return fixmu; }

    void SetFP(double fp) { f_psf = fp; }
    double GetFP() const { return f_psf; }

    void doTimings() { _shouldDoTimings = true; }
    void resetTimes() { _times.reset(); }
    const EllipseTimes& getTimes() { return _times; }

  private :

    bool DoMeasure(const std::vector<PixelList>& pix, 
	const std::vector<BVec>* psf, int order, double sigma,
	bool use_integ, long& flag, tmv::Matrix<double>* cov=0, 
	BVec* bret=0, tmv::Matrix<double>* bcov=0);
    void DoMeasureShapelet(const std::vector<PixelList>& pix, 
	const std::vector<BVec>* psf, BVec& bret,
	tmv::Matrix<double>* bcov=0) const;

    std::complex<double> cen,gamma;
    double mu; 

    bool fixcen,fixgam,fixmu;
    double f_psf;

    bool _shouldDoTimings;

    EllipseTimes _times;
};

inline std::ostream& operator<<(std::ostream& os, const Ellipse& s)
{ s.Write(os); return os; }

#endif
