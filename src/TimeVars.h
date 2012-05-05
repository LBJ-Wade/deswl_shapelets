#ifndef TIME_VARS_H
#define TIME_VARS_H

#include <sys/time.h>

struct EllipseTimes
{
    EllipseTimes() : 
        _tinteg(0.0), _tcentroid(0.0), _tgamma(0.0),
        _tmu(0.0), _tfixflux(0.0), _tfinal(0.0) {}

    EllipseTimes& operator+=(const EllipseTimes& rhs)
    {
        _tinteg += rhs._tinteg;
        _tcentroid += rhs._tcentroid;
        _tgamma += rhs._tgamma;
        _tmu += rhs._tmu;
        _tfixflux += rhs._tfixflux;
        _tfinal += rhs._tfinal;
        return *this;
    }

    void reset()
    { 
        _tinteg = _tcentroid = _tgamma = _tmu = _tfixflux = _tfinal = 0.0;
    }

    double _tinteg, _tcentroid, _tgamma, _tmu, _tfixflux, _tfinal;
};

struct OverallFitTimes 
{
    OverallFitTimes() :
        _nf_range1(0), _nf_range2(0),
        _nf_pixFlag1(0), _nf_npix1(0), _nf_small(0),
        _nf_pixFlag2(0), _nf_npix2(0),
        _ns_centroid(0), _nf_centroid(0),
        _ns_native(0), _nf_native(0),
        _ns_mu(0), _nf_mu(0),
        _ns_gamma(0), _nf_gamma(0) {}

    void successCentroid(const EllipseTimes& et)
    {
        _ts_centroid += et;
        ++_ns_centroid;
    }

    void failCentroid(const EllipseTimes& et)
    {
        _tf_centroid += et;
        ++_nf_centroid;
    }

    void success_native(const EllipseTimes& et)
    {
        _ts_native += et;
        ++_ns_native;
    }

    void fail_native(const EllipseTimes& et)
    {
        _tf_native += et;
        ++_nf_native;
    }

    void successMu(const EllipseTimes& et)
    {
        _ts_mu += et;
        ++_ns_mu;
    }

    void failMu(const EllipseTimes& et)
    {
        _tf_mu += et;
        ++_nf_mu;
    }

    void successGamma(const EllipseTimes& et)
    {
        _ts_gamma += et;
        ++_ns_gamma;
    }

    void failGamma(const EllipseTimes& et)
    {
        _tf_gamma += et;
        ++_nf_gamma;
    }

    OverallFitTimes& operator+=(const OverallFitTimes& rhs)
    {
        _ts_centroid += rhs._ts_centroid;
        _tf_centroid += rhs._tf_centroid;
        _ts_native += rhs._ts_native;
        _tf_native += rhs._tf_native;
        _ts_mu += rhs._ts_mu;
        _tf_mu += rhs._tf_mu;
        _ts_gamma += rhs._ts_gamma;
        _tf_gamma += rhs._tf_gamma;

        _nf_range1 += rhs._nf_range1;
        _nf_range2 += rhs._nf_range2;
        _nf_pixflag1 += rhs._nf_pixflag1;
        _nf_npix1 += rhs._nf_npix1;
        _nf_small += rhs._nf_small;
        _nf_pixflag2 += rhs._nf_pixflag2;
        _nf_npix2 += rhs._nf_npix2;
        _ns_centroid += rhs._ns_centroid;
        _nf_centroid += rhs._nf_centroid;
        _ns_native += rhs._ns_native;
        _nf_native += rhs._nf_native;
        _ns_mu += rhs._ns_mu;
        _nf_mu += rhs._nf_mu;
        _ns_gamma += rhs._ns_gamma;
        _nf_gamma += rhs._nf_gamma;
        return *this;
    }

    void write(std::ostream& os) const
    {
        os<<"N_Rejected for range error - distortion = "<<_nf_range1<<std::endl;
        os<<"N_Rejected for range error - psf interp = "<<_nf_range2<<std::endl;
        os<<"N_Rejected for pixel flag (#1) = "<<_nf_pixflag1<<std::endl;
        os<<"N_Rejected for too few pixels (#1) = "<<_nf_npix1<<std::endl;

        os<<"Centroid step:\n";
        os<<"  N_Success = "<<_ns_centroid<<std::endl;
        os<<"  Times = "<<_ts_centroid._tInteg<<"  "<<_ts_centroid._tcentroid<<"  "<<
            _ts_centroid._tgamma<<"  "<<_ts_centroid._tmu<<"  "<<
            _ts_centroid._tfixflux<<"  "<<_ts_centroid._tfinal<<std::endl;
        os<<"  N_Fail = "<<_nf_centroid<<std::endl;
        os<<"  Times = "<<_tf_centroid._tInteg<<"  "<<_tf_centroid._tcentroid<<"  "<<
            _tf_centroid._tgamma<<"  "<<_tf_centroid._tmu<<"  "<<
            _tf_centroid._tfixflux<<"  "<<_tf_centroid._tfinal<<std::endl;

        os<<"Native fits:\n";
        os<<"  N_Success = "<<_ns_native<<std::endl;
        os<<"  Times = "<<_ts_native._tInteg<<"  "<<_ts_native._tcentroid<<"  "<<
            _ts_native._tgamma<<"  "<<_ts_native._tmu<<"  "<<
            _ts_native._tfixflux<<"  "<<_ts_native._tfinal<<std::endl;
        os<<"  N_Fail = "<<_nf_native<<std::endl;
        os<<"  Times = "<<_tf_native._tInteg<<"  "<<_tf_native._tcentroid<<"  "<<
            _tf_native._tgamma<<"  "<<_tf_native._tmu<<"  "<<
            _tf_native._tfixflux<<"  "<<_tf_native._tfinal<<std::endl;

        os<<"N_Rejected for being too small = "<<_nf_small<<std::endl;
        os<<"N_Rejected for pixel flag (#2) = "<<_nf_pixflag2<<std::endl;
        os<<"N_Rejected for too few pixels (#2) = "<<_nf_npix2<<std::endl;

        os<<"Mu fits:\n";
        os<<"  N_Success = "<<_ns_mu<<std::endl;
        os<<"  Times = "<<_ts_mu._tInteg<<"  "<<_ts_mu._tcentroid<<"  "<<
            _ts_mu._tgamma<<"  "<<_ts_mu._tmu<<"  "<<
            _ts_mu._tfixflux<<"  "<<_ts_mu._tfinal<<std::endl;
        os<<"  N_Fail = "<<_nf_mu<<std::endl;
        os<<"  Times = "<<_tf_mu._tInteg<<"  "<<_tf_mu._tcentroid<<"  "<<
            _tf_mu._tgamma<<"  "<<_tf_mu._tmu<<"  "<<
            _tf_mu._tfixflux<<"  "<<_tf_mu._tfinal<<std::endl;

        os<<"Gamma fits:\n";
        os<<"  N_Success = "<<_ns_gamma<<std::endl;
        os<<"  Times = "<<_ts_gamma._tInteg<<"  "<<_ts_gamma._tcentroid<<"  "<<
            _ts_gamma._tgamma<<"  "<<_ts_gamma._tmu<<"  "<<
            _ts_gamma._tfixflux<<"  "<<_ts_gamma._tfinal<<std::endl;
        os<<"  N_Fail = "<<_nf_gamma<<std::endl;
        os<<"  Times = "<<_tf_gamma._tInteg<<"  "<<_tf_gamma._tcentroid<<"  "<<
            _tf_gamma._tgamma<<"  "<<_tf_gamma._tmu<<"  "<<
            _tf_gamma._tfixflux<<"  "<<_tf_gamma._tfinal<<std::endl;
    }

    EllipseTimes _ts_centroid;
    EllipseTimes _tf_centroid;
    EllipseTimes _ts_native;
    EllipseTimes _tf_native;
    EllipseTimes _ts_mu;
    EllipseTimes _tf_mu;
    EllipseTimes _ts_gamma;
    EllipseTimes _tf_gamma;

    int _nf_range1;
    int _nf_range2;
    int _nf_pixflag1;
    int _nf_npix1;
    int _nf_small;
    int _nf_pixflag2;
    int _nf_npix2;
    int _ns_centroid;
    int _nf_centroid;
    int _ns_native;
    int _nf_native;
    int _ns_mu;
    int _nf_mu;
    int _ns_gamma;
    int _nf_gamma;
};

inline std::ostream& operator<<(std::ostream& os, const OverallFitTimes& t)
{ t.write(os); return os; }

#endif
