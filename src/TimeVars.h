#ifndef TIME_VARS_H
#define TIME_VARS_H

#include <sys/time.h>

struct EllipseTimes
{
    EllipseTimes() : 
        _tInteg(0.0), _tCentroid(0.0), _tGamma(0.0),
        _tMu(0.0), _tFixFlux(0.0), _tFinal(0.0) {}

    EllipseTimes& operator+=(const EllipseTimes& rhs)
    {
        _tInteg += rhs._tInteg;
        _tCentroid += rhs._tCentroid;
        _tGamma += rhs._tGamma;
        _tMu += rhs._tMu;
        _tFixFlux += rhs._tFixFlux;
        _tFinal += rhs._tFinal;
        return *this;
    }

    void reset()
    { 
        _tInteg = _tCentroid = _tGamma = _tMu = _tFixFlux = _tFinal = 0.0;
    }

    double _tInteg, _tCentroid, _tGamma, _tMu, _tFixFlux, _tFinal;
};

struct OverallFitTimes 
{
    OverallFitTimes() :
        _nfRange1(0), _nfRange2(0),
        _nfPixFlag1(0), _nfNPix1(0), _nfSmall(0),
        _nfPixFlag2(0), _nfNPix2(0),
        _nsCentroid(0), _nfCentroid(0),
        _nsNative(0), _nfNative(0),
        _nsMu(0), _nfMu(0),
        _nsGamma(0), _nfGamma(0) {}

    void successCentroid(const EllipseTimes& et)
    {
        _tsCentroid += et;
        ++_nsCentroid;
    }

    void failCentroid(const EllipseTimes& et)
    {
        _tfCentroid += et;
        ++_nfCentroid;
    }

    void successNative(const EllipseTimes& et)
    {
        _tsNative += et;
        ++_nsNative;
    }

    void failNative(const EllipseTimes& et)
    {
        _tfNative += et;
        ++_nfNative;
    }

    void successMu(const EllipseTimes& et)
    {
        _tsMu += et;
        ++_nsMu;
    }

    void failMu(const EllipseTimes& et)
    {
        _tfMu += et;
        ++_nfMu;
    }

    void successGamma(const EllipseTimes& et)
    {
        _tsGamma += et;
        ++_nsGamma;
    }

    void failGamma(const EllipseTimes& et)
    {
        _tfGamma += et;
        ++_nfGamma;
    }

    OverallFitTimes& operator+=(const OverallFitTimes& rhs)
    {
        _tsCentroid += rhs._tsCentroid;
        _tfCentroid += rhs._tfCentroid;
        _tsNative += rhs._tsNative;
        _tfNative += rhs._tfNative;
        _tsMu += rhs._tsMu;
        _tfMu += rhs._tfMu;
        _tsGamma += rhs._tsGamma;
        _tfGamma += rhs._tfGamma;

        _nfRange1 += rhs._nfRange1;
        _nfRange2 += rhs._nfRange2;
        _nfPixFlag1 += rhs._nfPixFlag1;
        _nfNPix1 += rhs._nfNPix1;
        _nfSmall += rhs._nfSmall;
        _nfPixFlag2 += rhs._nfPixFlag2;
        _nfNPix2 += rhs._nfNPix2;
        _nsCentroid += rhs._nsCentroid;
        _nfCentroid += rhs._nfCentroid;
        _nsNative += rhs._nsNative;
        _nfNative += rhs._nfNative;
        _nsMu += rhs._nsMu;
        _nfMu += rhs._nfMu;
        _nsGamma += rhs._nsGamma;
        _nfGamma += rhs._nfGamma;
        return *this;
    }

    void write(std::ostream& os) const
    {
        os<<"N_Rejected for range error - distortion = "<<_nfRange1<<std::endl;
        os<<"N_Rejected for range error - psf interp = "<<_nfRange2<<std::endl;
        os<<"N_Rejected for pixel flag (#1) = "<<_nfPixFlag1<<std::endl;
        os<<"N_Rejected for too few pixels (#1) = "<<_nfNPix1<<std::endl;

        os<<"Centroid step:\n";
        os<<"  N_Success = "<<_nsCentroid<<std::endl;
        os<<"  Times = "<<_tsCentroid._tInteg<<"  "<<_tsCentroid._tCentroid<<"  "<<
            _tsCentroid._tGamma<<"  "<<_tsCentroid._tMu<<"  "<<
            _tsCentroid._tFixFlux<<"  "<<_tsCentroid._tFinal<<std::endl;
        os<<"  N_Fail = "<<_nfCentroid<<std::endl;
        os<<"  Times = "<<_tfCentroid._tInteg<<"  "<<_tfCentroid._tCentroid<<"  "<<
            _tfCentroid._tGamma<<"  "<<_tfCentroid._tMu<<"  "<<
            _tfCentroid._tFixFlux<<"  "<<_tfCentroid._tFinal<<std::endl;

        os<<"Native fits:\n";
        os<<"  N_Success = "<<_nsNative<<std::endl;
        os<<"  Times = "<<_tsNative._tInteg<<"  "<<_tsNative._tCentroid<<"  "<<
            _tsNative._tGamma<<"  "<<_tsNative._tMu<<"  "<<
            _tsNative._tFixFlux<<"  "<<_tsNative._tFinal<<std::endl;
        os<<"  N_Fail = "<<_nfNative<<std::endl;
        os<<"  Times = "<<_tfNative._tInteg<<"  "<<_tfNative._tCentroid<<"  "<<
            _tfNative._tGamma<<"  "<<_tfNative._tMu<<"  "<<
            _tfNative._tFixFlux<<"  "<<_tfNative._tFinal<<std::endl;

        os<<"N_Rejected for being too small = "<<_nfSmall<<std::endl;
        os<<"N_Rejected for pixel flag (#2) = "<<_nfPixFlag2<<std::endl;
        os<<"N_Rejected for too few pixels (#2) = "<<_nfNPix2<<std::endl;

        os<<"Mu fits:\n";
        os<<"  N_Success = "<<_nsMu<<std::endl;
        os<<"  Times = "<<_tsMu._tInteg<<"  "<<_tsMu._tCentroid<<"  "<<
            _tsMu._tGamma<<"  "<<_tsMu._tMu<<"  "<<
            _tsMu._tFixFlux<<"  "<<_tsMu._tFinal<<std::endl;
        os<<"  N_Fail = "<<_nfMu<<std::endl;
        os<<"  Times = "<<_tfMu._tInteg<<"  "<<_tfMu._tCentroid<<"  "<<
            _tfMu._tGamma<<"  "<<_tfMu._tMu<<"  "<<
            _tfMu._tFixFlux<<"  "<<_tfMu._tFinal<<std::endl;

        os<<"Gamma fits:\n";
        os<<"  N_Success = "<<_nsGamma<<std::endl;
        os<<"  Times = "<<_tsGamma._tInteg<<"  "<<_tsGamma._tCentroid<<"  "<<
            _tsGamma._tGamma<<"  "<<_tsGamma._tMu<<"  "<<
            _tsGamma._tFixFlux<<"  "<<_tsGamma._tFinal<<std::endl;
        os<<"  N_Fail = "<<_nfGamma<<std::endl;
        os<<"  Times = "<<_tfGamma._tInteg<<"  "<<_tfGamma._tCentroid<<"  "<<
            _tfGamma._tGamma<<"  "<<_tfGamma._tMu<<"  "<<
            _tfGamma._tFixFlux<<"  "<<_tfGamma._tFinal<<std::endl;
    }

    EllipseTimes _tsCentroid;
    EllipseTimes _tfCentroid;
    EllipseTimes _tsNative;
    EllipseTimes _tfNative;
    EllipseTimes _tsMu;
    EllipseTimes _tfMu;
    EllipseTimes _tsGamma;
    EllipseTimes _tfGamma;

    int _nfRange1;
    int _nfRange2;
    int _nfPixFlag1;
    int _nfNPix1;
    int _nfSmall;
    int _nfPixFlag2;
    int _nfNPix2;
    int _nsCentroid;
    int _nfCentroid;
    int _nsNative;
    int _nfNative;
    int _nsMu;
    int _nfMu;
    int _nsGamma;
    int _nfGamma;
};

inline std::ostream& operator<<(std::ostream& os, const OverallFitTimes& t)
{ t.write(os); return os; }

#endif
