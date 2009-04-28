#ifndef FitPSF_H
#define FitPSF_H

#include <vector>
#include <iostream>
#include "BVec.h"
#include "Bounds.h"
#include "ConfigFile.h"
#include "dbg.h"
#include "PSFCatalog.h"

class FittedPSFAtXY;

class FittedPSF {

  public :

    // Make from PSFCatalog
    FittedPSF(const PSFCatalog& psfcat, const ConfigFile& params);

    // Read from file
    FittedPSF(const ConfigFile& params);

    int GetPSFOrder() const { return psforder; }
    int GetPSFSize() const { return (psforder+1)*(psforder+2)/2; }
    int GetFitOrder() const { return fitorder; }
    int GetFitSize() const { return fitsize; }
    int GetNpca() const { return npca; }
    double GetSigma() const { return sigma; }

    double GetXMin() const {return bounds.GetXMin();}
    double GetXMax() const {return bounds.GetXMax();}
    double GetYMin() const {return bounds.GetYMin();}
    double GetYMax() const {return bounds.GetYMax();}

    void Write() const;
    void WriteAscii(std::string file) const;
    void WriteFits(std::string file) const;
    void WriteFitsOld(std::string file) const;

    void Read();
    void ReadAscii(std::string file);
    void ReadFits(std::string file);

    void Interpolate(Position pos, BVec& b) const
    {
      Assert(avepsf.get());
      Assert(V.get());
      Assert(b.GetOrder() == psforder);
      Assert(b.size() == avepsf->size());
      b.SetSigma(sigma);
      InterpolateVector(pos,b.View());
    }

    // This next construct with FittedPSFAtXY allows you to write:
    // b = psf(pos);
    // instead of:
    // psf.Interpolate(pos,b);
    // Both do the same thing.  I just like the first notation better.
    friend class FittedPSFAtXY;
    inline FittedPSFAtXY operator()(Position pos) const; // below...

  private :

    void InterpolateVector(Position pos,
	const tmv::VectorView<double>& b) const;

    const ConfigFile& params;

    int psforder;
    double sigma;
    int fitorder;
    int fitsize;
    int npca;
    Bounds bounds;
    std::auto_ptr<tmv::Vector<double> > avepsf;
    std::auto_ptr<tmv::Matrix<double,tmv::RowMajor> > V;
    std::auto_ptr<tmv::Matrix<double> > f;
};

class FittedPSFAtXY : public tmv::AssignableToVector<double> {

  public :

    FittedPSFAtXY(const FittedPSF& _psf, Position _pos) :
      psf(_psf), pos(_pos) {}

    size_t size() const { return psf.avepsf->size(); }
    void AssignToV(const tmv::VectorView<double>& b) const
    { psf.InterpolateVector(pos,b); }

    void AssignToV(const tmv::VectorView<std::complex<double> >& ) const
    { Assert(false); }

  private :

    const FittedPSF& psf;
    Position pos;
};

inline FittedPSFAtXY FittedPSF::operator()(Position pos) const
{ return FittedPSFAtXY(*this,pos); }

#endif
