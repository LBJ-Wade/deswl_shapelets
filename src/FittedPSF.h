#ifndef FitPSF_H
#define FitPSF_H

#include "BVec.h"
#include "Bounds.h"
#include "ConfigFile.h"
#include "dbg.h"
#include <vector>
#include <iostream>

class FittedPSFAtXY;

class FittedPSF {

  public :

    FittedPSF() {};
    FittedPSF(const std::vector<BVec>& psf,
	const std::vector<long>& flagvec,
	const std::vector<Position>& pos,
	const std::vector<double>& nu,
	double sigma_p, const ConfigFile& params);

    int GetOrder() const { return psforder; }
    int GetFitOrder() const { return fitorder; }
    int GetNpca() const { return npca; }
    double GetSigma() const { return sigma; }

    double GetXMin() const {return bounds.GetXMin();}
    double GetXMax() const {return bounds.GetXMax();}
    double GetYMin() const {return bounds.GetYMin();}
    double GetYMax() const {return bounds.GetYMax();}

    void Write(const ConfigFile& params) const;
    void Write(std::ostream& os) const;
    void WriteFits(std::string file, const ConfigFile& params) const;

    void Read(const ConfigFile& params);
    void Read(std::istream& os);
    void ReadFits(std::string file, int hdu, const ConfigFile& params);

    void Interpolate(Position pos, BVec& b) const
    {
      Assert(avepsf.get());
      Assert(V.get());
      Assert(b.GetOrder() == psforder);
      Assert(b.size() == avepsf->size());
      Assert(b.GetSigma() == sigma);
      InterpolateVector(pos,b.View());
    }

    // Instead of:
    // psf.Interpolate(pos,b);
    // This next construct allows you to instead write:
    // b = psf(pos);
    // Both do the same thing.  I just like the second notation better.
    friend class FittedPSFAtXY;
    inline FittedPSFAtXY operator()(Position pos) const; // below...

  private :

    void InterpolateVector(Position pos,
	const tmv::VectorView<double>& b) const;

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

inline std::ostream& operator<<(std::ostream& os, const FittedPSF& p)
{ p.Write(os); return os; }


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
