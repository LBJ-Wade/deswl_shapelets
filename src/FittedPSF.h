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

    FittedPSF(std::vector<BVec*>& psf, const std::vector<Position>& pos,
	const std::vector<double>& nu,
	double sigma_p, ConfigFile& params);
    FittedPSF(std::istream& is);

    int GetOrder() const { return psforder; }
    double GetSigma() const { return sigma; }


    void Write(std::ostream& os) const;

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
    inline FittedPSFAtXY operator()(Position pos); // below...

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
    std::auto_ptr<tmv::Matrix<double> > V;
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

    void AssignToV(const tmv::VectorView<std::complex<double> >& b) const
    { Assert(false); }

  private :

    const FittedPSF& psf;
    Position pos;
};

inline FittedPSFAtXY FittedPSF::operator()(Position pos)
{ return FittedPSFAtXY(*this,pos); }

#endif