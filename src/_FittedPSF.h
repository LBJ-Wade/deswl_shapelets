#ifndef FitPsf_H
#define FitPsf_H

#include <vector>
#include <iostream>
#include "BVec.h"
#include "Bounds.h"
#include "ConfigFile.h"
#include "dbg.h"
#include "PsfCatalog.h"

class FittedPsfAtXY;

class FittedPsf {

public :

    // Make from PsfCatalog
    FittedPsf(PsfCatalog& psfcat, const ConfigFile& params);

    // Read from file
    FittedPsf(const ConfigFile& params);
    // With this one we don't have to use the whole root= thing
    FittedPsf(const ConfigFile& params, std::string file);

    int getPsfOrder() const { return _psfOrder; }
    int getPsfSize() const { return (_psfOrder+1)*(_psfOrder+2)/2; }
    int getFitOrder() const { return _fitOrder; }
    int getFitSize() const { return _fitSize; }
    int getNpca() const { return _nPca; }
    double getSigma() const { return _sigma; }

    double getXMin() const { return _bounds.getXMin(); }
    double getXMax() const { return _bounds.getXMax(); }
    double getYMin() const { return _bounds.getYMin(); }
    double getYMax() const { return _bounds.getYMax(); }
    const Bounds& getBounds() const { return _bounds; }

    void write() const;
    void writeAscii(std::string file) const;
    void writeFits(std::string file) const;
    void writeFitsOld(std::string file) const;

    void read();
    void read(std::string file);
    void readAscii(std::string file);
    void readFits(std::string file);

    void interpolate(Position pos, BVec& b) const
    {
        Assert(_avePsf.get());
        Assert(_mV.get());
        Assert(b.GetOrder() == _psfOrder);
        Assert(b.size() == _avePsf->size());
        b.SetSigma(_sigma);
        interpolateVector(pos,b.View());
    }

    // This next construct with FittedPsfAtXY allows you to write:
    // b = psf(pos);
    // instead of:
    // psf.interpolate(pos,b);
    // Both do the same thing.  I just like the first notation better.
    friend class FittedPsfAtXY;
    inline FittedPsfAtXY operator()(Position pos) const; // below...

private :

    void interpolateVector(
        Position pos, const tmv::VectorView<double>& b) const;

    const ConfigFile& _params;

    int _psfOrder;
    double _sigma;
    int _fitOrder;
    int _fitSize;
    int _nPca;
    Bounds _bounds;
    std::auto_ptr<tmv::Vector<double> > _avePsf;
    std::auto_ptr<tmv::Matrix<double,tmv::RowMajor> > _mV;
    std::auto_ptr<tmv::Matrix<double> > _f;
};

class FittedPsfAtXY : public tmv::AssignableToVector<double> {

public :

    FittedPsfAtXY(const FittedPsf& _psf, Position _pos) :
        psf(_psf), pos(_pos) {}

    size_t size() const { return psf._avePsf->size(); }
    void AssignToV(const tmv::VectorView<double>& b) const
    { psf.interpolateVector(pos,b); }

    void AssignToV(const tmv::VectorView<std::complex<double> >& ) const
    { Assert(false); }

private :

    const FittedPsf& psf;
    Position pos;
};

inline FittedPsfAtXY FittedPsf::operator()(Position pos) const
{ return FittedPsfAtXY(*this,pos); }

#endif
