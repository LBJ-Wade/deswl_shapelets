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
        Assert(b.getOrder() == _psfOrder);
        Assert(b.size() == _avePsf->size());
        b.setSigma(_sigma);
        interpolateVector(pos,b.vec().View());
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

class FittedPsfAtXY : public AssignableToBVec
{

public :

    FittedPsfAtXY(const FittedPsf& psf, Position pos) :
        _psf(psf), _pos(pos) 
    {}

    size_t size() const { return _psf._avePsf->size(); }
    int getOrder() const { return _psf.getPsfOrder(); }
    double getSigma() const { return _psf.getSigma(); }

    void AssignTo(BVec& b) const
    { _psf.interpolate(_pos,b); }

private :

    const FittedPsf& _psf;
    Position _pos;
};

inline FittedPsfAtXY FittedPsf::operator()(Position pos) const
{ return FittedPsfAtXY(*this,pos); }

#endif
