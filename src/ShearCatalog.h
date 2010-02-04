#ifndef ShearCatalog_H
#define ShearCatalog_H

#include <vector>
#include "BVec.h"
#include "Bounds.h"
#include "ConfigFile.h"
#include "Name.h"
#include "MyMatrix.h"
#include "Transformation.h"
#include "InputCatalog.h"
#include "Image.h"
#include "Log.h"
#include "TimeVars.h"
#include "FittedPsf.h"

class ShearCatalog 
{
public :

    // Make from incat, trans
    ShearCatalog(
        const InputCatalog& inCat, const Transformation& trans,
        const ConfigFile& params);

    // Read from file
    ShearCatalog(const ConfigFile& params);

    size_t size() const { return _id.size(); }

    void write() const;
    void writeFits(std::string file) const;
    void writeAscii(std::string file, std::string delim = "  ") const;

    void read();
    void readFits(std::string file);
    void readAscii(std::string file, std::string delim = "  ");

    int measureShears(
        const Image<double>& im,
        const Image<double>* weightIm, const Transformation& trans,
        const FittedPsf& fitPsf, ShearLog& log);

    const std::vector<long> getIdList() const { return _id; }
    const std::vector<Position> getPosList() const { return _pos; }
    const std::vector<double> getSkyList() const { return _sky; }
    const std::vector<double> getNoiseList() const { return _noise; }
    const std::vector<long> getFlagsList() const { return _flags; }
    const std::vector<Position> getSkyPosList() const { return _skyPos; }
    const std::vector<std::complex<double> > getShearList() const 
    { return _shear; }
    const std::vector<double> getNuList() const { return _nu; }
    const std::vector<DSmallMatrix22>& getCovList() const 
    { return _cov; }
    const std::vector<BVec>& getShapeList() const { return _shape; }

    long getId(int i) const { return _id[i]; }
    Position getPos(int i) const { return _pos[i]; }
    double getSky(int i) const { return _sky[i]; }
    double getNoise(int i) const { return _noise[i]; }
    long getFlags(int i) const { return _flags[i]; }
    Position getSkyPos(int i) const { return _skyPos[i]; }
    std::complex<double> getShear(int i) const { return _shear[i]; }
    double getNu(int i) const { return _nu[i]; }
    const DSmallMatrix22& getCov(int i) const { return _cov[i]; }
    const BVec& getShape(int i) const { return _shape[i]; }

    const Bounds& getBounds() const { return _bounds; }
    const Bounds& getSkyBounds() const { return _skyBounds; }

private :

    std::vector<long> _id;
    std::vector<Position> _pos;
    std::vector<double> _sky;
    std::vector<double> _noise;
    std::vector<long> _flags;
    std::vector<Position> _skyPos;
    std::vector<std::complex<double> > _shear;
    std::vector<double> _nu;
    std::vector<DSmallMatrix22 > _cov;
    EIGEN_mutable std::vector<BVec> _shape;

    Bounds _bounds;
    Bounds _skyBounds;

    const ConfigFile& _params;

};

#endif
