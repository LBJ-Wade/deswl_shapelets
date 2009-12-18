#ifndef MeasureShearAlgo_H
#define MeasureShearAlgo_H

#include <vector>
#include "BVec.h"
#include "TMV_Small.h"
#include "Transformation.h"
#include "Image.h"
#include "Log.h"
#include "FittedPsf.h"
#include "Pixel.h"

void measureSingleShear(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans, const FittedPsf& fitPsf,
    double noise, double gain, const Image<double>* weightIm, 
    double galAperture, double maxAperture,
    int galOrder, int galOrder2,
    double fPsf, double minGalSize,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    tmv::SmallMatrix<double,2,2>& shearCov, BVec& shapelet,
    double& nu, long& flag);

void measureSingleShear1(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans, const std::vector<BVec>& psf,
    double noise, double gain, const Image<double>* weightIm, 
    double galAperture, double maxAperture,
    int galOrder, int galOrder2,
    double fPsf, double minGalSize,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    tmv::SmallMatrix<double,2,2>& shearCov, BVec& shapelet,
    double& nu, long& flag);

void measureMultiShear(
    const Position& cen, 
    const std::vector<PixelList>& allPix,
    const std::vector<BVec>& psf,
    double galAperture, double maxAperture,
    int galOrder, int galOrder2,
    double fPsf, double minGalSize,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    tmv::SmallMatrix<double,2,2>& shearCov, BVec& shapelet,
    double& nu, long& flag);

#endif
#ifndef MultiShearCatalog_H
#define MultiShearCatalog_H

#include <vector>
#include <string>
#include "TMV.h"

#include "CoaddCatalog.h"
#include "Bounds.h"
#include "ConfigFile.h"
#include "Image.h"
#include "Pixel.h"
#include "Transformation.h"
#include "FittedPsf.h"

class MultiShearCatalog 
{

public :

    MultiShearCatalog(const CoaddCatalog& coaddcat, const ConfigFile& params);
    MultiShearCatalog(const ConfigFile& params);
    ~MultiShearCatalog();

    size_t size() const { return _skyPos.size(); }
    size_t getNImages() const { return _imageFileList.size(); }

    void resize(int n);

    // Read the srclist file
    void readFileLists();

    // Get a set of bounds with a maximum linear extent in either direction
    // of side arcminutes on a side.
    std::vector<Bounds> splitBounds(double side);

    // Get pixel lists for the component images/catalogs
    int getPixels(const Bounds& b);
    void getImagePixelLists(int fnum, const Bounds& b);

    // Measure the shears
    int measureMultiShears(const Bounds& b, ShearLog& log);

    // Write output
    void write() const;
    void writeFits(std::string file) const;
    void writeAscii(std::string file, std::string delim = "  ") const;

    // Read from file
    void read();
    void readFits(std::string file);
    void readAscii(std::string file, std::string delim = "  ");

    // Calculate the current memory footprint of the entire structure in MB.
    // (Optionally output some info to os.)
    double calculateMemoryFootprint(bool getmax=false) const;

    const std::vector<long> getIdList() const { return _id; }
    const std::vector<Position> getSkyPosList() const { return _skyPos; }
    const std::vector<long> getFlagsList() const { return _flags; }
    const std::vector<std::complex<double> > getShearList() const 
    { return _shear; }
    const std::vector<double> getNuList() const { return _nu; }
    const std::vector<tmv::SmallMatrix<double,2,2> >& getCovList() const 
    { return _cov; }
    const std::vector<BVec>& getShapeList() const { return _shape; }

    long getId(int i) const { return _id[i]; }
    Position getSkyPos(int i) const { return _skyPos[i]; }
    long getFlags(int i) const { return _flags[i]; }
    std::complex<double> getShear(int i) const { return _shear[i]; }
    double getNu(int i) const { return _nu[i]; }
    const tmv::SmallMatrix<double,2,2>& getCov(int i) const { return _cov[i]; }
    const BVec& getShape(int i) const { return _shape[i]; }

    const Bounds& getSkyBounds() const { return _skyBounds; }

private :

    // flags related to i/o and psf interpolation
    std::vector<long> _inputFlags;

    // number of images each object was found in
    std::vector<int> _nImagesFound;
    // number of images for which pixels were extracted
    std::vector<int> _nImagesGotPix;

    // These have an element for each coadd object.
    std::vector<long> _id;
    std::vector<Position> _skyPos;
    std::vector<long> _flags;
    std::vector<std::complex<double> > _shear;
    std::vector<double> _nu;
    std::vector<tmv::SmallMatrix<double,2,2> > _cov;
    std::vector<BVec> _shape;

    // There is no non-sky bounds of course, but to be consistent with 
    // the skybounds name in other catalogs, we keep that prefix here.
    Bounds _skyBounds; 

    // We copy the parameter info, since we mess around with the parameters
    // to set all the file names of the component images.
    // Hence, not a const reference the way I usually do it.
    ConfigFile _params;

    // For each coadd object, we have a vector with an element for each 
    // single-epoch (se) image it is found in, and each of those elements 
    // is PixelList (basically a vector of Pixels).
    std::vector<std::vector<PixelList> > _pixList;
    // The PSF for each single-epoch object.
    std::vector<std::vector<BVec> > _psfList;
    // The single-epoch shear if known.
    std::vector<std::vector<std::complex<double> > > _seShearList;
    // The single-epoch size if known.
    std::vector<std::vector<double> > _seSizeList;

    // These each have an element for each single-epoch image
    std::vector<std::string> _imageFileList;
    std::vector<std::string> _shearFileList;
    std::vector<std::string> _fitPsfFileList;
    std::vector<Bounds> _savedSeSkyBounds;
};

#endif