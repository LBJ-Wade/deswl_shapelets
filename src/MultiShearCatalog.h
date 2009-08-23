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
#include "FittedPSF.h"
#include "dbg.h"
#include "Name.h"
#include "WlVersion.h"

//#define USE_POOL

#ifdef USE_POOL
#define MULTI_BLOCK 10*1024*1024  
// 10 MB blocks for the pools other than PixelList

typedef pool_allocator<PixelList,MULTI_BLOCK> pool_alloc1;
typedef pool_allocator<std::vector<PixelList>,MULTI_BLOCK> pool_alloc2;
typedef pool_allocator<int,MULTI_BLOCK> pool_alloc3;
typedef pool_allocator<std::vector<int>,MULTI_BLOCK> pool_alloc4;
typedef pool_allocator<Position,MULTI_BLOCK> pool_alloc5;
typedef pool_allocator<std::vector<Position>,MULTI_BLOCK> pool_alloc6;

#define POOL1 ,pool_alloc1
#define POOL2 ,pool_alloc2
#define POOL3 ,pool_alloc3
#define POOL4 ,pool_alloc4
#define POOL5 ,pool_alloc5
#define POOL6 ,pool_alloc6
#else
#define POOL1
#define POOL2
#define POOL3
#define POOL4
#define POOL5
#define POOL6
#endif

class MultiShearCatalog 
{

  public :

    MultiShearCatalog(const CoaddCatalog& coaddcat, const ConfigFile& params);
    MultiShearCatalog(const ConfigFile& params);
    ~MultiShearCatalog();

    size_t size() const { return skypos.size(); }
    size_t NImages() const { return fitpsf.size(); }

    void Resize(int n);

    // Read the srclist file
    void ReadFileLists();

    // Get a set of bounds with a maximum linear extent in either direction
    // of side arcminutes on a side.
    std::vector<Bounds> SplitBounds(double side);

    // Get pixel lists for the component images/catalogs
    int GetPixels(const Bounds& b);
    void GetImagePixelLists(int fnum, const Bounds& b);

    // Measure the shears
    int MeasureMultiShears(const Bounds& b, ShearLog& log);

    // Write output
    void Write() const;
    void WriteFits(std::string file) const;
    void WriteAscii(std::string file, std::string delim = "  ") const;

    // Read from file
    void Read();
    void ReadFits(std::string file);
    void ReadAscii(std::string file, std::string delim = "  ");

    // These next two have length NImages(), rather than size()
    // (aka the number of objects)
    std::vector<const Transformation*> trans;
    std::vector<const FittedPSF*> fitpsf;

    // flags related to i/o and psf interpolation
    std::vector<long> input_flags;

    // number of images each object was found in
    std::vector<int> nimages_found;
    // number of images for which pixels were extracted
    std::vector<int> nimages_gotpix;

    std::vector<long> id;
    std::vector<Position> skypos;
    std::vector<double> sky;
    std::vector<double> noise;
    std::vector<long> flags;

    Bounds skybounds; 
    // There is no non-sky bounds of course, but to be consistent with 
    // the skybounds name in other catalogs, we keep that prefix here.

    std::vector<std::complex<double> > shear;
    std::vector<double> nu;
    std::vector<tmv::SmallMatrix<double,2,2> > cov;
    std::vector<BVec> shape;

  private :

    // We copy the parameter info, since we mess around with the parameters
    // to set all the file names of the component images.
    ConfigFile params;

    // For each object, we have a vector with an element for each image it
    // is found in, and each of those elements is a vector of Pixel
    std::vector<std::vector<PixelList POOL1> POOL2> pixlist;
    // Which image index corresponds to each pixel list?
    std::vector<std::vector<int POOL3> POOL4> image_indexlist;
    // The skypos projected back onto the source image
    std::vector<std::vector<Position POOL5> POOL6> image_cenlist;

    std::vector<std::string> image_file_list;
    std::vector<Bounds> image_skybounds;
    std::vector<std::string> shear_file_list;
    std::vector<std::string> fitpsf_file_list;
};

void MeasureMultiShear1(
    const Position& cen, 
    const std::vector<PixelList POOL1> allpix,
    const std::vector<BVec>& psf,
    double gal_aperture, double max_aperture,
    int gal_order, int gal_order2,
    double f_psf, double min_gal_size, bool desqa,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    tmv::SmallMatrix<double,2,2>& shearcov, BVec& shapelet,
    long& flag);

void MeasureMultiShear(
    const Position& cen, 
    const std::vector<PixelList POOL1>& pix,
    const std::vector<int POOL3>& image_index,
    const std::vector<Position POOL5>& image_cen,
    const std::vector<const FittedPSF*>& fitpsf,
    double gal_aperture, double max_aperture,
    int gal_order, int gal_order2,
    double f_psf, double min_gal_size, bool desqa,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    tmv::SmallMatrix<double,2,2>& shearcov, BVec& shapelet,
    long& flag);

#endif
