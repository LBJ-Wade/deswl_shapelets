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

class MultiShearCatalog 
{

  public :

    MultiShearCatalog(const CoaddCatalog& coaddcat, const ConfigFile& params);
    MultiShearCatalog(const ConfigFile& params);
    ~MultiShearCatalog();

    size_t size() const { return skypos.size(); }
    size_t NImages() const { return image_file_list.size(); }

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

    // Calculate the current memory footprint of the entire structure in MB.
    // (Optionally output some info to os.)
    double CalcMemoryFootprint(bool getmax=false) const;

    // flags related to i/o and psf interpolation
    std::vector<long> input_flags;

    // number of images each object was found in
    std::vector<int> nimages_found;
    // number of images for which pixels were extracted
    std::vector<int> nimages_gotpix;

    std::vector<long> id;
    std::vector<Position> skypos;
    //std::vector<double> sky;
    //std::vector<double> noise;
    std::vector<long> flags;

    Bounds skybounds; 
    // There is no non-sky bounds of course, but to be consistent with 
    // the skybounds name in other catalogs, we keep that prefix here.

    // These have an element for each coadd object.
    std::vector<std::complex<double> > shear;
    std::vector<double> nu;
    std::vector<tmv::SmallMatrix<double,2,2> > cov;
    std::vector<BVec> shape;

  private :

    // We copy the parameter info, since we mess around with the parameters
    // to set all the file names of the component images.
    ConfigFile params;

    // For each coadd object, we have a vector with an element for each 
    // single-epoch (se) image it is found in, and each of those elements 
    // is PixelList (basically a vector of Pixels).
    std::vector<std::vector<PixelList> > pixlist;
    // The PSF for each single-epoch object.
    std::vector<std::vector<BVec> > psflist;
    // The single-epoch shear if known.
    std::vector<std::vector<std::complex<double> > > se_shearlist;
    // The single-epoch size if known.
    std::vector<std::vector<double> > se_sizelist;

    // These each have an element for each single-epoch image
    std::vector<std::string> image_file_list;
    std::vector<std::string> shear_file_list;
    std::vector<std::string> fitpsf_file_list;
    std::vector<Bounds> saved_se_skybounds;
};

void MeasureMultiShear(
    const Position& cen, 
    const std::vector<PixelList>& allpix,
    const std::vector<BVec>& psf,
    double gal_aperture, double max_aperture,
    int gal_order, int gal_order2,
    double f_psf, double min_gal_size,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    tmv::SmallMatrix<double,2,2>& shearcov, BVec& shapelet,
    long& flag);

#endif
