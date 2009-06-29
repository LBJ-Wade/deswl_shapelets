#ifndef MultiShearCatalog_H
#define MultiShearCatalog_H

#include <vector>
#include <string>
#include "TMV.h"
#include "Bounds.h"
#include "ConfigFile.h"
#include "Image.h"
#include "Pixel.h"
#include "Transformation.h"
#include "FittedPSF.h"
#include "dbg.h"
#include <CCfits/CCfits>
#include "Name.h"
#include "WlVersion.h"

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

    // Get pixel lists for the component images/catalogs
    void GetImagePixelLists();

    // Measure the shears
    int MeasureMultiShears(ShearLog& log);

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

    std::vector<std::complex<double> > shear;
    std::vector<double> nu;
    std::vector<tmv::SmallMatrix<double,2,2> > cov;
    std::vector<BVec> shape;

  private :

    // We copy the parameter info, since we mess around with the parameters
    // to set all the file names of the component images.
    ConfigFile params;

    // For each object, we have a vector with an element for each image it
    // is found in, and each of those elements is a vector if Pixel
    std::vector<std::vector<std::vector<Pixel> > > pixlist;
    // Which image index corresponds to each pixel list?
    std::vector<std::vector<int> > image_indexlist;
    // The skypos projected back onto the source image
    std::vector<std::vector<Position> > image_cenlist;

    std::vector<string> image_file_list;
    std::vector<string> fitpsf_file_list;


};

#endif
