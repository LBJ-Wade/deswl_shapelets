#ifndef CoaddCatalog_H
#define CoaddCatalog_H

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

class CoaddCatalog 
{

  public :

    CoaddCatalog(ConfigFile& _params);
    ~CoaddCatalog();

    size_t size() const { return pos.size(); }
    size_t NImages() const { return fitpsf.size(); }

    void ReadCatalog();
    void Resize(int n);

    void ReadPixelLists();
    void ReadFileLists();

    void GetImagePixelLists();

    int MeasureMultiShears(ShearLog& log);

    void WriteFits() const;

    // Leave these public, rather than use Get and Set methods.
    std::vector<long> id;
    std::vector<Position> pos;
    std::vector<Position> skypos;

    std::vector<double> sky;
    std::vector<double> noise;

    std::vector<long> flags;

    std::vector<float> ra;
    std::vector<float> dec;

    std::vector<float> mag;
    std::vector<float> mag_err;

    // These next two have length NImages(), rather than size()
    // (aka the number of objects)
    std::vector<const Transformation*> trans;
    std::vector<const FittedPSF*> fitpsf;

    //
    // shear data
    // 

    // flags related to i/o and psf interpolation
    std::vector<long> input_flags;

    // number of images each object was found in
    std::vector<int> nimages_found;
    // number of images for which pixels were extracted
    std::vector<int> nimages_gotpix;

    std::vector<std::complex<double> > shear;
    std::vector<double> nu;
    std::vector<tmv::SmallMatrix<double,2,2> > cov;
    std::vector<BVec> shape;


  private :

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
