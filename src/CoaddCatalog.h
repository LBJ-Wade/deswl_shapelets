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

    size_t size() const { return pos.size(); }

    void ReadCatalog();
    void Resize(int n);

    void ReadPixelLists();
    void ReadFileLists();

    void GetImagePixelLists();

    void MeasureMultiShears();

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

    std::vector<string> image_file_list;
    std::vector<string> fitpsf_file_list;


};

#endif
