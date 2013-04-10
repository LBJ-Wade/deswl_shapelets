#ifndef MultiShearCatalog_H
#define MultiShearCatalog_H

#include <vector>
#include <string>
#include "MyMatrix.h"

#include "dbg.h"
#include "CoaddCatalog.h"
#include "Bounds.h"
#include "ConfigFile.h"
#include "Image.h"
#include "Pixel.h"
#include "Transformation.h"
#include "FittedPsf.h"
#include "MEDSFile.h"

class MultiShearCatalog 
{

public :

    MultiShearCatalog(const CoaddCatalog& coaddcat, const ConfigFile& params);
    MultiShearCatalog(const ConfigFile& params);
    ~MultiShearCatalog();

    int size() const { return _skypos.size(); }
    int getNImages() const { return _image_file_list.size(); }
    int getNGalsWithPixels() const;

    // Read the srclist file
    void readFileLists();

    // Or can add images by hand.
    // The last two file names are optional.
    void addImage( 
        const std::string& image_filename,
        const std::string& fitpsf_filename,
        const std::string& shear_filename="",
        const std::string& skymap_filename="");

    // Get a set of bounds with a maximum linear extent in either direction
    // of params["multishear_section_size"] arcminutes on a side.
    std::vector<Bounds> splitBounds();

    // Get pixel lists for the component images/catalogs
    bool getPixels(const Bounds& b);
    bool getImagePixelLists(int fnum, const Bounds& b);

    // Measure the shears
    int measureMultiShears(const Bounds& b, ShearLog& log);
    int measureMEDS(const MEDSFile& meds, ShearLog& log);

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
    double calculateMemoryFootprint() const;

    const std::vector<long> getIdList() const { return _id; }
    const std::vector<Position> getSkyPosList() const { return _skypos; }
    const std::vector<long> getFlagsList() const { return _flags; }
    const std::vector<std::complex<double> > getShearList() const 
    { return _shear; }
    const std::vector<double> getNuList() const { return _nu; }
    const std::vector<DSmallMatrix22 >& getCovList() const 
    { return _cov; }
    const std::vector<int> getMeasGalOrderList() const 
    { return _meas_galorder; }
    const std::vector<BVec>& getShapeList() const { return _shape; }

    long getId(int i) const { return _id[i]; }
    Position getSkyPos(int i) const { return _skypos[i]; }
    long getFlags(int i) const { return _flags[i]; }
    std::complex<double> getShear(int i) const { return _shear[i]; }
    double getNu(int i) const { return _nu[i]; }
    const DSmallMatrix22& getCov(int i) const { return _cov[i]; }
    int getMeasGalOrder(int i) const { return _meas_galorder[i]; }
    const BVec& getShape(int i) const { return _shape[i]; }

    const Bounds& getSkyBounds() const { return _skybounds; }

private :

    // flags related to i/o and psf interpolation
    std::vector<long> _input_flags;

    // number of images each object was found in
    std::vector<int> _nimages_found;
    // number of images for which pixels were extracted
    std::vector<int> _nimages_gotpix;

    // These have an element for each coadd object.
    std::vector<long> _id;
    std::vector<Position> _chippos;
    std::vector<Position> _skypos;
    std::vector<long> _flags;
    std::vector<std::complex<double> > _shear;
    std::vector<double> _nu;
    std::vector<DSmallMatrix22> _cov;
    std::vector<int> _meas_galorder;
    EIGEN_mutable std::vector<BVec> _shape;

    // There is no non-sky bounds of course, but to be consistent with 
    // the skybounds name in other catalogs, we keep that prefix here.
    Bounds _skybounds; 

    // We copy the parameter info, since we mess around with the parameters
    // to set all the file names of the component images.
    // Hence, not a const reference the way I usually do it.
    ConfigFile _params;

    // For each coadd object, we have a vector with an element for each 
    // single-epoch (se) image it is found in, and each of those elements 
    // is PixelList (basically a vector of Pixels).
    std::vector<std::vector<PixelList> > _pix_list;
    // The PSF for each single-epoch object.
    std::vector<std::vector<BVec> > _psf_list;
    // The index number of each single-epoch image
    std::vector<std::vector<int> > _se_num;
    // The position in each single-epoch image
    std::vector<std::vector<Position> > _se_pos;

    // These each have an element for each single-epoch image
    std::vector<std::string> _image_file_list;
    std::vector<std::string> _shear_file_list;
    std::vector<std::string> _fitpsf_file_list;
    std::vector<std::string> _skymap_file_list;
    std::vector<Bounds> _saved_se_skybounds;
};

#endif
