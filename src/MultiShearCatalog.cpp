
#include <valarray>
#include "TMV.h"
#include <CCfits/CCfits>

#include "CoaddCatalog.h"
#include "ConfigFile.h"
#include "dbg.h"
#include "Params.h"
#include "BVec.h"
#include "Pixel.h"
#include "Image.h"
#include "FittedPsf.h"
#include "Log.h"
#include "TimeVars.h"
#include "MultiShearCatalog.h"
#include "ShearCatalog.h"
#include "Form.h"
#include "WriteParam.h"
#include "WlVersion.h"
#include "ShearCatalogTree.h"

//#define ONLY_N_IMAGES 30

MultiShearCatalog::MultiShearCatalog(
    const CoaddCatalog& coaddCat, const ConfigFile& params) :
    _id(coaddCat.getIdList()), _skyPos(coaddCat.getSkyPosList()),
    _flags(coaddCat.getFlagsList()), _skyBounds(coaddCat.getSkyBounds()),
    _params(params)
{
    resize(coaddCat.size());

    // Read the names of the component image and catalog files from
    // the srclist file (given as params.coadd_srclist)
    readFileLists();
}

std::vector<Bounds> MultiShearCatalog::splitBounds(double side)
{
    double xRange = _skyBounds.getXMax() - _skyBounds.getXMin();
    double yRange = _skyBounds.getYMax() - _skyBounds.getYMin();

    int nx = int(floor(xRange/side))+1;
    int ny = int(floor(yRange/side))+1;
    return _skyBounds.divide(nx,ny);
}

int MultiShearCatalog::getPixels(const Bounds& bounds)
{
    // The pixlist object takes up a lot of memory, so at the start 
    // of this function, I clear it out, along with psflist, etc.
    // Then we loop over sections of the coaddCat and only load the information
    // for each single-epoch observation that actually falls within the 
    // bounds for the section we are working on.

    const int nPix = _pixList.size();

    for (int i=0;i<nPix;++i) {
        _pixList[i].clear();
        _psfList[i].clear();
        _seShearList[i].clear();
        _seSizeList[i].clear();
    }

    dbg<<"Start GetPixels for b = "<<bounds<<std::endl;
    // Loop over the files and read pixel lists for each object.
    // The Transformation and FittedPSF constructors get the name
    // information from the parameter file, so we use that to set the 
    // names of each component image here.
    const int nFiles = _imageFileList.size();
    for (int iFile=0; iFile<nFiles; ++iFile) {
#ifdef ONLY_N_IMAGES
        if (iFile >= ONLY_N_IMAGES) break;
#endif

        // Get the file names
        std::string imageFile = _imageFileList[iFile];
        std::string shearFile = _shearFileList[iFile];
        std::string fitPsfFile = _fitPsfFileList[iFile];

        dbg<<"Reading image file: "<<imageFile<<"\n";
        // Set the appropriate parameters
        _params["image_file"] = imageFile;
        setRoot(_params,imageFile);
        _params["shear_file"] = shearFile;
        _params["fitpsf_file"] = fitPsfFile;

        if (_skyMapFileList.size() > 0) {
            std::string skyMapFile = _skyMapFileList[iFile];
            _params["skymap_file"] = skyMapFile;
        }

        // Load the pixels
        getImagePixelLists(iFile,bounds);
        dbg<<"\n";
    }

    _params["maxmem"] = calculateMemoryFootprint(true);

    const int nGals = size();
    int nGalsWithPix=0;
    for (int i=0;i<nGals;++i) if (_pixList[i].size() > 0) ++nGalsWithPix;
    return nGalsWithPix;
}

MultiShearCatalog::MultiShearCatalog(const ConfigFile& params) :
    _params(params)
{
    read();
}

MultiShearCatalog::~MultiShearCatalog()
{
}

// readFileLists reads the srclist file specified in params
// and reads the names of the images and fitpsf 
void MultiShearCatalog::readFileLists()
{
    std::string file = _params.get("coadd_srclist");
    if (!doesFileExist(file)) {
        throw FileNotFoundException(file);
    }

    try {
        dbg<<"Opening coadd srclist\n";
        std::ifstream flist(file.c_str(), std::ios::in);
        if (!flist) {
            throw ReadException("Unable to open source list file " + file);
        }

        _imageFileList.clear();
        _shearFileList.clear();
        _fitPsfFileList.clear();
        _skyMapFileList.clear();

        std::string imageFilename;
        std::string shearFilename;
        std::string fitPsfFilename;
        std::string skyMapFilename;
        bool isSkyMapInList = _params.read("multishear_skymap_in_list",false);

        while (flist >> imageFilename >> shearFilename >> fitPsfFilename) {
            _imageFileList.push_back(imageFilename);
            _shearFileList.push_back(shearFilename);
            _fitPsfFileList.push_back(fitPsfFilename);
            dbg<<"Files are :\n"<<imageFilename<<std::endl;
            dbg<<shearFilename<<std::endl;
            dbg<<fitPsfFilename<<std::endl;
            if (isSkyMapInList) {
                flist >> skyMapFilename;
                if (!flist) {
                    throw ReadException(
                        "Unable to read skyMapFilename in list file " + file);
                }
                _skyMapFileList.push_back(skyMapFilename);
            }
        }
        if (isSkyMapInList) {
            Assert(_skyMapFileList.size() == _imageFileList.size());
        }
    } catch (std::exception& e) {
        throw ReadException(
            "Error reading from "+file+" -- caught error\n" + e.what());
    } catch (...) {
        throw ReadException(
            "Error reading from "+file+" -- caught unknown error");
    }

    dbg<<"Done reading file lists\n";
    Assert(_shearFileList.size() == _imageFileList.size());
    Assert(_fitPsfFileList.size() == _imageFileList.size());
}

void MultiShearCatalog::resize(int n)
{
    _pixList.clear();
    _pixList.resize(n);

    _psfList.clear();
    _psfList.resize(n);

    _seShearList.clear();
    _seShearList.resize(n);

    _seSizeList.clear();
    _seSizeList.resize(n);

    _inputFlags.resize(n,0);
    _nImagesFound.resize(n, 0);
    _nImagesGotPix.resize(n, 0);

    _shear.resize(n);
    _nu.resize(n);
    _cov.resize(n);

    int galOrder = _params.get("shear_gal_order");
    BVec shapeDefault(galOrder,1.);
    shapeDefault.vec().SetAllTo(DEFVALNEG);
    _shape.resize(n,shapeDefault);

}

template <typename T>
inline long long getMemoryFootprint(const T& x)
{
    return sizeof(x); 
}

template <typename T>
inline long long getMemoryFootprint(T*const x)
{
    long long res=sizeof(x);
    res += getMemoryFootprint(*x);
    return res;
}

template <typename T>
inline long long getMemoryFootprint(const T*const x)
{
    long long res=sizeof(x);
    res += getMemoryFootprint(*x);
    return res;
}

template <typename T>
inline long long getMemoryFootprint(const std::auto_ptr<T>& x)
{
    long long res=sizeof(x);
    res += getMemoryFootprint(*x);
    return res;
}

template <typename T, class Alloc>
inline long long getMemoryFootprint(const std::vector<T,Alloc>& x)
{
    long long res=sizeof(x);
    const int xsize = x.size();
    for(int i=0;i<xsize;++i) res += getMemoryFootprint(x[i]);
    return res;
}

inline long long getMemoryFootprint(const PixelList& x)
{
    long long res=sizeof(x);
    res += x.size() * sizeof(Pixel);
    return res;
}

template <typename T>
inline long long getMemoryFootprint(const tmv::Vector<T>& x)
{
    long long res=sizeof(x);
    res += x.size() * sizeof(T);
    return res;
}

inline long long getMemoryFootprint(const BVec& x)
{
    long long res=sizeof(x);
    res += x.size() * sizeof(double);
    return res;
}

template <typename T>
inline long long getMemoryFootprint(const tmv::Matrix<T>& x)
{
    long long res=sizeof(x);
    res += x.colsize() * x.rowsize() * sizeof(T);
    return res;
}

double MultiShearCatalog::calculateMemoryFootprint(bool shouldGetMax) const
{

    static double maxMem=0.;

    dbg<<"Memory usage:\n";
    dbg<<"input_flags: "<<
        getMemoryFootprint(_inputFlags)/1024./1024.<<" MB\n";
    dbg<<"nimages_found: "<<
        getMemoryFootprint(_nImagesFound)/1024./1024.<<" MB\n";
    dbg<<"pixlist: "<<
        getMemoryFootprint(_pixList)/1024./1024.<<" MB\n";
    dbg<<"nimages_gotpix: "<<
        getMemoryFootprint(_nImagesGotPix)/1024./1024.<<" MB\n";
    dbg<<"id: "<<
        getMemoryFootprint(_id)/1024./1024.<<" MB\n";
    dbg<<"skypos: "<<
        getMemoryFootprint(_skyPos)/1024./1024.<<" MB\n";
    dbg<<"flags: "<<
        getMemoryFootprint(_flags)/1024./1024.<<" MB\n";
    dbg<<"shear: "<<
        getMemoryFootprint(_shear)/1024./1024.<<" MB\n";
    dbg<<"nu: "<<
        getMemoryFootprint(_nu)/1024./1024.<<" MB\n";
    dbg<<"cov: "<<
        getMemoryFootprint(_cov)/1024./1024.<<" MB\n";
    dbg<<"shape: "<<
        getMemoryFootprint(_shape)/1024./1024.<<" MB\n";
    dbg<<"psflist: "<<
        getMemoryFootprint(_psfList)/1024./1024.<<" MB\n";
    dbg<<"se_shearlist: "<<
        getMemoryFootprint(_seShearList)/1024./1024.<<" MB\n";
    dbg<<"se_sizelist: "<<
        getMemoryFootprint(_seSizeList)/1024./1024.<<" MB\n";
    dbg<<"image_file_list: "<<
        getMemoryFootprint(_imageFileList)/1024./1024.<<" MB\n";
    dbg<<"shear_file_list: "<<
        getMemoryFootprint(_shearFileList)/1024./1024.<<" MB\n";
    dbg<<"fitpsf_file_list: "<<
        getMemoryFootprint(_fitPsfFileList)/1024./1024.<<" MB\n";
    dbg<<"skymap_file_list: "<<
        getMemoryFootprint(_skyMapFileList)/1024./1024.<<" MB\n";
    dbg<<"saved_se_skybounds: "<<
        getMemoryFootprint(_savedSeSkyBounds)/1024./1024.<<" MB\n";

    double totMem = 
        getMemoryFootprint(_inputFlags) +
        getMemoryFootprint(_nImagesFound) +
        getMemoryFootprint(_pixList) +
        getMemoryFootprint(_nImagesGotPix) +
        getMemoryFootprint(_id) +
        getMemoryFootprint(_skyPos) +
        getMemoryFootprint(_flags) +
        getMemoryFootprint(_shear) +
        getMemoryFootprint(_nu) +
        getMemoryFootprint(_cov) +
        getMemoryFootprint(_shape) +
        getMemoryFootprint(_psfList) +
        getMemoryFootprint(_seShearList) +
        getMemoryFootprint(_seSizeList) +
        getMemoryFootprint(_imageFileList) +
        getMemoryFootprint(_shearFileList) +
        getMemoryFootprint(_fitPsfFileList) +
        getMemoryFootprint(_skyMapFileList) +
        getMemoryFootprint(_savedSeSkyBounds);
    totMem /= 1024.*1024.;  // B -> MB
    dbg<<"totmem = "<<totMem<<std::endl;

    if (totMem > maxMem) maxMem = totMem;

    return shouldGetMax ? maxMem : totMem;
}

void MultiShearCatalog::write() const
{
    std::vector<std::string> fileList = makeMultiName(_params, "multishear");

    const int nFiles = fileList.size();
    for(int i=0; i<nFiles; ++i) {
        const std::string& file = fileList[i];
        dbg<<"Writing multishear catalog to file: "<<file<<std::endl;

        bool isFitsIo = false;
        if (_params.keyExists("multishear_io")) {
            std::vector<std::string> ioList = _params["multishear_io"];
            Assert(ioList.size() == fileList.size());
            isFitsIo = (ioList[i] == "FITS");
        } else if (file.find("fits") != std::string::npos) {
            isFitsIo = true;
        }

        try {
            if (isFitsIo) {
                writeFits(file);
            } else {
                std::string delim = "  ";
                if (_params.keyExists("multishear_delim")) {
                    std::vector<std::string> delimList = 
                        _params["multishear_delim"];
                    Assert(delimList.size() == fileList.size());
                    delim = delimList[i];
                } else if (file.find("csv") != std::string::npos) {
                    delim = ",";
                }
                writeAscii(file,delim);
            }
        } catch (CCfits::FitsException& e) {
            throw WriteException(
                "Error writing to "+file+" -- caught error\n" + e.message());
        } catch (std::exception& e) {
            throw WriteException(
                "Error writing to "+file+" -- caught error\n" + e.what());
        } catch (...) {
            throw WriteException(
                "Error writing to "+file+" -- caught unknown error");
        }
    }
    dbg<<"Done Write ShearCatalog\n";
}

void MultiShearCatalog::writeFits(std::string file) const
{
    // ! means overwrite existing file
    CCfits::FITS fits("!"+file, CCfits::Write);

    std::vector<string> colNames;
    std::vector<string> colFmts;

    colNames.push_back(_params.get("multishear_id_col")); 
    colFmts.push_back("1J");
    int idColNum = colNames.size()-1;

    colNames.push_back(_params.get("multishear_flags_col"));
    colFmts.push_back("1J"); // flags
    int flagsColNum = colNames.size()-1;

    colNames.push_back(_params.get("multishear_ra_col"));
    colFmts.push_back("1D"); // ra
    int raColNum = colNames.size()-1;

    colNames.push_back(_params.get("multishear_dec_col"));
    colFmts.push_back("1D"); // decl
    int declColNum = colNames.size()-1;

    colNames.push_back(_params.get("multishear_shear1_col"));
    colFmts.push_back("1D"); // shear1
    int shear1ColNum = colNames.size()-1;

    colNames.push_back(_params.get("multishear_shear2_col"));
    colFmts.push_back("1D"); // shear2
    int shear2ColNum = colNames.size()-1;

    colNames.push_back(_params.get("multishear_nu_col"));
    colFmts.push_back("1D"); // nu
    int nuColNum = colNames.size()-1;

    colNames.push_back(_params.get("multishear_cov00_col"));
    colFmts.push_back("1D"); // cov00
    int cov00ColNum = colNames.size()-1;

    colNames.push_back(_params.get("multishear_cov01_col"));
    colFmts.push_back("1D"); // cov01
    int cov01ColNum = colNames.size()-1;

    colNames.push_back(_params.get("multishear_cov11_col"));
    colFmts.push_back("1D"); // cov11
    int cov11ColNum = colNames.size()-1;

    colNames.push_back(_params.get("multishear_order_col"));
    colFmts.push_back("1J"); // order
    int orderColNum = colNames.size()-1;

    colNames.push_back(_params.get("multishear_sigma_col"));
    colFmts.push_back("1D"); // sigma
    int sigmaColNum = colNames.size()-1;

    colNames.push_back(_params.get("multishear_coeffs_col"));

    int nCoeff = _shape[0].size();
    dbg<<"ncoeff = "<<nCoeff<<std::endl;
    std::stringstream coeffForm;
    coeffForm << nCoeff << "D";
    colFmts.push_back(coeffForm.str()); // shapelet coeffs

    int coeffsColNum = colNames.size()-1;

    colNames.push_back(_params.get("multishear_nimages_found_col"));
    colFmts.push_back("1J"); // nImages_found
    int nImagesFoundColNum = colNames.size()-1;

    colNames.push_back(_params.get("multishear_nimages_gotpix_col"));
    colFmts.push_back("1J"); // nImages_gotpix
    int nImagesGotPixColNum = colNames.size()-1;

    colNames.push_back(_params.get("multishear_input_flags_col"));
    colFmts.push_back("1J"); // inputFlags
    int inputFlagsColNum = colNames.size()-1;

    // We don't output units yet, just make it same size as others
    std::vector<string> colunits(colFmts.size());

    dbg<<"Before Create table"<<std::endl;
    CCfits::Table* table;
    table = fits.addTable("coadd_shearcat",size(),colNames,colFmts,colunits);

    // Header keywords
    std::string tmvvers = tmv::TMV_Version();
    std::string wlvers = getWlVersion();

    table->addKey("tmvvers", tmvvers, "version of TMV code");
    table->addKey("wlvers", wlvers, "version of weak lensing code");

    std::string str;
    double dbl;
    int intgr;
    writeParamToTable(_params, table, "noise_method", str);
    writeParamToTable(_params, table, "dist_method", str);

    writeParamToTable(_params, table, "shear_aperture", dbl);
    writeParamToTable(_params, table, "shear_max_aperture", dbl);
    writeParamToTable(_params, table, "shear_gal_order", intgr);
    writeParamToTable(_params, table, "shear_gal_order2", intgr);
    writeParamToTable(_params, table, "shear_min_gal_size", dbl);
    writeParamToTable(_params, table, "shear_f_psf", dbl);

    writeParamToTable(_params, table, "maxmem", dbl);

    // if merun= is sent we'll put it in the header.  This allows us to 
    // associate some more, possibly complicated, metadata with this file
    if (_params.keyExists("merun")) {
        writeParamToTable(_params, table, "merun", str);
    }

    const int nGals = size();
    // data
    // make vector copies for writing
    std::vector<double> ra(nGals);
    std::vector<double> decl(nGals);
    std::vector<double> shear1(nGals);
    std::vector<double> shear2(nGals);
    std::vector<double> cov00(nGals);
    std::vector<double> cov01(nGals);
    std::vector<double> cov11(nGals);

    for(int i=0;i<nGals;++i) {
        // internally we use arcseconds
        ra[i] = _skyPos[i].getX()/3600.;
        decl[i] = _skyPos[i].getY()/3600.;
        shear1[i] = real(_shear[i]);
        shear2[i] = imag(_shear[i]);
        cov00[i] = _cov[i](0,0);
        cov01[i] = _cov[i](0,1);
        cov11[i] = _cov[i](1,1);
    }

    int startRow=1;

    table->column(colNames[idColNum]).write(_id,startRow);
    table->column(colNames[flagsColNum]).write(_flags,startRow);
    table->column(colNames[raColNum]).write(ra,startRow);
    table->column(colNames[declColNum]).write(decl,startRow);
    table->column(colNames[shear1ColNum]).write(shear1,startRow);
    table->column(colNames[shear2ColNum]).write(shear2,startRow);
    table->column(colNames[nuColNum]).write(_nu,startRow);
    table->column(colNames[cov00ColNum]).write(cov00,startRow);
    table->column(colNames[cov01ColNum]).write(cov01,startRow);
    table->column(colNames[cov11ColNum]).write(cov11,startRow);

    for (int i=0; i<nGals; ++i) {
        int row = i+1;
        long bOrder = _shape[i].getOrder();
        double bSigma = _shape[i].getSigma();

        table->column(colNames[orderColNum]).write(&bOrder,1,row);
        table->column(colNames[sigmaColNum]).write(&bSigma,1,row);
        double* cptr = (double *) _shape[i].vec().cptr();
        table->column(colNames[coeffsColNum]).write(cptr, nCoeff, 1, row);

    }

    table->column(colNames[nImagesFoundColNum]).write(_nImagesFound,startRow);
    table->column(colNames[nImagesGotPixColNum]).write(_nImagesGotPix,startRow);
    table->column(colNames[inputFlagsColNum]).write(_inputFlags,startRow);
}

void MultiShearCatalog::writeAscii(std::string file, std::string delim) const
{
    Assert(_id.size() == size());
    Assert(_skyPos.size() == size());
    Assert(_flags.size() == size());
    Assert(_shear.size() == size());
    Assert(_nu.size() == size());
    Assert(_cov.size() == size());
    Assert(_shape.size() == size());
    Assert(_nImagesFound.size() == size());
    Assert(_nImagesGotPix.size() == size());
    Assert(_inputFlags.size() == size());

    std::ofstream fout(file.c_str());
    if (!fout) {
        throw WriteException("Error opening shear file"+file);
    }

    Form sci8; sci8.sci().trail(0).prec(8);
    Form fix3; fix3.fix().trail(0).prec(3);
    Form fix8; fix8.fix().trail(0).prec(8);
    Form fix6; fix6.fix().trail(0).prec(6);

    const int nGals = size();
    for(int i=0;i<nGals;++i) {
        fout
            << _id[i] << delim
            << fix8(_skyPos[i].getX()/3600.) << delim
            << fix8(_skyPos[i].getY()/3600.) << delim
            << _flags[i] << delim
            << sci8(real(_shear[i])) << delim
            << sci8(imag(_shear[i])) << delim
            << fix3(_nu[i]) << delim
            << sci8(_cov[i](0,0)) << delim
            << sci8(_cov[i](0,1)) << delim
            << sci8(_cov[i](1,1)) << delim
            << _nImagesFound[i] << delim
            << _nImagesGotPix[i] << delim
            << _inputFlags[i] << delim
            << _shape[i].getOrder() << delim
            << fix6(_shape[i].getSigma()) << delim;
        const int nCoeffs = _shape[i].size();
        for(int j=0;j<nCoeffs;++j)
            fout << delim << sci8(_shape[i](j));
        fout << std::endl;
    }
}

void MultiShearCatalog::read()
{
    std::string file = makeName(_params,"multishear",false,true);
    // false,true = input_prefix=false, mustexist=true.
    // It is an input here, but it is in the output_prefix directory.
    dbg<< "Reading Shear cat from file: " << file << std::endl;

    bool isFitsIo = false;
    if (_params.keyExists("multishear_io"))
        isFitsIo = (_params["multishear_io"] == "FITS");
    else if (file.find("fits") != std::string::npos)
        isFitsIo = true;

    if (!doesFileExist(file)) {
        throw FileNotFoundException(file);
    }
    try {
        if (isFitsIo) {
            readFits(file);
        } else {
            std::string delim = "  ";
            if (_params.keyExists("multishear_delim")) delim = 
                _params["multishear_delim"];
            else if (file.find("csv") != std::string::npos) delim = ",";
            readAscii(file,delim);
        }
    } catch (CCfits::FitsException& e) {
        throw ReadException(
            "Error reading from "+file+" -- caught error\n" + e.message());
    } catch (std::exception& e) {
        throw ReadException(
            "Error reading from "+file+" -- caught error\n" + e.what());
    } catch (...) {
        throw ReadException(
            "Error reading from "+file+" -- caught unknown error");
    }
    dbg<<"Done Read ShearCatalog\n";
}

void MultiShearCatalog::readFits(std::string file)
{
    int hdu = _params.read("multishear_hdu",2);

    dbg<<"Opening FITS file "<<file<<" at hdu "<<hdu<<std::endl;
    CCfits::FITS fits(file, CCfits::Read);
    if (hdu > 1) fits.read(hdu-1);

    CCfits::ExtHDU& table=fits.extension(hdu-1);

    long nRows=table.rows();

    dbg<<"  nrows = "<<nRows<<std::endl;
    if (nRows <= 0) {
        throw ReadException(
            "ShearCatalog found to have 0 rows.  Must have > 0 rows.");
    }

    std::string idCol=_params.get("multishear_id_col");
    std::string raCol=_params.get("multishear_ra_col");
    std::string declCol=_params.get("multishear_dec_col");
    std::string skyCol=_params.get("multishear_sky_col");
    std::string noiseCol=_params.get("multishear_noise_col");
    std::string flagsCol=_params.get("multishear_flags_col");
    std::string shear1Col=_params.get("multishear_shear1_col");
    std::string shear2Col=_params.get("multishear_shear2_col");
    std::string nuCol=_params.get("multishear_nu_col");
    std::string cov00Col=_params.get("multishear_cov00_col");
    std::string cov01Col=_params.get("multishear_cov01_col");
    std::string cov11Col=_params.get("multishear_cov11_col");
    std::string orderCol=_params.get("multishear_order_col");
    std::string sigmaCol=_params.get("multishear_sigma_col");
    std::string coeffsCol=_params.get("multishear_coeffs_col");
    std::string nImagesFoundCol=_params.get("multishear_nimages_found_col");
    std::string nImagesGotPixCol=_params.get("multishear_nimages_gotpix_col");
    std::string inputFlagsCol=_params.get("multishear_input_flags_col");

    long start=1;
    long end=nRows;

    dbg<<"Reading columns"<<std::endl;
    dbg<<"  "<<idCol<<std::endl;
    table.column(idCol).read(_id, start, end);

    dbg<<"  "<<raCol<<"  "<<declCol<<std::endl;
    _skyPos.resize(nRows);
    std::vector<double> ra;
    std::vector<double> decl;
    table.column(raCol).read(ra, start, end);
    table.column(declCol).read(decl, start, end);
    for(long i=0;i<nRows;++i) _skyPos[i] = Position(ra[i],decl[i]);

    dbg<<"  "<<flagsCol<<std::endl;
    table.column(flagsCol).read(_flags, start, end);

    dbg<<"  "<<shear1Col<<"  "<<shear2Col<<std::endl;
    _shear.resize(nRows);
    std::vector<double> shear1;
    std::vector<double> shear2;
    table.column(shear1Col).read(shear1, start, end);
    table.column(shear2Col).read(shear2, start, end);
    for(long i=0;i<nRows;++i) {
        _shear[i] = std::complex<double>(shear1[i],shear2[i]);
    }

    dbg<<"  "<<nuCol<<std::endl;
    table.column(nuCol).read(_nu, start, end);

    dbg<<"  "<<cov00Col<<"  "<<cov01Col<<"  "<<cov11Col<<std::endl;
    _cov.resize(nRows);
    std::vector<double> cov00;
    std::vector<double> cov01;
    std::vector<double> cov11;
    table.column(cov00Col).read(cov00, start, end);
    table.column(cov01Col).read(cov01, start, end);
    table.column(cov11Col).read(cov11, start, end);
    for(long i=0;i<nRows;++i) {
        _cov[i] = tmv::ListInit, cov00[i], cov01[i], cov01[i], cov11[i];
    }

    dbg<<"  "<<sigmaCol<<"  "<<orderCol<<std::endl;
    // temporary
    std::vector<double> sigma;
    std::vector<int> order;
    table.column(sigmaCol).read(sigma, start, end);
    table.column(orderCol).read(order, start, end);

    _shape.reserve(nRows);
    for (int i=0; i<nRows; ++i) {
        int row=i+1;

        _shape.push_back(BVec(order[i],sigma[i]));
        int nCoeff=(order[i]+1)*(order[i]+2)/2;

        std::valarray<double> coeffs;
        table.column(coeffsCol).read(coeffs, row);

        double* ptri = (double* ) _shape[i].vec().cptr();
        for (int j=0; j<nCoeff; ++j) {
            ptri[i] = coeffs[i];
        }
    }

    dbg<<"  "<<nImagesFoundCol<<std::endl;
    table.column(nImagesFoundCol).read(_nImagesFound, start, end);

    dbg<<"  "<<nImagesGotPixCol<<std::endl;
    table.column(nImagesGotPixCol).read(_nImagesGotPix, start, end);

    dbg<<"  "<<inputFlagsCol<<std::endl;
    table.column(inputFlagsCol).read(_inputFlags, start, end);
}

void MultiShearCatalog::readAscii(std::string file, std::string delim)
{
    std::ifstream fin(file.c_str());
    if (!fin) {
        throw ReadException("Error opening stars file"+file);
    }

    _id.clear(); _skyPos.clear(); _flags.clear();
    _shear.clear(); _nu.clear(); _cov.clear(); _shape.clear();
    _nImagesFound.clear(); _nImagesGotPix.clear(); _inputFlags.clear();

    if (delim == "  ") {
        ConvertibleString flag;
        long id1,bOrder,nfound,ngotpix,inflag;
        double ra,decl,s1,s2,nu1,c00,c01,c11,bSigma;
        while ( fin >> id1 >> ra >> decl >> 
                flag >> s1 >> s2 >> nu1 >>
                c00 >> c01 >> c11 >> 
                nfound >> ngotpix >> inflag >>
                bOrder >> bSigma) {
            _id.push_back(id1);
            _skyPos.push_back(Position(ra,decl));
            _flags.push_back(flag);
            _shear.push_back(std::complex<double>(s1,s2));
            _nu.push_back(nu1);
            _cov.push_back(tmv::SmallMatrix<double,2,2>());
            _cov.back() = tmv::ListInit, c00, c01, c01, c11;
            _nImagesFound.push_back(nfound);
            _nImagesGotPix.push_back(ngotpix);
            _inputFlags.push_back(inflag);
            _shape.push_back(BVec(bOrder,bSigma));
            const int nCoeffs = _shape.back().size();
            for(int j=0;j<nCoeffs;++j) 
                fin >> _shape.back()(j);
        }
    } else {
        if (delim.size() > 1) {
            // getline only works with single character delimiters.
            // Since I don't really expect a multicharacter delimiter to
            // be used ever, I'm just going to throw an exception here 
            // if we do need it, and I can write the workaround then.
            throw ParameterException(
                "ReadAscii delimiter must be a single character");
        }
        char d = delim[0];
        long bOrder;
        double ra,decl,s1,s2,bSigma,c00,c01,c11;
        ConvertibleString temp;
        while (getline(fin,temp,d)) {
            _id.push_back(temp);
            getline(fin,temp,d); ra = temp;
            getline(fin,temp,d); decl = temp;
            _skyPos.push_back(Position(ra,decl));
            getline(fin,temp,d); _flags.push_back(temp);
            getline(fin,temp,d); s1 = temp;
            getline(fin,temp,d); s2 = temp;
            _shear.push_back(std::complex<double>(s1,s2));
            getline(fin,temp,d); _nu.push_back(temp);
            getline(fin,temp,d); c00 = temp;
            getline(fin,temp,d); c01 = temp;
            getline(fin,temp,d); c11 = temp;
            _cov.push_back(tmv::SmallMatrix<double,2,2>());
            _cov.back() = tmv::ListInit, c00, c01, c01, c11;
            getline(fin,temp,d); _nImagesFound.push_back(temp);
            getline(fin,temp,d); _nImagesGotPix.push_back(temp);
            getline(fin,temp,d); _inputFlags.push_back(temp);
            getline(fin,temp,d); bOrder = temp;
            getline(fin,temp,d); bSigma = temp;
            _shape.push_back(BVec(bOrder,bSigma));
            const int nCoeffs = _shape.back().size();
            for(int j=0;j<nCoeffs-1;++j) {
                getline(fin,temp,d); _shape.back()(j) = temp;
            }
            // Last one doesn't have a following delimiter
            getline(fin,temp); _shape.back()(nCoeffs-1) = temp;
        }
    }
}

