
#include <valarray>
#include <CCfits/CCfits>
#include <sstream>

#include "dbg.h"
#include "CoaddCatalog.h"
#include "ConfigFile.h"
#include "Params.h"
#include "BVec.h"
#include "Pixel.h"
#include "Image.h"
#include "FittedPsf.h"
#include "Log.h"
#include "MultiShearCatalog.h"
#include "ShearCatalog.h"
#include "Form.h"
#include "WriteParam.h"
#include "WlVersion.h"
#include "ShearCatalogTree.h"

//#define ONLY_N_IMAGES 30

MultiShearCatalog::MultiShearCatalog(
    const CoaddCatalog& coaddCat, const ConfigFile& params) :
    _id(coaddCat.getIdList()), _chippos(coaddCat.getPosList()),
    _skypos(coaddCat.getSkyPosList()),
    _flags(coaddCat.getFlagsList()), _skybounds(coaddCat.getSkyBounds()),
    _params(params)
{
    dbg<<"Start MultiShearCatalog constructor\n";
    xdbg<<"memory_usage = "<<memory_usage()<<std::endl;

    std::complex<double> shear_default(DEFVALPOS,DEFVALPOS);
    DSmallMatrix22 cov_default;
    cov_default << DEFVALPOS, 0, 0, DEFVALPOS;

    const int n = coaddCat.size();
    _pix_list.resize(n);
    _psf_list.resize(n);
    _se_num.resize(n);
    _se_pos.resize(n);
    _input_flags.resize(n,0);
    _nimages_found.resize(n, 0);
    _nimages_gotpix.resize(n, 0);

    _meas_galorder.resize(n,DEFVALNEG);
    _shear.resize(n,shear_default);
    _nu.resize(n,DEFVALNEG);
    _cov.resize(n,cov_default);

    int galorder = _params.get("shear_gal_order");
    BVec shape_default(galorder,1.);
    shape_default.vec().TMV_setAllTo(DEFVALNEG);
    _shape.resize(n,shape_default);

    xdbg<<"after resize, memory_usage = "<<memory_usage()<<std::endl;

    // Read the names of the component image and catalog files from
    // the srclist file (given as params.coadd_srclist)
    readFileLists();
    xdbg<<"after readfilelists, memory_usage = "<<memory_usage()<<std::endl;
}

int MultiShearCatalog::getNGalsWithPixels() const
{
    const int ngals = size();
    int ngals_withpix=0;
    for (int i=0;i<ngals;++i) if (_pix_list[i].size() > 0) ++ngals_withpix;
    return ngals_withpix;
}

std::vector<Bounds> MultiShearCatalog::splitBounds()
{
    const double ARCSEC_PER_RAD = 206264.806247;

    double sectionSize = _params.get("multishear_section_size");
    xdbg<<"sectionSize = "<<sectionSize<<std::endl;
    sectionSize *= 60.; // arcmin -> arcsec

    double xRange = _skybounds.getXMax() - _skybounds.getXMin();
    double yRange = _skybounds.getYMax() - _skybounds.getYMin();
    double dec = _skybounds.getCenter().getY();
    double cosdec = cos(dec / ARCSEC_PER_RAD);
    xRange *= cosdec;

    int nx = int(floor(xRange/sectionSize))+1;
    int ny = int(floor(yRange/sectionSize))+1;
    return _skybounds.divide(nx,ny);
}

bool MultiShearCatalog::getPixels(const Bounds& bounds)
{
    // The pixlist object takes up a lot of memory, so at the start 
    // of this function, I clear it out, along with psflist, etc.
    // Then we loop over sections of the coaddCat and only load the information
    // for each single-epoch observation that actually falls within the 
    // bounds for the section we are working on.

    const int nPix = _pix_list.size();

    dbg<<"Start getPixels: memory_usage = "<<memory_usage()<<std::endl;
    for (int i=0;i<nPix;++i) {
        _pix_list[i].clear();
        _psf_list[i].clear();
    }
    dbg<<"After clear: memory_usage = "<<memory_usage()<<std::endl;

    bool des_qa = _params.read("des_qa",false); 

    try {
        dbg<<"Start GetPixels for b = "<<bounds<<std::endl;
        memory_usage(dbgout);
        // Loop over the files and read pixel lists for each object.
        // The Transformation and FittedPSF constructors get the name
        // information from the parameter file, so we use that to set the 
        // names of each component image here.
        const int nfiles = _image_file_list.size();
        for (int ifile=0; ifile<nfiles; ++ifile) {
#ifdef ONLY_N_IMAGES
            if (ifile >= ONLY_N_IMAGES) break;
#endif
            dbg<<"ifile = "<<ifile<<std::endl;

            // Get the file names
            Assert(ifile < int(_image_file_list.size()));
            Assert(ifile < int(_fitpsf_file_list.size()));
            std::string image_file = _image_file_list[ifile];
            std::string fitpsf_file = _fitpsf_file_list[ifile];

            dbg<<"Reading image file: "<<image_file<<"\n";
            // Set the appropriate parameters
            _params["image_file"] = image_file;
            SetRoot(_params,image_file);
            _params["fitpsf_file"] = fitpsf_file;

            if (_shear_file_list.size() > 0) {
                Assert(ifile < int(_shear_file_list.size()));
                std::string shear_file = _shear_file_list[ifile];
                _params["shear_file"] = shear_file;
            }

            if (_skymap_file_list.size() > 0) {
                Assert(ifile < int(_skymap_file_list.size()));
                std::string skymap_file = _skymap_file_list[ifile];
                _params["skymap_file"] = skymap_file;
            }

            // Load the pixels
            dbg<<"Before load pixels for file "<<ifile<<
                ": memory_usage = "<<memory_usage()<<std::endl;
            if (!getImagePixelLists(ifile,bounds)) {
                for (int i=0;i<nPix;++i) {
                    _pix_list[i].clear();
                    _psf_list[i].clear();
                }
                PixelList::reclaimMemory();
                return false;
            }
            dbg<<"After load pixels for file "<<ifile<<
                ": memory_usage = "<<memory_usage()<<std::endl;
        }
    } catch (std::bad_alloc) {
        dbg<<"Caught bad_alloc\n";
        double mem = memory_usage(dbgout);
        double peak_mem = peak_memory_usage();
        dbg<<"memory usage = "<<mem<<" MB\n";
        dbg<<"peak memory usage = "<<peak_mem<<" MB\n";
        if (des_qa) std::cerr<<"STATUS5BEG ";
        std::cerr
            << "Memory exhausted in MultShearCatalog.\n"
            << "Memory Usage in MultiShearCatalog = "
            << calculateMemoryFootprint()<<" MB \n"
            << "Actual Virtual Memory Usage = "
            << mem<<" MB \n"
            << "Try reducing multishear_section_size or "
            << "reducing mam_vmem\n"
            << "(Current values = "
            << _params["multishear_section_size"]
            << " , "<<_params["max_vmem"]<<")";
        if (des_qa) std::cerr<<" STATUS5END";
        std::cerr<<std::endl;
        dbg << "Memory exhausted in MultShearCatalog.\n"
            << "Memory Usage in MultiShearCatalog = "
            << calculateMemoryFootprint()<<" MB \n"
            << "Actual Virtual Memory Usage = "
            << mem<<" MB \n"
            << "Try reducing multishear_section_size or "
            << "reducing mam_vmem\n"
            << "(Current values = "
            << _params["multishear_section_size"]
            << " , "<<_params["max_vmem"]<<")"
            << std::endl;
        exit(1);
    }
    dbg <<"Done getPixels\n";
    dbg << "Memory Usage in MultiShearCatalog = "
        << calculateMemoryFootprint()<<" MB \n";
    double mem = memory_usage(dbgout);
    double peak_mem = peak_memory_usage();
    double max_mem = double(_params["max_vmem"])*1024.;
    _params["peak_mem"] = peak_mem; // Keep track...
    dbg<<"Actual memory usage = "<<mem<<" MB\n";
    dbg<<"Peak memory usage = "<<peak_mem<<" MB\n";
    dbg<<"Max allowed memory usage = "<<max_mem<<" MB\n";
    return true;
}

MultiShearCatalog::MultiShearCatalog(const ConfigFile& params) :
    _params(params)
{}

MultiShearCatalog::~MultiShearCatalog() 
{}

void MultiShearCatalog::addImage(
    const std::string& image_filename, const std::string& fitpsf_filename,
    const std::string& shear_filename, const std::string& skymap_filename)
{
    _image_file_list.push_back(image_filename);
    _fitpsf_file_list.push_back(fitpsf_filename);
    if (shear_filename != "") _shear_file_list.push_back(shear_filename);
    if (skymap_filename != "") _skymap_file_list.push_back(skymap_filename);
}

// readFileLists reads the srclist file specified in params
// and reads the names of the images and fitpsf 
void MultiShearCatalog::readFileLists()
{
    std::string file = _params.get("coadd_srclist");
    if (!DoesFileExist(file)) {
        throw FileNotFoundException(file);
    }

    try {
        dbg<<"Opening coadd srclist\n";
        std::ifstream flist(file.c_str(), std::ios::in);
        if (!flist) {
            throw ReadException("Unable to open source list file " + file);
        }

        _image_file_list.clear();
        _shear_file_list.clear();
        _fitpsf_file_list.clear();
        _skymap_file_list.clear();

        std::string image_filename;
        std::string shear_filename;
        std::string fitpsf_filename;
        std::string skymap_filename;
        bool isSkyMapIn_list = _params.read("multishear_skymap_in_list",false);

        while (flist >> image_filename >> shear_filename >> fitpsf_filename) {
            dbg<<"Files are :\n"<<image_filename<<std::endl;
            dbg<<shear_filename<<std::endl;
            dbg<<fitpsf_filename<<std::endl;
            if (isSkyMapIn_list) {
                flist >> skymap_filename;
                if (!flist) {
                    throw ReadException(
                        "Unable to read skymap_filename in list file " + file);
                }
                addImage(image_filename,fitpsf_filename,
                         shear_filename,skymap_filename);
            } else {
                addImage(image_filename,fitpsf_filename,shear_filename);
            }
        }
        if (isSkyMapIn_list) {
            Assert(_skymap_file_list.size() == _image_file_list.size());
        }
    } catch (std::exception& e) {
        xdbg<<"Caught std::exception: \n"<<e.what()<<std::endl;
        throw ReadException(
            "Error reading from "+file+" -- caught error\n" + e.what());
    } catch (...) {
        xdbg<<"Caught unknown exception"<<std::endl;
        throw ReadException(
            "Error reading from "+file+" -- caught unknown error");
    }

    dbg<<"Done reading file lists\n";
    Assert(_shear_file_list.size() == _image_file_list.size());
    Assert(_fitpsf_file_list.size() == _image_file_list.size());
}

void MultiShearCatalog::write() const
{
    xdbg<<"Start MultiShearCatalog.write\n";
    std::vector<std::string> file_list = MakeMultiName(_params, "multishear");

    const int nfiles = file_list.size();
    for(int i=0; i<nfiles; ++i) {
        const std::string& file = file_list[i];
        dbg<<"Writing multishear catalog to file: "<<file<<std::endl;

        bool isFitsIo = false;
        if (_params.keyExists("multishear_io")) {
            std::vector<std::string> io_list = _params["multishear_io"];
            Assert(io_list.size() == file_list.size());
            isFitsIo = (io_list[i] == "FITS");
        } else if (file.find("fits") != std::string::npos) {
            isFitsIo = true;
        }

        try {
            if (isFitsIo) {
                writeFits(file);
            } else {
                std::string delim = "  ";
                if (_params.keyExists("multishear_delim")) {
                    std::vector<std::string> delim_list = 
                        _params["multishear_delim"];
                    Assert(delim_list.size() == file_list.size());
                    delim = delim_list[i];
                } else if (file.find("csv") != std::string::npos) {
                    delim = ",";
                }
                writeAscii(file,delim);
            }
        } catch (CCfits::FitsException& e) {
            xdbg<<"Caught FitsException: \n"<<e.message()<<std::endl;
            throw WriteException(
                "Error writing to "+file+" -- caught error\n" + e.message());
        } catch (std::exception& e) {
            xdbg<<"Caught std::exception: \n"<<e.what()<<std::endl;
            throw WriteException(
                "Error writing to "+file+" -- caught error\n" + e.what());
        } catch (...) {
            xdbg<<"Caught unknown exception"<<std::endl;
            throw WriteException(
                "Error writing to "+file+" -- caught unknown error");
        }
    }
    dbg<<"Done Write ShearCatalog\n";
}

void MultiShearCatalog::writeFits(std::string file) const
{
    xdbg<<"Start MultiShearCatalog.writeFits "<<file<<std::endl;
    // ! means overwrite existing file
    CCfits::FITS fits("!"+file, CCfits::Write);
    xdbg<<"Opened file\n";

    std::vector<string> col_names;
    std::vector<string> col_fmts;

    col_names.push_back(_params.get("multishear_id_col")); 
    col_fmts.push_back("1J");
    int id_col_num = col_names.size()-1;

    col_names.push_back(_params.get("multishear_flags_col"));
    col_fmts.push_back("1J"); // flags
    int flags_col_num = col_names.size()-1;

    col_names.push_back(_params.get("multishear_ra_col"));
    col_fmts.push_back("1D"); // ra
    int ra_col_num = col_names.size()-1;

    col_names.push_back(_params.get("multishear_dec_col"));
    col_fmts.push_back("1D"); // decl
    int decl_col_num = col_names.size()-1;

    col_names.push_back(_params.get("multishear_shear1_col"));
    col_fmts.push_back("1D"); // shear1
    int shear1_col_num = col_names.size()-1;

    col_names.push_back(_params.get("multishear_shear2_col"));
    col_fmts.push_back("1D"); // shear2
    int shear2_col_num = col_names.size()-1;

    col_names.push_back(_params.get("multishear_nu_col"));
    col_fmts.push_back("1D"); // nu
    int nu_col_num = col_names.size()-1;

    col_names.push_back(_params.get("multishear_cov00_col"));
    col_fmts.push_back("1D"); // cov00
    int cov00_col_num = col_names.size()-1;

    col_names.push_back(_params.get("multishear_cov01_col"));
    col_fmts.push_back("1D"); // cov01
    int cov01_col_num = col_names.size()-1;

    col_names.push_back(_params.get("multishear_cov11_col"));
    col_fmts.push_back("1D"); // cov11
    int cov11_col_num = col_names.size()-1;

    col_names.push_back(_params.get("multishear_order_col"));
    col_fmts.push_back("1J"); // order
    int order_col_num = col_names.size()-1;

    col_names.push_back(_params.get("multishear_sigma_col"));
    col_fmts.push_back("1D"); // sigma
    int sigma_col_num = col_names.size()-1;

    col_names.push_back(_params.get("multishear_coeffs_col"));

    int ncoeff = _shape[0].size();
    dbg<<"ncoeff = "<<ncoeff<<std::endl;
    std::stringstream coeff_form;
    coeff_form << ncoeff << "D";
    col_fmts.push_back(coeff_form.str()); // shapelet coeffs

    int coeffs_col_num = col_names.size()-1;

    col_names.push_back(_params.get("multishear_nimages_found_col"));
    col_fmts.push_back("1J"); // nimages_found
    int nimages_found_col_num = col_names.size()-1;

    col_names.push_back(_params.get("multishear_nimages_gotpix_col"));
    col_fmts.push_back("1J"); // nimages_gotpix
    int nimages_gotpix_col_num = col_names.size()-1;

    col_names.push_back(_params.get("multishear_input_flags_col"));
    col_fmts.push_back("1J"); // input_flags
    int input_flags_col_num = col_names.size()-1;

    // We don't output units yet, just make it same size as others
    std::vector<string> col_units(col_fmts.size());

    xdbg<<"Before Create table"<<std::endl;
    CCfits::Table* table;
    table = fits.addTable("coadd_shearcat",size(),col_names,col_fmts,col_units);
    xdbg<<"Made table"<<std::endl;

    // Header keywords
#ifdef USE_TMV
    std::string tmvVers = tmv::TMV_Version();
#else
    std::string tmvVers = "Eigen";
#endif
    std::string wlVers = GetWlVersion();

    table->addKey("tmvvers", tmvVers, "version of TMV code");
    table->addKey("wlvers", wlVers, "version of weak lensing code");

    std::string str;
    double dbl;
    int intgr;
    WriteParamToTable(_params, table, "noise_method", str);
    WriteParamToTable(_params, table, "dist_method", str);

    WriteParamToTable(_params, table, "shear_aperture", dbl);
    WriteParamToTable(_params, table, "shear_max_aperture", dbl);
    WriteParamToTable(_params, table, "shear_gal_order", intgr);
    WriteParamToTable(_params, table, "shear_gal_order2", intgr);
    WriteParamToTable(_params, table, "shear_min_gal_size", dbl);
    WriteParamToTable(_params, table, "shear_f_psf", dbl);
    WriteParamToTable(_params, table, "peak_mem", dbl);

    // if wlmerun= is set we'll put it in the header.  This allows us to 
    // associate some more, possibly complicated, metadata with this file
    if (_params.keyExists("wlmerun")) {
        WriteParamToTable(_params, table, "wlmerun", str);
    }

    const int ngals = size();
    // data
    // make vector copies for writing
    std::vector<double> ra(ngals);
    std::vector<double> decl(ngals);
    std::vector<double> shear1(ngals);
    std::vector<double> shear2(ngals);
    std::vector<double> cov00(ngals);
    std::vector<double> cov01(ngals);
    std::vector<double> cov11(ngals);

    for(int i=0;i<ngals;++i) {
        // internally we use arcseconds
        ra[i] = _skypos[i].getX()/3600.;
        decl[i] = _skypos[i].getY()/3600.;
        shear1[i] = real(_shear[i]);
        shear2[i] = imag(_shear[i]);
        cov00[i] = _cov[i](0,0);
        cov01[i] = _cov[i](0,1);
        cov11[i] = _cov[i](1,1);
    }

    int start_row=1;

    table->column(col_names[id_col_num]).write(_id,start_row);
    table->column(col_names[flags_col_num]).write(_flags,start_row);
    table->column(col_names[ra_col_num]).write(ra,start_row);
    table->column(col_names[decl_col_num]).write(decl,start_row);
    table->column(col_names[shear1_col_num]).write(shear1,start_row);
    table->column(col_names[shear2_col_num]).write(shear2,start_row);
    table->column(col_names[nu_col_num]).write(_nu,start_row);
    table->column(col_names[cov00_col_num]).write(cov00,start_row);
    table->column(col_names[cov01_col_num]).write(cov01,start_row);
    table->column(col_names[cov11_col_num]).write(cov11,start_row);
    table->column(col_names[order_col_num]).write(_meas_galorder,start_row);

    for (int i=0; i<ngals; ++i) {
        int row = i+1;
        // Note: meas_galorder keeps track of the order of the shapelet that
        // was actually measured.  It is <= the full order of the shapelet
        // vector, but the higher order terms are set to zero.
        // Since we can't have different numbers of columns for each row, 
        // we have to write out the full shapelet order for all rows
        // including the zeros.

        double bSigma = _shape[i].getSigma();
        table->column(col_names[sigma_col_num]).write(&bSigma,1,row);

        // There is a bug in the CCfits Column::write function that the 
        // first argument has to be S* rather than const S*, even though
        // it is only read from. 
        // Hence the const_cast.
        double* ptr = const_cast<double*>(TMV_cptr(_shape[i].vec()));
        table->column(col_names[coeffs_col_num]).write(ptr, ncoeff, 1, row);
    }

    table->column(col_names[nimages_found_col_num]).write(_nimages_found,start_row);
    table->column(col_names[nimages_gotpix_col_num]).write(_nimages_gotpix,start_row);
    table->column(col_names[input_flags_col_num]).write(_input_flags,start_row);
}

void MultiShearCatalog::writeAscii(std::string file, std::string delim) const
{
    Assert(int(_id.size()) == size());
    Assert(int(_skypos.size()) == size());
    Assert(int(_flags.size()) == size());
    Assert(int(_shear.size()) == size());
    Assert(int(_nu.size()) == size());
    Assert(int(_cov.size()) == size());
    Assert(int(_meas_galorder.size()) == size());
    Assert(int(_shape.size()) == size());
    Assert(int(_nimages_found.size()) == size());
    Assert(int(_nimages_gotpix.size()) == size());
    Assert(int(_input_flags.size()) == size());

    std::ofstream fout(file.c_str());
    if (!fout) {
        throw WriteException("Error opening shear file"+file);
    }

    Form sci8; sci8.sci().trail(0).prec(8);
    Form fix3; fix3.fix().trail(0).prec(3);
    Form fix8; fix8.fix().trail(0).prec(8);
    Form fix6; fix6.fix().trail(0).prec(6);

    const int ngals = size();
    for(int i=0;i<ngals;++i) {
        fout
            << _id[i] << delim
            << fix8(_skypos[i].getX()/3600.) << delim
            << fix8(_skypos[i].getY()/3600.) << delim
            << _flags[i] << delim
            << sci8(real(_shear[i])) << delim
            << sci8(imag(_shear[i])) << delim
            << fix3(_nu[i]) << delim
            << sci8(_cov[i](0,0)) << delim
            << sci8(_cov[i](0,1)) << delim
            << sci8(_cov[i](1,1)) << delim
            << _nimages_found[i] << delim
            << _nimages_gotpix[i] << delim
            << _input_flags[i] << delim
            << _meas_galorder[i] << delim
            << fix6(_shape[i].getSigma()) << delim;
        const int ncoeffs = _shape[i].size();
        for(int j=0;j<ncoeffs;++j)
            fout << delim << sci8(_shape[i](j));
        fout << std::endl;
    }
}

void MultiShearCatalog::read()
{
    std::string file = MakeName(_params,"multishear",false,true);
    // false,true = input_prefix=false, mustexist=true.
    // It is an input here, but it is in the output_prefix directory.
    dbg<< "Reading Shear cat from file: " << file << std::endl;

    bool isFitsIo = false;
    if (_params.keyExists("multishear_io"))
        isFitsIo = (_params["multishear_io"] == "FITS");
    else if (file.find("fits") != std::string::npos)
        isFitsIo = true;

    if (!DoesFileExist(file)) {
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
        xdbg<<"Caught FitsException: \n"<<e.message()<<std::endl;
        throw ReadException(
            "Error reading from "+file+" -- caught error\n" + e.message());
    } catch (std::exception& e) {
        xdbg<<"Caught std::exception: \n"<<e.what()<<std::endl;
        throw ReadException(
            "Error reading from "+file+" -- caught error\n" + e.what());
    } catch (...) {
        xdbg<<"Caught unknown exception"<<std::endl;
        throw ReadException(
            "Error reading from "+file+" -- caught unknown error");
    }
    dbg<<"Done Read ShearCatalog\n";
}

void MultiShearCatalog::readFits(std::string file)
{
    int hdu = GetHdu(_params,"multishear",file,2);

    dbg<<"Opening MultiShearCatalog file "<<file<<" at hdu "<<hdu<<std::endl;
    CCfits::FITS fits(file, CCfits::Read);
    if (hdu > 1) fits.read(hdu-1);

    CCfits::ExtHDU& table=fits.extension(hdu-1);

    long nrows=table.rows();

    dbg<<"  nrows = "<<nrows<<std::endl;
    if (nrows <= 0) {
        throw ReadException(
            "ShearCatalog found to have 0 rows.  Must have > 0 rows.");
    }

    std::string id_col=_params.get("multishear_id_col");
    std::string ra_col=_params.get("multishear_ra_col");
    std::string decl_col=_params.get("multishear_dec_col");
    std::string sky_col=_params.get("multishear_sky_col");
    std::string noise_col=_params.get("multishear_noise_col");
    std::string flags_col=_params.get("multishear_flags_col");
    std::string shear1_col=_params.get("multishear_shear1_col");
    std::string shear2_col=_params.get("multishear_shear2_col");
    std::string nu_col=_params.get("multishear_nu_col");
    std::string cov00_col=_params.get("multishear_cov00_col");
    std::string cov01_col=_params.get("multishear_cov01_col");
    std::string cov11_col=_params.get("multishear_cov11_col");
    std::string order_col=_params.get("multishear_order_col");
    std::string sigma_col=_params.get("multishear_sigma_col");
    std::string coeffs_col=_params.get("multishear_coeffs_col");
    std::string nimages_found_col=_params.get("multishear_nimages_found_col");
    std::string nimages_gotpix_col=_params.get("multishear_nimages_gotpix_col");
    std::string input_flags_col=_params.get("multishear_input_flags_col");

    // MJ: This is not really optimal.
    // We should get this value from the fits header file, since it 
    // is stored there.  What CCfits command does this?
    int fullOrder = _params.get("shear_gal_order");

    long start=1;
    long end=nrows;

    dbg<<"Reading columns"<<std::endl;
    dbg<<"  "<<id_col<<std::endl;
    table.column(id_col).read(_id, start, end);

    dbg<<"  "<<ra_col<<"  "<<decl_col<<std::endl;
    _skypos.resize(nrows);
    std::vector<double> ra;
    std::vector<double> decl;
    table.column(ra_col).read(ra, start, end);
    table.column(decl_col).read(decl, start, end);
    for(long i=0;i<nrows;++i) _skypos[i] = Position(ra[i],decl[i]);

    dbg<<"  "<<flags_col<<std::endl;
    table.column(flags_col).read(_flags, start, end);

    dbg<<"  "<<shear1_col<<"  "<<shear2_col<<std::endl;
    _shear.resize(nrows);
    std::vector<double> shear1;
    std::vector<double> shear2;
    table.column(shear1_col).read(shear1, start, end);
    table.column(shear2_col).read(shear2, start, end);
    for(long i=0;i<nrows;++i) {
        _shear[i] = std::complex<double>(shear1[i],shear2[i]);
    }

    dbg<<"  "<<nu_col<<std::endl;
    table.column(nu_col).read(_nu, start, end);

    dbg<<"  "<<cov00_col<<"  "<<cov01_col<<"  "<<cov11_col<<std::endl;
    _cov.resize(nrows);
    std::vector<double> cov00;
    std::vector<double> cov01;
    std::vector<double> cov11;
    table.column(cov00_col).read(cov00, start, end);
    table.column(cov01_col).read(cov01, start, end);
    table.column(cov11_col).read(cov11, start, end);
    for(long i=0;i<nrows;++i) {
        _cov[i] << cov00[i], cov01[i], cov01[i], cov11[i];
    }

    dbg<<"  "<<sigma_col<<"  "<<order_col<<std::endl;
    _meas_galorder.resize(nrows);
    std::vector<double> sigma; // temporary
    table.column(sigma_col).read(sigma, start, end);
    table.column(order_col).read(_meas_galorder, start, end);

    _shape.reserve(nrows);
    for (int i=0; i<nrows; ++i) {
        int row=i+1;

        _shape.push_back(BVec(fullOrder,sigma[i]));
        int ncoeff=(fullOrder+1)*(fullOrder+2)/2;

        std::valarray<double> coeffs;
        table.column(coeffs_col).read(coeffs, row);

        double* ptri = TMV_ptr(_shape[i].vec());
        for (int j=0; j<ncoeff; ++j) {
            ptri[i] = coeffs[i];
        }
    }

    dbg<<"  "<<nimages_found_col<<std::endl;
    table.column(nimages_found_col).read(_nimages_found, start, end);

    dbg<<"  "<<nimages_gotpix_col<<std::endl;
    table.column(nimages_gotpix_col).read(_nimages_gotpix, start, end);

    dbg<<"  "<<input_flags_col<<std::endl;
    table.column(input_flags_col).read(_input_flags, start, end);
}

void MultiShearCatalog::readAscii(std::string file, std::string delim)
{
    std::ifstream fin(file.c_str());
    if (!fin) {
        throw ReadException("Error opening stars file"+file);
    }

    int fullOrder = _params.get("shear_gal_order");

    _id.clear(); _skypos.clear(); _flags.clear();
    _shear.clear(); _nu.clear(); _cov.clear();
    _meas_galorder.clear(); _shape.clear();
    _nimages_found.clear(); _nimages_gotpix.clear(); _input_flags.clear();

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
            _skypos.push_back(Position(ra,decl));
            _flags.push_back(flag);
            _shear.push_back(std::complex<double>(s1,s2));
            _nu.push_back(nu1);
            _cov.push_back(DSmallMatrix22());
            _cov.back() << c00, c01, c01, c11;
            _nimages_found.push_back(nfound);
            _nimages_gotpix.push_back(ngotpix);
            _input_flags.push_back(inflag);
            _meas_galorder.push_back(bOrder);
            _shape.push_back(BVec(fullOrder,bSigma));
            const int ncoeffs = _shape.back().size();
            for(int j=0;j<ncoeffs;++j) 
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
            _skypos.push_back(Position(ra,decl));
            getline(fin,temp,d); _flags.push_back(temp);
            getline(fin,temp,d); s1 = temp;
            getline(fin,temp,d); s2 = temp;
            _shear.push_back(std::complex<double>(s1,s2));
            getline(fin,temp,d); _nu.push_back(temp);
            getline(fin,temp,d); c00 = temp;
            getline(fin,temp,d); c01 = temp;
            getline(fin,temp,d); c11 = temp;
            _cov.push_back(DSmallMatrix22());
            _cov.back() << c00, c01, c01, c11;
            getline(fin,temp,d); _nimages_found.push_back(temp);
            getline(fin,temp,d); _nimages_gotpix.push_back(temp);
            getline(fin,temp,d); _input_flags.push_back(temp);
            getline(fin,temp,d); bOrder = temp;
            getline(fin,temp,d); bSigma = temp;
            _meas_galorder.push_back(bOrder);
            _shape.push_back(BVec(fullOrder,bSigma));
            const int ncoeffs = _shape.back().size();
            for(int j=0;j<ncoeffs-1;++j) {
                getline(fin,temp,d); _shape.back()(j) = temp;
            }
            // Last one doesn't have a following delimiter
            getline(fin,temp); _shape.back()(ncoeffs-1) = temp;
        }
    }
}

