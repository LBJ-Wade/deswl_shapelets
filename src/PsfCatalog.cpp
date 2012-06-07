
#include <sstream>
#include <valarray>
#include <CCfits/CCfits>

#include "dbg.h"
#include "PsfCatalog.h"
#include "Params.h"
#include "Pixel.h"
#include "Ellipse.h"
#include "Params.h"
#include "Name.h"
#include "Form.h"
#include "WlVersion.h"
#include "WriteParam.h"

void MeasureSinglePsf1(
    Position& cen, const Image<double>& im, double sky,
    const Transformation& trans,
    double noise, const Image<double>* weight_image,
    double sigma_p, const ConfigFile& params,
    PsfLog& log, BVec& psf, double& nu, long& flag)
{
    int psf_order = params.read<int>("psf_order");
    bool fixcen = params.read("psf_fix_centroid",false);
    double psf_ap = params.read<double>("psf_aperture");
    int maxm = params.read("psf_maxm",psf_order);

    std::vector<PixelList> pix(1);
    GetPixList(im,pix[0],cen,sky,noise,weight_image,trans,psf_ap,params,flag);

    int npix = pix[0].size();
    xdbg<<"npix = "<<npix<<std::endl;

    Ellipse ell;
    ell.fixGam();
    ell.fixMu();
    if (fixcen) ell.fixCen();
    else {
        ell.peakCentroid(pix[0],psf_ap/3.);
        ell.crudeMeasure(pix[0],sigma_p);
    }
    DMatrix cov(int(psf.size()),int(psf.size()));

    // First make sure it is centered.
    if (!(ell.measure(pix,psf_order,psf_order+4,maxm,sigma_p,flag,1.e-4))) {
        xdbg<<"Initial measurement failed\n";
        ++log._nf_psf;
        dbg<<"FLAG MEASURE_PSF_FAILED\n";
        flag |= MEASURE_PSF_FAILED;
        return;
    } 
    
    // Then measure it with the accurate center.
    psf.setSigma(sigma_p);
    if (ell.measureShapelet(pix,psf,psf_order,psf_order+4,maxm,&cov)) {
        ++log._ns_psf;
    } else {
        xdbg<<"Shapelet Measurement failed\n";
        ++log._nf_psf;
        dbg<<"FLAG MEASURE_PSF_FAILED\n";
        flag |= MEASURE_PSF_FAILED;
        return;
    }

    // Calculate the 0th order shapelet to make sure the flux is consistent 
    // with what we get at the full order.  It should be within a factor of 3
    // of the same value.  (Within 10s of percent usually, but we only call it
    // an error if it is more than a factor of 3 different.)
    BVec flux(0,sigma_p);
    DMatrix flux_cov(1,1);
    if (!ell.measureShapelet(pix,flux,0,0,0,&flux_cov)) {
        xdbg<<"Measurement of flux failed.\n";
        ++log._nf_psf;
        dbg<<"FLAG PSF_BAD_FLUX\n";
        flag |= PSF_BAD_FLUX;
    }
    nu = flux(0) / std::sqrt(flux_cov(0,0));
    dbg<<"nu = "<<flux(0)<<" / sqrt("<<flux_cov(0,0)<<") = "<<nu<<std::endl;
    dbg<<" or  "<<psf(0)<<" / sqrt("<<cov(0,0)<<") = "<<
        psf(0)/std::sqrt(cov(0,0))<<std::endl;
    if (!(flux(0) > 0.0 &&
          psf(0) >= flux(0)/3. &&
          psf(0) <= flux(0)*3.)) {
        dbg<<"Bad flux value: \n";
        dbg<<"flux = "<<flux(0)<<std::endl;
        dbg<<"psf = "<<psf.vec()<<std::endl;
        dbg<<"FLAG PSF_BAD_FLUX\n";
        flag |= PSF_BAD_FLUX;
    }

    xdbg<<"psf = "<<psf.vec()<<std::endl;
    psf.normalize();  // Divide by (0,0) element
    xdbg<<"Normalized psf: "<<psf.vec()<<std::endl;
    if (!fixcen) cen += ell.getCen();
}

void MeasureSinglePsf(
    Position& cen, const Image<double>& im, double sky,
    const Transformation& trans,
    double noise, const Image<double>* weight_image,
    double sigma_p, const ConfigFile& params,
    PsfLog& log, BVec& psf, double& nu, long& flag)
{
    try {
        // We don't need to save skyPos.  We just want to catch the range
        // error here, so we don't need to worry about it for dudx, etc.
        Position skyPos;
        trans.transform(cen,skyPos);
        dbg<<"skypos = "<<skyPos<<std::endl;
    } catch (RangeException& e) {
        xdbg<<"skip: transformation range error: \n";
        xdbg<<"p = "<<cen<<", b = "<<e.getBounds()<<std::endl;
        ++log._nf_range;
        dbg<<"FLAG TRANSFORM_EXCEPTION\n";
        flag |= TRANSFORM_EXCEPTION;
        return;
    }

    try {
        MeasureSinglePsf1(
            cen,im,sky,trans,noise,weight_image,
            sigma_p,params,log,psf,nu,flag);
#ifdef USE_TMV
    } catch (tmv::Error& e) {
        dbg<<"TMV Error thrown in MeasureSinglePSF\n";
        dbg<<e<<std::endl;
        ++log._nf_tmv_error;
        dbg<<"FLAG TMV_EXCEPTION\n";
        flag |= TMV_EXCEPTION;
#endif
    } catch (...) {
        dbg<<"unkown exception in MeasureSinglePSF\n";
        ++log._nf_other_error;
        dbg<<"FLAG UNKNOWN_EXCEPTION\n";
        flag |= UNKNOWN_EXCEPTION;
    }
}

PsfCatalog::PsfCatalog(
    const StarCatalog& starcat, const ConfigFile& params) :
    _params(params)
{
    dbg<<"Create PSFCatalog\n";
    const int ntot = starcat.size();
    dbg<<"ntot = "<<ntot<<std::endl;
    const std::vector<bool>& is_star = starcat.getIsStarList();
    const int nstars = std::count(is_star.begin(),is_star.end(),1);
    dbg<<"count stars = "<<nstars<<std::endl;
    _id.reserve(nstars);
    _pos.reserve(nstars);
    _flags.reserve(nstars);
    for (int i=0; i<ntot; ++i) {
        if (is_star[i]) {
            _id.push_back( starcat.getId(i) );
            _pos.push_back( starcat.getPos(i) );
            _sky.push_back( starcat.getSky(i) );
            _noise.push_back( starcat.getNoise(i) );
            _flags.push_back( starcat.getFlags(i) );
        }
    }
    Assert(int(_id.size()) == nstars);
    dbg<<"nstars = "<<nstars<<std::endl;

    // Fix flags to only have INPUT_FLAG set.
    // i.e. Ignore any flags set by StarCatalog.
    for (int i=0; i<nstars; ++i) {
        if (_flags[i]) _flags[i] = INPUT_FLAG;
    }

    // Set up a default psf vector for output when an object measurement
    // fails
    long psf_order = _params.read<long>("psf_order");
    BVec psfDefault(psf_order,1.);
    const int psfSize = psfDefault.size();
    for (int i=0; i<psfSize; ++i) {
        psfDefault(i) = DEFVALNEG;
    }
    double nuDefault = DEFVALNEG;

    _nu.resize(_id.size(),nuDefault);
    _psf.resize(_id.size(),psfDefault);

    Assert(int(_id.size()) == size());
    Assert(int(_pos.size()) == size());
    Assert(int(_sky.size()) == size());
    Assert(int(_noise.size()) == size());
    Assert(int(_flags.size()) == size());
    Assert(int(_nu.size()) == size());
    Assert(int(_psf.size()) == size());
}

PsfCatalog::PsfCatalog(const ConfigFile& params) : _params(params)
{
    Assert(int(_id.size()) == size());
    Assert(int(_pos.size()) == size());
    Assert(int(_sky.size()) == size());
    Assert(int(_noise.size()) == size());
    Assert(int(_flags.size()) == size());
    Assert(int(_nu.size()) == size());
    Assert(int(_psf.size()) == size());
}

void PsfCatalog::writeFits(std::string file) const
{
    Assert(int(_id.size()) == size());
    Assert(int(_pos.size()) == size());
    Assert(int(_nu.size()) == size());
    Assert(int(_flags.size()) == size());
    Assert(int(_psf.size()) == size());
    const int npsf = size();

    // ! means overwrite existing file
    CCfits::FITS fits("!"+file, CCfits::Write);

    const int nFields = 10;

    std::vector<string> col_names(nFields);
    std::vector<string> col_fmts(nFields);
    std::vector<string> col_units(nFields);

    col_names[0] = _params["psf_id_col"];
    col_names[1] = _params["psf_x_col"];
    col_names[2] = _params["psf_y_col"];
    col_names[3] = _params["psf_sky_col"];
    col_names[4] = _params["psf_noise_col"];
    col_names[5] = _params["psf_flags_col"];
    col_names[6] = _params["psf_nu_col"];
    col_names[7] = _params["psf_order_col"];
    col_names[8] = _params["psf_sigma_col"];
    col_names[9] = _params["psf_coeffs_col"];

    int ncoeff = _psf[0].size();
    dbg<<"ncoeff = "<<ncoeff<<std::endl;
    std::stringstream coeff_form;
    coeff_form << ncoeff << "d";

    col_fmts[0] = "1J"; //id
    col_fmts[1] = "1D"; //x
    col_fmts[2] = "1D"; //y
    col_fmts[3] = "1D"; //sky
    col_fmts[4] = "1D"; //noise
    col_fmts[5] = "1J"; //flags
    col_fmts[6] = "1D"; //nu
    col_fmts[7] = "1J"; //order
    col_fmts[8] = "1D"; //sigma
    col_fmts[9] = coeff_form.str(); //coeffs

    col_units[0] = "None";     // id
    col_units[1] = "pixels";   // x
    col_units[2] = "pixels";   // y
    col_units[3] = "ADU";      // sky
    col_units[4] = "ADU^2";    // noise
    col_units[5] = "None";     // flags
    col_units[6] = "None";     // nu
    col_units[7] = "None";     // order
    col_units[8] = "Arcsec";   // sigma
    col_units[9] = "None";     // coeffs


    dbg<<"Before Create table"<<std::endl;
    CCfits::Table* table;
    table = fits.addTable("psfcat",npsf,col_names,col_fmts,col_units);


    // Header keywords
    std::string str;
    double dbl;
    int intgr;

#ifdef USE_TMV
    std::string tmvVers = tmv::TMV_Version();
#else
    std::string tmvVers = "Eigen";
#endif
    std::string wlVers = GetWlVersion();

    table->addKey("tmvvers", tmvVers, "version of TMV code");
    table->addKey("wlvers", wlVers, "version of weak lensing code");

    // if wlserun= is sent we'll put it in the header.  This allows us to 
    // associate some more, possibly complicated, metadata with this file
    if ( _params.keyExists("wlserun") ) {
        WriteParamToTable(_params, table, "wlserun", str);
    }

    WriteParamToTable(_params, table, "noise_method", str);
    WriteParamToTable(_params, table, "dist_method", str);

    WriteParamToTable(_params, table, "psf_aperture", dbl);
    WriteParamToTable(_params, table, "psf_order", intgr);
    WriteParamToTable(_params, table, "psf_seeing_est", dbl);


    // Data columns

    // make vector copies for writing
    std::vector<double> x(npsf);
    std::vector<double> y(npsf);
    for(int i=0;i<npsf;++i) {
        x[i] = _pos[i].getX();
        y[i] = _pos[i].getY();
    }

    int startrow=1;

    table->column(col_names[0]).write(_id,startrow);
    table->column(col_names[1]).write(x,startrow);
    table->column(col_names[2]).write(y,startrow);
    table->column(col_names[3]).write(_sky,startrow);
    table->column(col_names[4]).write(_noise,startrow);
    table->column(col_names[5]).write(_flags,startrow);

    table->column(col_names[6]).write(_nu,startrow);

    for (int i=0; i<npsf; ++i) {
        int row = i+1;
        long border = _psf[i].getOrder();
        double bsigma = _psf[i].getSigma();

        table->column(col_names[7]).write(&border,1,row);
        table->column(col_names[8]).write(&bsigma,1,row);
        // There is a bug in the CCfits Column::write function that the 
        // first argument has to be S* rather than const S*, even though
        // it is only read from. 
        // Hence the const_cast.
        double* ptr = const_cast<double*>(TMV_cptr(_psf[i].vec()));
        table->column(col_names[9]).write(ptr, ncoeff, 1, row);
    }
}

void PsfCatalog::writeAscii(std::string file, std::string delim) const
{
    Assert(int(_id.size()) == size());
    Assert(int(_pos.size()) == size());
    Assert(int(_sky.size()) == size());
    Assert(int(_noise.size()) == size());
    Assert(int(_flags.size()) == size());
    Assert(int(_nu.size()) == size());
    Assert(int(_psf.size()) == size());
    const int npsf = size();

    std::ofstream fout(file.c_str());
    if (!fout) {
        throw WriteException("Error opening psf file"+file);
    }

    Form sci8; sci8.sci().trail(0).prec(8);
    Form fix3; fix3.fix().trail(0).prec(3);
    Form fix6; fix6.fix().trail(0).prec(3);

    for(int i=0;i<npsf;++i) {
        fout
            << _id[i] << delim
            << fix3(_pos[i].getX()) << delim
            << fix3(_pos[i].getY()) << delim
            << fix3(_sky[i]) << delim
            << sci8(_noise[i]) << delim
            << _flags[i] << delim
            << fix3(_nu[i]) << delim
            << _psf[i].getOrder() << delim
            << fix6(_psf[i].getSigma());
        const int ncoeff = _psf[i].size();
        for(int j=0;j<ncoeff;++j) {
            fout << delim << sci8(_psf[i](j));
        }
        fout << std::endl;
    }
}

void PsfCatalog::write() const
{
    std::vector<std::string> file_list = MakeMultiName(_params, "psf");  

    const int nFiles = file_list.size();
    for(int i=0; i<nFiles; ++i) {
        const std::string& file = file_list[i];
        dbg<<"Writing psf catalog to file: "<<file<<std::endl;

        bool isFitsIo = false;
        if (_params.keyExists("psf_io")) {
            std::vector<std::string> io_list = _params["psf_io"];
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
                if (_params.keyExists("psf_delim")) {
                    std::vector<std::string> delim_list = _params["psf_delim"];
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
    dbg<<"Done Write PSFCatalog\n";
}

void PsfCatalog::readFits(std::string file)
{
    int hdu = GetHdu(_params,"psf",file,2);

    dbg<<"Opening PsfCatalog file "<<file<<" at hdu "<<hdu<<std::endl;
    CCfits::FITS fits(file, CCfits::Read);
    if (hdu > 1) fits.read(hdu-1);

    CCfits::ExtHDU& table=fits.extension(hdu-1);

    long nrows=table.rows();

    dbg<<"  nrows = "<<nrows<<std::endl;
    if (nrows <= 0) {
        throw ReadException(
            "PSFCatalog found to have 0 rows.  Must have > 0 rows.");
    }

    std::string id_col=_params.get("psf_id_col");
    std::string x_col=_params.get("psf_x_col");
    std::string y_col=_params.get("psf_y_col");
    std::string sky_col=_params.get("psf_sky_col");
    std::string noise_col=_params.get("psf_noise_col");
    std::string flags_col=_params.get("psf_flags_col");
    std::string nu_col=_params.get("psf_nu_col");
    std::string order_col=_params.get("psf_order_col");
    std::string sigma_col=_params.get("psf_sigma_col");
    std::string coeffs_col=_params.get("psf_coeffs_col");

    long start=1;
    long end=nrows;


    dbg<<"Reading columns"<<std::endl;
    dbg<<"  "<<id_col<<std::endl;
    table.column(id_col).read(_id, start, end);

    dbg<<"  "<<x_col<<"  "<<y_col<<std::endl;
    _pos.resize(nrows);
    std::vector<double> x;
    std::vector<double> y;
    table.column(x_col).read(x, start, end);
    table.column(y_col).read(y, start, end);
    for(long i=0;i<nrows;++i) _pos[i] = Position(x[i],y[i]);

    dbg<<"  "<<sky_col<<std::endl;
    table.column(sky_col).read(_sky, start, end);

    dbg<<"  "<<noise_col<<std::endl;
    table.column(noise_col).read(_noise, start, end);

    dbg<<"  "<<flags_col<<std::endl;
    table.column(flags_col).read(_flags, start, end);

    dbg<<"  "<<nu_col<<std::endl;
    table.column(nu_col).read(_nu, start, end);

    // these are temporary
    std::vector<double> sigma;
    std::vector<long> order;

    dbg<<"  "<<sigma_col<<std::endl;
    table.column(sigma_col).read(sigma, start, end);

    dbg<<"  "<<order_col<<std::endl;
    table.column(order_col).read(order, start, end);


    // gotta loop for this one
    _psf.reserve(nrows);
    for (int i=0; i<nrows; ++i) {
        int row=i+1;

        _psf.push_back(BVec(order[i],sigma[i]));
        int ncoeff=(order[i]+1)*(order[i]+2)/2;
        // although we are allowed to write lots of different ways, the
        // reading is less flexible.  We can *only* read a vector
        // column into a valarray, period
        std::valarray<double> coeffs;
        table.column(coeffs_col).read(coeffs, row);
        for (int j=0; j<ncoeff; ++j) _psf[i](j) = coeffs[j];
    }
}

void PsfCatalog::readAscii(std::string file, std::string delim)
{
    std::ifstream fin(file.c_str());
    if (!fin) {
        throw ReadException("Error opening psf file"+file);
    }

    _id.clear(); _pos.clear(); _sky.clear(); _noise.clear(); _flags.clear();
    _nu.clear(); _psf.clear();

    if (delim == "  ") {
        ConvertibleString flag;
        long id1,border;
        double x,y,sky1,noise1,nu1,bsigma;
        while ( fin
                >> id1
                >> x >> y
                >> sky1
                >> noise1
                >> flag
                >> nu1
                >> border >> bsigma) {
            _id.push_back(id1);
            _pos.push_back(Position(x,y));
            _sky.push_back(sky1);
            _noise.push_back(noise1);
            _flags.push_back(flag);
            _nu.push_back(nu1);
            _psf.push_back(BVec(border,bsigma));
            const int ncoeffs = _psf.back().size();
            for(int j=0;j<ncoeffs;++j) fin >> _psf.back()(j);
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
        ConvertibleString temp;
        long border;
        double x,y,bsigma;
        while (getline(fin,temp,d)) {
            _id.push_back(temp);
            getline(fin,temp,d); x = temp;
            getline(fin,temp,d); y = temp;
            _pos.push_back(Position(x,y));
            getline(fin,temp,d); _sky.push_back(temp);
            getline(fin,temp,d); _noise.push_back(temp);
            getline(fin,temp,d); _flags.push_back(temp);
            getline(fin,temp,d); _nu.push_back(temp);
            getline(fin,temp,d); border = temp;
            getline(fin,temp,d); bsigma = temp;
            _psf.push_back(BVec(border,bsigma));
            const int ncoeffs = _psf.back().size();
            for(int j=0;j<ncoeffs-1;++j) {
                getline(fin,temp,d); _psf.back()(j) = temp;
            }
            // Last one doesn't have a following delimiter
            getline(fin,temp); _psf.back()(ncoeffs-1) = temp;
        }
    }
}

void PsfCatalog::read()
{
    std::string file = MakeName(_params,"psf",false,true);
    read(file);
}

void PsfCatalog::read(std::string file)
{
    // false,true = input_prefix=false, mustexist=true.
    // It is an input here, but it is in the output_prefix directory.
    dbg<< "Reading PSF cat from file: " << file << std::endl;

    bool isFitsIo = false;
    if (_params.keyExists("psf_io")) 
        isFitsIo = (_params["psf_io"] == "FITS");
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
            if (_params.keyExists("psf_delim")) delim = _params["psf_delim"];
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
    dbg<<"Done Read PSFCatalog\n";
}

