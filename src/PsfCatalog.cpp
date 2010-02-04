
#include <sstream>
#include <valarray>
#include <CCfits/CCfits>

#include "PsfCatalog.h"
#include "dbg.h"
#include "Params.h"
#include "Pixel.h"
#include "Ellipse.h"
#include "Params.h"
#include "Name.h"
#include "Form.h"
#include "WlVersion.h"
#include "WriteParam.h"

void measureSinglePsf1(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans,
    double noise, double gain, const Image<double>* weightIm,
    double sigmaP, double psfAp, int psfOrder,
    PsfLog& log, BVec& psf, double& nu, long& flag)
{
    std::vector<PixelList> pix(1);
    getPixList(im,pix[0],cen,sky,noise,gain,weightIm,trans,psfAp,flag);

    int nPix = pix[0].size();
    xdbg<<"npix = "<<nPix<<std::endl;

    Ellipse ell;
    ell.fixGam();
    ell.fixMu();
    ell.peakCentroid(pix[0],psfAp/3.);
    ell.crudeMeasure(pix[0],sigmaP);
    DMatrix cov(psf.size(),psf.size());

    long flag1=0;
    if (ell.measure(pix,psfOrder,sigmaP,false,flag1,0,&psf,&cov)) {
        ++log._nsPsf;
    } else {
        xdbg<<"Measurement failed\n";
        ++log._nfPsf;
        flag |= MEASURE_PSF_FAILED;
    }

    // Calculate the 0th order shapelet to make sure the flux is consistent 
    // with what we get at the full order.  It should be within a factor of 3
    // of the same value.  (Within 10s of percent usually, but we only call it
    // an error if it is more than a factor of 3 different.)
    BVec flux(0,sigmaP);
    DMatrix fluxCov(1,1);
    ell.measureShapelet(pix,flux,0,&fluxCov);
    nu = flux(0) / std::sqrt(fluxCov(0,0));
    dbg<<"nu = "<<flux(0)<<" / sqrt("<<fluxCov(0,0)<<") = "<<nu<<std::endl;
    dbg<<" or  "<<psf(0)<<" / sqrt("<<cov(0,0)<<") = "<<
        psf(0)/std::sqrt(cov(0,0))<<std::endl;
    if (!(flux(0) > 0.0 &&
          psf(0) >= flux(0)/3. &&
          psf(0) <= flux(0)*3.)) {
        dbg<<"Bad flux value: \n";
        dbg<<"flux = "<<flux(0)<<std::endl;
        dbg<<"psf = "<<psf.vec()<<std::endl;
        flag |= PSF_BAD_FLUX;
    }

    xdbg<<"psf = "<<psf.vec()<<std::endl;
    psf.normalize();  // Divide by (0,0) element
    xdbg<<"Normalized psf: "<<psf.vec()<<std::endl;
}

void measureSinglePsf(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans,
    double noise, double gain, const Image<double>* weightIm,
    double sigmaP, double psfAp, int psfOrder,
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
        ++log._nfRange;
        flag |= TRANSFORM_EXCEPTION;
        return;
    }

    try {
        measureSinglePsf1(
	    cen,im,sky,trans,noise,gain,weightIm,
	    sigmaP,psfAp,psfOrder,log,psf,nu,flag);
#ifdef USE_TMV
    } catch (tmv::Error& e) {
        dbg<<"TMV Error thrown in MeasureSinglePSF\n";
        dbg<<e<<std::endl;
        ++log._nfTmvError;
        flag |= TMV_EXCEPTION;
#endif
    } catch (...) {
        dbg<<"unkown exception in MeasureSinglePSF\n";
        ++log._nfOtherError;
        flag |= UNKNOWN_EXCEPTION;
    }
}

PsfCatalog::PsfCatalog(
    const StarCatalog& starCat, const ConfigFile& params) :
    _params(params)
{
    dbg<<"Create PSFCatalog\n";
    const int nTot = starCat.size();
    dbg<<"ntot = "<<nTot<<std::endl;
    const std::vector<bool>& isAStar = starCat.getIsAStarList();
    const int nStars = std::count(isAStar.begin(),isAStar.end(),1);
    dbg<<"count stars = "<<nStars<<std::endl;
    _id.reserve(nStars);
    _pos.reserve(nStars);
    _flags.reserve(nStars);
    for (int i=0; i<nTot; ++i) {
        if (starCat.isAStar(i)) {
            _id.push_back( starCat.getId(i) );
            _pos.push_back( starCat.getPos(i) );
            _sky.push_back( starCat.getSky(i) );
            _noise.push_back( starCat.getNoise(i) );
            _flags.push_back( starCat.getFlags(i) );
        }
    }
    Assert(int(_id.size()) == nStars);
    dbg<<"nstars = "<<nStars<<std::endl;

    // Fix flags to only have INPUT_FLAG set.
    // i.e. Ignore any flags set by StarCatalog.
    for (int i=0; i<nStars; ++i) {
        _flags[i] &= INPUT_FLAG;
    }

    // Set up a default psf vector for output when an object measurement
    // fails
    long psfOrder = _params.read<long>("psf_order");
    BVec psfDefault(psfOrder,1.);
    const int psfSize = psfDefault.size();
    for (int i=0; i<psfSize; ++i) {
        psfDefault(i) = DEFVALNEG;
    }
    double nuDefault = DEFVALNEG;

    _nu.resize(_id.size(),nuDefault);
    _psf.resize(_id.size(),psfDefault);

    Assert(_id.size() == size());
    Assert(_pos.size() == size());
    Assert(_sky.size() == size());
    Assert(_noise.size() == size());
    Assert(_flags.size() == size());
    Assert(_nu.size() == size());
    Assert(_psf.size() == size());
}

PsfCatalog::PsfCatalog(const ConfigFile& params) : _params(params)
{
    read();

    Assert(_id.size() == size());
    Assert(_pos.size() == size());
    Assert(_sky.size() == size());
    Assert(_noise.size() == size());
    Assert(_flags.size() == size());
    Assert(_nu.size() == size());
    Assert(_psf.size() == size());
}

// this one we don't have to deal with root= stuff
PsfCatalog::PsfCatalog(const ConfigFile& params, std::string file) : 
    _params(params)
{
    read(file);

    Assert(_id.size() == size());
    Assert(_pos.size() == size());
    Assert(_sky.size() == size());
    Assert(_noise.size() == size());
    Assert(_flags.size() == size());
    Assert(_nu.size() == size());
    Assert(_psf.size() == size());
}

void PsfCatalog::writeFits(std::string file) const
{
    Assert(_id.size() == size());
    Assert(_pos.size() == size());
    Assert(_nu.size() == size());
    Assert(_flags.size() == size());
    Assert(_psf.size() == size());
    const int nPsf = size();

    // ! means overwrite existing file
    CCfits::FITS fits("!"+file, CCfits::Write);

    const int nFields = 10;

    std::vector<string> colNames(nFields);
    std::vector<string> colFmts(nFields);
    std::vector<string> colUnits(nFields);

    colNames[0] = _params["psf_id_col"];
    colNames[1] = _params["psf_x_col"];
    colNames[2] = _params["psf_y_col"];
    colNames[3] = _params["psf_sky_col"];
    colNames[4] = _params["psf_noise_col"];
    colNames[5] = _params["psf_flags_col"];
    colNames[6] = _params["psf_nu_col"];
    colNames[7] = _params["psf_order_col"];
    colNames[8] = _params["psf_sigma_col"];
    colNames[9] = _params["psf_coeffs_col"];

    int nCoeff = _psf[0].size();
    dbg<<"ncoeff = "<<nCoeff<<std::endl;
    std::stringstream coeffForm;
    coeffForm << nCoeff << "d";

    colFmts[0] = "1J"; //id
    colFmts[1] = "1D"; //x
    colFmts[2] = "1D"; //y
    colFmts[3] = "1D"; //sky
    colFmts[4] = "1D"; //noise
    colFmts[5] = "1J"; //flags
    colFmts[6] = "1D"; //nu
    colFmts[7] = "1J"; //order
    colFmts[8] = "1D"; //sigma
    colFmts[9] = coeffForm.str(); //coeffs

    colUnits[0] = "None";     // id
    colUnits[1] = "pixels";   // x
    colUnits[2] = "pixels";   // y
    colUnits[3] = "ADU";      // sky
    colUnits[4] = "ADU^2";    // noise
    colUnits[5] = "None";     // flags
    colUnits[6] = "None";     // nu
    colUnits[7] = "None";     // order
    colUnits[8] = "Arcsec";   // sigma
    colUnits[9] = "None";     // coeffs


    dbg<<"Before Create table"<<std::endl;
    CCfits::Table* table;
    table = fits.addTable("psfcat",nPsf,colNames,colFmts,colUnits);


    // Header keywords
    std::string str;
    double dbl;
    int intgr;

#ifdef USE_TMV
    std::string tmvVers = tmv::TMV_Version();
#else
    std::string tmvVers = "Eigen";
#endif
    std::string wlVers = getWlVersion();

    table->addKey("tmvvers", tmvVers, "version of TMV code");
    table->addKey("wlvers", wlVers, "version of weak lensing code");


    // if serun= is sent we'll put it in the header.  This allows us to 
    // associate some more, possibly complicated, metadata with this file
    if ( _params.keyExists("serun") ) {
        writeParamToTable(_params, table, "serun", str);
    }

    writeParamToTable(_params, table, "noise_method", str);
    writeParamToTable(_params, table, "dist_method", str);

    writeParamToTable(_params, table, "psf_aperture", dbl);
    writeParamToTable(_params, table, "psf_order", intgr);
    writeParamToTable(_params, table, "psf_seeing_est", dbl);



    // Data columns

    // make vector copies for writing
    std::vector<double> x(nPsf);
    std::vector<double> y(nPsf);
    for(int i=0;i<nPsf;++i) {
        x[i] = _pos[i].getX();
        y[i] = _pos[i].getY();
    }

    int startrow=1;

    table->column(colNames[0]).write(_id,startrow);
    table->column(colNames[1]).write(x,startrow);
    table->column(colNames[2]).write(y,startrow);
    table->column(colNames[3]).write(_sky,startrow);
    table->column(colNames[4]).write(_noise,startrow);
    table->column(colNames[5]).write(_flags,startrow);

    table->column(colNames[6]).write(_nu,startrow);

    for (int i=0; i<nPsf; ++i) {
        int row = i+1;
        long bOrder = _psf[i].getOrder();
        double bSigma = _psf[i].getSigma();

        table->column(colNames[7]).write(&bOrder,1,row);
        table->column(colNames[8]).write(&bSigma,1,row);
        double* cptr = const_cast<double *>(TMV_cptr(_psf[i].vec()));
        table->column(colNames[9]).write(cptr, nCoeff, 1, row);
    }
}

void PsfCatalog::writeAscii(std::string file, std::string delim) const
{
    Assert(_id.size() == size());
    Assert(_pos.size() == size());
    Assert(_sky.size() == size());
    Assert(_noise.size() == size());
    Assert(_flags.size() == size());
    Assert(_nu.size() == size());
    Assert(_psf.size() == size());
    const int nPsf = size();

    std::ofstream fout(file.c_str());
    if (!fout) {
        throw WriteException("Error opening psf file"+file);
    }

    Form sci8; sci8.sci().trail(0).prec(8);
    Form fix3; fix3.fix().trail(0).prec(3);
    Form fix6; fix6.fix().trail(0).prec(3);

    for(int i=0;i<nPsf;++i) {
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
        const int nCoeff = _psf[i].size();
        for(int j=0;j<nCoeff;++j) {
            fout << delim << sci8(_psf[i](j));
        }
        fout << std::endl;
    }
}

void PsfCatalog::write() const
{
    std::vector<std::string> fileList = makeMultiName(_params, "psf");  

    const int nFiles = fileList.size();
    for(int i=0; i<nFiles; ++i) {
        const std::string& file = fileList[i];
        dbg<<"Writing psf catalog to file: "<<file<<std::endl;

        bool isFitsIo = false;
        if (_params.keyExists("psf_io")) {
            std::vector<std::string> ioList = _params["psf_io"];
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
                if (_params.keyExists("psf_delim")) {
                    std::vector<std::string> delimList = _params["psf_delim"];
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
    dbg<<"Done Write PSFCatalog\n";
}

void PsfCatalog::readFits(std::string file)
{
    int hdu = getHdu(_params,"psf",file,2);

    dbg<<"Opening PsfCatalog file "<<file<<" at hdu "<<hdu<<std::endl;
    CCfits::FITS fits(file, CCfits::Read);
    if (hdu > 1) fits.read(hdu-1);

    CCfits::ExtHDU& table=fits.extension(hdu-1);

    long nRows=table.rows();

    dbg<<"  nrows = "<<nRows<<std::endl;
    if (nRows <= 0) {
        throw ReadException(
            "PSFCatalog found to have 0 rows.  Must have > 0 rows.");
    }

    std::string idCol=_params.get("psf_id_col");
    std::string xCol=_params.get("psf_x_col");
    std::string yCol=_params.get("psf_y_col");
    std::string skyCol=_params.get("psf_sky_col");
    std::string noiseCol=_params.get("psf_noise_col");
    std::string flagsCol=_params.get("psf_flags_col");
    std::string nuCol=_params.get("psf_nu_col");
    std::string orderCol=_params.get("psf_order_col");
    std::string sigmaCol=_params.get("psf_sigma_col");
    std::string coeffsCol=_params.get("psf_coeffs_col");

    long start=1;
    long end=nRows;


    dbg<<"Reading columns"<<std::endl;
    dbg<<"  "<<idCol<<std::endl;
    table.column(idCol).read(_id, start, end);

    dbg<<"  "<<xCol<<"  "<<yCol<<std::endl;
    _pos.resize(nRows);
    std::vector<double> x;
    std::vector<double> y;
    table.column(xCol).read(x, start, end);
    table.column(yCol).read(y, start, end);
    for(long i=0;i<nRows;++i) _pos[i] = Position(x[i],y[i]);

    dbg<<"  "<<skyCol<<std::endl;
    table.column(skyCol).read(_sky, start, end);

    dbg<<"  "<<noiseCol<<std::endl;
    table.column(noiseCol).read(_noise, start, end);

    dbg<<"  "<<flagsCol<<std::endl;
    table.column(flagsCol).read(_flags, start, end);

    dbg<<"  "<<nuCol<<std::endl;
    table.column(nuCol).read(_nu, start, end);

    // these are temporary
    std::vector<double> sigma;
    std::vector<long> order;

    dbg<<"  "<<sigmaCol<<std::endl;
    table.column(sigmaCol).read(sigma, start, end);

    dbg<<"  "<<orderCol<<std::endl;
    table.column(orderCol).read(order, start, end);


    // gotta loop for this one
    _psf.reserve(nRows);
    for (int i=0; i<nRows; ++i) {
        int row=i+1;

        _psf.push_back(BVec(order[i],sigma[i]));
        int nCoeff=(order[i]+1)*(order[i]+2)/2;
        // although we are allowed to write lots of different ways, the
        // reading is less flexible.  We can *only* read a vector
        // column into a valarray, period
        std::valarray<double> coeffs;
        table.column(coeffsCol).read(coeffs, row);
        for (int j=0; j<nCoeff; ++j) _psf[i](j) = coeffs[j];
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
        long id1,bOrder;
        double x,y,sky1,noise1,nu1,bSigma;
        while ( fin
                >> id1
                >> x >> y
                >> sky1
                >> noise1
                >> flag
                >> nu1
                >> bOrder >> bSigma) {
            _id.push_back(id1);
            _pos.push_back(Position(x,y));
            _sky.push_back(sky1);
            _noise.push_back(noise1);
            _flags.push_back(flag);
            _nu.push_back(nu1);
            _psf.push_back(BVec(bOrder,bSigma));
            const int nCoeffs = _psf.back().size();
            for(int j=0;j<nCoeffs;++j) fin >> _psf.back()(j);
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
        long bOrder;
        double x,y,bSigma;
        while (getline(fin,temp,d)) {
            _id.push_back(temp);
            getline(fin,temp,d); x = temp;
            getline(fin,temp,d); y = temp;
            _pos.push_back(Position(x,y));
            getline(fin,temp,d); _sky.push_back(temp);
            getline(fin,temp,d); _noise.push_back(temp);
            getline(fin,temp,d); _flags.push_back(temp);
            getline(fin,temp,d); _nu.push_back(temp);
            getline(fin,temp,d); bOrder = temp;
            getline(fin,temp,d); bSigma = temp;
            _psf.push_back(BVec(bOrder,bSigma));
            const int nCoeffs = _psf.back().size();
            for(int j=0;j<nCoeffs-1;++j) {
                getline(fin,temp,d); _psf.back()(j) = temp;
            }
            // Last one doesn't have a following delimiter
            getline(fin,temp); _psf.back()(nCoeffs-1) = temp;
        }
    }
}

void PsfCatalog::read()
{
    std::string file = makeName(_params,"psf",false,true);
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

    if (!doesFileExist(file)) {
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
        throw ReadException(
            "Error reading from "+file+" -- caught error\n" + e.message());
    } catch (std::exception& e) {
        throw ReadException(
            "Error reading from "+file+" -- caught error\n" + e.what());
    } catch (...) {
        throw ReadException(
            "Error reading from "+file+" -- caught unknown error");
    }
    dbg<<"Done Read PSFCatalog\n";
}

