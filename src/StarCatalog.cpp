
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <valarray>
#include <CCfits/CCfits>

#include "StarCatalog.h"
#include "StarFinder.h"
#include "Bounds.h"
#include "ConfigFile.h"
#include "Image.h"
#include "Name.h"
#include "Transformation.h"
#include "Pixel.h"
#include "Params.h"
#include "Ellipse.h"
#include "Log.h"
#include "Form.h"
#include "WlVersion.h"
#include "WriteParam.h"

static void calculateSigma1(
    double& sigma,
    const Image<double>& im, const Position& pos, double sky,
    double noise, double gain, const Image<double>* weightIm, 
    const Transformation& trans, double psfAp, double xOffset, double yOffset,
    long& flag, bool shouldUseShapeletSigma)
{
    std::vector<PixelList> pix(1);
    long flag1 = 0;
    try {
        getPixList(im, pix[0], pos, sky, noise, gain, weightIm, trans, 
                   psfAp, xOffset, yOffset, flag1);
    } catch (RangeException) {
        flag1 |= TRANSFORM_EXCEPTION;
    }
    if (flag1) {
        flag |= flag1;
        sigma = DEFVALNEG;
        return;
    }

    Ellipse ell;
    ell.fixGam();
    ell.peakCentroid(pix[0],psfAp/3.);
    ell.crudeMeasure(pix[0],sigma);
    xdbg<<"Crude Measure: centroid = "<<ell.getCen();
    xdbg<<", mu = "<<ell.getMu()<<std::endl;
    if (shouldUseShapeletSigma) {
        if (ell.measure(pix,2,sigma,true,flag1)) { // true means use integ first
            xdbg<<"Successful 2nd order measure.\n";
            xdbg<<"mu = "<<ell.getMu()<<std::endl;
        } else {
            flag |= NATIVE_FAILED;
            xdbg<<"Ellipse measure returned flag "<<flag1<<std::endl;
            sigma = DEFVALNEG;
            return;
        }
    }

    double mu = ell.getMu();
    sigma *= exp(mu);
    dbg<<"sigma = "<<sigma<<std::endl;
    Assert(sigma > 0);
}

void calculateSigma(
    double& sigma,
    const Image<double>& im, const Position& pos, double sky,
    double noise, double gain, const Image<double>* weightIm, 
    const Transformation& trans, double psfAp, double xOffset, double yOffset,
    long& flag, bool shouldUseShapeletSigma)
{
    try {
        calculateSigma1(
            sigma,
            im, pos, sky, noise, gain, weightIm,
            trans, psfAp, xOffset, yOffset, flag, shouldUseShapeletSigma);
        dbg<<"objsize: "<<sigma<<std::endl;
        dbg<<"flags: "<<flag<<std::endl;
#ifdef USE_TMV
    } catch (tmv::Error& e) {
        dbg<<"Caught: "<<e<<std::endl;
        sigma = DEFVALNEG;
        flag |= TMV_EXCEPTION;
#endif
    } catch (std::exception& e) {
        dbg<<"Caught: "<<e.what()<<std::endl;
        sigma = DEFVALNEG;
        flag |= STD_EXCEPTION;
    } catch (...) {
        dbg<<"Caught unknown exception"<<std::endl;
        sigma = DEFVALNEG;
        flag |= UNKNOWN_EXCEPTION;
    }
}

StarCatalog::StarCatalog(
        const StarCatalog& inStarCat) :
    _id(inStarCat.getIdList()), 
    _pos(inStarCat.getPosList()), 
    _sky(inStarCat.getSkyList()), 
    _noise(inStarCat.getNoiseList()),
    _flags(inStarCat.getFlagsList()), 
    _mag(inStarCat.getMagList()), 
    _objSize(inStarCat.getObjSizeList()),
    _isAStar(inStarCat.getIsAStarList()),
    _params(inStarCat.getParams())
{
  Assert(_id.size() == size());
  Assert(_pos.size() == size());
  Assert(_sky.size() == size());
  Assert(_noise.size() == size());
  Assert(_flags.size() == size());
  Assert(_mag.size() == size());
  Assert(_objSize.size() == size());
  Assert(_isAStar.size() == size());


}


StarCatalog::StarCatalog(
    const StarCatalog& inStarCat, const std::vector<long> indices) :
    _params(inStarCat.getParams())
{

  size_t nind = indices.size();
  _id.resize(nind);
  _pos.resize(nind);
  _sky.resize(nind);
  _noise.resize(nind);
  _flags.resize(nind);

  _mag.resize(nind);
  _objSize.resize(nind);
  _isAStar.resize(nind);

  for (size_t i=0; i<nind; i++) {
    long ind = indices[i];
    _id[i]      = inStarCat.getId(ind);
    _pos[i]     = inStarCat.getPos(ind);
    _sky[i]     = inStarCat.getSky(ind);
    _noise[i]   = inStarCat.getNoise(ind);
    _flags[i]   = inStarCat.getFlags(ind);

    _mag[i]     = inStarCat.getMag(ind);
    _objSize[i] = inStarCat.getObjSize(ind);
    _isAStar[i] = inStarCat.getIsAStar(ind);
  }

}



StarCatalog::StarCatalog(
    const InputCatalog& inCat,
    const ConfigFile& params, std::string fsPrefix) :
    _id(inCat.getIdList()), _pos(inCat.getPosList()), 
    _sky(inCat.getSkyList()), _noise(inCat.getNoiseList()),
    _flags(inCat.getFlagsList()), _mag(inCat.getMagList()), 
    _objSize(inCat.getObjSizeList()),
    _isAStar(_id.size(),0), _params(params), _prefix(fsPrefix)
{
    Assert(_id.size() == size());
    Assert(_pos.size() == size());
    Assert(_sky.size() == size());
    Assert(_noise.size() == size());
    Assert(_flags.size() == size());
    Assert(_mag.size() == size());
    if (_objSize.size() == 0) _objSize.resize(size(),DEFVALNEG);
    Assert(_objSize.size() == size());
    Assert(_isAStar.size() == size());
}

StarCatalog::StarCatalog(const ConfigFile& params, std::string fsPrefix) :
    _params(params), _prefix(fsPrefix)
{ 
    Assert(_id.size() == size());
    Assert(_pos.size() == size());
    Assert(_sky.size() == size());
    Assert(_noise.size() == size());
    Assert(_flags.size() == size());
    Assert(_objSize.size() == size());
    Assert(_mag.size() == size());
    Assert(_isAStar.size() == size());
}


void StarCatalog::splitInTwo(
    const std::string f1, const std::string f2) const
{

    // use a deterministic seed: number of objects in catalog plus number
    // that are stars
    unsigned int nstar=0;
    for (unsigned int i=0;i<this->size();i++) {
        if (this->isAStar(i)) {
            nstar++;
        }
    }
    unsigned int seed = this->size() + nstar*100000;
    this->splitInTwo(f1,f2,seed);
}

void StarCatalog::splitInTwo(
    const std::string f1, const std::string f2,
    const unsigned int seed) const
{
    srand(seed);

    std::vector<long> id1;
    std::vector<long> id2;

    double split_point = 0.5;


    for (size_t i=0; i<this->size(); i++) {

        // number between [0,1)
        double r = ( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        if (r > split_point) {
            id1.push_back(i);
        } else {
            id2.push_back(i);
        }

    }

    StarCatalog cat1(*this, id1);
    StarCatalog cat2(*this, id2);

    cat1.writeFits(f1);
    cat2.writeFits(f2);
}

#if 0
std::string StarCatalog::getFitsName(std::string extra) 
{
    std::string file = makeName(_params, "stars",false,false);  

    size_t fits_pos = file.rfind(".fits");
    size_t dat_pos = file.rfind(".dat");

    if (fits_pos != std::string::npos) {
        if (extra != "") {
            file.insert(fits_pos, extra);
        }
    } else if (dat_pos != std::string::npos) {
        std::string e = extra+".fits";
        file.replace(dat_pos, 4, e);
    } else {

        // didn't find the expected ending .fits or .dat;  we will just tac it on
        // the end
        file += extra+".fits";

    }

    return file;
}
#endif

int StarCatalog::findStars(FindStarsLog& log)
{
    StarFinder finder(_params,_prefix);

    std::vector<PotentialStar*> maybestars;

    // First get a list of potential stars
    dbg<<"Finding stars"<<std::endl;
    long count=0;
    const int nObj = size();
    log._nTot = nObj;

    for (int i=0; i<nObj; ++i) {
        // A series of checks
        // Only objects with no flags in SExtractor or updated size calculation
        if (_flags[i]) {
            xdbg<<"Reject "<<i<<" for input flags: "<<_flags[i]<<"\n";
            ++log._nrFlag;
            continue;
        }
        // Range checking
        if (!finder.isOkSize(_objSize[i])) {
            xdbg<<"Reject "<<i<<" for size "<<_objSize[i]<<" outside range "<<
                finder.getMinSize()<<" -- "<<finder.getMaxSize()<<std::endl;
            ++log._nrSize;
            continue;
        }
        if (!finder.isOkMag(_mag[i])) {
            xdbg<<"Reject "<<i<<" for mag "<<_mag[i]<<" outside range "<<
                finder.getMinMag()<<" -- "<<finder.getMaxMag()<<std::endl;
            ++log._nrMag;
            continue;
        }
        xdbg<<"OK: "<<_objSize[i]<<"  "<<_mag[i]<<std::endl;

        ++count;
        double logObjSize = finder.convertToLogSize(_objSize[i]);
        maybestars.push_back(
            new PotentialStar(_pos[i],_mag[i],logObjSize,i,""));
    }
    log._nObj = count;
    dbg<<"  Possible Stars: "<<count<<"/"<<size()<<"\n";

    dbg<<"  Running FindStars\n";
    std::vector<PotentialStar*> stars = finder.findStars(maybestars);
    const int nStars = stars.size();
    dbg<<"  Found "<<nStars<<"\n";
    log._nAllStars = nStars;

    _isAStar.resize(size(),0);
    count = 0;
    for (int k=0; k<nStars;++k) {
        int i=stars[k]->getIndex();
        if (finder.isOkOutputMag(_mag[i])) {
            _isAStar[i] = 1;
            ++count;
        }
    }
    dbg<<"  Cut to "<<count<<" by maxoutmag cut\n";

    log._nStars = count;

    if (_params.read("des_qa",false)) {
        if (count < 100) {
            try {
                std::cout
                    << "STATUS3BEG Warning: Only "<<count
                    << " stars found for Name="
                    << makeName(_params,"stars",false,false)
                    << ". STATUS3END"<<std::endl;
            } catch (AssertFailureException& e) {
                std::cout
                    << "STATUS3BEG Warning: Only "<<count
                    << " stars found. STATUS3END"<<std::endl;
            }

        }
    }

    return count;
}

void StarCatalog::writeFits(std::string file) const 
{
    Assert(_id.size() == size());
    Assert(_pos.size() == size());
    Assert(_sky.size() == size());
    Assert(_noise.size() == size());
    Assert(_flags.size() == size());
    Assert(_mag.size() == size());
    Assert(_objSize.size() == size());
    Assert(_isAStar.size() == size());

    const int nTot = size();

    // ! means overwrite existing file
    CCfits::FITS fits("!"+file, CCfits::Write);

    const int nFields = 9;

    std::vector<string> colNames(nFields);
    std::vector<string> colFmts(nFields);
    std::vector<string> colUnits(nFields);

    colNames[0] = _params["stars_id_col"];
    colNames[1] = _params["stars_x_col"];
    colNames[2] = _params["stars_y_col"];
    colNames[3] = _params["stars_sky_col"];
    colNames[4] = _params["stars_noise_col"];
    colNames[5] = _params["stars_flags_col"];
    colNames[6] = _params["stars_mag_col"];
    colNames[7] = _params["stars_objsize_col"];
    colNames[8] = _params["stars_isastar_col"];

    colFmts[0] = "1J"; // id
    colFmts[1] = "1D"; // x
    colFmts[2] = "1D"; // y
    colFmts[3] = "1D"; // sky
    colFmts[4] = "1D"; // noise
    colFmts[5] = "1J"; // flags
    colFmts[6] = "1E"; // mag
    colFmts[7] = "1D"; // sigma
    colFmts[8] = "1J"; // star flag

    colUnits[0] = "None";     // id
    colUnits[1] = "pixels";   // x
    colUnits[2] = "pixels";   // y
    colUnits[3] = "ADU";      // sky
    colUnits[4] = "ADU^2";    // noise
    colUnits[5] = "None";     // flags
    colUnits[6] = "mags";     // mag
    colUnits[7] = "Arcsec";   // sigma0
    colUnits[8] = "None";     //star flag

    CCfits::Table* table;
    table = fits.addTable("findstars",nTot,colNames,colFmts,colUnits);

    // Header Keywords
#ifdef USE_TMV
    std::string tmvVers = tmv::TMV_Version();
#else
    std::string tmvVers = "Eigen";
#endif
    std::string wlVers = getWlVersion();

    table->addKey("tmvvers", tmvVers, "version of TMV code");
    table->addKey("wlvers", wlVers, "version of weak lensing code");

    // kind of kludgy but is more flexible than using type numbers or strings
    double dbl;
    int intgr;
    std::string str;


    // if wlserun= is sent we'll put it in the header.  This allows us to 
    // associate some more, possibly complicated, metadata with this file
    if ( _params.keyExists("wlserun") ) {
        writeParamToTable(_params, table, "wlserun", str);
    }


    writeParamToTable(_params, table, "noise_method", str);
    writeParamToTable(_params, table, "dist_method", str);

    if (_params.keyExists("stars_minsize")) {
        writeParamToTable(_params, table, "stars_minsize", dbl);
        writeParamToTable(_params, table, "stars_maxsize", dbl);
        writeParamToTable(_params, table, "stars_minmag", dbl);
        writeParamToTable(_params, table, "stars_maxmag", dbl);
        writeParamToTable(_params, table, "stars_ndivx", intgr);
        writeParamToTable(_params, table, "stars_ndivy", intgr);

        writeParamToTable(_params, table, "stars_startn1", dbl);
        writeParamToTable(_params, table, "stars_starfrac", dbl);
        writeParamToTable(_params, table, "stars_magstep1", dbl);
        writeParamToTable(_params, table, "stars_miniter1", intgr);
        writeParamToTable(_params, table, "stars_reject1", dbl);
        writeParamToTable(_params, table, "stars_binsize1", dbl);
        writeParamToTable(_params, table, "stars_maxratio1", dbl);
        writeParamToTable(_params, table, "stars_okvalcount", intgr);
        writeParamToTable(_params, table, "stars_maxrms", dbl);
        writeParamToTable(_params, table, "stars_starsperbin", intgr);

        writeParamToTable(_params, table, "stars_fitorder", intgr);
        writeParamToTable(_params, table, "stars_fitsigclip", dbl);
        writeParamToTable(_params, table, "stars_startn2", dbl);
        writeParamToTable(_params, table, "stars_magstep2", dbl);
        writeParamToTable(_params, table, "stars_miniter2", intgr);
        writeParamToTable(_params, table, "stars_minbinsize", dbl);
        writeParamToTable(_params, table, "stars_reject2", dbl);

        writeParamToTable(_params, table, "stars_purityratio", dbl);
        writeParamToTable(_params, table, "stars_maxrefititer", intgr);
    }


    // Now the data columns

    // make vector copies for writing
    std::vector<double> x(nTot);
    std::vector<double> y(nTot);
    for(int i=0;i<nTot;++i) {
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
    table->column(colNames[6]).write(_mag,startrow);
    table->column(colNames[7]).write(_objSize,startrow);
    table->column(colNames[8]).write(_isAStar,startrow);

}

void StarCatalog::writeAscii(std::string file, std::string delim) const 
{
    Assert(_id.size() == size());
    Assert(_pos.size() == size());
    Assert(_sky.size() == size());
    Assert(_noise.size() == size());
    Assert(_flags.size() == size());
    Assert(_mag.size() == size());
    Assert(_objSize.size() == size());
    Assert(_isAStar.size() == size());

    std::ofstream fout(file.c_str());
    if (!fout) {
        throw WriteException("Error opening stars file"+file);
    }

    Form sci8; sci8.sci().trail(0).prec(8);
    Form fix3; fix3.fix().trail(0).prec(3);
    Form fix6; fix6.fix().trail(0).prec(6);

    const int nTot = size();
    for (int i=0; i<nTot; ++i) {
        fout 
            << _id[i] << delim
            << fix3(_pos[i].getX()) << delim
            << fix3(_pos[i].getY()) << delim
            << fix3(_sky[i]) << delim
            << sci8(_noise[i]) << delim
            << _flags[i] << delim
            << fix3(_mag[i]) << delim
            << fix6(_objSize[i]) << delim
            << _isAStar[i] << std::endl;
    }

}

void StarCatalog::write() const 
{
    std::vector<std::string> files = makeMultiName(_params, "stars");  

    const int nFiles = files.size();
    for(int i=0; i<nFiles; ++i) {
        const std::string& file = files[i];
        dbg<<"Writing star catalog to file: "<<file<<std::endl;

        bool isFitsIo = false;
        if (_params.keyExists("stars_io")) {
            std::vector<std::string> ios = _params["stars_io"];
            Assert(ios.size() == files.size());
            isFitsIo = (ios[i] == "FITS");
        } else if (file.find("fits") != std::string::npos) {
            isFitsIo = true;
        }

        try {
            if (isFitsIo) {
                writeFits(file);
            } else {
                std::string delim = "  ";
                if (_params.keyExists("stars_delim")) {
                    std::vector<std::string> delims = _params["stars_delim"];
                    Assert(delims.size() == files.size());
                    delim = delims[i];
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
    dbg<<"Done Write StarCatalog\n";
}

void StarCatalog::readFits(std::string file) 
{
    int hdu = getHdu(_params,"stars",file,2);

    dbg<<"Opening StarCatalog file "<<file<<" at hdu "<<hdu<<std::endl;
    CCfits::FITS fits(file, CCfits::Read);
    if (hdu > 1) fits.read(hdu-1);

    CCfits::ExtHDU& table=fits.extension(hdu-1);

    long nRows=table.rows();

    dbg<<"  nrows = "<<nRows<<std::endl;
    if (nRows <= 0) {
        throw ReadException(
            "StarCatalog found to have 0 rows.  Must have > 0 rows.");
    }


    std::string idCol=_params.get("stars_id_col");
    std::string xCol=_params.get("stars_x_col");
    std::string yCol=_params.get("stars_y_col");
    std::string skyCol=_params.get("stars_sky_col");
    std::string noiseCol=_params.get("stars_noise_col");
    std::string flagsCol=_params.get("stars_flags_col");
    std::string magCol=_params.get("stars_mag_col");
    std::string objSizeCol=_params.get("stars_objsize_col");
    std::string isAStarCol=_params.get("stars_isastar_col");

    long start=1;
    long end=nRows;

    dbg<<"Reading columns"<<std::endl;
    // will copy these out to positions types
    std::vector<double> x(nRows);
    std::vector<double> y(nRows);

    dbg<<"  "<<idCol<<std::endl;
    table.column(idCol).read(_id, start, end);
    dbg<<"  "<<xCol<<"  "<<yCol<<std::endl;
    table.column(xCol).read(x, start, end);
    table.column(yCol).read(y, start, end);
    dbg<<"  "<<skyCol<<std::endl;
    table.column(skyCol).read(_sky, start, end);
    dbg<<"  "<<noiseCol<<std::endl;
    table.column(noiseCol).read(_noise, start, end);
    dbg<<"  "<<flagsCol<<std::endl;
    table.column(flagsCol).read(_flags, start, end);
    dbg<<"  "<<magCol<<std::endl;
    table.column(magCol).read(_mag, start, end);
    dbg<<"  "<<objSizeCol<<std::endl;
    table.column(objSizeCol).read(_objSize, start, end);
    dbg<<"  "<<isAStarCol<<std::endl;
    table.column(isAStarCol).read(_isAStar, start, end);

    _pos.resize(nRows);
    for(long i=0;i<nRows;++i) {
        _pos[i] = Position(x[i],y[i]);
    }

}

void StarCatalog::readAscii(std::string file, std::string delim)
{
    std::ifstream fin(file.c_str());
    if (!fin) {
        throw ReadException("Error opening stars file"+file);
    }

    _id.clear(); _pos.clear(); _sky.clear(); _noise.clear(); _flags.clear();
    _mag.clear(); _objSize.clear(); _isAStar.clear();

    if (delim == "  ") {
        ConvertibleString flag;
        long id1,star;
        double x,y,sky1,n,s,m;
        while ( fin >> id1 >> x >> y >> sky1 >> n 
                >> flag >> m >> s >> star) {
            _id.push_back(id1);
            _pos.push_back(Position(x,y));
            _sky.push_back(sky1);
            _noise.push_back(n);
            _flags.push_back(flag);
            _mag.push_back(m);
            _objSize.push_back(s);
            _isAStar.push_back(star);
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
        dbg<<"Reading with delimeter "<<d<<std::endl;
        double x,y;
        ConvertibleString temp;
        while (getline(fin,temp,d)) {
            dbg<<"First elemnt in line = "<<temp<<std::endl;
            _id.push_back(temp);
            getline(fin,temp,d); x = temp;
            getline(fin,temp,d); y = temp;
            _pos.push_back(Position(x,y));
            getline(fin,temp,d); _sky.push_back(temp);
            getline(fin,temp,d); _noise.push_back(temp);
            getline(fin,temp,d); _flags.push_back(temp);
            getline(fin,temp,d); _mag.push_back(temp);
            getline(fin,temp,d); _objSize.push_back(temp);
            getline(fin,temp); _isAStar.push_back(temp);
            dbg<<"Last elemnt in line = "<<temp<<std::endl;
        }
        dbg<<"nlines = "<<size()<<std::endl;
    }
}

void StarCatalog::read()
{
    std::string file = makeName(_params,"stars",false,true);
    // false,true = input_prefix=false, mustexist=true.
    // It is an input here, but it is in the output_prefix directory.
    dbg<< "Reading star catalog from file: " << file << std::endl;

    bool isFitsIo = false;
    if (_params.keyExists("stars_io")) {
        isFitsIo = (_params["stars_io"] == "FITS");
    } else if (file.find("fits") != std::string::npos) {
        isFitsIo = true;
    }

    if (!doesFileExist(file)) {
        throw FileNotFoundException(file);
    }
    try {
        if (isFitsIo) {
            readFits(file);
        } else {
            std::string delim = "  ";
            if (_params.keyExists("stars_delim")) {
                delim = _params["stars_delim"];
            } else if (file.find("csv") != std::string::npos) {
                delim = ",";
            }
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
    dbg<<"Done Read StarCatalog\n";
}


bool StarCatalog::isAStar(long id) const
{
    std::vector<long>::const_iterator p = find(_id.begin(),_id.end(),id);
    if (p == _id.end()) return false;
    int i = p - _id.begin();
    Assert(i < _isAStar.size());
    return _isAStar[i];
}

void StarCatalog::printall(int i) 
{
    if (int(_id.size()) > i) {
        std::cout<<"  StarCatalog::printall id["<<i<<"]: "<<_id[i]<<"\n";
    } else {
        std::cout<<"  StarCatalog::printall index "<<i<<" is larger than id.size\n";
    }

    if (int(_pos.size()) > i) {
        std::cout<<"  StarCatalog::printall pos["<<i<<"]: "<<_pos[i]<<"\n";
    } else {
        std::cout<<"  StarCatalog::printall index "<<i<<" is larger than pos.size\n";
    }

    if (int(_sky.size()) > i) {
        std::cout<<"  StarCatalog::printall sky["<<i<<"]: "<<_sky[i]<<"\n";
    } else {
        std::cout<<"  StarCatalog::printall index "<<i<<" is larger than sky.size\n";
    }

    if (int(_mag.size()) > i) {
        std::cout<<"  StarCatalog::printall mag["<<i<<"]: "<<_mag[i]<<"\n";
    } else {
        std::cout<<"  StarCatalog::printall index "<<i<<" is larger than mag.size\n";
    }

    if (int(_objSize.size()) > i) {
        std::cout<<"  StarCatalog::printall objSize["<<i<<"]: "<<_objSize[i]<<"\n";
    } else {
        std::cout<<"  StarCatalog::printall index "<<i<<" is larger than objSize.size\n";
    }

    if (int(_flags.size()) > i) {
        std::cout<<"  StarCatalog::printall flags["<<i<<"]: "<<this->_flags[i]<<"\n";
    } else {
        std::cout<<"  StarCatalog::printall index "<<i<<" is larger than flags.size\n";
    }

    if (int(_noise.size()) > i) {
        std::cout<<"  StarCatalog::printall noise["<<i<<"]: "<<this->_noise[i]<<"\n";
    } else {
        std::cout<<"  StarCatalog::printall index "<<i<<" is larger than noise.size\n";
    }
    if (int(_isAStar.size()) > i) {
        std::cout<<"  StarCatalog::printall isAStar["<<i<<"]: "<<this->_isAStar[i]<<"\n";
    } else {
        std::cout<<"  StarCatalog::printall index "<<i<<" is larger than isAStar.size\n";
    }


}
