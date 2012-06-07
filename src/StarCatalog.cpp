
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
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

static void CalculateSigma1(
    double& sigma,
    const Image<double>& im, const Position& pos, double sky,
    double noise, const Image<double>* weight_image, 
    const Transformation& trans, const ConfigFile& params, long& flag,
    bool use_shapelet_sigma)
{
    double psf_ap = params.read<double>("psf_aperture");

    std::vector<PixelList> pix(1);
    long flag1 = 0;
    try {
        GetPixList(im, pix[0], pos, sky, noise, weight_image, trans, 
                   psf_ap, params, flag1);
    } catch (RangeException& e) {
        dbg<<"distortion range error: \n";
        xdbg<<"center = "<<pos<<", b = "<<e.getBounds()<<std::endl;
        dbg<<"FLAG TRANSFORM_EXCEPTION\n";
        flag1 |= TRANSFORM_EXCEPTION;
    }
    if (flag1) {
        dbg<<"getPixList returned flag "<<FlagText(flag1)<<std::endl;
        flag |= flag1;
        sigma = DEFVALNEG;
        return;
    }

    Ellipse ell;
    ell.fixGam();
    ell.peakCentroid(pix[0],psf_ap/3.);
    ell.crudeMeasure(pix[0],sigma);
    xdbg<<"Crude Measure: centroid = "<<ell.getCen();
    xdbg<<", mu = "<<ell.getMu()<<std::endl;
    if (use_shapelet_sigma) {
        if (ell.measure(pix,2,6,2,sigma,flag1,0.01)) {
            xdbg<<"Successful 2nd order measure.\n";
            xdbg<<"mu = "<<ell.getMu()<<std::endl;
        } else {
            dbg<<"FLAG NATIVE_FAILED\n";
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

void CalculateSigma(
    double& sigma,
    const Image<double>& im, const Position& pos, double sky,
    double noise, const Image<double>* weight_image, 
    const Transformation& trans, const ConfigFile& params,
    long& flag, bool use_shapelet_sigma)
{
    try {
        CalculateSigma1(
            sigma,
            im, pos, sky, noise, weight_image,
            trans, params, flag, use_shapelet_sigma);
        dbg<<"objsize: "<<sigma<<std::endl;
        dbg<<"flags: "<<flag<<std::endl;
#ifdef USE_TMV
    } catch (tmv::Error& e) {
        dbg<<"Caught: "<<e<<std::endl;
        sigma = DEFVALNEG;
        dbg<<"FLAG TMV_EXCEPTION\n";
        flag |= TMV_EXCEPTION;
#endif
    } catch (std::exception& e) {
        dbg<<"Caught: "<<e.what()<<std::endl;
        sigma = DEFVALNEG;
        dbg<<"FLAG STD_EXCEPTION\n";
        flag |= STD_EXCEPTION;
    } catch (...) {
        dbg<<"Caught unknown exception"<<std::endl;
        sigma = DEFVALNEG;
        dbg<<"FLAG UNKNOWN_EXCEPTION\n";
        flag |= UNKNOWN_EXCEPTION;
    }
}

StarCatalog::StarCatalog(const StarCatalog& instarcat) :
    _id(instarcat.getIdList()), 
    _pos(instarcat.getPosList()), 
    _sky(instarcat.getSkyList()), 
    _noise(instarcat.getNoiseList()),
    _flags(instarcat.getFlagsList()), 
    _mag(instarcat.getMagList()), 
    _sg(instarcat.getSgList()),
    _objsize(instarcat.getObjSizeList()),
    _is_star(instarcat.getIsStarList()),
    _params(instarcat.getParams())
{
    Assert(int(_id.size()) == size());
    Assert(int(_pos.size()) == size());
    Assert(int(_sky.size()) == size());
    Assert(int(_noise.size()) == size());
    Assert(int(_flags.size()) == size());
    Assert(int(_mag.size()) == size());
    Assert(int(_sg.size()) == size());
    Assert(int(_objsize.size()) == size());
    Assert(int(_is_star.size()) == size());


}


StarCatalog::StarCatalog(
    const StarCatalog& instarcat, const std::vector<long> indices) :
    _params(instarcat.getParams())
{
    int nind = indices.size();
    _id.resize(nind);
    _pos.resize(nind);
    _sky.resize(nind);
    _noise.resize(nind);
    _flags.resize(nind);

    _mag.resize(nind);
    _sg.resize(nind);
    _objsize.resize(nind);
    _is_star.resize(nind);

    for (int i=0; i<nind; i++) {
        long ind = indices[i];
        _id[i]      = instarcat.getId(ind);
        _pos[i]     = instarcat.getPos(ind);
        _sky[i]     = instarcat.getSky(ind);
        _noise[i]   = instarcat.getNoise(ind);
        _flags[i]   = instarcat.getFlags(ind);

        _mag[i]     = instarcat.getMag(ind);
        _sg[i]      = instarcat.getSg(ind);
        _objsize[i] = instarcat.getObjSize(ind);
        _is_star[i] = instarcat.getIsStar(ind);
    }
}


StarCatalog::StarCatalog(
    const InputCatalog& incat,
    const ConfigFile& params, std::string fs_prefix) :
    _id(incat.getIdList()), _pos(incat.getPosList()), 
    _sky(incat.getSkyList()), _noise(incat.getNoiseList()),
    _flags(incat.getFlagsList()), _mag(incat.getMagList()), 
    _sg(incat.getSgList()), _objsize(incat.getObjSizeList()),
    _is_star(_id.size(),0), _params(params), _prefix(fs_prefix)
{
    Assert(int(_id.size()) == size());
    Assert(int(_pos.size()) == size());
    Assert(int(_sky.size()) == size());
    Assert(int(_noise.size()) == size());
    Assert(int(_flags.size()) == size());
    if (_sg.size() == 0) _sg.resize(size(),DEFVALNEG);
    Assert(int(_sg.size()) == size());
    if (_objsize.size() == 0) _objsize.resize(size(),DEFVALNEG);
    Assert(int(_objsize.size()) == size());
    Assert(int(_is_star.size()) == size());
    if (params.read("cat_all_stars",false)) {
        for(int i=0;i<size();++i) _is_star[i] = 1;
    } else if (params.read("stars_trust_sg",false)) {
        double minsg = params.get("stars_minsg");
        double maxsg = params.get("stars_maxsg");
        Assert(int(_sg.size()) == size());
        for(int i=0;i<size();++i) 
            _is_star[i] = (_sg[i] >= minsg && _sg[i] <= maxsg);
    } else {
        Assert(int(_mag.size()) == size());
    }
}

StarCatalog::StarCatalog(const ConfigFile& params, std::string fs_prefix) :
    _params(params), _prefix(fs_prefix)
{ 
    Assert(int(_id.size()) == size());
    Assert(int(_pos.size()) == size());
    Assert(int(_sky.size()) == size());
    Assert(int(_noise.size()) == size());
    Assert(int(_flags.size()) == size());
    Assert(int(_objsize.size()) == size());
    Assert(int(_mag.size()) == size());
    Assert(int(_sg.size()) == size());
    Assert(int(_is_star.size()) == size());
}


void StarCatalog::splitInTwo(
    const std::string f1, const std::string f2) const
{

    // use a deterministic seed: number of objects in catalog plus number
    // that are stars
    int nstar=0;
    for (int i=0;i<this->size();i++) {
        if (this->isStar(i)) {
            nstar++;
        }
    }
    int seed = this->size() + nstar*100000;
    this->splitInTwo(f1,f2,seed);
}

void StarCatalog::splitInTwo(
    const std::string f1, const std::string f2,
    const int seed) const
{
    srand(seed);

    std::vector<long> id1;
    std::vector<long> id2;

    double split_point = 0.5;


    for (int i=0; i<this->size(); i++) {

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


int StarCatalog::findStars(FindStarsLog& log)
{
    StarFinder finder(_params,_prefix);
    std::vector<PotentialStar*> maybestars;

    // First get a list of potential stars
    dbg<<"Finding stars"<<std::endl;
    long count=0;
    const int nobj = size();
    log._ntot = nobj;

    for (int i=0; i<nobj; ++i) {
        // A series of checks
        // Only objects with no flags in SExtractor or updated size calculation
        if (_flags[i]) {
            xdbg<<"Reject "<<i<<" for input flags: "<<_flags[i]<<"\n";
            ++log._nr_flag;
            continue;
        }
        // Range checking
        if (!finder.isOkSize(_objsize[i])) {
            xdbg<<"Reject "<<i<<" for size "<<_objsize[i]<<" outside range "<<
                finder.getMinSize()<<" -- "<<finder.getMaxSize()<<std::endl;
            ++log._nr_size;
            continue;
        }
        if (!finder.isOkMag(_mag[i])) {
            xdbg<<"Reject "<<i<<" for mag "<<_mag[i]<<" outside range "<<
                finder.getMinMag()<<" -- "<<finder.getMaxMag()<<std::endl;
            ++log._nr_mag;
            continue;
        }
        // Count how many are initially selected by star-galaxy cut
        if (finder.isOkSg(_sg[i]) && finder.isOkSgMag(_mag[i])) {
          ++log._nsg;
        }
        
        ++count;
        double logObjSize = finder.convertToLogSize(_objsize[i]);
        maybestars.push_back(
            new PotentialStar(_pos[i],_mag[i],logObjSize,_sg[i],i,""));
    }
    log._nobj = count;
    dbg<<"  Possible Stars: "<<count<<"/"<<size()<<"\n";

    dbg<<"  Running FindStars\n";
    std::vector<PotentialStar*> stars = finder.findStars(maybestars);
    const int nstars = stars.size();
    dbg<<"  Found "<<nstars<<"\n";
    log._nallstars = nstars;

    _is_star.resize(size(),0);
    count = 0;

    for (int k=0; k<nstars;++k) {
        int i=stars[k]->getIndex();
        if (finder.isOkOutputMag(_mag[i])) {
            _is_star[i] = 1;
            ++count;
        }
    }
    dbg<<"  Cut to "<<count<<" by maxoutmag cut\n";

    log._nstars = count;

    if (_params.read("des_qa",false)) {
        if (count < 100) {
            try {
                std::string name = MakeName(_params,"stars",false,false);
                std::cout
                    << "STATUS3BEG Warning: Only "<<count
                    << " stars found for Name="<<name
                    << ". STATUS3END"<<std::endl;
            } catch (AssertFailureException& ) {
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
    Assert(int(_id.size()) == size());
    Assert(int(_pos.size()) == size());
    Assert(int(_sky.size()) == size());
    Assert(int(_noise.size()) == size());
    Assert(int(_flags.size()) == size());
    Assert(int(_mag.size()) == size());
    Assert(int(_sg.size()) == size());
    Assert(int(_objsize.size()) == size());
    Assert(int(_is_star.size()) == size());

    const int ntot = size();

    // ! means overwrite existing file
    CCfits::FITS fits("!"+file, CCfits::Write);

    const int nFields = 10;

    std::vector<string> col_names(nFields);
    std::vector<string> col_fmts(nFields);
    std::vector<string> col_units(nFields);

    col_names[0] = _params["stars_id_col"];
    col_names[1] = _params["stars_x_col"];
    col_names[2] = _params["stars_y_col"];
    col_names[3] = _params["stars_sky_col"];
    col_names[4] = _params["stars_noise_col"];
    col_names[5] = _params["stars_flags_col"];
    col_names[6] = _params["stars_mag_col"];
    col_names[7] = _params["stars_sg_col"];
    col_names[8] = _params["stars_objsize_col"];
    col_names[9] = _params["stars_isastar_col"];

    col_fmts[0] = "1J"; // id
    col_fmts[1] = "1D"; // x
    col_fmts[2] = "1D"; // y
    col_fmts[3] = "1D"; // sky
    col_fmts[4] = "1D"; // noise
    col_fmts[5] = "1J"; // flags
    col_fmts[6] = "1E"; // mag
    col_fmts[7] = "1E"; // star-galaxy
    col_fmts[8] = "1D"; // sigma
    col_fmts[9] = "1J"; // star flag



    col_units[0] = "None";     // id
    col_units[1] = "pixels";   // x
    col_units[2] = "pixels";   // y
    col_units[3] = "ADU";      // sky
    col_units[4] = "ADU^2";    // noise
    col_units[5] = "None";     // flags
    col_units[6] = "mags";     // mag
    col_units[7] = "None";     // star-galaxy
    col_units[8] = "Arcsec";   // sigma0
    col_units[9] = "None";     // star flag


    CCfits::Table* table;
    table = fits.addTable("findstars",ntot,col_names,col_fmts,col_units);

    // Header Keywords
#ifdef USE_TMV
    std::string tmvVers = tmv::TMV_Version();
#else
    std::string tmvVers = "Eigen";
#endif
    std::string wlVers = GetWlVersion();

    table->addKey("tmvvers", tmvVers, "version of TMV code");
    table->addKey("wlvers", wlVers, "version of weak lensing code");

    // kind of kludgy but is more flexible than using type numbers or strings
    double dbl;
    int intgr;
    std::string str;


    // if wlserun= is sent we'll put it in the header.  This allows us to 
    // associate some more, possibly complicated, metadata with this file
    if ( _params.keyExists("wlserun") ) {
        WriteParamToTable(_params, table, "wlserun", str);
    }


    WriteParamToTable(_params, table, "noise_method", str);
    WriteParamToTable(_params, table, "dist_method", str);

    if (_params.keyExists("stars_minsize")) {
        WriteParamToTable(_params, table, "stars_minsize", dbl);
        WriteParamToTable(_params, table, "stars_maxsize", dbl);
        WriteParamToTable(_params, table, "stars_minmag", dbl);
        WriteParamToTable(_params, table, "stars_maxmag", dbl);
        WriteParamToTable(_params, table, "stars_ndivx", intgr);
        WriteParamToTable(_params, table, "stars_ndivy", intgr);

        WriteParamToTable(_params, table, "stars_startn1", dbl);
        WriteParamToTable(_params, table, "stars_starfrac", dbl);
        WriteParamToTable(_params, table, "stars_magstep1", dbl);
        WriteParamToTable(_params, table, "stars_miniter1", intgr);
        WriteParamToTable(_params, table, "stars_reject1", dbl);
        WriteParamToTable(_params, table, "stars_binsize1", dbl);
        WriteParamToTable(_params, table, "stars_maxratio1", dbl);
        WriteParamToTable(_params, table, "stars_okvalcount", intgr);
        WriteParamToTable(_params, table, "stars_maxrms", dbl);
        WriteParamToTable(_params, table, "stars_starsperbin", intgr);

        WriteParamToTable(_params, table, "stars_fitorder", intgr);
        WriteParamToTable(_params, table, "stars_fitsigclip", dbl);
        WriteParamToTable(_params, table, "stars_startn2", dbl);
        WriteParamToTable(_params, table, "stars_magstep2", dbl);
        WriteParamToTable(_params, table, "stars_miniter2", intgr);
        WriteParamToTable(_params, table, "stars_minbinsize", dbl);
        WriteParamToTable(_params, table, "stars_reject2", dbl);

        WriteParamToTable(_params, table, "stars_purityratio", dbl);
        WriteParamToTable(_params, table, "stars_maxrefititer", intgr);
        WriteParamToTable(_params, table, "stars_minsg", dbl);
        WriteParamToTable(_params, table, "stars_maxsg", dbl);
        WriteParamToTable(_params, table, "stars_minsgmag", dbl);
        WriteParamToTable(_params, table, "stars_maxsgmag", dbl);
    }
  

    // Now the data columns

    // make vector copies for writing
    std::vector<double> x(ntot);
    std::vector<double> y(ntot);
    for(int i=0;i<ntot;++i) {
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
    table->column(col_names[6]).write(_mag,startrow);
    table->column(col_names[7]).write(_sg,startrow);
    table->column(col_names[8]).write(_objsize,startrow);
    table->column(col_names[9]).write(_is_star,startrow);


}

void StarCatalog::writeAscii(std::string file, std::string delim) const 
{
    Assert(int(_id.size()) == size());
    Assert(int(_pos.size()) == size());
    Assert(int(_sky.size()) == size());
    Assert(int(_noise.size()) == size());
    Assert(int(_flags.size()) == size());
    Assert(int(_mag.size()) == size());
    Assert(int(_sg.size()) == size());
    Assert(int(_objsize.size()) == size());
    Assert(int(_is_star.size()) == size());

    std::ofstream fout(file.c_str());
    if (!fout) {
        throw WriteException("Error opening stars file"+file);
    }

    Form sci8; sci8.sci().trail(0).prec(8);
    Form fix3; fix3.fix().trail(0).prec(3);
    Form fix6; fix6.fix().trail(0).prec(6);

    const int ntot = size();
    for (int i=0; i<ntot; ++i) {
        fout 
            << _id[i] << delim
            << fix3(_pos[i].getX()) << delim
            << fix3(_pos[i].getY()) << delim
            << fix3(_sky[i]) << delim
            << sci8(_noise[i]) << delim
            << _flags[i] << delim
            << fix3(_mag[i]) << delim
            << fix3(_sg[i]) << delim
            << fix6(_objsize[i]) << delim
            << _is_star[i] << std::endl;
    }

}

void StarCatalog::write() const 
{
    std::vector<std::string> files = MakeMultiName(_params, "stars");  

    const int nfiles = files.size();
    for(int i=0; i<nfiles; ++i) {
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
    dbg<<"Done Write StarCatalog\n";
}

void StarCatalog::readFits(std::string file) 
{
    int hdu = GetHdu(_params,"stars",file,2);

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


    std::string id_col=_params.get("stars_id_col");
    std::string x_col=_params.get("stars_x_col");
    std::string y_col=_params.get("stars_y_col");
    std::string sky_col=_params.get("stars_sky_col");
    std::string noise_col=_params.get("stars_noise_col");
    std::string flags_col=_params.get("stars_flags_col");
    std::string mag_col=_params.get("stars_mag_col");
    std::string sg_col=_params.get("stars_sg_col");
    std::string objsize_col=_params.get("stars_objsize_col");
    std::string is_star_col=_params.get("stars_isastar_col");

    long start=1;
    long end=nRows;

    dbg<<"Reading columns"<<std::endl;
    // will copy these out to positions types
    std::vector<double> x(nRows);
    std::vector<double> y(nRows);

    dbg<<"  "<<id_col<<std::endl;
    table.column(id_col).read(_id, start, end);
    dbg<<"  "<<x_col<<"  "<<y_col<<std::endl;
    table.column(x_col).read(x, start, end);
    table.column(y_col).read(y, start, end);
    dbg<<"  "<<sky_col<<std::endl;
    table.column(sky_col).read(_sky, start, end);
    dbg<<"  "<<noise_col<<std::endl;
    table.column(noise_col).read(_noise, start, end);
    dbg<<"  "<<flags_col<<std::endl;
    table.column(flags_col).read(_flags, start, end);
    dbg<<"  "<<mag_col<<std::endl;
    table.column(mag_col).read(_mag, start, end);
    dbg<<"  "<<sg_col<<std::endl;
    table.column(sg_col).read(_sg, start, end);
    dbg<<"  "<<objsize_col<<std::endl;
    table.column(objsize_col).read(_objsize, start, end);
    dbg<<"  "<<is_star_col<<std::endl;
    table.column(is_star_col).read(_is_star, start, end);
    
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
    _mag.clear(); _sg.clear(); _objsize.clear(); _is_star.clear();

    if (delim == "  ") {
        ConvertibleString flag;
        long id1,star;
        double x,y,sky1,n,s,m,sg;
        while ( fin >> id1 >> x >> y >> sky1 >> n 
                >> flag >> m >> sg >> s >> star) {
            _id.push_back(id1);
            _pos.push_back(Position(x,y));
            _sky.push_back(sky1);
            _noise.push_back(n);
            _flags.push_back(flag);
            _mag.push_back(m);
            _sg.push_back(sg);
            _objsize.push_back(s);
            _is_star.push_back(star);
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
            getline(fin,temp,d); _sg.push_back(temp);
            getline(fin,temp,d); _objsize.push_back(temp);
            getline(fin,temp); _is_star.push_back(temp);
            dbg<<"Last elemnt in line = "<<temp<<std::endl;
        }
        dbg<<"nlines = "<<size()<<std::endl;
    }
}

void StarCatalog::read()
{
    std::string file = MakeName(_params,"stars",false,true);
    // false,true = input_prefix=false, mustexist=true.
    // It is an input here, but it is in the output_prefix directory.
    dbg<< "Reading star catalog from file: " << file << std::endl;

    bool isFitsIo = false;
    if (_params.keyExists("stars_io")) {
        isFitsIo = (_params["stars_io"] == "FITS");
    } else if (file.find("fits") != std::string::npos) {
        isFitsIo = true;
    }

    if (!DoesFileExist(file)) {
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
    dbg<<"Done Read StarCatalog\n";
}


bool StarCatalog::isStar(long id) const
{
    std::vector<long>::const_iterator p = find(_id.begin(),_id.end(),id);
    if (p == _id.end()) return false;
    int i = p - _id.begin();
    Assert(i < int(_is_star.size()));
    return _is_star[i];
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

    if (int(_sg.size()) > i) {
        std::cout<<"  StarCatalog::printall sg["<<i<<"]: "<<_sg[i]<<"\n";
    } else {
        std::cout<<"  StarCatalog::printall index "<<i<<" is larger than sg.size\n";
    }

    if (int(_objsize.size()) > i) {
        std::cout<<"  StarCatalog::printall objsize["<<i<<"]: "<<_objsize[i]<<"\n";
    } else {
        std::cout<<"  StarCatalog::printall index "<<i<<" is larger than objsize.size\n";
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
    if (int(_is_star.size()) > i) {
        std::cout<<"  StarCatalog::printall is_star["<<i<<"]: "<<this->_is_star[i]<<"\n";
    } else {
        std::cout<<"  StarCatalog::printall index "<<i<<" is larger than is_star.size\n";
    }


}
