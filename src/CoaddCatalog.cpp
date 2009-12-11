
#include <valarray>
#include <CCfits/CCfits>

#include "CoaddCatalog.h"
#include "dbg.h"
#include "Params.h"
#include "Name.h"

CoaddCatalog::CoaddCatalog(ConfigFile& params) :
    _params(params)
{
    readCatalog();

    Assert(_pos.size() == _id.size());
    Assert(_skyPos.size() == _id.size());
    Assert(_sky.size() == _id.size());
    Assert(_mag.size() == _id.size());
    Assert(_magErr.size() == _id.size());
    Assert(_flags.size() == _id.size());

    std::vector<int> flagCount(9,0);
    std::vector<int> magCount(10,0);
    const int nObj = size();
    for (int i=0; i<nObj; ++i) {
        if (_flags[i] & 1) ++flagCount[0]; // neighbor or badpix
        if (_flags[i] & 2) ++flagCount[1]; // deblended
        if (_flags[i] & 4) ++flagCount[2]; // saturated
        if (_flags[i] & 8) ++flagCount[3]; // edge
        if (_flags[i] & 16) ++flagCount[4]; // incomplete aperture
        if (_flags[i] & 32) ++flagCount[5]; // incomplete isophotal
        if (_flags[i] & 64) ++flagCount[6]; // memory overflow in deblend
        if (_flags[i] & 128) ++flagCount[7]; // memory overflow in extract
        if (!_flags[i]) ++flagCount[8]; // no flags

        if (_mag[i] < 19) ++magCount[0];
        else if (_mag[i] < 20) ++magCount[1];
        else if (_mag[i] < 21) ++magCount[2];
        else if (_mag[i] < 22) ++magCount[3];
        else if (_mag[i] < 23) ++magCount[4];
        else if (_mag[i] < 24) ++magCount[5];
        else if (_mag[i] < 25) ++magCount[6];
        else if (_mag[i] < 26) ++magCount[7];
        else if (_mag[i] < 27) ++magCount[8];
        else ++magCount[9];
    }

    bool shouldOutputDots = _params.read("output_dots",false);
    if (shouldOutputDots) {
        std::cerr<<"total ojbects = "<<size()<<std::endl;
        std::cerr<<"flag counts = ";
        for(int i=0;i<8;++i) std::cerr<<flagCount[i]<<" ";
        std::cerr<<"  no flag: "<<flagCount[8]<<std::endl;
        std::cerr<<"mag counts = ";
        for(int i=0;i<10;++i) std::cerr<<magCount[i]<<" ";
        std::cerr<<std::endl;
    }

    // Convert input flags into our flag schema
    if (_flags.size() == 0) {
        dbg<<"No flags read in -- starting all with 0\n";
        _flags.resize(size(),0);
    } else {
        long ignoreFlags = ~0L;
        dbg<<std::hex<<std::showbase;
        if (_params.keyExists("coaddcat_ignore_flags")) {
            ignoreFlags = _params["coaddcat_ignore_flags"];
            dbg<<"Using ignore flag parameter = "<<ignoreFlags<<std::endl;
        } else if (_params.keyExists("coaddcat_ok_flags")) {
            ignoreFlags = _params["coaddcat_ok_flags"];
            dbg<<"Using ok flag parameter = "<<ignoreFlags<<std::endl;
            ignoreFlags = ~ignoreFlags;
            dbg<<"ignore flag = "<<ignoreFlags<<std::endl;
        } else {
            dbg<<"No ok_flags or ignore_flags parameter: use ignore_flags = "<<
                ignoreFlags<<std::endl;
        }
        Assert(_flags.size() == size());
        for(int i=0;i<nObj;++i) {
            _flags[i] = (_flags[i] & ignoreFlags) ? INPUT_FLAG : 0;
        }
        dbg<<std::dec<<std::noshowbase;
    }
    if (_params.keyExists("coaddcat_minnu_mag")) {
        double minNu = _params["coaddcat_minnu_mag"];
        dbg<<"Limiting signal to noise in input magnitude to >= "<<
            minNu<<std::endl;
        dbg<<"Marking fainter objects with flag INPUT_FLAG = "<<
            INPUT_FLAG<<std::endl;
        // S/N = 1.086 / magErr (= 2.5*log(e) / magErr)
        double maxMagErr = 1.086/minNu;
        dbg<<"(Corresponds to max mag_err = "<<minNu<<")\n";
        for(int i=0;i<nObj;++i) {
            if (_magErr[i] > maxMagErr) _flags[i] = INPUT_FLAG;
        }
    }
    int goodCount = std::count(_flags.begin(),_flags.end(),0);
    if (shouldOutputDots) {
        std::cerr<<"# good objects = "<<goodCount<<std::endl;
    }
    Bounds skyBounds2; // In degrees, rather than arcsec
    for(int i=0;i<nObj;++i) {
        _skyBounds += _skyPos[i];
        Position temp = _skyPos[i];
        temp /= 3600.;
        skyBounds2 += temp;
    }
    dbg<<"skyBounds = "<<_skyBounds<<std::endl;
    dbg<<"in degrees: "<<skyBounds2<<std::endl;
}

CoaddCatalog::~CoaddCatalog()
{
}


void CoaddCatalog::readCatalog()
{
    std::string file=_params.get("coaddcat_file");
    int hdu = _params.read<int>("coaddcat_hdu");

    if (!doesFileExist(file)) {
        throw FileNotFoundException(file);
    }
    try {
        dbg<<"Opening FITS file "<<file<<" at hdu "<<hdu<<std::endl;
        CCfits::FITS fits(file, CCfits::Read, hdu-1);

        CCfits::ExtHDU& table=fits.extension(hdu-1);

        long nRows=table.rows();

        dbg<<"  nrows = "<<nRows<<std::endl;
        if (nRows <= 0) {
            throw ReadException(
                "CoaddCatalog found to have 0 rows.  Must have > 0 rows.");
        }

        std::string idCol=_params.get("coaddcat_id_col");
        std::string xCol=_params.get("coaddcat_x_col");
        std::string yCol=_params.get("coaddcat_y_col");
        std::string skyCol=_params.get("coaddcat_sky_col");

        std::string magCol=_params.get("coaddcat_mag_col");
        std::string magErrCol=_params.get("coaddcat_mag_err_col");

        std::string flagsCol=_params.get("coaddcat_flag_col");
        std::string raCol=_params.get("coaddcat_ra_col");
        std::string declCol=_params.get("coaddcat_dec_col");

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

        dbg<<"  "<<skyCol<<std::endl;
        table.column(skyCol).read(_sky, start, end);

        dbg<<"  "<<magCol<<std::endl;
        table.column(magCol).read(_mag, start, end);

        dbg<<"  "<<magErrCol<<std::endl;
        table.column(magErrCol).read(_magErr, start, end);

        dbg<<"  "<<flagsCol<<std::endl;
        table.column(flagsCol).read(_flags, start, end);

        dbg<<"  "<<raCol<<"  "<<declCol<<std::endl;
        _skyPos.resize(nRows);
        std::vector<float> ra;
        std::vector<float> decl;
        table.column(raCol).read(ra, start, end);
        table.column(declCol).read(decl, start, end);

        xdbg<<"list of coaddcatalog: (id, pos, ra, dec, mag, flag)\n";
        for(long i=0;i<nRows;++i) {
            _pos[i] = Position(x[i],y[i]);
            _skyPos[i] = Position(ra[i],decl[i]);
            // The convention for Position is to use arcsec for everything.
            // ra and dec come in as degrees.  So wee need to convert to arcsec.
            _skyPos[i] *= 3600.;  // deg -> arcsec
            xdbg<<_id[i]<<"  "<<x[i]<<"  "<<y[i]<<"  "<<ra[i]/15.<<"  "<<
                decl[i]<<"  "<<_mag[i]<<"  "<<_flags[i]<<std::endl;
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
}

