
#include <valarray>
#include <CCfits/CCfits>

#include "CoaddCatalog.h"
#include "dbg.h"
#include "Params.h"
#include "Name.h"

CoaddCatalog::CoaddCatalog(const ConfigFile& params) :
    _params(params)
{}

CoaddCatalog::~CoaddCatalog() 
{}


void CoaddCatalog::read()
{
    std::string file=_params.get("coaddcat_file");

    if (!DoesFileExist(file)) {
        throw FileNotFoundException(file);
    }

    try {
        // No ASCII version currently.
        readFits(file);
    } catch (CCfits::FitsException& e) {
        xdbg<<"Caught FitsException: \n"<<e.message()<<std::endl;
        throw ReadException(
            "Error reading from "+file+" -- caught error\n" + e.message());
    } catch (std::exception& e) {
        xdbg<<"Caught std::exception: \n"<<e.what()<<std::endl;
        throw ReadException(
            "Error reading from "+file+" -- caught error\n" + e.what());
    } catch (...) {
        xdbg<<"Caught unknown exception: "<<std::endl;
        throw ReadException(
            "Error reading from "+file+" -- caught unknown error");
    }

    Assert(_pos.size() == _id.size());
    Assert(_skypos.size() == _id.size());
    Assert(_sky.size() == _id.size());
    Assert(_mag.size() == _id.size());
    Assert(_mag_err.size() == _id.size());
    Assert(_flags.size() == _id.size());

    std::vector<int> flag_count(9,0);
    std::vector<int> mag_count(10,0);
    const int nobj = size();
    for (int i=0; i<nobj; ++i) {
        if (_flags[i] & 1) ++flag_count[0]; // neighbor or badpix
        if (_flags[i] & 2) ++flag_count[1]; // deblended
        if (_flags[i] & 4) ++flag_count[2]; // saturated
        if (_flags[i] & 8) ++flag_count[3]; // edge
        if (_flags[i] & 16) ++flag_count[4]; // incomplete aperture
        if (_flags[i] & 32) ++flag_count[5]; // incomplete isophotal
        if (_flags[i] & 64) ++flag_count[6]; // memory overflow in deblend
        if (_flags[i] & 128) ++flag_count[7]; // memory overflow in extract
        if (!_flags[i]) ++flag_count[8]; // no flags

        if (_mag[i] < 19) ++mag_count[0];
        else if (_mag[i] < 20) ++mag_count[1];
        else if (_mag[i] < 21) ++mag_count[2];
        else if (_mag[i] < 22) ++mag_count[3];
        else if (_mag[i] < 23) ++mag_count[4];
        else if (_mag[i] < 24) ++mag_count[5];
        else if (_mag[i] < 25) ++mag_count[6];
        else if (_mag[i] < 26) ++mag_count[7];
        else if (_mag[i] < 27) ++mag_count[8];
        else ++mag_count[9];
    }

    bool output_dots = _params.read("output_dots",false);
    dbg<<"total ojbects = "<<size()<<std::endl;
    dbg<<"flag counts = ";
    for(int i=0;i<8;++i) dbg<<flag_count[i]<<" ";
    dbg<<"  no flag: "<<flag_count[8]<<std::endl;
    dbg<<"mag counts = ";
    for(int i=0;i<10;++i) dbg<<mag_count[i]<<" ";
    dbg<<std::endl;

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
        Assert(int(_flags.size()) == size());
        for(int i=0;i<nobj;++i) {
            _flags[i] = (_flags[i] & ignoreFlags) ? INPUT_FLAG : 0;
        }
        dbg<<std::dec<<std::noshowbase;
    }
    if (_params.keyExists("coaddcat_minnu_mag")) {
        double min_nu = _params["coaddcat_minnu_mag"];
        dbg<<"Limiting signal to noise in input magnitude to >= "<<
            min_nu<<std::endl;
        dbg<<"Marking fainter objects with flag INPUT_FLAG = "<<
            INPUT_FLAG<<std::endl;
        // S/N = 1.086 / mag_err (= 2.5*log(e) / mag_err)
        double max_mag_err = 1.086/min_nu;
        dbg<<"(Corresponds to max mag_err = "<<max_mag_err<<")\n";
        for(int i=0;i<nobj;++i) {
            if (_mag_err[i] > max_mag_err) _flags[i] = INPUT_FLAG;
        }
    }
    int good_count = std::count(_flags.begin(),_flags.end(),0);
    if (output_dots) {
        std::cerr<<"Total # good objects = "<<good_count<<std::endl;
    }
    Bounds skybounds2; // In degrees, rather than arcsec
    for(int i=0;i<nobj;++i) if (!_flags[i]) {
        _skybounds += _skypos[i];
        Position temp = _skypos[i];
        temp /= 3600.;
        skybounds2 += temp;
    }
    dbg<<"skybounds = "<<_skybounds<<std::endl;
    dbg<<"in degrees: "<<skybounds2<<std::endl;
}

void CoaddCatalog::readFits(const std::string& file)
{
    dbg<<"Reading coadd catalog from FITS file: "<<file<<std::endl;

    int hdu = GetHdu(_params,"coaddcat",file,1);

    dbg<<"Opening CoaddCatalog file "<<file<<" at hdu "<<hdu<<std::endl;
    CCfits::FITS fits(file, CCfits::Read);
    if (hdu > 1) fits.read(hdu-1);

    CCfits::ExtHDU& table=fits.extension(hdu-1);

    long nrows=table.rows();

    dbg<<"  nrows = "<<nrows<<std::endl;
    if (nrows <= 0) {
        throw ReadException(
            "CoaddCatalog found to have 0 rows.  Must have > 0 rows.");
    }

    std::string id_col=_params.get("coaddcat_id_col");
    std::string x_col=_params.get("coaddcat_x_col");
    std::string y_col=_params.get("coaddcat_y_col");
    std::string sky_col=_params.get("coaddcat_sky_col");

    std::string mag_col=_params.get("coaddcat_mag_col");
    std::string mag_err_col=_params.get("coaddcat_mag_err_col");

    std::string flags_col=_params.get("coaddcat_flag_col");
    std::string ra_col=_params.get("coaddcat_ra_col");
    std::string decl_col=_params.get("coaddcat_dec_col");

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

    dbg<<"  "<<sky_col<<std::endl;
    table.column(sky_col).read(_sky, start, end);

    dbg<<"  "<<mag_col<<std::endl;
    table.column(mag_col).read(_mag, start, end);

    dbg<<"  "<<mag_err_col<<std::endl;
    table.column(mag_err_col).read(_mag_err, start, end);

    dbg<<"  "<<flags_col<<std::endl;
    table.column(flags_col).read(_flags, start, end);

    dbg<<"  "<<ra_col<<"  "<<decl_col<<std::endl;
    _skypos.resize(nrows);
    std::vector<double> ra;
    std::vector<double> decl;
    table.column(ra_col).read(ra, start, end);
    table.column(decl_col).read(decl, start, end);

    xdbg<<"list of coaddcatalog: (id, pos, ra, dec, mag, flag)\n";
    for(long i=0;i<nrows;++i) {
        _pos[i] = Position(x[i],y[i]);
        _skypos[i] = Position(ra[i],decl[i]);
        // The convention for Position is to use arcsec for everything.
        // ra and dec come in as degrees.  So wee need to convert to arcsec.
        _skypos[i] *= 3600.;  // deg -> arcsec
        xdbg<<_id[i]<<"  "<<x[i]<<"  "<<y[i]<<"  "<<ra[i]/15.<<"  "<<
            decl[i]<<"  "<<_mag[i]<<"  "<<_flags[i]<<std::endl;
    }
}

