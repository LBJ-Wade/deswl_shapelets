#include "CoaddCatalog.h"
#include "Ellipse.h"
#include "ConfigFile.h"
#include "dbg.h"
#include "Params.h"
#include "BVec.h"
#include "Ellipse.h"
#include "TMV.h"
#include "Pixel.h"
#include "Image.h"
#include "FittedPSF.h"
#include "Log.h"
#include "TimeVars.h"

#define UseInverseTransform

//#define OnlyNImages 10

CoaddCatalog::CoaddCatalog(ConfigFile& _params):
  params(_params)
{
  ReadCatalog();
}

CoaddCatalog::~CoaddCatalog()
{
}


void CoaddCatalog::ReadCatalog()
{
  std::string file=params.get("coaddcat_file");
  // I use read rather than get here to make sure we turn any
  // ConvertibleStringError into a ConfigFile_ParameterError
  int hdu = params.read<int>("coaddcat_hdu");

  if (!FileExists(file))
  {
    throw FileNotFound(file);
  }
  try
  {
    dbg<<"Opening FITS file at hdu "<<hdu<<std::endl;
    // true means read all as part of the construction
    CCfits::FITS fits(file, CCfits::Read, hdu-1, true);

    CCfits::ExtHDU& table=fits.extension(hdu-1);

    long nrows=table.rows();

    dbg<<"  nrows = "<<nrows<<std::endl;
    if (nrows <= 0) {
      throw ReadError("CoaddCatalog found to have 0 rows.  Must have > 0 rows.");
    }

    std::string id_col=params.get("coaddcat_id_col");
    std::string x_col=params.get("coaddcat_x_col");
    std::string y_col=params.get("coaddcat_y_col");
    std::string sky_col=params.get("coaddcat_sky_col");

    std::string mag_col=params.get("coaddcat_mag_col");
    std::string mag_err_col=params.get("coaddcat_mag_err_col");

    //std::string noise_col=params.get("coaddcat_noise_col");
    std::string flags_col=params.get("coaddcat_flag_col");
    std::string ra_col=params.get("coaddcat_ra_col");
    std::string dec_col=params.get("coaddcat_dec_col");

    long start=1;
    long end=nrows;

    dbg<<"Reading columns"<<std::endl;
    dbg<<"  "<<id_col<<std::endl;
    table.column(id_col).read(id, start, end);

    dbg<<"  "<<x_col<<"  "<<y_col<<std::endl;
    pos.resize(nrows);
    std::vector<double> x;
    std::vector<double> y;
    table.column(x_col).read(x, start, end);
    table.column(y_col).read(y, start, end);
    for(long i=0;i<nrows;++i) pos[i] = Position(x[i],y[i]);

    dbg<<"  "<<sky_col<<std::endl;
    table.column(sky_col).read(sky, start, end);

    //dbg<<"  "<<noise_col<<std::endl;
    //table.column(noise_col).read(noise, start, end);
    noise.resize(nrows,0);

    dbg<<"  "<<flags_col<<std::endl;
    table.column(flags_col).read(flags, start, end);

    dbg<<"  "<<ra_col<<"  "<<dec_col<<std::endl;

    std::vector<double> ra;
    std::vector<double> dec;
    table.column(ra_col).read(ra, start, end);
    table.column(dec_col).read(dec, start, end);

    skypos.resize(nrows);
    for(long i=0;i<nrows;++i) {
      skypos[i] = Position(ra[i],dec[i]);
      // The convention for Position is to use arcsec for everything.
      // ra and dec come in as degrees.  So wee need to convert to arcsec.
      skypos[i] *= 3600.;  // deg -> arcsec
    }
  }
  catch (std::runtime_error& e)
  {
    throw ReadError("Error reading from "+file+" -- caught error\n" +
	e.what());
  }
  catch (...)
  {
    throw ReadError("Error reading from "+file+" -- caught unknown error");
  }
}

