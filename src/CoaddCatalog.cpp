
#include <valarray>
#include "TMV.h"
#include <CCfits/CCfits>

#include "CoaddCatalog.h"
#include "Ellipse.h"
#include "ConfigFile.h"
#include "dbg.h"
#include "Params.h"
#include "BVec.h"
#include "Ellipse.h"
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

  Assert(pos.size() == id.size());
  Assert(skypos.size() == id.size());
  Assert(ra.size() == id.size());
  Assert(dec.size() == id.size());
  Assert(sky.size() == id.size());
  Assert(mag.size() == id.size());
  Assert(mag_err.size() == id.size());
  Assert(flags.size() == id.size());

  std::vector<int> flagcount(9,0);
  std::vector<int> magcount(10,0);
  for (size_t i=0; i<size(); i++) 
  {
    if (flags[i] & 1) flagcount[0]++; // neighbor or badpix
    if (flags[i] & 2) flagcount[1]++; // deblended
    if (flags[i] & 4) flagcount[2]++; // saturated
    if (flags[i] & 8) flagcount[3]++; // edge
    if (flags[i] & 16) flagcount[4]++; // incomplete aperture
    if (flags[i] & 32) flagcount[5]++; // incomplete isophotal
    if (flags[i] & 64) flagcount[6]++; // memory overflow in deblend
    if (flags[i] & 128) flagcount[7]++; // memory overflow in extract
    if (!flags[i]) flagcount[8]++; // no flags

    if (mag[i] < 19) magcount[0]++;
    else if (mag[i] < 20) magcount[1]++;
    else if (mag[i] < 21) magcount[2]++;
    else if (mag[i] < 22) magcount[3]++;
    else if (mag[i] < 23) magcount[4]++;
    else if (mag[i] < 24) magcount[5]++;
    else if (mag[i] < 25) magcount[6]++;
    else if (mag[i] < 26) magcount[7]++;
    else if (mag[i] < 27) magcount[8]++;
    else magcount[9]++;
  }

  std::cout<<"total ojbects = "<<size()<<std::endl;
  std::cout<<"flag counts = ";
  for(int i=0;i<8;++i) std::cout<<flagcount[i]<<" ";
  std::cout<<"  no flag: "<<flagcount[8]<<std::endl;
  std::cout<<"mag counts = ";
  for(int i=0;i<10;++i) std::cout<<magcount[i]<<" ";
  std::cout<<std::endl;

  // Convert input flags into our flag schema
  if (flags.size() == 0) 
  {
    dbg<<"No flags read in -- starting all with 0\n";
    flags.resize(size(),0);
  }
  else 
  {
    long ignore_flags = ~0L;
    dbg<<std::hex<<std::showbase;
    if (params.keyExists("coaddcat_ignore_flags")) 
    {
      ignore_flags = params["coaddcat_ignore_flags"];
      dbg<<"Using ignore flag parameter = "<<ignore_flags<<std::endl;
    }
    else if (params.keyExists("coaddcat_ok_flags")) 
    {
      ignore_flags = params["coaddcat_ok_flags"];
      dbg<<"Using ok flag parameter = "<<ignore_flags<<std::endl;
      ignore_flags = ~ignore_flags;
      dbg<<"ignore flag = "<<ignore_flags<<std::endl;
    } 
    else 
    {
      dbg<<"No ok_flags or ignore_flags parameter: use ignore_flags = "<<
	ignore_flags<<std::endl;
    }
    Assert(flags.size() == size());
    for(size_t i=0;i<size();++i) 
    {
      flags[i] = (flags[i] & ignore_flags) ? INPUT_FLAG : 0;
    }
    if (params.keyExists("coaddcat_max_mag")) 
    {
      double max_mag = params["coaddcat_max_mag"];
      dbg<<"Limiting input magnitude to mag <= "<<max_mag<<std::endl;
      dbg<<"Marking fainter objects with flag INPUT_FLAG = "<<INPUT_FLAG<<std::endl;
      for(size_t i=0;i<size();++i) 
      {
	if (mag[i] > max_mag) flags[i] = INPUT_FLAG;
      }
    }
    int goodcount = std::count(flags.begin(),flags.end(),0);
    std::cout<<"# good objects = "<<goodcount<<std::endl;

    dbg<<std::dec<<std::noshowbase;
  }
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

    dbg<<"  "<<sky_col<<std::endl;
    table.column(sky_col).read(sky, start, end);

    //dbg<<"  "<<noise_col<<std::endl;
    //table.column(noise_col).read(noise, start, end);
    noise.resize(nrows,0);

    dbg<<"  "<<mag_col<<std::endl;
    table.column(mag_col).read(mag, start, end);

    dbg<<"  "<<mag_err_col<<std::endl;
    table.column(mag_err_col).read(mag_err, start, end);

    dbg<<"  "<<flags_col<<std::endl;
    table.column(flags_col).read(flags, start, end);

    dbg<<"  "<<ra_col<<"  "<<dec_col<<std::endl;

    table.column(ra_col).read(ra, start, end);
    table.column(dec_col).read(dec, start, end);

    xdbg<<"list of coaddcatalog: (id, pos, ra, dec, mag, flag)\n";
    skypos.resize(nrows);
    for(long i=0;i<nrows;++i) {
      pos[i] = Position(x[i],y[i]);
      skypos[i] = Position(ra[i],dec[i]);
      // The convention for Position is to use arcsec for everything.
      // ra and dec come in as degrees.  So wee need to convert to arcsec.
      skypos[i] *= 3600.;  // deg -> arcsec
      skybounds += skypos[i];
      xdbg<<id[i]<<"  "<<x[i]<<"  "<<y[i]<<"  "<<ra[i]/15.<<"  "<<dec[i]<<"  "<<mag[i]<<"  "<<flags[i]<<std::endl;
    }
  }
  catch (CCfits::FitsException& e)
  {
    throw ReadError("Error reading from "+file+" -- caught error\n" +
	e.message());
  }
  catch (std::exception& e)
  {
    throw ReadError("Error reading from "+file+" -- caught error\n" +
	e.what());
  }
  catch (...)
  {
    throw ReadError("Error reading from "+file+" -- caught unknown error");
  }
}

