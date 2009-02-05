
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "StarCatalog.h"
#include "StarFinder.h"
#include "Bounds.h"
#include "ConfigFile.h"
#include "FitsFile.h"
#include "Image.h"
#include "Name.h"
#include "Transformation.h"
#include "Pixel.h"
#include "Params.h"
#include "Ellipse.h"
#include "Log.h"

void CalcSigma(
    double& sigma,
    const Image<double>& im, const Position& pos, double sky,
    double noise, double gain, const Image<double>* weight_im, 
    const Transformation& trans, double psfap, long& flag)
{
  std::vector<Pixel> pix;
  long flag1 = 0;
  try {
    GetPixList(im, pix, pos, sky, noise, gain, weight_im, trans, psfap, flag1);
  } catch (Range_error& e) {
    xxdbg<<"transformation range error: \n";
    xxdbg<<"p = "<<pos<<", b = "<<e.b<<std::endl;
    flag1 |= TRANSFORM_EXCEPTION;
  }
  if (flag1) {
    flag |= flag1;
    sigma = DEFVALNEG;
    return;
  }
  xxdbg<<"npix = "<<pix.size()<<std::endl;

  Ellipse ell;
  ell.PeakCentroid(pix,psfap/3.);
  ell.CrudeMeasure(pix,sigma);
  xdbg<<"Crude Measure: centroid = "<<ell.GetCen();
  xdbg<<", mu = "<<ell.GetMu()<<std::endl;
  double mu = ell.GetMu();

  sigma *= exp(mu);
}


StarCatalog::StarCatalog(const InputCatalog& incat,
    ConfigFile& _params, std::string key_prefix) :
  id(incat.id), pos(incat.pos), sky(incat.sky), noise(incat.noise),
  flags(incat.flags), mag(incat.mag), objsize(incat.objsize),
  isastar(id.size(),0), params(_params), prefix(key_prefix)
{
  Assert(id.size() == size());
  Assert(pos.size() == size());
  Assert(sky.size() == size());
  Assert(noise.size() == size());
  Assert(flags.size() == size());
  Assert(mag.size() == size());
  if (objsize.size() == 0) objsize.resize(size(),DEFVALNEG);
  Assert(objsize.size() == size());
  Assert(isastar.size() == size());
}

StarCatalog::StarCatalog(ConfigFile& _params, std::string key_prefix) :
  params(_params), prefix(key_prefix)
{ 
  Read();

  Assert(id.size() == size());
  Assert(pos.size() == size());
  Assert(sky.size() == size());
  Assert(noise.size() == size());
  Assert(flags.size() == size());
  Assert(objsize.size() == size());
  Assert(mag.size() == size());
  Assert(isastar.size() == size());
}

void StarCatalog::CalcSizes(const Image<double>& im, 
    const Image<double>* weight_im, const Transformation& trans)
{
  double psfap = params.get("psf_aperture"); 
  double gain = 0.;
  if (params.keyExists("gain")) gain = params["gain"];
  Assert(pos.size() == size());
  Assert(sky.size() == size());
  Assert(noise.size() == size());
  Assert(objsize.size() == size());
  Assert(flags.size() == size());
  int n = pos.size();
  dbg<<"n = "<<n<<std::endl;
#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
  for (int i=0; i<n; i++) if (!flags[i]) {
    dbg<<"use i = "<<i<<std::endl;

    try {
      // Negative value indicates not set yet.  Start with 1 then.
      if (objsize[i] <= 0.) objsize[i] = 1.;
      CalcSigma(
	  objsize[i],
	  im, pos[i], sky[i], noise[i], gain, weight_im, 
	  trans, psfap, flags[i]);
      dbg<<"objsize["<<i<<"]: "<<objsize[i]<<std::endl;
      dbg<<"flags["<<i<<"]: "<<flags[i]<<std::endl;
    } catch (tmv::Error& e) {
      dbg<<"Caught: "<<e<<std::endl;
      objsize[i] = DEFVALNEG;
      flags[i] |= TMV_EXCEPTION;
    } catch (std::exception& e) {
      dbg<<"Caught: "<<e.what()<<std::endl;
      objsize[i] = DEFVALNEG;
      flags[i] |= STD_EXCEPTION;
    } catch (...) {
      dbg<<"Caught unknown exception"<<std::endl;
      objsize[i] = DEFVALNEG;
      flags[i] |= UNKNOWN_EXCEPTION;
    }
  } // End omp parallel for
  dbg<<"Done MeasureSigmas\n";
}

int StarCatalog::FindStars(FindStarsLog& log)
{
  StarFinder finder(params,prefix);

  std::vector<PotentialStar*> maybestars;

  // First get a list of potential stars
  dbg<<"Finding stars"<<std::endl;
  long count=0;
  log.nobj = size();
  for (size_t i=0; i<pos.size(); i++)
  {
    // A series of checks
    // Only objects with no flags in SExtractor or updated size calculation
    if (flags[i]) {
      xdbg<<"Reject "<<i<<" for input flags\n";
      ++log.nr_flag;
      continue;
    }
    // Range checking
    if (objsize[i] < finder.minsize || objsize[i] > finder.maxsize) {
      xdbg<<"Reject "<<i<<" for size "<<objsize[i]<<" outside range "<<
	finder.minsize<<" -- "<<finder.maxsize<<std::endl;
      ++log.nr_size;
      continue;
    }
    if (mag[i] < finder.minmag || mag[i] > finder.maxmag) {
      xdbg<<"Reject "<<i<<" for mag "<<mag[i]<<" outside range "<<
	finder.minmag<<" -- "<<finder.maxmag<<std::endl;
      ++log.nr_mag;
      continue;
    }
    xdbg<<"OK: "<<objsize[i]<<"  "<<mag[i]<<std::endl;

    count++;
    maybestars.push_back(
	new PotentialStar(pos[i],mag[i],objsize[i],i,""));
  }
  log.nobj = count;
  dbg<<"  Possible Stars: "<<count<<"/"<<pos.size()<<"\n";

  dbg<<"  Running FindStars\n";
  std::vector<PotentialStar*> stars = finder.FindStars(maybestars);
  dbg<<"  Found "<<stars.size()<<"\n";
  log.nallstars = stars.size();

  isastar.resize(pos.size(),0);
  count = 0;
  for (size_t k=0; k<stars.size();k++) {
    int i=stars[k]->GetIndex();
    if (mag[i] <= finder.maxoutmag) {
      isastar[i] = 1;
      count++;
    }
  }
  dbg<<"  Cut to "<<count<<" by maxoutmag cut\n";

  log.nstars = count;

  if (count < 100) {
    std::cout<<"STATUS3BEG Warning: Only "<<count
      <<" stars found for Name="<<Name(params,"stars")
      <<". STATUS3END"<<std::endl;
  }

  return count;
}

void StarCatalog::WriteFits(std::string file) const 
{
  Assert(id.size() == size());
  Assert(pos.size() == size());
  Assert(sky.size() == size());
  Assert(noise.size() == size());
  Assert(flags.size() == size());
  Assert(mag.size() == size());
  Assert(objsize.size() == size());
  Assert(isastar.size() == size());

  FitsFile fits(file, READWRITE, true);

  // Setup the binary table
  std::string id_col=params.get(prefix + "id_col");
  std::string x_col=params.get(prefix + "x_col");
  std::string y_col=params.get(prefix + "y_col");
  std::string sky_col=params.get(prefix + "sky_col");
  std::string noise_col=params.get(prefix + "noise_col");
  std::string flags_col=params.get(prefix + "flags_col");
  std::string mag_col=params.get(prefix + "mag_col");
  std::string objsize_col=params.get(prefix + "objsize_col");
  std::string isastar_col=params.get(prefix + "isastar_col");

  const int nfields = 9;
  std::string table_cols[nfields] = {
    id_col,
    x_col,
    y_col,
    sky_col,
    noise_col,
    flags_col,
    mag_col,
    objsize_col,
    isastar_col
  };
  std::string table_types[nfields] = {
    "1j",
    "1d",
    "1d",
    "1d",
    "1d",
    "1i",
    "1d",
    "1d",
    "1b"
  };
  std::string table_units[nfields] =  {
    "None",
    "pixels",
    "pixels",
    "ADU",
    "ADU^2",
    "None",
    "mags",
    "Arcsec",
    "None"
  };
  fits.CreateBinaryTable(size(),nfields,table_cols,table_types,table_units);

  // Write the header keywords
  fits.WriteParKey(params, "version", TSTRING);
  fits.WriteParKey(params, "noise_method", TSTRING);
  fits.WriteParKey(params, "dist_method", TSTRING);

  fits.WriteParKey(params, prefix + "minsize", TDOUBLE);
  fits.WriteParKey(params, prefix + "maxsize", TDOUBLE);
  fits.WriteParKey(params, prefix + "minmag", TDOUBLE);
  fits.WriteParKey(params, prefix + "maxmag", TDOUBLE);
  fits.WriteParKey(params, prefix + "ndivx", TLONG);
  fits.WriteParKey(params, prefix + "ndivy", TLONG);

  fits.WriteParKey(params, prefix + "startn1", TDOUBLE);
  fits.WriteParKey(params, prefix + "starfrac", TDOUBLE);
  fits.WriteParKey(params, prefix + "magstep1", TDOUBLE);
  fits.WriteParKey(params, prefix + "miniter1", TLONG);
  fits.WriteParKey(params, prefix + "reject1", TDOUBLE);
  fits.WriteParKey(params, prefix + "binsize1", TDOUBLE);
  fits.WriteParKey(params, prefix + "maxratio1", TDOUBLE);
  fits.WriteParKey(params, prefix + "okvalcount", TLONG);
  fits.WriteParKey(params, prefix + "maxrms", TDOUBLE);
  fits.WriteParKey(params, prefix + "starsperbin", TLONG);

  fits.WriteParKey(params, prefix + "fitorder", TLONG);
  fits.WriteParKey(params, prefix + "fitsigclip", TDOUBLE);
  fits.WriteParKey(params, prefix + "startn2", TDOUBLE);
  fits.WriteParKey(params, prefix + "magstep2", TDOUBLE);
  fits.WriteParKey(params, prefix + "miniter2", TLONG);
  fits.WriteParKey(params, prefix + "minbinsize", TDOUBLE);
  fits.WriteParKey(params, prefix + "reject2", TDOUBLE);

  fits.WriteParKey(params, prefix + "purityratio", TDOUBLE);
  fits.WriteParKey(params, prefix + "maxrefititer", TLONG);

  // Write the data
  fits.WriteColumn(TLONG, 1, 1, 1, size(), &id[0]);
  std::vector<double> x(pos.size());
  std::vector<double> y(pos.size());
  for(size_t i=0;i<pos.size();i++) {
    x[i] = pos[i].GetX();
    y[i] = pos[i].GetY();
  }
  fits.WriteColumn(TDOUBLE, 2, 1, 1, size(), &x[0]);
  fits.WriteColumn(TDOUBLE, 3, 1, 1, size(), &y[0]);
  fits.WriteColumn(TDOUBLE, 4, 1, 1, size(), &sky[0]);
  fits.WriteColumn(TDOUBLE, 5, 1, 1, size(), &noise[0]);
  fits.WriteColumn(TLONG, 6, 1, 1, size(), &flags[0]);
  fits.WriteColumn(TFLOAT, 7, 1, 1, size(), &mag[0]);
  fits.WriteColumn(TDOUBLE, 8, 1, 1, size(), &objsize[0]);
  fits.WriteColumn(TINT, 9, 1, 1, size(), &isastar[0]);
}

void StarCatalog::WriteAscii(std::string file, std::string delim) const 
{
  Assert(id.size() == size());
  Assert(pos.size() == size());
  Assert(sky.size() == size());
  Assert(noise.size() == size());
  Assert(flags.size() == size());
  Assert(mag.size() == size());
  Assert(objsize.size() == size());
  Assert(isastar.size() == size());

  std::ofstream fout(file.c_str());
  if (!fout) {
    throw std::runtime_error("Error opening stars file");
  }

  for (size_t i=0; i<size(); i++) {
    fout 
      << id[i] << delim
      << pos[i].GetX() << delim
      << pos[i].GetY() << delim
      << sky[i] << delim
      << noise[i] << delim
      << flags[i] << delim
      << mag[i] << delim
      << objsize[i] << delim
      << isastar[i] << std::endl;
  }
}

void StarCatalog::Write() const 
{
  std::string file = Name(params, "stars");  
  dbg<<"Writing to stars file: "<<file<<std::endl;

  if (file.find("fits") != std::string::npos) {
    WriteFits(file);
  } else {
    std::string delim = "  ";
    if (params.keyExists(prefix + "delim")) delim = params[prefix + "delim"];
    WriteAscii(file,delim);
  }
  dbg<<"Done Write\n";
}

void StarCatalog::ReadFits(std::string file) 
{
  int hdu = 2;
  if (params.keyExists(prefix + "hdu")) hdu = params[prefix + "hdu"];

  FitsFile fits(file);

  dbg<<"Moving to HDU #"<<hdu<<std::endl;
  fits.GotoHDU(hdu);

  long nrows=0;
  fits.ReadKey("NAXIS2", TLONG, &nrows);

  dbg<<"  nrows = "<<nrows<<std::endl;
  if (nrows <= 0) {
    throw std::runtime_error("nrows must be >= 0");
  }

  std::string id_col=params.get(prefix + "id_col");
  std::string x_col=params.get(prefix + "x_col");
  std::string y_col=params.get(prefix + "y_col");
  std::string sky_col=params.get(prefix + "sky_col");
  std::string noise_col=params.get(prefix + "noise_col");
  std::string flags_col=params.get(prefix + "flags_col");
  std::string mag_col=params.get(prefix + "mag_col");
  std::string objsize_col=params.get(prefix + "objsize_col");
  std::string isastar_col=params.get(prefix + "isastar_col");

  dbg<<"Reading columns"<<std::endl;
  dbg<<"  "<<id_col<<std::endl;
  id.resize(nrows);
  fits.ReadScalarCol(id_col,TLONG,&id[0], nrows);

  dbg<<"  "<<x_col<<"  "<<y_col<<std::endl;
  pos.resize(nrows);
  std::vector<double> x(nrows);
  std::vector<double> y(nrows);
  fits.ReadScalarCol(x_col,TDOUBLE,&x[0], nrows);
  fits.ReadScalarCol(y_col,TDOUBLE,&y[0], nrows);
  for(long i=0;i<nrows;++i) pos[i] = Position(x[i],y[i]);

  dbg<<"  "<<sky_col<<std::endl;
  sky.resize(nrows);
  fits.ReadScalarCol(sky_col,TDOUBLE,&sky[0], nrows);

  dbg<<"  "<<noise_col<<std::endl;
  noise.resize(nrows);
  fits.ReadScalarCol(noise_col,TDOUBLE,&noise[0], nrows);

  dbg<<"  "<<flags_col<<std::endl;
  flags.resize(nrows);
  fits.ReadScalarCol(flags_col,TLONG,&flags[0], nrows);

  dbg<<"  "<<mag_col<<std::endl;
  mag.resize(nrows);
  fits.ReadScalarCol(mag_col,TDOUBLE,&mag[0], nrows);

  dbg<<"  "<<objsize_col<<std::endl;
  objsize.resize(nrows);
  fits.ReadScalarCol(objsize_col,TDOUBLE,&objsize[0], nrows);

  dbg<<"  "<<isastar_col<<std::endl;
  isastar.resize(nrows);
  fits.ReadScalarCol(isastar_col,TLONG,&isastar[0], nrows);
}

void StarCatalog::ReadAscii(std::string file, std::string delim)
{
  std::ifstream fin(file.c_str());
  if (!fin) {
    throw std::runtime_error("Error opening stars file");
  }

  if (delim == "  ") {
    long id1,flag,star;
    double x,y,sky1,n,s,m;
    while ( fin >> id1 >> x >> y >> sky1 >> n >> flag >> m >> s >> star)
    {
      id.push_back(id1);
      pos.push_back(Position(x,y));
      sky.push_back(sky1);
      noise.push_back(n);
      flags.push_back(flag);
      mag.push_back(m);
      objsize.push_back(s);
      isastar.push_back(star);
    } 
  } else {
    if (delim.size() > 1) {
      // getline only works with single character delimiters.
      // Since I don't really expect a multicharacter delimiter to
      // be used ever, I'm just going to throw an exception here 
      // if we do need it, and I can write the workaround then.
      throw std::runtime_error(
	  "ReadAscii delimiter must be a single character");
    }
    char d = delim[0];
    double x,y;
    ConvertibleString temp;
    while (getline(fin,temp,d))
    {
      id.push_back(temp);
      getline(fin,temp,d); x = temp;
      getline(fin,temp,d); y = temp;
      pos.push_back(Position(x,y));
      getline(fin,temp,d); sky.push_back(temp);
      getline(fin,temp,d); noise.push_back(temp);
      getline(fin,temp,d); flags.push_back(temp);
      getline(fin,temp,d); mag.push_back(temp);
      getline(fin,temp,d); objsize.push_back(temp);
      getline(fin,temp); isastar.push_back(temp);
    }
  }
}

void StarCatalog::Read()
{
  std::string file = Name(params,"stars",false,true);
  // false,true = input_prefix=false, mustexist=true.
  // It is an input here, but it is in the output_prefix directory.
  dbg<< "Reading StarCatalog from file: " << file << std::endl;

  if (file.find("fits") != std::string::npos) {
    ReadFits(file);
  } else {
    std::string delim = "  ";
    if (params.keyExists(prefix + "delim")) delim = params[prefix + "delim"];
    ReadAscii(file,delim);
  }
  dbg<<"Done Read StarCatalog\n";
}

