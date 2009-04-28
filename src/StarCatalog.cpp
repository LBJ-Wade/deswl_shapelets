
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
#include "Form.h"

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
  Assert(sigma > 0);
}


StarCatalog::StarCatalog(const InputCatalog& incat,
    ConfigFile& _params, std::string fs_prefix) :
  id(incat.id), pos(incat.pos), sky(incat.sky), noise(incat.noise),
  flags(incat.flags), mag(incat.mag), objsize(incat.objsize),
  isastar(id.size(),0), params(_params), prefix(fs_prefix)
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

StarCatalog::StarCatalog(ConfigFile& _params, std::string fs_prefix) :
  params(_params), prefix(fs_prefix)
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
    const Image<double>*const weight_im, const Transformation& trans)
{
  Assert(pos.size() == size());
  Assert(sky.size() == size());
  Assert(noise.size() == size());
  Assert(objsize.size() == size());
  Assert(flags.size() == size());

  const int n = pos.size();
  dbg<<"n = "<<n<<std::endl;
  double psfap = params.get("psf_aperture"); 
  double gain = params.read("image_gain",0.);

#ifdef _OPENMP
#pragma omp parallel 
  {
#pragma omp for schedule(guided)
#endif
    for (int i=0; i<n; i++) if (!flags[i]) {
      ompdbg<<"use i = "<<i<<std::endl;

      try {
	// Negative value indicates not set yet.  Start with 1 then.
	if (objsize[i] <= 0.) objsize[i] = 1.;
	CalcSigma(
	    objsize[i],
	    im, pos[i], sky[i], noise[i], gain, weight_im, 
	    trans, psfap, flags[i]);
	ompdbg<<"objsize["<<i<<"]: "<<objsize[i]<<std::endl;
	ompdbg<<"flags["<<i<<"]: "<<flags[i]<<std::endl;
      } catch (tmv::Error& e) {
	ompdbg<<"Caught: "<<e<<std::endl;
	objsize[i] = DEFVALNEG;
	flags[i] |= TMV_EXCEPTION;
      } catch (std::exception& e) {
	ompdbg<<"Caught: "<<e.what()<<std::endl;
	objsize[i] = DEFVALNEG;
	flags[i] |= STD_EXCEPTION;
      } catch (...) {
	ompdbg<<"Caught unknown exception"<<std::endl;
	objsize[i] = DEFVALNEG;
	flags[i] |= UNKNOWN_EXCEPTION;
      }
    }
#ifdef _OPENMP
  } // End omp parallel 
#endif
  dbg<<"Done MeasureSigmas\n";
}

int StarCatalog::FindStars(FindStarsLog& log)
{
  StarFinder finder(params,prefix);

  std::vector<PotentialStar*> maybestars;

  // First get a list of potential stars
  dbg<<"Finding stars"<<std::endl;
  long count=0;
  log.ntot = size();
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
    double logobjsize = finder.logsize ? objsize[i] : std::log(objsize[i]);
    maybestars.push_back(
	new PotentialStar(pos[i],mag[i],logobjsize,i,""));
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

  if (params.read("des_qa",false)) {
    if (count < 100) {
      std::cout<<"STATUS3BEG Warning: Only "<<count
	<<" stars found for Name="<<Name(params,"stars")
	<<". STATUS3END"<<std::endl;
    }
  }

  return count;
}

#include <CCfits/CCfits>
void StarCatalog::WriteFitsCCfits(std::string file) const 
{
  Assert(id.size() == size());
  Assert(pos.size() == size());
  Assert(sky.size() == size());
  Assert(noise.size() == size());
  Assert(flags.size() == size());
  Assert(mag.size() == size());
  Assert(objsize.size() == size());
  Assert(isastar.size() == size());

  // ! means overwrite existing file
  CCfits::FITS fits("!"+file, CCfits::Write);

  const int nfields = 9;

  std::vector<string> colnames(nfields);
  std::vector<string> colfmts(nfields);
  std::vector<string> colunits(nfields);

  colnames[0] = params["stars_id_col"];
  colnames[1] = params["stars_x_col"];
  colnames[2] = params["stars_y_col"];
  colnames[3] = params["stars_sky_col"];
  colnames[4] = params["stars_noise_col"];
  colnames[5] = params["stars_flags_col"];
  colnames[6] = params["stars_mag_col"];
  colnames[7] = params["stars_objsize_col"];
  colnames[8] = params["stars_isastar_col"];

  colfmts[0] = "1J"; // id
  colfmts[1] = "1D"; // x
  colfmts[2] = "1D"; // y
  colfmts[3] = "1D"; // sky
  colfmts[4] = "1D"; // noise
  colfmts[5] = "1J"; // flags
  colfmts[6] = "1E"; // mag
  colfmts[7] = "1D"; // sigma
  colfmts[8] = "1J"; // star flag

  colunits[0] = "None";     // id
  colunits[1] = "pixels";   // x
  colunits[2] = "pixels";   // y
  colunits[3] = "ADU";      // sky
  colunits[4] = "ADU^2";    // noise
  colunits[5] = "None";     // flags
  colunits[6] = "mags";     // mag
  colunits[7] = "Arcsec";   // sigma0
  colunits[8] = "None";     //star flag

  CCfits::Table* table;
  table = fits.addTable("findstars",size(),colnames,colfmts,colunits);

  // make vector copies for writing
  std::vector<double> x(pos.size());
  std::vector<double> y(pos.size());
  for(size_t i=0;i<pos.size();i++) {
    x[i] = pos[i].GetX();
    y[i] = pos[i].GetY();
  }

  int startrow=1;
  table->column(colnames[0]).write(id,startrow);
  table->column(colnames[1]).write(x,startrow);
  table->column(colnames[2]).write(y,startrow);
  table->column(colnames[3]).write(sky,startrow);
  table->column(colnames[4]).write(noise,startrow);
  table->column(colnames[5]).write(flags,startrow);
  table->column(colnames[6]).write(mag,startrow);
  table->column(colnames[7]).write(objsize,startrow);
  table->column(colnames[8]).write(isastar,startrow);

  // kind of kludgy but is more flexible than using type numbers or strings
  double dbl;
  int intgr;
  std::string str;

  WriteParKeyCCfits(params, table, "version", str);
  WriteParKeyCCfits(params, table, "noise_method", str);
  WriteParKeyCCfits(params, table, "dist_method", str);

  WriteParKeyCCfits(params, table, "stars_minsize", dbl);
  WriteParKeyCCfits(params, table, "stars_maxsize", dbl);
  WriteParKeyCCfits(params, table, "stars_minmag", dbl);
  WriteParKeyCCfits(params, table, "stars_maxmag", dbl);
  WriteParKeyCCfits(params, table, "stars_ndivx", intgr);
  WriteParKeyCCfits(params, table, "stars_ndivy", intgr);

  WriteParKeyCCfits(params, table, "stars_startn1", dbl);
  WriteParKeyCCfits(params, table, "stars_starfrac", dbl);
  WriteParKeyCCfits(params, table, "stars_magstep1", dbl);
  WriteParKeyCCfits(params, table, "stars_miniter1", intgr);
  WriteParKeyCCfits(params, table, "stars_reject1", dbl);
  WriteParKeyCCfits(params, table, "stars_binsize1", dbl);
  WriteParKeyCCfits(params, table, "stars_maxratio1", dbl);
  WriteParKeyCCfits(params, table, "stars_okvalcount", intgr);
  WriteParKeyCCfits(params, table, "stars_maxrms", dbl);
  WriteParKeyCCfits(params, table, "stars_starsperbin", intgr);

  WriteParKeyCCfits(params, table, "stars_fitorder", intgr);
  WriteParKeyCCfits(params, table, "stars_fitsigclip", dbl);
  WriteParKeyCCfits(params, table, "stars_startn2", dbl);
  WriteParKeyCCfits(params, table, "stars_magstep2", dbl);
  WriteParKeyCCfits(params, table, "stars_miniter2", intgr);
  WriteParKeyCCfits(params, table, "stars_minbinsize", dbl);
  WriteParKeyCCfits(params, table, "stars_reject2", dbl);

  WriteParKeyCCfits(params, table, "stars_purityratio", dbl);
  WriteParKeyCCfits(params, table, "stars_maxrefititer", intgr);

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
  std::string id_col=params.get("stars_id_col");
  std::string x_col=params.get("stars_x_col");
  std::string y_col=params.get("stars_y_col");
  std::string sky_col=params.get("stars_sky_col");
  std::string noise_col=params.get("stars_noise_col");
  std::string flags_col=params.get("stars_flags_col");
  std::string mag_col=params.get("stars_mag_col");
  std::string objsize_col=params.get("stars_objsize_col");
  std::string isastar_col=params.get("stars_isastar_col");

  const int nfields = 9;
  std::string table_cols[nfields] = {
    id_col,
    x_col, y_col,
    sky_col,
    noise_col,
    flags_col,
    mag_col,
    objsize_col,
    isastar_col
  };
  int table_nelem[nfields] = {
    1,1,1,1,1,1,1,1,1
  };
  Type_FITS table_types[nfields] = {
    XLONG,
    XDOUBLE, XDOUBLE,
    XDOUBLE,
    XDOUBLE,
    XLONG,
    XFLOAT,
    XDOUBLE,
    XINT
  };
  std::string table_units[nfields] =  {
    "None",
    "pixels", "pixels",
    "ADU",
    "ADU^2",
    "None",
    "mags",
    "Arcsec",
    "None"
  };
  fits.CreateBinaryTable(size(),nfields,
      table_cols,table_nelem,table_types,table_units);

  // Write the header keywords
  WriteParKey("version");
  WriteParKey("noise_method");
  WriteParKey("dist_method");

  WriteParKey("stars_minsize");
  WriteParKey("stars_maxsize");
  WriteParKey("stars_minmag");
  WriteParKey("stars_maxmag");
  WriteParKey("stars_ndivx");
  WriteParKey("stars_ndivy");

  WriteParKey("stars_startn1");
  WriteParKey("stars_starfrac");
  WriteParKey("stars_magstep1");
  WriteParKey("stars_miniter1");
  WriteParKey("stars_reject1");
  WriteParKey("stars_binsize1");
  WriteParKey("stars_maxratio1");
  WriteParKey("stars_okvalcount");
  WriteParKey("stars_maxrms");
  WriteParKey("stars_starsperbin");

  WriteParKey("stars_fitorder");
  WriteParKey("stars_fitsigclip");
  WriteParKey("stars_startn2");
  WriteParKey("stars_magstep2");
  WriteParKey("stars_miniter2");
  WriteParKey("stars_minbinsize");
  WriteParKey("stars_reject2");

  WriteParKey("stars_purityratio");
  WriteParKey("stars_maxrefititer");

  // Write the data
  fits.WriteColumn(XLONG, 1, 1, 1, size(), &id[0]);
  std::vector<double> x(pos.size());
  std::vector<double> y(pos.size());
  for(size_t i=0;i<pos.size();i++) {
    x[i] = pos[i].GetX();
    y[i] = pos[i].GetY();
  }
  fits.WriteColumn(XDOUBLE, 2, 1, 1, size(), &x[0]);
  fits.WriteColumn(XDOUBLE, 3, 1, 1, size(), &y[0]);
  fits.WriteColumn(XDOUBLE, 4, 1, 1, size(), &sky[0]);
  fits.WriteColumn(XDOUBLE, 5, 1, 1, size(), &noise[0]);
  fits.WriteColumn(XLONG, 6, 1, 1, size(), &flags[0]);
  fits.WriteColumn(XFLOAT, 7, 1, 1, size(), &mag[0]);
  fits.WriteColumn(XDOUBLE, 8, 1, 1, size(), &objsize[0]);
  fits.WriteColumn(XINT, 9, 1, 1, size(), &isastar[0]);
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

  Form hexform; hexform.hex().trail(0);

  for (size_t i=0; i<size(); i++) {
    fout 
      << id[i] << delim
      << pos[i].GetX() << delim
      << pos[i].GetY() << delim
      << sky[i] << delim
      << noise[i] << delim
      << hexform(flags[i]) << delim
      << mag[i] << delim
      << objsize[i] << delim
      << isastar[i] << std::endl;
  }
}

void StarCatalog::Write() const 
{
  std::vector<std::string> files = MultiName(params, "stars");  

  for(size_t i=0; i<files.size(); ++i) {
    const std::string& file = files[i];
    dbg<<"Writing star catalog to file: "<<file<<std::endl;

    bool fitsio = false;
    if (params.keyExists("stars_io")) {
      std::vector<std::string> ios = params["stars_io"];
      Assert(ios.size() == files.size());
      fitsio = (ios[i] == "FITS");
    }
    else if (file.find("fits") != std::string::npos) 
      fitsio = true;

    if (fitsio) {
      //WriteFits(file);
      WriteFitsCCfits(file);
    } else {
      std::string delim = "  ";
      if (params.keyExists("stars_delim")) {
	std::vector<std::string> delims = params["stars_delim"];
	Assert(delims.size() == files.size());
	delim = delims[i];
      }
      else if (file.find("csv") != std::string::npos) delim = ",";
      WriteAscii(file,delim);
    }
  }
  dbg<<"Done Write StarCatalog\n";
}

void StarCatalog::ReadFits(std::string file) 
{
  int hdu = params.read("stars_hdu",2);

  FitsFile fits(file);

  dbg<<"Moving to HDU #"<<hdu<<std::endl;
  fits.GotoHDU(hdu);

  long nrows=0;
  fits.ReadKey("NAXIS2", XLONG, &nrows);

  dbg<<"  nrows = "<<nrows<<std::endl;
  if (nrows <= 0) {
    throw std::runtime_error("nrows must be >= 0");
  }

  std::string id_col=params.get("stars_id_col");
  std::string x_col=params.get("stars_x_col");
  std::string y_col=params.get("stars_y_col");
  std::string sky_col=params.get("stars_sky_col");
  std::string noise_col=params.get("stars_noise_col");
  std::string flags_col=params.get("stars_flags_col");
  std::string mag_col=params.get("stars_mag_col");
  std::string objsize_col=params.get("stars_objsize_col");
  std::string isastar_col=params.get("stars_isastar_col");

  dbg<<"Reading columns"<<std::endl;
  dbg<<"  "<<id_col<<std::endl;
  id.resize(nrows);
  fits.ReadScalarCol(id_col,XLONG,&id[0], nrows);

  dbg<<"  "<<x_col<<"  "<<y_col<<std::endl;
  pos.resize(nrows);
  std::vector<double> x(nrows);
  std::vector<double> y(nrows);
  fits.ReadScalarCol(x_col,XDOUBLE,&x[0], nrows);
  fits.ReadScalarCol(y_col,XDOUBLE,&y[0], nrows);
  for(long i=0;i<nrows;++i) pos[i] = Position(x[i],y[i]);

  dbg<<"  "<<sky_col<<std::endl;
  sky.resize(nrows);
  fits.ReadScalarCol(sky_col,XDOUBLE,&sky[0], nrows);

  dbg<<"  "<<noise_col<<std::endl;
  noise.resize(nrows);
  fits.ReadScalarCol(noise_col,XDOUBLE,&noise[0], nrows);

  dbg<<"  "<<flags_col<<std::endl;
  flags.resize(nrows);
  fits.ReadScalarCol(flags_col,XLONG,&flags[0], nrows);

  dbg<<"  "<<mag_col<<std::endl;
  mag.resize(nrows);
  fits.ReadScalarCol(mag_col,XFLOAT,&mag[0], nrows);

  dbg<<"  "<<objsize_col<<std::endl;
  objsize.resize(nrows);
  fits.ReadScalarCol(objsize_col,XDOUBLE,&objsize[0], nrows);

  dbg<<"  "<<isastar_col<<std::endl;
  isastar.resize(nrows);
  fits.ReadScalarCol(isastar_col,XINT,&isastar[0], nrows);
}

void StarCatalog::ReadAscii(std::string file, std::string delim)
{
  std::ifstream fin(file.c_str());
  if (!fin) {
    throw std::runtime_error("Error opening stars file");
  }

  id.clear(); pos.clear(); sky.clear(); noise.clear(); flags.clear();
  mag.clear(); objsize.clear(); isastar.clear();

  if (delim == "  ") {
    ConvertibleString flag;
    long id1,star;
    double x,y,sky1,n,s,m;
    while ( fin >> id1 >> x >> y >> sky1 >> n 
	>> flag >> m >> s >> star)
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
    dbg<<"Reading with delimeter "<<d<<std::endl;
    double x,y;
    ConvertibleString temp;
    while (getline(fin,temp,d))
    {
      dbg<<"First elemnt in line = "<<temp<<std::endl;
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
      dbg<<"Last elemnt in line = "<<temp<<std::endl;
    }
    dbg<<"nlines = "<<id.size()<<std::endl;
  }
}

void StarCatalog::Read()
{
  std::string file = Name(params,"stars",false,true);
  // false,true = input_prefix=false, mustexist=true.
  // It is an input here, but it is in the output_prefix directory.
  dbg<< "Reading star catalog from file: " << file << std::endl;

  bool fitsio = false;
  if (params.keyExists("stars_io")) 
    fitsio = (params["stars_io"] == "FITS");
  else if (file.find("fits") != std::string::npos) 
    fitsio = true;

  if (fitsio) {
    ReadFits(file);
  } else {
    std::string delim = "  ";
    if (params.keyExists("stars_delim")) delim = params["stars_delim"];
    else if (file.find("csv") != std::string::npos) delim = ",";
    ReadAscii(file,delim);
  }
  dbg<<"Done Read StarCatalog\n";
}

