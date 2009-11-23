
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <valarray>
#include "TMV.h"
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
#include "WriteParKey.h"

static void CalcSigma1(
    double& sigma,
    const Image<double>& im, const Position& pos, double sky,
    double noise, double gain, const Image<double>* weight_im, 
    const Transformation& trans, double psfap, long& flag)
{
  std::vector<PixelList> pix(1);
  long flag1 = 0;
  try {
    GetPixList(im, pix[0], pos, sky, noise, gain, weight_im, trans, psfap, 
	flag1);
  } catch (Range_error& e) {
    flag1 |= TRANSFORM_EXCEPTION;
  }
  if (flag1) {
    flag |= flag1;
    sigma = DEFVALNEG;
    return;
  }
  xxdbg<<"npix = "<<pix[0].size()<<std::endl;

  Ellipse ell;
  ell.PeakCentroid(pix[0],psfap/3.);
  ell.CrudeMeasure(pix[0],sigma);
  xdbg<<"Crude Measure: centroid = "<<ell.GetCen();
  xdbg<<", mu = "<<ell.GetMu()<<std::endl;
  if (ell.Measure(pix,2,sigma,true,flag1)) { // true means use integ first
    xdbg<<"Successful 2nd order measure.\n";
    xdbg<<"mu = "<<ell.GetMu()<<std::endl;
  } else {
    flag |= flag1;
    xdbg<<"Ellipse measure returned flag "<<flag1<<std::endl;
    sigma = DEFVALNEG;
    return;
  }

  double mu = ell.GetMu();
  sigma *= exp(mu);
  dbg<<"sigma = "<<sigma<<std::endl;
  Assert(sigma > 0);
}

void CalcSigma(
    double& sigma,
    const Image<double>& im, const Position& pos, double sky,
    double noise, double gain, const Image<double>* weight_im, 
    const Transformation& trans, double psfap, long& flag)
{
  try {
    CalcSigma1(
	sigma,
	im, pos, sky, noise, gain, weight_im,
	trans, psfap, flag);
    dbg<<"objsize: "<<sigma<<std::endl;
    dbg<<"flags: "<<flag<<std::endl;
  } catch (tmv::Error& e) {
    dbg<<"Caught: "<<e<<std::endl;
    sigma = DEFVALNEG;
    flag |= TMV_EXCEPTION;
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
		std::cout
			<<"STATUS3BEG Warning: Only "<<count
			<<" stars found for Name="<<Name(params,"stars")
			<<". STATUS3END"<<std::endl;
    }
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

  // Header Keywords
  std::string tmvvers = tmv::TMV_Version();
  std::string wlvers = WlVersion();

  table->addKey("tmvvers", tmvvers, "version of TMV code");
  table->addKey("wlvers", wlvers, "version of weak lensing code");


  // kind of kludgy but is more flexible than using type numbers or strings
  double dbl;
  int intgr;
  std::string str;


  // if serun= is sent we'll put it in the header.  This allows us to 
  // associate some more, possibly complicated, metadata with this file
  if ( params.keyExists("serun") ) {
	  CCfitsWriteParKey(params, table, "serun", str);
  }


  //CCfitsWriteParKey(params, table, "version", str);
  CCfitsWriteParKey(params, table, "noise_method", str);
  CCfitsWriteParKey(params, table, "dist_method", str);

  if (params.keyExists("stars_minsize")) 
  {
    CCfitsWriteParKey(params, table, "stars_minsize", dbl);
    CCfitsWriteParKey(params, table, "stars_maxsize", dbl);
    CCfitsWriteParKey(params, table, "stars_minmag", dbl);
    CCfitsWriteParKey(params, table, "stars_maxmag", dbl);
    CCfitsWriteParKey(params, table, "stars_ndivx", intgr);
    CCfitsWriteParKey(params, table, "stars_ndivy", intgr);

    CCfitsWriteParKey(params, table, "stars_startn1", dbl);
    CCfitsWriteParKey(params, table, "stars_starfrac", dbl);
    CCfitsWriteParKey(params, table, "stars_magstep1", dbl);
    CCfitsWriteParKey(params, table, "stars_miniter1", intgr);
    CCfitsWriteParKey(params, table, "stars_reject1", dbl);
    CCfitsWriteParKey(params, table, "stars_binsize1", dbl);
    CCfitsWriteParKey(params, table, "stars_maxratio1", dbl);
    CCfitsWriteParKey(params, table, "stars_okvalcount", intgr);
    CCfitsWriteParKey(params, table, "stars_maxrms", dbl);
    CCfitsWriteParKey(params, table, "stars_starsperbin", intgr);

    CCfitsWriteParKey(params, table, "stars_fitorder", intgr);
    CCfitsWriteParKey(params, table, "stars_fitsigclip", dbl);
    CCfitsWriteParKey(params, table, "stars_startn2", dbl);
    CCfitsWriteParKey(params, table, "stars_magstep2", dbl);
    CCfitsWriteParKey(params, table, "stars_miniter2", intgr);
    CCfitsWriteParKey(params, table, "stars_minbinsize", dbl);
    CCfitsWriteParKey(params, table, "stars_reject2", dbl);

    CCfitsWriteParKey(params, table, "stars_purityratio", dbl);
    CCfitsWriteParKey(params, table, "stars_maxrefititer", intgr);
  }


  // Now the data columns

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
    throw WriteError("Error opening stars file"+file);
  }

  //Form hexform; hexform.hex().trail(0);

  for (size_t i=0; i<size(); i++) {
    fout 
      << id[i] << delim
      << pos[i].GetX() << delim
      << pos[i].GetY() << delim
      << sky[i] << delim
      << noise[i] << delim
      //<< hexform(flags[i]) << delim
      << flags[i] << delim
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

    try 
    {
      if (fitsio) {
	WriteFits(file);
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
    catch (CCfits::FitsException& e)
    {
      throw WriteError("Error writing to "+file+" -- caught error\n" +
	  e.message());
    }
    catch (std::exception& e)
    { 
      throw WriteError("Error writing to "+file+" -- caught error\n" + 
	  e.what());
    }
    catch (...)
    { 
      throw WriteError("Error writing to "+file+" -- caught unknown error");
    }
  }
  dbg<<"Done Write StarCatalog\n";
}

void StarCatalog::ReadFits(std::string file) 
{
  int hdu = params.read("stars_hdu",2);

  dbg<<"Opening FITS file at hdu "<<hdu<<std::endl;
  // true means read all as part of the construction
  CCfits::FITS fits(file, CCfits::Read, hdu-1, true);

  CCfits::ExtHDU& table=fits.extension(hdu-1);

  long nrows=table.rows();

  dbg<<"  nrows = "<<nrows<<std::endl;
  if (nrows <= 0) {
    throw ReadError("StarCatalog found to have 0 rows.  Must have > 0 rows.");
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

  long start=1;
  long end=nrows;

  dbg<<"Reading columns"<<std::endl;
  // will copy these out to positions types
  std::vector<double> x(nrows);
  std::vector<double> y(nrows);

  dbg<<"  "<<id_col<<std::endl;
  table.column(id_col).read(id, start, end);
  dbg<<"  "<<x_col<<"  "<<y_col<<std::endl;
  table.column(x_col).read(x, start, end);
  table.column(y_col).read(y, start, end);
  dbg<<"  "<<sky_col<<std::endl;
  table.column(sky_col).read(sky, start, end);
  dbg<<"  "<<noise_col<<std::endl;
  table.column(noise_col).read(noise, start, end);
  dbg<<"  "<<flags_col<<std::endl;
  table.column(flags_col).read(flags, start, end);
  dbg<<"  "<<mag_col<<std::endl;
  table.column(mag_col).read(mag, start, end);
  dbg<<"  "<<objsize_col<<std::endl;
  table.column(objsize_col).read(objsize, start, end);
  dbg<<"  "<<isastar_col<<std::endl;
  table.column(isastar_col).read(isastar, start, end);

  pos.resize(nrows);
  for(long i=0;i<nrows;++i) {
    pos[i] = Position(x[i],y[i]);
  }

}

void StarCatalog::ReadAscii(std::string file, std::string delim)
{
  std::ifstream fin(file.c_str());
  if (!fin) {
    throw ReadError("Error opening stars file"+file);
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
      throw ParameterError("ReadAscii delimiter must be a single character");
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

  if (!FileExists(file))
  {
    throw FileNotFound(file);
  }
  try 
  {
    if (fitsio) {
      ReadFits(file);
    } else {
      std::string delim = "  ";
      if (params.keyExists("stars_delim")) delim = params["stars_delim"];
      else if (file.find("csv") != std::string::npos) delim = ",";
      ReadAscii(file,delim);
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
  dbg<<"Done Read StarCatalog\n";
}

