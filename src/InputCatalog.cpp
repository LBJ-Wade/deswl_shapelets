
#include "InputCatalog.h"
#include "dbg.h"
#include "Name.h"
#include "Params.h"
#include "Image.h"
#include "FitsFile.h"

enum NOISE_METHOD {VALUE, CATALOG, CATALOG_SIGMA, GAIN_VALUE, GAIN_FITS,
    WEIGHT_IMAGE};

static void ReadGain(const std::string& file, ConfigFile& params)
{
  xdbg<<"ReadGain: from fits file "<<file<<std::endl;
  FitsFile fits(file);
  float gain, rdnoise;

  std::vector<std::string> gain_key = params.get("image_gain_key");
  std::vector<std::string> rdn_key = params.get("image_readnoise_key");

  for(size_t k=0;k<gain_key.size();k++) {
    xdbg<<"try "<<gain_key[k]<<std::endl;
    try {
      fits.ReadKey(gain_key[k], XFLOAT, &gain);
      break;
    } 
    catch (FitsException)
    { if (k == gain_key.size()-1) throw; }
  }

  for(size_t k=0;k<rdn_key.size();k++) {
    xdbg<<"try "<<rdn_key[k]<<std::endl;
    try {
      fits.ReadKey(rdn_key[k], XFLOAT, &rdnoise);
      break;
    } 
    catch (FitsException)
    { if (k == rdn_key.size()-1) throw; }
  }

  params["image_gain"] = gain;
  params["image_readnoise"] = rdnoise;
}

InputCatalog::InputCatalog(ConfigFile& _params, const Image<double>* im) : 
  params(_params)
{
  // Setup noise calculation:
  Assert(params.keyExists("noise_method"));
  NOISE_METHOD nm;
  double noise_value=0.;
  double gain=0.;
  double readnoise=0.;
  
  if (params["noise_method"] == "VALUE") {
    nm = VALUE;
    noise_value = params.get("noise");
  } 
  else if (params["noise_method"] == "CATALOG") {
    nm = CATALOG;
    Assert(params.keyExists("cat_noise_col"));
  }
  else if (params["noise_method"] == "CATALOG_SIGMA") {
    nm = CATALOG_SIGMA;
    Assert(params.keyExists("cat_noise_col"));
  }
  else if (params["noise_method"] == "GAIN_VALUE") {
    nm = GAIN_VALUE;
    gain = params.get("image_gain");
    readnoise = params.get("image_readnoise");
    dbg<<"gain, readnoise = "<<gain<<"  "<<readnoise<<std::endl;
  } 
  else if (params["noise_method"] == "GAIN_FITS") {
    ReadGain(Name(params,"image",true),_params);
    xdbg<<"Read gain = "<<params["image_gain"]<<
      ", rdn = "<<params["image_readnoise"]<<std::endl;
    gain = params.get("image_gain");
    readnoise = params.get("image_readnoise");
    nm = GAIN_VALUE;
    dbg<<"gain, readnoise = "<<gain<<"  "<<readnoise<<std::endl;
  } 
  else if (params["noise_method"] == "WEIGHT_IMAGE") {
    nm = WEIGHT_IMAGE;
    Assert(params.keyExists("weight_ext") || params.keyExists("weight_file"));
  }
  else {
    throw std::runtime_error("Unknown noise method");
  }

  // Read catalog from fits or ascii file as appropriate
  Read();

  // These are the only two fields guaranteed to be set.
  Assert(id.size() == pos.size());

  // Fix sky if necessary
  double bad_sky_val = params.read("cat_bad_sky",-999.);
  if (sky.size() == 0) sky.resize(id.size(),bad_sky_val);
  if (std::find(sky.begin(), sky.end(), bad_sky_val) != sky.end()) {
    double glob_sky = 0.;
    if (params.keyExists("cat_globalsky")) {
      glob_sky = params["cat_globalsky"];
    }
    else {
      std::auto_ptr<Image<double> > im1;
      if (!im) {
	im1.reset(new Image<double>(params));
	im = im1.get();
      }
      glob_sky = im->Median();
      dbg<<"Found global sky from image median value.\n";
    }
    dbg<<"Set global value of sky to "<<glob_sky<<std::endl;
    for(size_t i=0;i<sky.size();++i) 
      if (sky[i] == bad_sky_val) sky[i] = glob_sky;
  }
  Assert(sky.size() == pos.size());

  // MJ: <100 have basically no chance to find the stars
  int nrows = id.size();
  if (params.read("des_qa",false)) {
    int minrows = params.read("cat_nrows",0);
    if (nrows <= minrows) {
      std::cout<<"STATUS3BEG Warning: Input catalog only has "
	<<nrows<<" rows for Name="<<Name(params,"cat",true)
	<<". STATUS3END"<<std::endl;
    }
  }

  // Update noise calculation if necessary
  Assert(nm == VALUE || nm == CATALOG || nm == CATALOG_SIGMA ||
      nm == GAIN_VALUE || nm == WEIGHT_IMAGE);
  if (nm == VALUE) {
    noise.resize(nrows,0);
    for(size_t i=0;i<noise.size();++i) {
      noise[i] = noise_value;
    }
    xdbg<<"Set all noise to "<<noise_value<<std::endl;
  }
  else if (nm == CATALOG) {
    Assert(noise.size() == id.size());
  }
  else if (nm == CATALOG_SIGMA) {
    Assert(noise.size() == id.size());
    for(size_t i=0;i<noise.size();++i) {
      noise[i] = noise[i]*noise[i];
    }
    xdbg<<"Squared noise values from catalog\n";
  }
  else if (nm == GAIN_VALUE) {
    noise.resize(nrows,0);
    double extrasky=params.read("image_extra_sky",0.);
    xdbg<<"Calculate noise from sky, gain, readnoise\n";
    for(size_t i=0;i<noise.size();++i) {
      noise[i] = (sky[i]+extrasky)/gain + readnoise;
      xdbg<<"("<<sky[i]<<" + "<<extrasky<<")/"<<gain<<" + "<<readnoise<<" = "<<noise[i]<<std::endl;
      xdbg<<"ID="<<id[i]<<" NOISE="<<noise[i]<<std::endl;
    }
  }
  else if (nm == WEIGHT_IMAGE) {
    // Then we don't need the noise vector, but it is easier
    // to just fill it with 0's.
    noise.resize(nrows,0);
    xdbg<<"Set noise values to 0, since using weight_im, not noise vector.\n";
  }

  // Convert input flags into our flag schema
  if (flags.size() == 0) {
    dbg<<"No flags read in -- starting all with 0\n";
    flags.resize(id.size(),0);
  } else {
    long ignore_flags = ~0L;
    dbg<<std::hex<<std::showbase;
    if (params.keyExists("cat_ignore_flags")) {
      ignore_flags = params["cat_ignore_flags"];
      dbg<<"Using ignore flag parameter = "<<ignore_flags<<std::endl;
    }
    else if (params.keyExists("cat_ok_flags")) {
      ignore_flags = params["cat_ok_flags"];
      dbg<<"Using ok flag parameter = "<<ignore_flags<<std::endl;
      ignore_flags = ~ignore_flags;
      dbg<<"ignore flag = "<<ignore_flags<<std::endl;
    } else {
      dbg<<"No ok or ignore parameter: use ignore flag = "<<
	ignore_flags<<std::endl;
    }
    Assert(flags.size() == id.size());
    for(size_t i=0;i<flags.size();++i) {
      xdbg<<"flags[i] = "<<flags[i];
      flags[i] = (flags[i] & ignore_flags) ? INPUT_FLAG : 0;
      xdbg<<" => "<<flags[i]<<std::endl;
    }
    dbg<<std::dec<<std::noshowbase;
  }

  // At this point the vectors that are guaranteed to be filled are:
  // id, pos, sky, noise, flags
  Assert(pos.size() == id.size());
  Assert(sky.size() == id.size());
  Assert(noise.size() == id.size());
  Assert(flags.size() == id.size());
}

void InputCatalog::Read()
{
  std::string file = Name(params,"cat",true);
  // true = mustexist=true.
  dbg<< "Reading input cat from file: " << file << std::endl;

  bool fitsio = false;
  if (params.keyExists("cat_io")) 
    fitsio = (params["cat_io"] == "FITS");
  else if (file.find("fits") != std::string::npos) 
    fitsio = true;

  if (fitsio) {
    ReadFits(file);
  } else {
    std::string delim = "  ";
    if (params.keyExists("cat_delim")) delim = params["cat_delim"];
    else if (file.find("csv") != std::string::npos) delim = ",";
    ReadAscii(file,delim);
  }
  dbg<<"Done Read InputCatalog\n";
}

void InputCatalog::ReadFits(std::string file)
{
  dbg<< "Reading cat from FITS file: " << file << std::endl;

  FitsFile fits(file);

  int hdu = params.read("cat_hdu" , 1);
  dbg<<"Moving to HDU #"<<hdu<<std::endl;
  fits.GotoHDU(hdu);

  // These are not very portable.
  long nrows=0;
  fits.ReadKey("NAXIS2", XLONG, &nrows);
  dbg<<"  nrows = "<<nrows<<std::endl;
  if (nrows <= 0) {
    std::cerr<<"nrows must be > 0"<<std::endl;
    throw FAILURE_FORMAT_ERROR;
  }

  // Read each column in turn:
  dbg<<"Reading columns"<<std::endl;

  // ID
  id.resize(nrows,0);
  if (params.keyExists("cat_id_col")) {
    std::string id_col=params["cat_id_col"];
    dbg<<"  "<<id_col<<std::endl;
    fits.ReadScalarCol(id_col,XLONG,&id[0], nrows);
  } else {
    for (size_t i=0;i<id.size();i++) id[i] = i+1;
  }

  // Position (on chip)
  pos.resize(nrows);
  std::string x_col=params.get("cat_x_col");
  std::string y_col=params.get("cat_y_col");
  std::vector<float> pos_x(nrows);
  std::vector<float> pos_y(nrows);
  dbg<<"  "<<x_col<<std::endl;
  fits.ReadScalarCol(x_col,XFLOAT,&pos_x[0], nrows);
  dbg<<"  "<<y_col<<std::endl;
  fits.ReadScalarCol(y_col,XFLOAT,&pos_y[0], nrows);
  double x_offset = params.read("cat_x_offset",0.);
  double y_offset = params.read("cat_y_offset",0.);
  for (long i=0; i< nrows; i++) {
    pos[i] = Position(pos_x[i]-x_offset, pos_y[i]-y_offset);
    xdbg<<"pos["<<i<<"] = "<<pos[i]<<std::endl;
  }

  // Local sky
  if (params.keyExists("cat_sky_col")) {
    sky.resize(nrows,0);
    std::string sky_col=params["cat_sky_col"];
    dbg<<"  "<<sky_col<<std::endl;
    fits.ReadScalarCol(sky_col,XDOUBLE,&sky[0], nrows);
  }

  // Magnitude
  if (params.keyExists("cat_mag_col")) {
    mag.resize(nrows,0);
    std::string mag_col=params["cat_mag_col"];
    dbg<<"  "<<mag_col<<std::endl;
    fits.ReadScalarCol(mag_col,XFLOAT,&mag[0], nrows);
  }

  // Magnitude Error
  if (params.keyExists("cat_mag_err_col")) {
    mag_err.resize(nrows,0);
    std::string mag_err_col=params["cat_mag_err_col"];
    dbg<<"  "<<mag_err_col<<std::endl;
    fits.ReadScalarCol(mag_err_col,XFLOAT,&mag_err[0], nrows);
  }

  // Size
  if (params.keyExists("cat_size_col")) {
    objsize.resize(nrows,0);
    std::string size_col=params["cat_size_col"];
    dbg<<"  "<<size_col<<std::endl;
    fits.ReadScalarCol(size_col,XDOUBLE,&objsize[0], nrows);
    if (params.keyExists("cat_size2_col")) {
      std::string size2_col=params["cat_size2_col"];
      dbg<<"  "<<size2_col<<std::endl;
      std::vector<double> objsize2(objsize.size());
      fits.ReadScalarCol(size2_col,XDOUBLE,&objsize2[0], nrows);
      for(size_t i=0;i<objsize.size();++i) objsize[i] += objsize2[i];
    }
  }

  // Error flags
  if (params.keyExists("cat_flag_col")) {
    flags.resize(nrows,0);
    std::string flag_col=params["cat_flag_col"];
    dbg<<"  "<<flag_col<<std::endl;
    fits.ReadScalarCol(flag_col,XLONG,&flags[0], nrows);
  }

  // RA
  if (params.keyExists("cat_ra_col")) {
    ra.resize(nrows,0);
    std::string ra_col=params["cat_ra_col"];
    dbg<<"  "<<ra_col<<std::endl;
    fits.ReadScalarCol(ra_col,XFLOAT,&ra[0], nrows);
  }

  // Declination
  if (params.keyExists("cat_dec_col")) {
    dec.resize(nrows,0);
    std::string dec_col=params["cat_dec_col"];
    dbg<<"  "<<dec_col<<std::endl;
    fits.ReadScalarCol(dec_col,XFLOAT,&dec[0], nrows);
  }

  // Noise
  if (params.keyExists("cat_noise_col")) {
    noise.resize(nrows,0);
    std::string noise_col=params["cat_noise_col"];
    dbg<<"  "<<noise_col<<std::endl;
    fits.ReadScalarCol(noise_col,XDOUBLE,&noise[0], nrows);
  }
}

// Helper functino to convert ASCII line into vector of tokens
static void GetTokens(std::string delim, 
    std::string line, std::vector<ConvertibleString>& tokens)
{
  std::istringstream linein(line);
  if (delim == "  ") {
    std::string temp;
    while (linein >> temp) tokens.push_back(temp);
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
    std::string temp;
    while (getline(linein,temp,d)) {
      //std::cerr<<"Found token: "<<temp<<std::endl;
      tokens.push_back(temp);
    }
  }
}

void InputCatalog::ReadAscii(std::string file, std::string delim)
{
  std::ifstream catin(file.c_str());
  if (!catin) {
    throw std::runtime_error("Error opening input catalog file");
  }
  xdbg<<"Opened catalog "<<file<<std::endl;

  // x,y is required
  std::string line;
  size_t x_col = params.get("cat_x_col");
  size_t y_col = params.get("cat_y_col");

  // Set column numbers for optional columns
  size_t id_col = params.read("cat_id_col",0);

  size_t mag_col = params.read("cat_mag_col",0);
  size_t mag_err_col = params.read("cat_mag_err_col",0);

  size_t size_col = params.read("cat_size_col",0);
  size_t size2_col = params.read("cat_size2_col",0);

  size_t flag_col = params.read("cat_flag_col",0);

  size_t sky_col = params.read("cat_sky_col",0);

  size_t ra_col = params.read("cat_ra_col",0);
  size_t dec_col = params.read("cat_dec_col",0);

  size_t noise_col = params.read("cat_noise_col",0);

  // Set up allowed comment markers
  std::vector<std::string> comment_marker = 
    params.read("cat_comment_marker",std::vector<std::string>(1,"#"));

  // Allow for offset from given x,y values
  double x_offset = params.read("cat_x_offset",0.);
  double y_offset = params.read("cat_y_offset",0.);

  // Keep running id value when id_col = 0
  int id_val = 0;

  // Read each line from catalog file
  while (getline(catin,line)) {

    // Skip if this is a comment.
    bool skip = false;
    for(size_t k=0;k<comment_marker.size();++k) 
      if (std::string(line,0,comment_marker[k].size()) == comment_marker[k])
	skip = true;
    if (skip) continue;

    // Convert line into vector of tokens
    std::vector<ConvertibleString> tokens;
    GetTokens(delim,line,tokens);

    // ID
    if (id_col) {
      Assert(id_col <= tokens.size());
      id_val = tokens[id_col-1];
    } else {
      // if not reading id, then just increment to get sequential values
      ++id_val;
    }
    id.push_back(id_val);
    xdbg<<"ID="<<id_val<<"  ";

    // Position
    Assert(x_col <= tokens.size());
    Assert(y_col <= tokens.size());
    double x = tokens[x_col-1];
    double y = tokens[y_col-1];
    pos.push_back(Position(x-x_offset,y-y_offset));
    xdbg<<"POS=("<<pos.back()<<")  ";

    // Sky
    if (sky_col) {
      // Note: if sky not read in, then sky.size() is still 0
      // This is indicator to update with global given or median value later.
      double sky_val = 0.;
      Assert(sky_col <= tokens.size());
      sky_val = tokens[sky_col-1];
      sky.push_back(sky_val);
      xdbg<<"SKY="<<sky_val<<"  ";
    } 

    // Magnitude
    if (mag_col) {
      double mag_val = 0.;
      Assert(mag_col <= tokens.size());
      mag_val = tokens[mag_col-1];
      mag.push_back(mag_val);
      xdbg<<"MAG="<<mag_val<<"  ";
    } 

    // Magnitude error
    if (mag_err_col) {
      double mag_err_val = 0.;
      Assert(mag_err_col <= tokens.size());
      mag_err_val = tokens[mag_err_col-1];
      mag_err.push_back(mag_err_val);
      xdbg<<"MAG_ERR="<<mag_err_val<<"  ";
    } 

    // Size
    if (size_col) {
      double size_val = 0.;
      Assert(size_col <= tokens.size());
      size_val = tokens[size_col-1];
      if (size2_col) {
	Assert(size2_col <= tokens.size());
	size_val += double(tokens[size2_col-1]);
      } 
      objsize.push_back(size_val);
      xdbg<<"SIZE="<<size_val<<"  ";
    } 

    // Flags
    if (flag_col) {
      long flag_val = 0;
      Assert(flag_col <= tokens.size());
      flag_val = tokens[flag_col-1];
      flags.push_back(flag_val);
      xdbg<<"FLAG="<<flag_val<<"  ";
    } 

    // RA
    if (ra_col) {
      double ra_val = 0.;
      Assert(ra_col <= tokens.size());
      ra_val = tokens[ra_col-1];
      ra.push_back(ra_val);
      xdbg<<"RA="<<ra_val<<"  ";
    } 

    // Declination
    if (dec_col) {
      double dec_val = 0.;
      Assert(dec_col <= tokens.size());
      dec_val = tokens[dec_col-1];
      dec.push_back(dec_val);
      xdbg<<"DEC="<<dec_val<<"  ";
    }

    // Noise
    if (noise_col) {
      double noise_val = 0.;
      Assert(noise_col <= tokens.size());
      noise_val = tokens[noise_col-1];
      noise.push_back(noise_val);
      xdbg<<"NOISE="<<noise_val<<"  ";
    }
    xdbg<<std::endl;
  }
}
