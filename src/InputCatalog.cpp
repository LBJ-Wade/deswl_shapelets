
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
  xdbg<<"stars ReadGain: "<<file<<std::endl;
  FitsFile fits(file);
  float gain, rdnoise;

  std::vector<std::string> gain_key = params["image_gain_key"];
  std::vector<std::string> rdn_key = params["image_readnoise_key"];

  for(size_t k=0;k<gain_key.size();k++) {
    xdbg<<"try "<<gain_key[k]<<std::endl;
    try {
      fits.ReadKey(gain_key[k], TFLOAT, &gain);
      break;
    } 
    catch (FitsException)
    { if (k == gain_key.size()-1) throw; }
  }

  for(size_t k=0;k<rdn_key.size();k++) {
    xdbg<<"try "<<rdn_key[k]<<std::endl;
    try {
      fits.ReadKey(rdn_key[k], TFLOAT, &rdnoise);
      break;
    } 
    catch (FitsException)
    { if (k == rdn_key.size()-1) throw; }
  }

  params["image_gain"] = gain;
  params["image_readnoise"] = rdnoise;
}

InputCatalog::InputCatalog(ConfigFile& _params, std::string key_prefix) : 
  params(_params), prefix(key_prefix)
{
  // Setup noise calculation:
  Assert(params.keyExists("noise_method"));
  NOISE_METHOD nm;
  double noise_value=0.;
  double gain=0.;
  double readnoise=0.;
  int weight_hdu = 1;
  
  if (params["noise_method"] == "VALUE") {
    nm = VALUE;
    Assert(params.keyExists("noise"));
    noise_value = params["noise"];
  } 
  else if (params["noise_method"] == "CATALOG") {
    nm = CATALOG;
    Assert(params.keyExists(prefix + "noise_col"));
  }
  else if (params["noise_method"] == "CATALOG_SIGMA") {
    nm = CATALOG_SIGMA;
    Assert(params.keyExists(prefix + "noise_col"));
  }
  else if (params["noise_method"] == "GAIN_VALUE") {
    nm = GAIN_VALUE;
    Assert(params.keyExists("image_gain"));
    Assert(params.keyExists("image_readnoise"));
    gain = params["image_gain"];
    readnoise = params["image_readnoise"];
  } 
  else if (params["noise_method"] == "GAIN_FITS") {
    ReadGain(Name(params,"image",true),_params);
    xdbg<<"Read gain = "<<params["image_gain"]<<
      ", rdn = "<<params["image_readnoise"]<<std::endl;
    nm = GAIN_VALUE;
    gain = params["image_gain"];
    readnoise = params["image_readnoise"];
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
  if (sky.size() == 0) {
    double glob_sky = 0.;
    if (params.keyExists("image_sky")) glob_sky = params["image_sky"];
    else {
      // This is probably inefficient, since image is probably already
      // allocated, but I don't want to require the InputCatalog
      // constructor to have an Image open already, so I don't want to
      // pass it as a parameter.  
      // Anyway, hopefully this will be a rare thing to need to do, so
      // it won't be a problem.
      int image_hdu = 1;
      if (params.keyExists("image_hdu")) image_hdu = params["image_hdu"];
      Image<double> im(Name(params,"image",true),image_hdu);
      xdbg<<"Opened image "<<Name(params,"image",true)<<std::endl;
      glob_sky = im.Median();
      dbg<<"Found global sky from image median value.\n";
    }
    dbg<<"Set global value of sky to "<<glob_sky<<std::endl;
    sky.resize(id.size());
    std::fill(sky.begin(),sky.end(),glob_sky);
  }
  Assert(sky.size() == pos.size());

  // MJ: <100 have basically no chance to find the stars
  int nrows = id.size();
  int minrows = 0;
  if (params.keyExists(prefix + "minrows")) 
    minrows = params[prefix + "minrows"];
  if (nrows <= minrows) {
    std::cout<<"STATUS3BEG Warning: Input catalog only has "
      <<nrows<<" rows for Name="<<Name(params,"cat",true)
      <<". STATUS3END"<<std::endl;
  }

  // Update noise calculation if necessary
  if (nm == VALUE) {
    for(size_t i=0;i<noise.size();++i) {
      noise[i] = noise_value;
    }
  }
  else if (nm == CATALOG_SIGMA) {
    for(size_t i=0;i<noise.size();++i) {
      noise[i] = noise[i]*noise[i];
    }
  }
  else if (nm == GAIN_VALUE) {
    double extrasky=0.;
    if (params.keyExists("image_extra_sky")) 
      extrasky = params["image_extra_sky"];
    for(size_t i=0;i<noise.size();++i) {
      noise[i] = (sky[i]+extrasky)/gain + readnoise;
    }
  }
  else if (nm == WEIGHT_IMAGE) {
    // Then we don't need the noise vector, but it is easier
    // to just fill it with 0's.
    noise.resize(nrows,0);
  }

  // Convert input flags into our flag schema
  if (flags.size() == 0) {
    flags.resize(id.size(),0);
  } else {
    long ignore_flags = ~0;
    if (params.keyExists(prefix + "ignore_flags"))
      ignore_flags = params[prefix + "ignore_flags"];
    else if (params.keyExists(prefix + "ok_flags")) {
      ignore_flags = params[prefix + "ok_flags"];
      ignore_flags = ~ignore_flags;
    }
    Assert(flags.size() == id.size());
    for(size_t i=0;i<flags.size();++i) {
      flags[i] = (flags[i] & ignore_flags) ? INPUT_FLAG : 0;
    }
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
  // Read the file using a different routine depending on if it is
  // a fits file or ascii:
  std::string file = Name(params,"cat",true);
  if (file.find("fits") != std::string::npos) {
    ReadFits(file);
  } else {
    std::string delim = "  ";
    if (params.keyExists(prefix+"delim")) 
      delim = params[prefix+"delim"];
    ReadAscii(file,delim);
  }
}

void InputCatalog::ReadFits(std::string file)
{
  dbg<< "Reading cat from FITS file: " << file << std::endl;

  FitsFile fits(file);

  int hdu = 1;
  if (params.keyExists(prefix + "hdu")) hdu = params[prefix + "hdu"];
  dbg<<"Moving to HDU #"<<hdu<<std::endl;
  fits.GotoHDU(hdu);

  // These are not very portable.
  long nrows=0;
  fits.ReadKey("NAXIS2", TLONG, &nrows);
  dbg<<"  nrows = "<<nrows<<std::endl;
  if (nrows <= 0) {
    std::cerr<<"nrows must be > 0"<<std::endl;
    throw FAILURE_FORMAT_ERROR;
  }

  // Read each column in turn:
  dbg<<"Reading columns"<<std::endl;

  // ID
  id.resize(nrows,0);
  if (params.keyExists(prefix + "id_col")) {
    std::string id_col=params[prefix + "id_col"];
    dbg<<"  "<<id_col<<std::endl;
    fits.ReadScalarCol(id_col,TLONG,&id[0], nrows);
  } else {
    for (size_t i=0;i<id.size();i++) id[i] = i+1;
  }

  // Position (on chip)
  pos.resize(nrows);
  std::string x_col=params.get(prefix + "x_col");
  std::string y_col=params.get(prefix + "y_col");
  std::vector<float> pos_x(nrows);
  std::vector<float> pos_y(nrows);
  dbg<<"  "<<x_col<<std::endl;
  fits.ReadScalarCol(x_col,TFLOAT,&pos_x[0], nrows);
  dbg<<"  "<<y_col<<std::endl;
  fits.ReadScalarCol(y_col,TFLOAT,&pos_y[0], nrows);
  double x_offset = 0.;
  double y_offset = 0.;
  if (params.keyExists(prefix + "x_offset")) 
    x_offset = params[prefix + "x_offset"];
  if (params.keyExists(prefix + "y_offset")) 
    y_offset = params[prefix + "y_offset"];
  for (long i=0; i< nrows; i++) {
    pos[i] = Position(pos_x[i]-x_offset, pos_y[i]-y_offset);
    xdbg<<"pos["<<i<<"] = "<<pos[i]<<std::endl;
  }

  // Local sky
  if (params.keyExists(prefix + "sky_col")) {
    sky.resize(nrows,0);
    std::string sky_col=params[prefix + "sky_col"];
    dbg<<"  "<<sky_col<<std::endl;
    fits.ReadScalarCol(sky_col,TDOUBLE,&sky[0], nrows);
  }

  // Magnitude
  if (params.keyExists(prefix + "mag_col")) {
    mag.resize(nrows,0);
    std::string mag_col=params[prefix + "mag_col"];
    dbg<<"  "<<mag_col<<std::endl;
    fits.ReadScalarCol(mag_col,TFLOAT,&mag[0], nrows);
  }

  // Magnitude Error
  if (params.keyExists(prefix + "mag_err_col")) {
    mag_err.resize(nrows,0);
    std::string mag_err_col=params[prefix + "mag_err_col"];
    dbg<<"  "<<mag_err_col<<std::endl;
    fits.ReadScalarCol(mag_err_col,TFLOAT,&mag_err[0], nrows);
  }

  // Size
  if (params.keyExists(prefix + "size_col")) {
    objsize.resize(nrows,0);
    std::string size_col=params[prefix + "size_col"];
    dbg<<"  "<<size_col<<std::endl;
    fits.ReadScalarCol(size_col,TDOUBLE,&objsize[0], nrows);
    if (params.keyExists(prefix + "size2_col")) {
      std::string size2_col=params[prefix + "size2_col"];
      dbg<<"  "<<size2_col<<std::endl;
      std::vector<double> objsize2(objsize.size());
      fits.ReadScalarCol(size2_col,TDOUBLE,&objsize2[0], nrows);
      for(size_t i=0;i<objsize.size();++i) objsize[i] += objsize2[i];
    }
  }

  // Error flags
  if (params.keyExists(prefix + "flags_col")) {
    flags.resize(nrows,0);
    std::string flags_col=params[prefix + "flags_col"];
    dbg<<"  "<<flags_col<<std::endl;
    fits.ReadScalarCol(flags_col,TSHORT,&flags[0], nrows);
  }

  // RA
  if (params.keyExists(prefix + "ra_col")) {
    ra.resize(nrows,0);
    std::string ra_col=params[prefix + "ra_col"];
    dbg<<"  "<<ra_col<<std::endl;
    fits.ReadScalarCol(ra_col,TFLOAT,&ra[0], nrows);
  }

  // Declination
  if (params.keyExists(prefix + "dec_col")) {
    dec.resize(nrows,0);
    std::string dec_col=params[prefix + "dec_col"];
    dbg<<"  "<<dec_col<<std::endl;
    fits.ReadScalarCol(dec_col,TFLOAT,&dec[0], nrows);
  }

  // Noise
  if (params.keyExists(prefix + "noise_col")) {
    noise.resize(nrows,0);
    std::string noise_col=params[prefix + "noise_col"];
    dbg<<"  "<<noise_col<<std::endl;
    fits.ReadScalarCol(noise_col,TDOUBLE,&noise[0], nrows);
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
  Assert(params.keyExists(prefix + "x_col"));
  Assert(params.keyExists(prefix + "y_col"));
  size_t x_col = params[prefix + "x_col"];
  size_t y_col = params[prefix + "y_col"];

  // Set column numbers for optional columns
  size_t id_col = 0;
  if (params.keyExists(prefix + "id_col")) id_col = params[prefix + "id_col"];

  size_t mag_col = 0;
  if (params.keyExists(prefix + "mag_col")) 
    mag_col = params[prefix + "mag_col"];
  size_t mag_err_col = 0;
  if (params.keyExists(prefix + "mag_err_col")) 
    mag_err_col = params[prefix + "mag_err_col"];

  size_t size_col = 0;
  if (params.keyExists(prefix + "size_col")) 
    size_col = params[prefix + "size_col"];
  size_t size2_col = 0;
  if (params.keyExists(prefix + "size2_col")) 
    size2_col = params[prefix + "size2_col"];

  size_t flags_col = 0;
  if (params.keyExists(prefix + "flags_col")) 
    flags_col = params[prefix + "flags_col"];

  size_t sky_col = 0;
  if (params.keyExists(prefix + "sky_col")) 
    sky_col = params[prefix + "sky_col"];

  size_t ra_col = 0;
  if (params.keyExists(prefix + "ra_col")) 
    ra_col = params[prefix + "ra_col"];
  size_t dec_col = 0;
  if (params.keyExists(prefix + "dec_col")) 
    dec_col = params[prefix + "dec_col"];

  size_t noise_col = 0;
  if (params.keyExists(prefix + "noise_col")) 
    noise_col = params[prefix + "noise_col"];

  // Set up allowed comment markers
  std::vector<std::string> comment_marker;
  if (params.keyExists(prefix + "comment_marker")) 
    comment_marker = params[prefix + "comment_marker"];
  else comment_marker.push_back("#");

  // Allow for offset from given x,y values
  double x_offset = 0.;
  double y_offset = 0.;
  if (params.keyExists(prefix + "x_offset")) 
    x_offset = params[prefix + "x_offset"];
  if (params.keyExists(prefix + "y_offset")) 
    y_offset = params[prefix + "y_offset"];

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
    xdbg<<"tokens size = "<<tokens.size()<<std::endl;

    // ID
    if (id_col) {
      Assert(id_col <= tokens.size());
      id_val = tokens[id_col-1];
    } else {
      // if not reading id, then just increment to get sequential values
      ++id_val;
    }
    id.push_back(id_val);

    // Position
    Assert(x_col <= tokens.size());
    Assert(y_col <= tokens.size());
    double x = tokens[x_col-1];
    double y = tokens[y_col-1];
    pos.push_back(Position(x-x_offset,y-y_offset));

    // Sky
    if (sky_col) {
      // Note: if sky not read in, then sky.size() is still 0
      // This is indicator to update with global given or median value later.
      double sky_val = 0.;
      Assert(sky_col <= tokens.size());
      sky_val = tokens[sky_col-1];
      sky.push_back(sky_val);
    } 

    // Magnitude
    if (mag_col) {
      double mag_val = 0.;
      Assert(mag_col <= tokens.size());
      mag_val = tokens[mag_col-1];
      mag.push_back(mag_val);
    } 

    // Magnitude error
    if (mag_err_col) {
      double mag_err_val = 0.;
      Assert(mag_err_col <= tokens.size());
      mag_err_val = tokens[mag_err_col-1];
      mag_err.push_back(mag_err_val);
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
    } 

    // Flags
    if (flags_col) {
      long flags_val = 0;
      Assert(flags_col <= tokens.size());
      flags_val = tokens[flags_col-1];
      flags.push_back(flags_val);
    } 

    // RA
    if (ra_col) {
      double ra_val = 0.;
      Assert(ra_col <= tokens.size());
      ra_val = tokens[ra_col-1];
      ra.push_back(ra_val);
    } 

    // Declination
    if (dec_col) {
      double dec_val = 0.;
      Assert(dec_col <= tokens.size());
      dec_val = tokens[dec_col-1];
      dec.push_back(dec_val);
    }

    // Noise
    if (noise_col) {
      double noise_val = 0.;
      Assert(noise_col <= tokens.size());
      noise_val = tokens[noise_col-1];
      noise.push_back(noise_val);
    }
  }
}

