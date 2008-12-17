
#include "DoMeasure.h"
#include "Pixel.h"
#include "dbg.h"
#include "Name.h"
#include <fstream>
#include <iostream>
#include "fitsio.h"

enum NOISE_METHOD {VALUE, CATALOG, CATALOG_SIGMA, GAIN_VALUE, GAIN_FITS,
    WEIGHT_IMAGE};

static void ReadGain(const std::string& fitsname, ConfigFile& params)
{
  xdbg<<"stars ReadGain: "<<fitsname<<std::endl;
  fitsfile *fitsptr;
  int status=0;
  float gain, rdnoise;

  std::cerr<<"Opening file: "<<fitsname.c_str()<<std::endl;
  if (fits_open_file(&fitsptr,fitsname.c_str(),READONLY,&status))
    std::cerr<<"fits open file gave status #: "<<status<<std::endl;
  std::vector<std::string> gain_key = params["image_gain_key"];
  std::vector<std::string> rdn_key = params["image_readnoise_key"];

  for(size_t k=0;k<gain_key.size();k++) {
    xdbg<<"try "<<gain_key[k]<<std::endl;
    if (!fits_read_key(fitsptr,TFLOAT,(char*)(gain_key[k].c_str()),&gain,
	NULL,&(status=0)))
      break;
    xdbg<<"status = "<<status<<std::endl;
  }
  xdbg<<"done status = "<<status<<std::endl;
  if (status)
    std::cerr<<"fits read gain gave status #: "<<status<<std::endl;

  for(size_t k=0;k<rdn_key.size();k++) {
    xdbg<<"try "<<rdn_key[k]<<std::endl;
    if (!fits_read_key(fitsptr,TFLOAT,(char*)(rdn_key[k].c_str()),&rdnoise,
	  NULL,&(status=0)))
      break;
    xdbg<<"status = "<<status<<std::endl;
  }
  xdbg<<"done status = "<<status<<std::endl;
  if (status)
    std::cerr<<"fits read readnoise gave status #: "<<status<<std::endl;

  if (fits_close_file(fitsptr,&status))
    std::cerr<<"fits close file gave status #: "<<status<<std::endl;
  if (status) {
    char errmsg[80];
    fits_get_errstatus(status,errmsg);
    std::cerr<<"fits error status "<<status<<" = "<<errmsg<<std::endl;
    while (fits_read_errmsg(errmsg))
      std::cerr<<"fits error: "<<errmsg<<std::endl;
    throw(std::runtime_error("reading gain/rdnoise"));
  }

  params["image_gain"] = gain;
  params["image_readnoise"] = rdnoise;
}

#define DEF_BUFFER_SIZE 500

static void GetTokens(ConfigFile& params, std::string incat, std::string line,
    std::vector<ConvertibleString>& tokens)
{
  std::istringstream linein(line);
  if (params.keyExists(incat+"_delim")) { 
    //std::cerr<<"delim exists"<<std::endl;
    char delim = params[incat+"_delim"];
    int bufsize = DEF_BUFFER_SIZE;
    if (params.keyExists("csv_bufsize")) bufsize = params["csv_bufsize"];
    std::auto_ptr<char> temp(new char[bufsize]);
    while (linein.getline(temp.get(),bufsize,delim)) {
      //std::cerr<<"Found token: "<<temp<<std::endl;
      tokens.push_back(std::string(temp.get()));
    }
  } else {
    std::string temp;
    while (linein >> temp) tokens.push_back(temp);
  }
}

void ReadCatalog(ConfigFile& params, std::string incat,
    std::vector<Position>& all_pos, std::vector<double>& all_sky,
    std::vector<double>& all_noise, double& gain, Image<double>*& weight_im,
    std::vector<Position>& all_skypos)
{
  // Setup noise calculation:
  Assert(params.keyExists("noise_method"));
  NOISE_METHOD nm;
  double noise_value=0.;
  gain=0.;
  double readnoise=0.;
  size_t i_noise=0;
  int weight_hdu = 1;
  
  if (params["noise_method"] == "VALUE") {
    nm = VALUE;
    Assert(params.keyExists("noise"));
    noise_value = params["noise"];
  } 
  else if (params["noise_method"] == "CATALOG") {
    nm = CATALOG;
    Assert(params.keyExists("cat_noise_num"));
    i_noise = params["cat_noise_num"];
  }
  else if (params["noise_method"] == "CATALOG_SIGMA") {
    nm = CATALOG_SIGMA;
    Assert(params.keyExists("cat_noise_num"));
    i_noise = params["cat_noise_num"];
  }
  else if (params["noise_method"] == "GAIN_VALUE") {
    nm = GAIN_VALUE;
    Assert(params.keyExists("image_gain"));
    Assert(params.keyExists("image_readnoise"));
    gain = params["image_gain"];
    readnoise = params["image_readnoise"];
  } 
  else if (params["noise_method"] == "GAIN_FITS") {
    ReadGain(Name(params,"image",true),params);
    xdbg<<"Read gain = "<<params["image_gain"]<<
      ", rdn = "<<params["image_readnoise"]<<std::endl;
    nm = GAIN_VALUE;
    gain = params["image_gain"];
    readnoise = params["image_readnoise"];
  } 
  else if (params["noise_method"] == "WEIGHT_IMAGE") {
    nm = WEIGHT_IMAGE;
    Assert(params.keyExists("weight_ext") || params.keyExists("weight_file"));
    if (params.keyExists("weight_hdu")) weight_hdu = params["weight_hdu"];
  }
  else {
    throw std::runtime_error("Unknown noise method");
  }
  double extrasky=0.;
  if (params.keyExists("image_extra_sky")) 
    extrasky = params["image_extra_sky"];

  // Read input catalog:
  std::string incatfile = Name(params,incat,false,true);
  std::ifstream catin(incatfile.c_str());
  Assert(catin);
  xdbg<<"Opened catalog "<<incatfile<<std::endl;
  std::string line;
  Assert(params.keyExists("cat_x_num"));
  Assert(params.keyExists("cat_y_num"));
  size_t i_x = params["cat_x_num"];
  size_t i_y = params["cat_y_num"];
  size_t i_flags = 0;
  if (params.keyExists("cat_flags_num")) 
    i_flags = params["cat_flags_num"];
  long ignore_flags = -1;
  if (params.keyExists("cat_ignore_flags"))
    ignore_flags = params["cat_ignore_flags"];
  else if (params.keyExists("cat_ok_flags")) {
    ignore_flags = params["cat_ok_flags"];
    ignore_flags = ~ignore_flags;
  }

  std::vector<std::string> comment_marker;
  if (params.keyExists("cat_comment_marker")) 
    comment_marker = params["cat_comment_marker"];
  else comment_marker.push_back("#");
  size_t i_sky = 0;
  if (params.keyExists("cat_local_sky_num")) 
    i_sky = params["cat_local_sky_num"];
  size_t i_ra = 0;
  if (params.keyExists("cat_ra_num")) i_ra = params["cat_ra_num"];
  size_t i_dec = 0;
  if (params.keyExists("cat_dec_num")) i_dec = params["cat_dec_num"];
  double x_offset = 0.;
  double y_offset = 0.;
  if (params.keyExists("cat_x_offset")) x_offset = params["cat_x_offset"];
  if (params.keyExists("cat_y_offset")) y_offset = params["cat_y_offset"];
  while (getline(catin,line)) {
    bool skip = false;
    for(size_t k=0;k<comment_marker.size();++k) 
      if (std::string(line,0,comment_marker[k].size()) == comment_marker[k])
	skip = true;
    if (skip) continue;
    std::vector<ConvertibleString> tokens;
    GetTokens(params,incat,line,tokens);
    xdbg<<"tokens size = "<<tokens.size()<<std::endl;
    Assert(i_x <= tokens.size());
    Assert(i_y <= tokens.size());
    double x = tokens[i_x-1];
    x -= x_offset;
    double y = tokens[i_y-1];
    y -= y_offset;
    xdbg<<"x,y = "<<x<<','<<y<<std::endl;

    if (i_flags) {
      long flags = 0;
      Assert(i_flags <= tokens.size());
      flags = tokens[i_flags-1];
      if (flags & ignore_flags) continue;
    } 

    double sky = 0.;
    if (i_sky) {
      Assert(i_sky <= tokens.size());
      sky = tokens[i_sky-1];
    } 
    double ra = 0.;
    double dec = 0.;
    if (i_ra && i_dec) {
      Assert(i_ra <= tokens.size());
      Assert(i_dec <= tokens.size());
      ra = tokens[i_ra-1];
      dec = tokens[i_dec-1];
    } 

    double noise = 0.;
    if (nm == VALUE) {
      noise = noise_value;
    } else if (nm == CATALOG) {
      Assert(i_noise <= tokens.size());
      noise = tokens[i_noise-1];
    } else if (nm == CATALOG_SIGMA) {
      Assert(i_noise <= tokens.size());
      double sigma = tokens[i_noise-1];
      noise = sigma*sigma;
    } else if (nm == GAIN_VALUE) {
      noise = extrasky/gain + readnoise;
    }
    all_pos.push_back(Position(x,y));
    all_noise.push_back(noise);
    if (i_sky) {
      all_sky.push_back(sky);
    }
    if (i_ra && i_dec) {
      all_skypos.push_back(Position(ra,dec));
    }
  }

  // Read weight image if appropriate
  if (nm == WEIGHT_IMAGE) {
    dbg<<"Before open weight image: "<<Name(params,"weight",true)<<std::endl;
    weight_im = new Image<double>(Name(params,"weight",true),weight_hdu);
    dbg<<"Opened weight image.\n";
  }
}

