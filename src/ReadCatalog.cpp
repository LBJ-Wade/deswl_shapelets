
#include "DoMeasure.h"
#include "Pixel.h"
#include "dbg.h"
#include "Name.h"
#include <fstream>
#include <iostream>
#include "fitsio.h"

enum NOISE_METHOD {VALUE, CATALOG, CATALOG_SIGMA, GAIN_VALUE, GAIN_FITS,
    WEIGHT_IMAGE};

void ReadGain(const std::string& fitsname, ConfigFile& params)
{
  xdbg<<"stars ReadGain: "<<fitsname<<std::endl;
  fitsfile *fitsptr;
  int status=0;
  float gain, rdnoise;

  std::cout<<"Opening file: "<<fitsname.c_str()<<std::endl;
  if (fits_open_file(&fitsptr,fitsname.c_str(),READONLY,&status))
    std::cerr<<"fits open file gave status #: "<<status<<std::endl;
  std::vector<std::string> gain_key = params["gain_key"];
  std::vector<std::string> rdn_key = params["readnoise_key"];

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
#ifdef NOTHROW
    exit(1);
#else
    throw(std::runtime_error(("reading gain/rdnoise")));
#endif
  }

  params["gain"] = gain;
  params["readnoise"] = rdnoise;
}

#define DEF_BUFFER_SIZE 500

void GetTokens(ConfigFile& params, std::string line,
    std::vector<ConvertibleString>& tokens)
{
  std::istringstream linein(line);
  if (params.keyExists("delim")) { 
    //std::cout<<"delim exists"<<std::endl;
    char delim = params["delim"];
    int bufsize = DEF_BUFFER_SIZE;
    if (params.keyExists("bufsize")) bufsize = params["bufsize"];
    char temp[bufsize];
    while (linein.getline(temp,bufsize,delim)) {
      //std::cout<<"Found token: "<<temp<<std::endl;
      tokens.push_back(std::string(temp));
    }
  } else {
    std::string temp;
    while (linein >> temp) tokens.push_back(temp);
  }
}

void ReadCatalog(ConfigFile& params, std::string incat,
    std::vector<Position>& all_pos, std::vector<double>& all_sky,
    std::vector<double>& all_noise, double& gain, Image<double>*& weight_im)
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
    Assert(params.keyExists("i_noise"));
    i_noise = params["i_noise"];
  }
  else if (params["noise_method"] == "CATALOG_SIGMA") {
    nm = CATALOG_SIGMA;
    Assert(params.keyExists("i_noise"));
    i_noise = params["i_noise"];
  }
  else if (params["noise_method"] == "GAIN_VALUE") {
    nm = GAIN_VALUE;
    Assert(params.keyExists("gain"));
    Assert(params.keyExists("readnoise"));
    gain = params["gain"];
    readnoise = params["readnoise"];
  } 
  else if (params["noise_method"] == "GAIN_FITS") {
    ReadGain(Name(params,"image"),params);
    xdbg<<"Read gain = "<<params["gain"]<<", rdn = "<<params["readnoise"]<<std::endl;
    nm = GAIN_VALUE;
    gain = params["gain"];
    readnoise = params["readnoise"];
  } 
  else if (params["noise_method"] = "WEIGHT_IMAGE") {
    nm = WEIGHT_IMAGE;
    Assert(params.keyExists("weight_ext") || params.keyExists("weight_name"));
    if (params.keyExists("weight_hdu")) weight_hdu = params["weight_hdu"];
  }
  else {
#ifdef NOTHROW
    std::cerr<<"Unknown noise method.\n"; exit(1);
#else
    throw std::runtime_error("Unknown noise method");
#endif
  }
  double extrasky=0.;
  if (params.keyExists("extra_sky")) extrasky = params["extra_sky"];

  // Read input catalog:
  std::string incatfile = Name(params,incat);
  std::ifstream catin(incatfile.c_str());
  Assert(catin);
  std::string line;
  Assert(params.keyExists("i_x"));
  Assert(params.keyExists("i_y"));
  size_t i_x = params["i_x"];
  size_t i_y = params["i_y"];
  size_t i_errcode = 0;
  if (params.keyExists("i_errcode")) i_errcode = params["i_errcode"];
  std::vector<std::string> comment_marker;
  if (params.keyExists("comment_marker")) 
    comment_marker = params["comment_marker"];
  else comment_marker.push_back("#");
  size_t i_sky = 0;
  if (params.keyExists("i_sky")) i_sky = params["i_sky"];
  double x_offset = 0.;
  double y_offset = 0.;
  if (params.keyExists("x_offset")) x_offset = params["x_offset"];
  if (params.keyExists("y_offset")) y_offset = params["y_offset"];
  while (getline(catin,line)) {
    bool skip = false;
    for(size_t k=0;k<comment_marker.size();++k) 
      if (std::string(line,0,comment_marker[k].size()) == comment_marker[k])
	skip = true;
    if (skip) continue;
    std::vector<ConvertibleString> tokens;
    GetTokens(params,line,tokens);
    xdbg<<"tokens size = "<<tokens.size()<<std::endl;
    Assert(i_x <= tokens.size());
    Assert(i_y <= tokens.size());
    double x = tokens[i_x-1];
    x -= x_offset;
    double y = tokens[i_y-1];
    y -= y_offset;
    xdbg<<"x,y = "<<x<<','<<y<<std::endl;

    if (i_errcode) {
      int errcode = 0;
      Assert(i_errcode <= tokens.size());
      errcode = tokens[i_errcode-1];
      if (errcode != 0) continue;
    } 

    double sky = 0.;
    if (i_sky) {
      Assert(i_sky <= tokens.size());
      sky = tokens[i_sky-1];
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
    all_sky.push_back(sky);
  }

  // Read weight image if appropriate
  if (nm == WEIGHT_IMAGE) {
    // Do mask image too?
    weight_im = new Image<double>(Name(params,"weight"),weight_hdu);
  }
}

