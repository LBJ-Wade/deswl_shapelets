#include "SXCat.h"

void ResizeSXCat(SXCAT_STRUCT& cat, long n)
{
  cat.id.resize(n,0);
  cat.x.resize(n,0);
  cat.y.resize(n,0);
  cat.pos.resize(n);
  cat.mag.resize(n,0);
  cat.mag_err.resize(n,0);
  cat.flags.resize(n,0);

  cat.sigma.resize(n,0);
  cat.size_flags.resize(n,0);
  cat.star_flag.resize(n,0);

  // Not used
  cat.local_sky.resize(n,0);
  cat.noise.resize(n,0);
}

void ReadSXCat(ConfigFile& params, SXCAT_STRUCT& cat)
{

  int cat_hdu = 1;
  if (params.keyExists("cat_hdu")) cat_hdu = params["cat_hdu"];
  std::string file = Name(params,"cat",true);
  int minrows = 100;
  if (params.keyExists("sx_minrows")) minrows = params["sx_minrows"];

  std::cerr<< "Reading cat from file: " << file << std::endl;

  FitsFile fits(file);

  //int hdutype;

  dbg<<"Moving to HDU #"<<cat_hdu<<std::endl;
  fits.GotoHDU(cat_hdu);

  // These are not very portable.
  //int nfound;
  long nrows=0;
  fits.ReadKey("NAXIS2", TLONG, (char *)&nrows);

  std::cerr<<"  nrows = "<<nrows<<std::endl;
  if (nrows <= 0) {
    std::cerr<<"nrows must be >= 0"<<std::endl;
    throw FAILURE_FORMAT_ERROR;
  }

  // MJ: <100 have basically no chance to find the stars
  if (nrows <= minrows) {
    std::cerr<<"Too few rows in catalog"<<std::endl;
    throw FAILURE_FORMAT_ERROR;
  }

  // Allocate memory for the columns we will read
  ResizeSXCat(cat, nrows);

  std::string id_name=params["sx_id_name"];
  std::string x_name=params["sx_x_name"];
  std::string y_name=params["sx_y_name"];
  std::string local_sky_name=params["sx_local_sky_name"];
  std::string mag_name=params["sx_mag_name"];
  std::string mag_err_name=params["sx_mag_err_name"];
  std::string flags_name=params["sx_flags_name"];

  std::cerr<<"Reading columns"<<std::endl;
  std::cerr<<"  "<<id_name<<std::endl;
  fits.ReadScalarCol((char *)id_name.c_str(),TLONG,(char *)&cat.id[0], nrows);
  std::cerr<<"  "<<x_name<<std::endl;
  fits.ReadScalarCol((char *)x_name.c_str(),TFLOAT,(char *)&cat.x[0], nrows);
  std::cerr<<"  "<<y_name<<std::endl;
  fits.ReadScalarCol((char *)y_name.c_str(),TFLOAT,(char *)&cat.y[0], nrows);
  std::cerr<<"  "<<local_sky_name<<std::endl;
  fits.ReadScalarCol((char *)local_sky_name.c_str(),TDOUBLE,(char *)&cat.local_sky[0], nrows);
  std::cerr<<"  "<<mag_name<<std::endl;
  fits.ReadScalarCol((char *)mag_name.c_str(),TFLOAT,(char *)&cat.mag[0], nrows);
  std::cerr<<"  "<<mag_err_name<<std::endl;
  fits.ReadScalarCol((char *)mag_err_name.c_str(),TFLOAT,(char *)&cat.mag_err[0], nrows);
  std::cerr<<"  "<<flags_name<<std::endl;
  fits.ReadScalarCol((char *)flags_name.c_str(),TSHORT,(char *)&cat.flags[0], nrows);

  for (long i=0; i< nrows; i++) {
    cat.pos[i] = Position(cat.x[i], cat.y[i]);
  }

  fits.Close();
}


