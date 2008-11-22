
#include "ConfigFile.h"
#include "DoMeasure.h"
#include "TMV.h"
#include "dbg.h"
#include "Name.h"

#include "fitsio.h"
#include "FitsFile.h"
#include "StarFinder.h"

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

std::ostream* dbgout = 0;
bool XDEBUG = false;


using std::cout;
using std::cerr;
using std::endl;

typedef struct {
  vector<long> id;

  vector<Position> pos;
  vector<float> x;
  vector<float> y;
  vector<float> mag;
  vector<float> mag_err;
  vector<int> flags;

  // calculated here
  vector<double> sigma;
  vector<int> size_flags;
  vector<int> star_flag;

  // These unused
  vector<double> sky;
  vector<double> noise;
} SXCAT_STRUCT;

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
  cat.sky.resize(n,0);
  cat.noise.resize(n,0);
}

int ReadFitsCat(ConfigFile& params, SXCAT_STRUCT& cat)
{

  int cat_hdu = 1;
  if (params.keyExists("cat_hdu")) cat_hdu = params["cat_hdu"];
  std::string file = Name(params,"cat",true);
  int minrows = 100;
  if (params.keyExists("sx_minrows")) minrows = params["sx_minrows"];

  cout<< "Reading cat from file: " << file << endl;

  FitsFile fits(file);

  int hdutype;

  dbg<<"Moving to HDU #"<<cat_hdu<<endl;
  fits.GotoHDU(cat_hdu);

  // These are not very portable.
  int nfound;
  long nrows=0;
  fits.ReadKey("NAXIS2", TLONG, (char *)&nrows);

  cout<<"  nrows = "<<nrows<<endl;
  if (nrows <= 0) {
    std::cerr<<"nrows must be >= 0"<<std::endl;
    return FAILURE_FORMAT_ERROR;
  }

  // MJ: <100 have basically no chance to find the stars
  if (nrows <= minrows) {
    std::cerr<<"Too few rows in catalog"<<std::endl;
    return FAILURE_FORMAT_ERROR;
  }

  // Allocate memory for the columns we will read
  ResizeSXCat(cat, nrows);

  std::string id_name=params["sx_id_name"];
  std::string x_name=params["sx_x_name"];
  std::string y_name=params["sx_y_name"];
  std::string mag_name=params["sx_mag_name"];
  std::string mag_err_name=params["sx_mag_err_name"];
  std::string flags_name=params["sx_flags_name"];

  std::cout<<"Reading columns"<<std::endl;
  fits.ReadScalarCol((char *)id_name.c_str(),TLONG,(char *)&cat.id[0], nrows);
  std::cout<<"  "<<id_name<<std::endl;
  fits.ReadScalarCol((char *)x_name.c_str(),TFLOAT,(char *)&cat.x[0], nrows);
  std::cout<<"  "<<x_name<<std::endl;
  fits.ReadScalarCol((char *)y_name.c_str(),TFLOAT,(char *)&cat.y[0], nrows);
  std::cout<<"  "<<y_name<<std::endl;
  fits.ReadScalarCol((char *)mag_name.c_str(),TFLOAT,(char *)&cat.mag[0], nrows);
  std::cout<<"  "<<mag_name<<std::endl;
  fits.ReadScalarCol((char *)mag_err_name.c_str(),TFLOAT,(char *)&cat.mag_err[0], nrows);
  std::cout<<"  "<<mag_err_name<<std::endl;
  fits.ReadScalarCol((char *)flags_name.c_str(),TSHORT,(char *)&cat.flags[0], nrows);
  std::cout<<"  "<<flags_name<<std::endl;

  for (long i=0; i< nrows; i++) {
    cat.pos[i] = Position(cat.x[i], cat.y[i]);
  }

  fits.Close();

}



int DoMeasureSizes(ConfigFile& params)
{
  // Load image:
  int image_hdu = 1;
  if (params.keyExists("image_hdu")) image_hdu = params["image_hdu"];
  Image<double> im(Name(params,"image",true),image_hdu);

  // load weight image
  int weight_hdu = 1;
  if (params.keyExists("weight_hdu")) weight_hdu = params["weight_hdu"];
  Image<double> weight_im(Name(params,"weight",true),weight_hdu);

  SXCAT_STRUCT sxcat;
  ReadFitsCat(params, sxcat);

  double psfap = double(params["psf_aperture"]); 

  // Read distortion function
  Transformation trans;
  if (params.keyExists("dist_ext") || params.keyExists("dist_file")) {
    std::string distfile = Name(params,"dist",true);
    std::ifstream distin(distfile.c_str());
    Assert(distin);
    distin >> trans;
    xdbg<<"Done read distortion "<<Name(params,"dist",true)<<std::endl;
  } // else stay with identity transformation.

  // we are using weight image, so sk, noise, gain are dummy
  MeasureSigmas(
      im, 
      sxcat.pos, sxcat.sky, sxcat.noise, 0.0, 
      weight_im, 
      trans, 
      psfap, 
      sxcat.sigma,
      sxcat.size_flags);

  // Eventually separate this out
  std::cout<<"Loading Star Finder"<<std::endl;
  StarFinder sf(params);

  
  vector<int> starflags; 
  sf.RunFindStars(
      sxcat.flags,
      sxcat.size_flags,
      sxcat.x,
      sxcat.y,
      sxcat.sigma,
      sxcat.mag,
      sxcat.star_flag);


  std::string output_file = Name(params, "allcat");

  cout<<"Writing to allcat file: "<<output_file<<endl;
  std::ofstream output(output_file.c_str());

  std::string delim = params["allcat_delim"];
  for (long i=0; i<sxcat.pos.size(); i++) {
    output
      <<sxcat.id[i]<<delim
      <<sxcat.x[i]<<delim
      <<sxcat.y[i]<<delim
      <<sxcat.sigma[i]<<delim
      <<sxcat.mag[i]<<delim
      <<sxcat.size_flags[i]<<delim
      <<sxcat.star_flag[i]
      <<endl;
  }


}


int main(int argc, char **argv) 
#ifndef NOTHROW
  try 
#endif
{
  // Read parameters
  if (argc < 2) {
    std::cerr<<"Usage: measure_sizes configfile [param=value ...]\n";
    std::cerr<<"\tThe first parameter is the configuration file that has \n";
    std::cerr<<"\tall the parameters for this run. \n";
    std::cerr<<"\tThese values may be modified on the command line by \n";
    std::cerr<<"\tentering param/value pais as param=value. \n";
    std::cerr<<"Note: root is not usuallly given in the parameter file, \n";
    std::cerr<<"\tso the normal command line would be something like:\n";
    std::cerr<<"\tmeasurepsf measurepsf.config root=img123\n";
    return EXIT_FAILURE;
  }
  ConfigFile params(argv[1]);
  for(int k=2;k<argc;k++) params.Append(argv[k]);

  std::string logfile = "measure_sizes.log";
  if (params.keyExists("log_file") || params.keyExists("log_ext")) {
    logfile = Name(params,"log");
  }

  std::string logdelim = "  ";
  if (params.keyExists("log_delim")) logdelim = params["log_delim"];
  PSFLog log(logfile,logdelim); 
  // This automatically writes its output when it goes out of scope
  // whether that is naturally in after an exception is caught.
  // Log output is:  (all on one line)
  //    exitcode  nstars  nsuccess 
  //    nfail_range_distortion
  //    nfail_edge  nfail_npix<10
  //    nfail_tmv_error  nfail_other_error
  //    nfail_psf_measurement


#ifndef NOTHROW
  try 
#endif
  {
    // Setup debugging
    if (params.keyExists("verbose") && int(params["verbose"]) > 0) {
      if (params.keyExists("debug_file") || params.keyExists("debug_ext")) {
	dbgout = new std::ofstream(Name(params,"debug").c_str());
      }
      else dbgout = &std::cout;
      if (int(params["verbose"]) > 1) XDEBUG = true;
    }

#ifdef _OPENMP
    if (params.keyExists("num_threads")) {
      int num_threads = params["num_threads"];
      omp_set_num_threads(num_threads);
    }
#endif

    dbg<<"Config params = \n"<<params<<std::endl;

    DoMeasureSizes(params);

    if (dbgout && dbgout != &std::cout) {delete dbgout; dbgout=0;}
    return EXIT_SUCCESS;
  }
#ifndef NOTHROW
#if 0
  // Change to 1 to let gdb see where the program bombed out.
  catch(int) {}
#else
  catch (file_not_found& e)
  {
    std::cerr<<"Caught \n"<<e.what()<<std::endl;
    log.exitcode = FAILURE_FILE_NOT_FOUND;
    return EXIT_FAILURE; // = 1 typically
  }
  catch (ConfigFile::file_not_found& e)
  {
    std::cerr<<"Caught \n"<<e.what()<<std::endl;
    log.exitcode = FAILURE_CONFIGFILE_ERROR;
    return EXIT_FAILURE;
  }
  catch (ConfigFile::key_not_found& e)
  {
    std::cerr<<"Caught \n"<<e.what()<<std::endl;
    log.exitcode = FAILURE_CONFIGFILE_ERROR;
    return EXIT_FAILURE;
  }
  catch (tmv::Error& e)
  {
    std::cerr<<"Caught \n"<<e<<std::endl;
    log.exitcode = FAILURE_TMV_ERROR;
    return EXIT_FAILURE;
  }
  catch (std::exception& e)
  {
    std::cerr<<"Caught \n"<<e.what()<<std::endl;
    log.exitcode = FAILURE_STD_EXCEPTION;
    return EXIT_FAILURE;
  }
  catch (...)
  {
    log.exitcode = FAILURE;
    std::cerr<<"Caught Unknown error\n";
    return EXIT_FAILURE;
  }
#endif
#endif
}
#ifndef NOTHROW
catch (std::exception& e)
{
  std::cerr<<"Caught \n"<<e.what()<<std::endl;
  std::cerr<<"outside of the normal try..catch block.";
  std::cerr<<"Unable to write to the log file.\n";
  return EXIT_FAILURE;
}
catch (...)
{
  std::cerr<<"Cought an exception outside of the normal try..catch block.";
  std::cerr<<"Unable to write to the log file.\n";
  return EXIT_FAILURE;
}
#endif
