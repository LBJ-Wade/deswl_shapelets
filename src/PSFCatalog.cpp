
#include <sstream>

#include "PSFCatalog.h"
#include "FitsFile.h"
#include "dbg.h"
#include "Params.h"
#include "Pixel.h"
#include "Ellipse.h"
#include "Params.h"
#include "Name.h"

//#define SINGLESTAR 18
//#define NSTARS 10
//#define STARTAT 85
//#define ENDAT 95


static void MeasureSinglePSF(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans,
    double noise, double gain, const Image<double>* weight_im,
    double sigma_p, double psfap, int psforder, PSFLog& log,
    BVec& psf, double& nu, long& flag)
{
  try {
    // We don't need to save skypos.  We just want to catch the range
    // error here, so we don't need to worry about it for dudx, etc.
    Position skypos;
    trans.Transform(cen,skypos);
    dbg<<"skypos = "<<skypos<<std::endl;
  } catch (Range_error& e) {
    xdbg<<"skip: transformation range error: \n";
    xdbg<<"p = "<<cen<<", b = "<<e.b<<std::endl;
    log.nf_range++;
    flag |= TRANSFORM_EXCEPTION;
    return;
  }

  try {

    std::vector<std::vector<Pixel> > pix(1);
    GetPixList(im,pix[0],cen,sky,noise,gain,weight_im,trans,psfap,flag);

    int npix = pix[0].size();
    xdbg<<"npix = "<<npix<<std::endl;

    Ellipse ell;
    ell.FixGam();
    ell.FixMu();
    ell.PeakCentroid(pix[0],psfap/3.);
    ell.CrudeMeasure(pix[0],sigma_p);
    tmv::Matrix<double> cov(psf.size(),psf.size());

    if (ell.Measure(pix,psforder,sigma_p,false,0,&psf,&cov)) {
      xdbg<<"psf = "<<psf<<std::endl;
      nu = psf[0] / std::sqrt(cov(0,0));
      //nu = 1.;
      psf.Normalize();  // Divide by (0,0) element
      xdbg<<"Normalized psf: "<<psf<<std::endl;
      log.ns_psf++;
    }
    else {
      xdbg<<"Measurement failed\n";
      log.nf_psf++;
      nu = psf[0] / std::sqrt(cov(0,0));
      psf.Normalize();  // Divide by (0,0) element
      flag |= MEASURE_PSF_FAILED;
    }

  } catch (tmv::Error& e) {
    dbg<<"TMV Error thrown in MeasureSinglePSF\n";
    dbg<<e<<std::endl;
    log.nf_tmverror++;
    flag |= TMV_EXCEPTION;
  } catch (...) {
    dbg<<"unkown exception in MeasureSinglePSF\n";
    log.nf_othererror++;
    flag |= UNKNOWN_EXCEPTION;
  }
}

double PSFCatalog::EstimateSigma(const Image<double>& im,
    const Image<double>* weight_im, const Transformation& trans)
{
  // Initial sigma_p for shapelet measurements
  double sigma_p = 1.;
  if (params.keyExists("psf_seeing_est")) {
    double seeing = params["psf_seeing_est"];
    // seeing is given as FWHM
    // for a gaussian 0.5 = exp(-((FWHM/2)/sigma)^2/2)
    // FWHM/sigma = 2.*sqrt(2 ln(2)) = 2.35
    sigma_p = seeing / 2.35;
  }

  // Calculate a good value of sigma to use:
  // (For this calculation, psfap is psf_aperture * 1 arcsec.)
  double gain = 0.;
  if (params.keyExists("gain")) gain = params["gain"];
  double psfap = params["psf_aperture"];
  dbg<<"psfap = "<<psfap<<std::endl;

  // Current value of sigma_p is taken as initial guess
  std::vector<double> sigmas(pos.size(),sigma_p);

  // Make new flags vector, rather than propagate flags to catalog.
  std::vector<long> flags(pos.size(),0);

  int nstars = pos.size();
  double meanmu = 0.;
  int count = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(guided) reduction(+ : meanmu) reduction(+ : count)
#endif
  for (int i=0; i<nstars; i++) if (!flags[i]) {
    dbg<<"use i = "<<i<<std::endl;

    try {
      double sigma = sigma_p;
      long flag1 = 0; // Ignore flags set by CalcSigma
      CalcSigma(
	  sigma,
	  im, pos[i], sky[i], noise[i], gain, weight_im, 
	  trans, psfap, flag1);
      meanmu += log(sigma);
      count++;
    } catch (...) {
      // Ignore errors -- just don't add to meanmu
    }
  } // End omp parallel for

  if (count < nstars/3) {
    std::ostringstream msgout;
    msgout<<"Too many objects were rejected. \n";
    msgout<<"nstars = "<<nstars<<", but only "<<count<<" successful measurements.\n";
    std::string msg = msgout.str();
    dbg<<msg<<std::endl;
    throw std::runtime_error(msg);
  }
  meanmu /= count;
  xdbg<<"meanmu = "<<meanmu<<std::endl;
  xdbg<<"input sigma_p = "<<sigma_p<<std::endl;
  sigma_p = exp(meanmu);
  dbg<<"sigma_p = "<<sigma_p<<std::endl;
  return sigma_p;
}

int PSFCatalog::MeasurePSF(const Image<double>& im,
    const Image<double>* weight_im,
    const Transformation& trans, double sigma_p, PSFLog& log)
{
  // Read some needed parameters

  Assert(params.keyExists("psf_order"));
  int psforder = params["psf_order"];

  bool output_dots=false;
  if (params.keyExists("output_dots")) output_dots=true;

  double gain = 0.;
  if (params.keyExists("gain")) gain = params["gain"];

  Assert(params.keyExists("psf_aperture"));
  double psfap = params["psf_aperture"];
  dbg<<"psfap = "<<psfap<<std::endl;
  psfap *= sigma_p;  // arcsec
  dbg<<"psfap => "<<psfap<<std::endl;

  int nstars = size();
  dbg<<"nstars = "<<nstars<<std::endl;

#ifdef ENDAT
  nstars = ENDAT;
#endif

  log.nstars = nstars;
#ifdef STARTAT
  log.nstars -= STARTAT;
#endif
#ifdef SINGLESTAR
  log.nstars = 1;
#endif

  // Main loop to measure psf shapelets:
#ifdef _OPENMP
#pragma omp parallel 
  {
    try {
#endif
      PSFLog log1;  // Just for this thread
      log1.NoWriteLog();
#ifdef _OPENMP
#pragma omp for schedule(guided)
#endif
      for(int i=0;i<nstars;i++) {
#ifdef STARTAT
	if (i < STARTAT) continue;
#endif
#ifdef SINGLESTAR
	if (i < SINGLESTAR) continue;
	if (i > SINGLESTAR) break;
	XDEBUG = true;
#endif
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  if (output_dots) { std::cerr<<"."; std::cerr.flush(); }
	  dbg<<"star "<<i<<":\n";
	  dbg<<"pos["<<i<<"] = "<<pos[i]<<std::endl;
	}

	// Start with an error code of unknown failure, in case
	// something happens that I don't detect as an error.
	flags[i] = UNKNOWN_FAILURE;
	Assert(i<pos.size());
	Assert(i<sky.size());
	Assert(i<noise.size());
	Assert(i<psf.size());
	Assert(i<nu.size());
	Assert(i<flags.size());
	dbg<<"Before MeasureSinglePSF1"<<std::endl;
	long flag1 = 0;
	MeasureSinglePSF(
	    // Input data:
	    pos[i], im, sky[i], trans, 
	    // Noise values:
	    noise[i], gain, weight_im,
	    // Parameters:
	    sigma_p, psfap, psforder,
	    // Log information
	    log1,
	    // Ouput value:
	    psf[i], nu[i], flag1);
	dbg<<"After MeasureSinglePSF"<<std::endl;

	flags[i] = flag1;
	if (!flag1) {
	  dbg<<"Successful psf measurement: "<<psf[i]<<std::endl;
	}
	else {
	  dbg<<"Unsuccessful psf measurement\n"; 
	}
      }
#ifdef _OPENMP
#pragma omp critical
#endif
      {
	log += log1;
      }
#ifdef _OPENMP
    }
    catch (...)
    {
      // This isn't supposed to happen.
      std::cerr<<"STATUS5BEG Caught some error in parallel region STATUS5END\n";
      exit(1); 
    }
  }
#endif

  int nsuccess = log.ns_psf;
  dbg<<nsuccess<<" successful star measurements, ";
  dbg<<nstars-nsuccess<<" unsuccessful\n";

  if (output_dots) {
    std::cerr
      <<std::endl
      <<"Success rate: "<<nsuccess<<"/"<<nstars
      <<std::endl;
  }

  return nsuccess;
}

PSFCatalog::PSFCatalog(const StarCatalog& starcat,
    const ConfigFile& _params, std::string key_prefix) :
  params(_params), prefix(key_prefix)
{
  dbg<<"Create PSFCatalog\n";
  dbg<<"ntot = "<<starcat.id.size()<<std::endl;
  Assert(starcat.isastar.size() == starcat.id.size());
  size_t nstars = std::count(starcat.isastar.begin(),starcat.isastar.end(),1);
  dbg<<"count stars = "<<nstars<<std::endl;
  id.reserve(nstars);
  pos.reserve(nstars);
  flags.reserve(nstars);
  for (size_t i=0; i<starcat.id.size(); i++) {
    if (starcat.isastar[i]) {
      id.push_back( starcat.id[i] );
      pos.push_back( starcat.pos[i] );
      sky.push_back( starcat.sky[i] );
      noise.push_back( starcat.noise[i] );
      flags.push_back( starcat.flags[i] );
    }
  }
  Assert(id.size() == nstars);
  dbg<<"nstars = "<<id.size()<<std::endl;

  // Set up a default psf vector for output when an object measurement
  // fails
  long psf_order = params.get(prefix + "order");
  BVec psf_default(psf_order,1.);
  for (size_t i=0; i<psf_default.size(); i++) {
    psf_default[i] = DEFVALNEG;
  }
  double nu_default = DEFVALNEG;

  nu.resize(id.size(),nu_default);
  psf.resize(id.size(),psf_default);

  Assert(id.size() == size());
  Assert(pos.size() == size());
  Assert(sky.size() == size());
  Assert(noise.size() == size());
  Assert(flags.size() == size());
  Assert(nu.size() == size());
  Assert(psf.size() == size());
}

PSFCatalog::PSFCatalog(const ConfigFile& _params, std::string key_prefix) :
  params(_params), prefix(key_prefix)
{
  Read();

  Assert(id.size() == size());
  Assert(pos.size() == size());
  Assert(sky.size() == size());
  Assert(noise.size() == size());
  Assert(flags.size() == size());
  Assert(nu.size() == size());
  Assert(psf.size() == size());
}

void PSFCatalog::WriteFits(std::string file) const
{
  Assert(id.size() == size());
  Assert(pos.size() == size());
  Assert(nu.size() == size());
  Assert(flags.size() == size());
  Assert(psf.size() == size());

  FitsFile fits(file, READWRITE, true);

  int ncoeff = psf[0].size();
  dbg<<"ncoeff = "<<ncoeff<<std::endl;

  std::stringstream coeff_form;
  coeff_form << ncoeff << "d";

  std::string id_col=params.get(prefix + "id_col");
  std::string x_col=params.get(prefix + "x_col");
  std::string y_col=params.get(prefix + "y_col");
  std::string sky_col=params.get(prefix + "sky_col");
  std::string noise_col=params.get(prefix + "noise_col");
  std::string flags_col=params.get(prefix + "flags_col");
  std::string nu_col=params.get(prefix + "nu_col");
  std::string order_col=params.get(prefix + "order_col");
  std::string sigma_col=params.get(prefix + "sigma_col");
  std::string coeffs_col=params.get(prefix + "coeffs_col");

  const int nfields=10;
  std::string table_cols[nfields] = {
    id_col,
    x_col,
    y_col,
    sky_col,
    noise_col,
    flags_col,
    nu_col,
    order_col,
    sigma_col,
    coeffs_col,
  };
  std::string table_types[nfields] = {
    "1j", // id
    "1d", // x
    "1d", // y
    "1d", // sky
    "1d", // noise
    "1i", // flags
    "1d", // nu
    "1j", // order
    "1d", // sigma
    coeff_form.str()
  };
  std::string table_units[nfields] = {
    "None",
    "Pixels",
    "Pixels",
    "ADU",
    "ADU^2",
    "None",
    "None",
    "None",
    "Arcsec",
    "None"
  };

  // Create a binary table
  fits.CreateBinaryTable(size(),nfields,table_cols,table_types,table_units);

  // Write the header keywords
  fits.WriteParKey(params, "version", TSTRING);
  fits.WriteParKey(params, "noise_method", TSTRING);
  fits.WriteParKey(params, "dist_method", TSTRING);

  fits.WriteParKey(params, prefix + "aperture", TDOUBLE);
  fits.WriteParKey(params, prefix + "order", TLONG);
  fits.WriteParKey(params, prefix + "seeing_est", TDOUBLE);

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
  fits.WriteColumn(TDOUBLE, 7, 1, 1, size(), &nu[0]);
  
  // Now we have to loop through each psf decomposition and write
  // separately.  This is pretty dumb.
  for (size_t i=0; i<size(); i++) {
    size_t row = i+1;
    long b_order = psf[i].GetOrder();
    double b_sigma = psf[i].GetSigma();
    fits.WriteColumn(TLONG, 8, row, 1, 1, &b_order);
    fits.WriteColumn(TDOUBLE, 9, row, 1, 1, &b_sigma);
    fits.WriteColumn(TDOUBLE, 10, row, 1, ncoeff, psf[i].cptr());
  }
}

void PSFCatalog::WriteAscii(std::string file, std::string delim) const
{
  Assert(id.size() == size());
  Assert(pos.size() == size());
  Assert(sky.size() == size());
  Assert(noise.size() == size());
  Assert(flags.size() == size());
  Assert(nu.size() == size());
  Assert(psf.size() == size());

  std::ofstream fout(file.c_str());
  if (!fout) {
    throw std::runtime_error("Error opening psf file");
  }

  for(size_t i=0;i<size();i++) {
    fout
      << id[i] << delim
      << pos[i].GetX() << delim
      << pos[i].GetY() << delim
      << sky[i] << delim
      << noise[i] << delim
      << flags[i] << delim
      << nu[i] << delim
      << psf[i].GetOrder() << delim
      << psf[i].GetSigma();
    for(size_t j=0;j<psf[i].size();++j)
      fout << delim << psf[i][j];
    fout << std::endl;
  }
}

void PSFCatalog::Write() const
{
  std::string file = Name(params,"psf");
  dbg<<"Writing to psf file: "<<file<<std::endl;

  if (file.find("fits") != std::string::npos) {
    WriteFits(file);
  } else {
    std::string delim = "  ";
    if (params.keyExists(prefix + "delim")) delim = params[prefix + "delim"];
    WriteAscii(file,delim);
  }
  dbg<<"Done Write\n";
}

void PSFCatalog::ReadFits(std::string file)
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
    throw std::runtime_error("nrows must be > 0");
  }

  std::string id_col=params.get(prefix + "id_col");
  std::string x_col=params.get(prefix + "x_col");
  std::string y_col=params.get(prefix + "y_col");
  std::string sky_col=params.get(prefix + "sky_col");
  std::string noise_col=params.get(prefix + "noise_col");
  std::string flags_col=params.get(prefix + "flags_col");
  std::string nu_col=params.get(prefix + "nu_col");
  std::string order_col=params.get(prefix + "order_col");
  std::string sigma_col=params.get(prefix + "sigma_col");
  std::string coeffs_col=params.get(prefix + "coeffs_col");

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
  sky.resize(nrows,0);
  fits.ReadScalarCol(sky_col,TDOUBLE,&sky[0], nrows);

  dbg<<"  "<<noise_col<<std::endl;
  noise.resize(nrows,0);
  fits.ReadScalarCol(noise_col,TDOUBLE,&noise[0], nrows);

  dbg<<"  "<<flags_col<<std::endl;
  flags.resize(nrows,0);
  fits.ReadScalarCol(flags_col,TLONG,&flags[0], nrows);

  dbg<<"  "<<nu_col<<std::endl;
  nu.resize(nrows,0);
  fits.ReadScalarCol(nu_col,TDOUBLE,&nu[0], nrows);

  // gotta loop for this one
  psf.reserve(nrows);
  for (size_t i=0; i<size(); ++i) {
    size_t row=i+1;
    double b_sigma;
    long b_order;
    fits.ReadCell(order_col,TLONG,&b_order,row,1);
    fits.ReadCell(sigma_col,TDOUBLE,&b_sigma,row,1);
    psf.push_back(BVec(b_order,b_sigma));
    int ncoeff=(b_order+1)*(b_order+2)/2;
    fits.ReadCell(coeffs_col,TDOUBLE,psf[i].ptr(),row,ncoeff);
  }
}

void PSFCatalog::ReadAscii(std::string file, std::string delim)
{
  std::ifstream fin(file.c_str());
  if (!fin) {
    throw std::runtime_error("Error opening psf file");
  }

  if (delim == "  ") {
    long id1,flag,b_order;
    double x,y,sky1,noise1,nu1,b_sigma;
    while ( fin
	>> id1
	>> x >> y
	>> sky1
	>> noise1
	>> flag
	>> nu1
	>> b_order >> b_sigma
	) 
    {
      id.push_back(id1);
      pos.push_back(Position(x,y));
      sky.push_back(sky1);
      noise.push_back(noise1);
      flags.push_back(flag);
      nu.push_back(nu1);
      psf.push_back(BVec(b_order,b_sigma));
      for(size_t j=0;j<psf.back().size();++j)
	fin >> psf.back()[j];
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
    ConvertibleString temp;
    long b_order;
    double x,y,b_sigma;
    while (getline(fin,temp,d))
    { 
      id.push_back(temp);
      getline(fin,temp,d); x = temp;
      getline(fin,temp,d); y = temp;
      pos.push_back(Position(x,y));
      getline(fin,temp,d); sky.push_back(temp);
      getline(fin,temp,d); noise.push_back(temp);
      getline(fin,temp,d); flags.push_back(temp);
      getline(fin,temp,d); nu.push_back(temp);
      getline(fin,temp,d); b_order = temp;
      getline(fin,temp,d); b_sigma = temp;
      psf.push_back(BVec(b_order,b_sigma));
      for(size_t j=0;j<psf.back().size()-1;++j) {
	getline(fin,temp,d); psf.back()[j] = temp;
      }
      getline(fin,temp); psf.back()[psf.back().size()-1] = temp;
    }
  }
}


void PSFCatalog::Read()
{
  std::string file = Name(params,"psf",false,true);
  // false,true = input_prefix=false, mustexist=true.
  // It is an input here, but it is in the output_prefix directory.
  dbg<< "Reading PSF cat from file: " << file << std::endl;

  if (file.find("fits") != std::string::npos) {
    ReadFits(file);
  } else {
    std::string delim = "  ";
    if (params.keyExists(prefix + "delim")) delim = params[prefix + "delim"];
    ReadAscii(file,delim);
  }
  dbg<<"Done Read PSFCatalog\n";
}

