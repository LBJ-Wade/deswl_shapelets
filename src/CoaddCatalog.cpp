#include "CoaddCatalog.h"


CoaddCatalog::CoaddCatalog(ConfigFile& _params):
  params(_params)
{
  ReadCatalog();
}

void CoaddCatalog::ReadCatalog()
{
  std::string file=params.get("coaddcat_file");
  int hdu = params.get("coaddcat_hdu");

  dbg<<"Opening FITS file at hdu "<<hdu<<std::endl;
  // true means read all as part of the construction
  CCfits::FITS fits(file, CCfits::Read, hdu-1, true);

  CCfits::ExtHDU& table=fits.extension(hdu-1);

  long nrows=table.rows();

  dbg<<"  nrows = "<<nrows<<std::endl;
  if (nrows <= 0) {
    throw std::runtime_error("nrows must be > 0");
  }

  std::string id_col=params.get("coaddcat_id_col");
  std::string x_col=params.get("coaddcat_x_col");
  std::string y_col=params.get("coaddcat_y_col");
  std::string sky_col=params.get("coaddcat_sky_col");

  std::string mag_col=params.get("coaddcat_mag_col");
  std::string mag_err_col=params.get("coaddcat_mag_err_col");

  //std::string noise_col=params.get("coaddcat_noise_col");
  std::string flags_col=params.get("coaddcat_flag_col");
  std::string ra_col=params.get("coaddcat_ra_col");
  std::string dec_col=params.get("coaddcat_dec_col");

  long start=1;
  long end=nrows;

  dbg<<"Reading columns"<<std::endl;
  dbg<<"  "<<id_col<<std::endl;
  table.column(id_col).read(id, start, end);

  dbg<<"  "<<x_col<<"  "<<y_col<<std::endl;
  pos.resize(nrows);
  std::vector<double> x;
  std::vector<double> y;
  table.column(x_col).read(x, start, end);
  table.column(y_col).read(y, start, end);
  for(long i=0;i<nrows;++i) pos[i] = Position(x[i],y[i]);

  dbg<<"  "<<sky_col<<std::endl;
  table.column(sky_col).read(sky, start, end);

  //dbg<<"  "<<noise_col<<std::endl;
  //table.column(noise_col).read(noise, start, end);

  dbg<<"  "<<flags_col<<std::endl;
  table.column(flags_col).read(flags, start, end);

  dbg<<"  "<<ra_col<<"  "<<dec_col<<std::endl;
  skypos.resize(nrows);
  std::vector<double> ra;
  std::vector<double> dec;
  table.column(ra_col).read(ra, start, end);
  table.column(dec_col).read(dec, start, end);
  for(long i=0;i<nrows;++i) {
    skypos[i] = Position(ra[i],dec[i]);
    // The convention for Position is to use arcsec for everything.
    // ra and dec come in as degrees.  So wee need to convert to arcsec.
    skypos[i] *= 3600.;  // deg -> arcsec
  }

}

void CoaddCatalog::ReadFileLists()
{
  std::string srclist_file = params.get("coadd_srclist");

  std::string image_filename;
  std::string psf_filename;

  dbg<<"Opening coadd srclist\n";
  std::ifstream flist(srclist_file.c_str(), std::ios::in);

  image_file_list.clear();
  fitpsf_file_list.clear();

  if (flist.is_open() ) {
    while (flist >> image_filename) {
      if (flist >> psf_filename) {

	image_file_list.push_back(image_filename);
	fitpsf_file_list.push_back(psf_filename);

      } else {
	throw std::runtime_error("Failed to read from file");
      }
    }
  }

  Assert(image_file_list.size() == fitpsf_file_list.size());

}

void CoaddCatalog::ReadPixelLists()
{

  ReadFileLists();

  pixlist.clear();
  pixlist.resize(skypos.size());

  getpixlist_flags.resize(skypos.size(),0);
  nimages_found.resize(skypos.size(), 0);
  nimages_gotpix.resize(skypos.size(), 0);

  // Loop over the files and read pixel lists for each object
  for (int fnum=0; fnum<image_file_list.size(); fnum++) {
    std::string image_file = image_file_list[fnum];
    std::string psf_file = fitpsf_file_list[fnum];

    std::cout<<"  Reading image file: "<<image_file<<"\n";
    params["image_file"] = image_file;
    params["weight_file"] = image_file;
    params["dist_file"] = image_file;
    params["dist_hdu"] = 2;
    params["dist_method"] = "WCS";

    params["fitpsf_file"] = psf_file;

    GetImagePixelLists();
  }

}


// Get pixel lists from the file specified in params
void CoaddCatalog::GetImagePixelLists()
{

  std::auto_ptr<Image<double> > weight_im;
  Image<double> im(params,weight_im);

  int maxi=im.GetMaxI();
  int maxj=im.GetMaxJ();
  std::cout<<"MaxI: "<<maxi<<" MaxJ: "<<maxj<<"\n";

  std::cout<<"  Reading distortion\n";
  Transformation trans(params);

  // read the psf
  std::cout<<"Reading fitted psf\n";
  FittedPSF fitpsf(params);

  // loop over the the objects, if the object falls on the image get
  // the pixel list

  double noise = 0.0;
  double gain=0.0;
  double gal_aperture = params.get("shear_aperture");
  for (int i=0; i<skypos.size(); i++) {

    // convert ra/dec to x,y in this image

    Position pxy((double)maxi/2.,(double)maxj/2.);

    if (!trans.InverseTransform(skypos[i], pxy) ) {
      throw std::runtime_error("transformation failed");
    }

    double x=pxy.GetX();
    double y=pxy.GetY();
    if ( (x >= 0) && (x <= maxi) && (y >= 0) && (y <= maxj) ) {
      dbg<<"("<<skypos[i]<<")  x="<<x<<" y="<<y<<"\n";

      nimages_found[i]++;

      // Calculate the psf from the fitted-psf formula:
      std::vector<BVec> psf(1,
	  BVec(fitpsf.GetPSFOrder(), fitpsf.GetSigma()));
      try {
	dbg<<"for fittedpsf cen = "<<pxy<<std::endl;
	psf[0] = fitpsf(pxy);
      } catch (Range_error& e) {
	dbg<<"fittedpsf range error: \n";
	xdbg<<"p = "<<pxy<<", b = "<<e.b<<std::endl;
	getpixlist_flags[i] |= FITTEDPSF_EXCEPTION;
	continue;
      }

      double sigma_p=0;
      Assert(psf.size() > 0);
      for(size_t j=0;j<psf.size();j++) sigma_p += 1./psf[j].GetSigma();
      sigma_p = double(psf.size()) / sigma_p;

      // we are using the weight image so the noise and gain are 
      // dummy variables
      double sigma_obs = 2.*sigma_p;
      double galap = gal_aperture * sigma_obs;

      long flag = 0;
      std::vector<Pixel> tmp_pixlist;
      GetPixList(
	  im,tmp_pixlist,pxy,
	  sky[i],noise,gain,weight_im.get(),trans,galap,flag);

      // make sure not (edge or < 10 pixels) although edge is already
      // checked above
      if (flag == 0) {
	pixlist[i].push_back(tmp_pixlist);
	nimages_gotpix[i]++;
      } else {
	getpixlist_flags[i] |= flag;
      }
    }
  }
}
