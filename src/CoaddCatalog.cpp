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
  }

}

void CoaddCatalog::ReadPixelLists()
{
  // source images and psf files
  std::string srclist = params.get("coadd_srclist");

  dbg<<"Opening coadd srclist\n";
  std::ifstream flistfile(srclist.c_str(), std::ios::in);

  std::string image_file;
  std::string psf_file;

  // these may not be present so I'll define them
  if (!params.keyExists("image_file")) {
    params.Append("image_file=junk");
  }
  if (!params.keyExists("weight_file")) {
    params.Append("weight_file=junk");
  }

  if (!params.keyExists("dist_file")) {
    params.Append("dist_file=junk");
  }
  if (!params.keyExists("dist_hdu")) {
    params.Append("dist_hdu=junk");
  }
  if (!params.keyExists("dist_method")) {
    params.Append("dist_method=junk");
  }

  if (flistfile.is_open() ) {
    while (! flistfile.eof() ) {

      flistfile >> image_file;
      flistfile >> psf_file;

      if (image_file != "" && psf_file != "") {

	std::cout<<"  Reading image file: "<<image_file<<"\n";
	params["image_file"] = image_file;
	params["weight_file"] = image_file;
	params["dist_file"] = image_file;
	params["dist_hdu"] = 2;
	params["dist_method"] = "WCS";

	std::auto_ptr<Image<double> > weight_im;
	Image<double> im(params,weight_im);

	int maxi=im.GetMaxI();
	int maxj=im.GetMaxJ();
	std::cout<<"MaxI: "<<maxi<<" MaxJ: "<<maxj<<"\n";

	std::cout<<"  Reading distortion\n";
	Transformation trans(params);

	// loop over the the objects, if the object falls on the image get
	// the pixel list

	for (int i=0; i<skypos.size(); i++) {
	  // convert ra/dec to x,y in this image
	
	  //Position pxy((double)maxi/2.,(double)maxj/2.);
	  // This only works when we are testing the non-coadd catalog
	  Position pxy = pos[i];

	  std::cout<<"("<<skypos[i]<<")  initial: ("<<pxy<<")";
	  if (!trans.InverseTransform(skypos[i], pxy) ) {
	    std::cout<<"Transform failed ("<<skypos[i]<<")\n";
	  }
	  std::cout<<" final: ("<<pxy<<")\n";

	  Position skypos_back;
	  trans.Transform(pxy, skypos_back);

	  std::cout<<"skypos in: ("<<skypos[i]<<") out: ("<<skypos_back<<")\n";

	  double x=pxy.GetX();
	  double y=pxy.GetY();
	  if ( (x >= 0) && (x <= maxi) && (y >= 0) && (y <= maxj) ) {
	    std::cout<<"("<<skypos[i]<<")  x="<<x<<" y="<<y<<"\n";
	  }
	  return;
	}
      }
      // for now only first image
      //return;
    }
  }

}

