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

void CoaddCatalog::ReadPixelLists()
{

  pixlist.resize(skypos.size());
  getpixlist_flags.resize(skypos.size());

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

  if (!params.keyExists("fitpsf_file")) {
    params.Append("fitpsf_file=junk");
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

	params["fitpsf_file"] = psf_file;

	double gal_aperture = params.get("shear_aperture");

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

	for (int i=0; i<skypos.size(); i++) {

	  long flag=0;

	  // convert ra/dec to x,y in this image
	
	  Position pxy((double)maxi/2.,(double)maxj/2.);
	  // This only works when we are testing the non-coadd catalog
	  //Position pxy = pos[i];

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

	    // Calculate the psf from the fitted-psf formula:
	    std::vector<BVec> psf(1,
		BVec(fitpsf.GetPSFOrder(), fitpsf.GetSigma()));
	    try {
	      dbg<<"for fittedpsf cen = "<<pxy<<std::endl;
	      psf[0] = fitpsf(pxy);
	    } catch (Range_error& e) {
	      dbg<<"fittedpsf range error: \n";
	      xdbg<<"p = "<<pxy<<", b = "<<e.b<<std::endl;
	      flag |= FITTEDPSF_EXCEPTION;
	      return;
	    }


	    double sigma_p=0;
	    Assert(psf.size() > 0);
	    for(size_t j=0;j<psf.size();j++) sigma_p += 1./psf[j].GetSigma();
	    sigma_p = double(psf.size()) / sigma_p;

	    // we are using the weight image so the noise is a dummy variable
	    // same for gain
	    double noise = 0.0;
	    double gain=0.0;
	    double sigma_obs = 2.*sigma_p;
	    double galap = gal_aperture * sigma_obs;

	    flag =0;
	    std::vector<Pixel> tmp_pixlist;
	    GetPixList(
		im,
		tmp_pixlist,
		pxy,
		sky[i],noise,
		gain,weight_im.get(),trans,
		galap,
		flag);

	    // make sure not (edge or < 10 pixels) although edge is already
	    // checked above
	    if (flag == 0) {
	      pixlist[i].push_back(tmp_pixlist);
	    }
	  }
	  //return;
	}
      }
      // for now only first image
      //return;
    }
  }

}

