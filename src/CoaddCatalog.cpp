#include "CoaddCatalog.h"

#define UseInverseTransform

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
  noise.resize(nrows,0);

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
	throw std::runtime_error("Failed to read from source list file");
      }
    }
  }

  Assert(image_file_list.size() == fitpsf_file_list.size());

}

void CoaddCatalog::Resize(int n)
{
  pixlist.clear();

  pixlist.resize(n);

  input_flags.resize(n,0);
  nimages_found.resize(n, 0);
  nimages_gotpix.resize(n, 0);

  shear.resize(n);
  nu.resize(n);
  cov.resize(n);

  int gal_order = params.get("shear_gal_order");
  BVec shape_default(gal_order,1.);
  shape_default.SetAllTo(DEFVALNEG);
  shape.resize(n,shape_default);


}

void CoaddCatalog::ReadPixelLists()
{

  ReadFileLists();

  // Loop over the files and read pixel lists for each object
  for (int fnum=0; fnum<image_file_list.size(); fnum++) {
    std::string image_file = image_file_list[fnum];
    std::string psf_file = fitpsf_file_list[fnum];

    dbg<<"Reading image file: "<<image_file<<"\n";
    params["image_file"] = image_file;
    params["weight_file"] = image_file;
    params["dist_file"] = image_file;
    params["dist_hdu"] = 2;
    params["dist_method"] = "WCS";

    params["fitpsf_file"] = psf_file;

    GetImagePixelLists();
    dbg<<"\n";
  }

}


// Get pixel lists from the file specified in params
void CoaddCatalog::GetImagePixelLists()
{

  std::auto_ptr<Image<double> > weight_im;
  Image<double> im(params,weight_im);

  int maxi=im.GetMaxI();
  int maxj=im.GetMaxJ();
  xdbg<<"MaxI: "<<maxi<<" MaxJ: "<<maxj<<"\n";

  // read transformation between ra/dec and x/y
  Transformation trans(params);

  // read the psf
  FittedPSF fitpsf(params);

#ifdef UseInverseTransform
  Transformation invtrans;
  Bounds invb = invtrans.MakeInverseOf(trans,im.GetBounds(),4);
#endif

  // loop over the the objects, if the object falls on the image get
  // the pixel list

  double noise = 0.0;
  double gain=0.0;
  double gal_aperture = params.get("shear_aperture");
  dbg<<"Extracting pixel lists\n";
  for (int i=0; i<skypos.size(); i++) {

    // convert ra/dec to x,y in this image

#ifdef UseInverseTransform
    // Figure out a good starting point for the nonlinear solver:
    Position pxy;
    xdbg<<"skypos = "<<skypos[i]<<std::endl;
    if (!invb.Includes(skypos[i])) {
      xdbg<<"skypos "<<skypos[i]<<" not in "<<invb<<std::endl;
      continue;
    }
    invtrans.Transform(skypos[i],pxy);
    xdbg<<"invtrans(skypos) = "<<pxy<<std::endl;
#else
    Position pxy((double)maxi/2.,(double)maxj/2.);
#endif

    if (!trans.InverseTransform(skypos[i], pxy) ) {
      throw std::runtime_error("transformation failed");
    }

    double x=pxy.GetX();
    double y=pxy.GetY();
    if ( (x >= 0) && (x <= maxi) && (y >= 0) && (y <= maxj) ) {
      xdbg<<"("<<skypos[i]<<")  x="<<x<<" y="<<y<<"\n";

      nimages_found[i]++;

      // Calculate the psf from the fitted-psf formula:
      std::vector<BVec> psf(1,
	  BVec(fitpsf.GetPSFOrder(), fitpsf.GetSigma()));
      try {
	xdbg<<"for fittedpsf cen = "<<pxy<<std::endl;
	psf[0] = fitpsf(pxy);
      } catch (Range_error& e) {
	xdbg<<"fittedpsf range error: \n";
	xdbg<<"p = "<<pxy<<", b = "<<e.b<<std::endl;
	input_flags[i] |= FITTEDPSF_EXCEPTION;
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
	input_flags[i] |= flag;
      }
    }
  } // loop over objects
  dbg<<"Done extracting pixel lists\n";
}

void CoaddCatalog::MeasureMultiShears()
{

}


void CoaddCatalog::WriteFits() const
{
  std::string file=Name(params, "multishear");

  // ! means overwrite existing file
  CCfits::FITS fits("!"+file, CCfits::Write);

  const int nfields=20;
  std::vector<string> colnames(nfields);
  std::vector<string> colfmts(nfields);
  std::vector<string> colunits(nfields);

  colnames[0] = params.get("coaddshear_id_col");
  colnames[1] = params.get("coaddshear_x_col");
  colnames[2] = params.get("coaddshear_y_col");
  colnames[3] = params.get("coaddshear_sky_col");
  colnames[4] = params.get("coaddshear_noise_col");
  colnames[5] = params.get("coaddshear_flags_col");
  colnames[6] = params.get("coaddshear_ra_col");
  colnames[7] = params.get("coaddshear_dec_col");
  colnames[8] = params.get("coaddshear_shear1_col");
  colnames[9] = params.get("coaddshear_shear2_col");
  colnames[10] = params.get("coaddshear_nu_col");
  colnames[11] = params.get("coaddshear_cov00_col");
  colnames[12] = params.get("coaddshear_cov01_col");
  colnames[13] = params.get("coaddshear_cov11_col");
  colnames[14] = params.get("coaddshear_order_col");
  colnames[15] = params.get("coaddshear_sigma_col");
  colnames[16] = params.get("coaddshear_coeffs_col");

  colnames[17] = params.get("coaddshear_nimages_found_col");
  colnames[18] = params.get("coaddshear_nimages_gotpix_col");
  colnames[19] = params.get("coaddshear_input_flags_col");

  int ncoeff = shape[0].size();
  dbg<<"ncoeff = "<<ncoeff<<std::endl;
  std::stringstream coeff_form;
  coeff_form << ncoeff << "D";
  colfmts[16] = coeff_form.str(); // shapelet coeffs
 
  colfmts[0] = "1J"; // id
  colfmts[1] = "1D"; // x
  colfmts[2] = "1D"; // y
  colfmts[3] = "1D"; // sky
  colfmts[4] = "1D"; // noise
  colfmts[5] = "1J"; // flags
  colfmts[6] = "1D"; // ra
  colfmts[7] = "1D"; // dec
  colfmts[8] = "1D"; // shear1
  colfmts[9] = "1D"; // shear2
  colfmts[10] = "1D"; // nu
  colfmts[11] = "1D"; // cov00
  colfmts[12] = "1D"; // cov01
  colfmts[13] = "1D"; // cov11
  colfmts[14] = "1J"; // order
  colfmts[15] = "1D"; // sigma

  colfmts[17] = "1J"; // nimages_found
  colfmts[18] = "1J"; // nimages_gotpix
  colfmts[19] = "1J"; // input_flags


  dbg<<"Before Create table"<<std::endl;
  CCfits::Table* table;
  table = fits.addTable("coadd_shearcat",size(),colnames,colfmts,colunits);

  // Header keywords
  std::string tmvvers = tmv::TMV_Version();
  std::string wlvers = WlVersion();

  table->addKey("tmvvers", tmvvers, "version of TMV code");
  table->addKey("wlvers", wlvers, "version of weak lensing code");

  std::string str;
  double dbl;
  int intgr;
  //CCfitsWriteParKey(params, table, "version", str);
  CCfitsWriteParKey(params, table, "noise_method", str);
  CCfitsWriteParKey(params, table, "dist_method", str);

  CCfitsWriteParKey(params, table, "shear_aperture", dbl);
  CCfitsWriteParKey(params, table, "shear_max_aperture", dbl);
  CCfitsWriteParKey(params, table, "shear_gal_order", intgr);
  CCfitsWriteParKey(params, table, "shear_gal_order2", intgr);
  CCfitsWriteParKey(params, table, "shear_min_gal_size", dbl);
  CCfitsWriteParKey(params, table, "shear_f_psf", dbl);

  // data
  // make vector copies for writing
  std::vector<double> x(pos.size());
  std::vector<double> y(pos.size());
  std::vector<double> ra(size());
  std::vector<double> dec(size());
  std::vector<double> shear1(size());
  std::vector<double> shear2(size());
  std::vector<double> cov00(size());
  std::vector<double> cov01(size());

  for(size_t i=0;i<pos.size();i++) {
    x[i] = pos[i].GetX();
    y[i] = pos[i].GetY();
  }
  for(size_t i=0;i<size();i++) { 
    ra[i] = skypos[i].GetX();
    dec[i] = skypos[i].GetY();
  }
  for(size_t i=0;i<size();i++) { 
    shear1[i] = real(shear[i]);
    shear2[i] = imag(shear[i]);
  }
  std::vector<double> cov11(size());
  for(size_t i=0;i<size();i++) { 
    cov00[i] = cov[i](0,0);
    cov01[i] = cov[i](0,1);
    cov11[i] = cov[i](1,1);
  }

  int startrow=1;

  table->column(colnames[0]).write(id,startrow);
  table->column(colnames[1]).write(x,startrow);
  table->column(colnames[2]).write(y,startrow);
  table->column(colnames[3]).write(sky,startrow);
  table->column(colnames[4]).write(noise,startrow);
  table->column(colnames[5]).write(flags,startrow);
  table->column(colnames[6]).write(ra,startrow);
  table->column(colnames[7]).write(dec,startrow);
  table->column(colnames[8]).write(shear1,startrow);
  table->column(colnames[9]).write(shear2,startrow);
  table->column(colnames[10]).write(nu,startrow);
  table->column(colnames[11]).write(cov00,startrow);
  table->column(colnames[12]).write(cov01,startrow);
  table->column(colnames[13]).write(cov11,startrow);

  for (size_t i=0; i<size(); i++) {
    size_t row = i+1;
    long b_order = shape[i].GetOrder();
    double b_sigma = shape[i].GetSigma();

    table->column(colnames[14]).write(&b_order,1,row);
    table->column(colnames[15]).write(&b_sigma,1,row);
    double* cptr = (double *) shape[i].cptr();
    table->column(colnames[16]).write(cptr, ncoeff, 1, row);

  }

  table->column(colnames[17]).write(nimages_found,startrow);
  table->column(colnames[18]).write(nimages_gotpix,startrow);
  table->column(colnames[19]).write(input_flags,startrow);


}
