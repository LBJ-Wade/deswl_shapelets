#include "SXCat.h"

void ResizeSXCat(SXCAT_STRUCT& cat, long n)
{
  cat.id.resize(n,0);
  cat.x.resize(n,0);
  cat.y.resize(n,0);
  cat.pos.resize(n);
  cat.mag.resize(n,0);
  cat.mag_err.resize(n,0);
  cat.local_sky.resize(n,0);
  cat.flags.resize(n,0);

  // generated
  //cat.sigma.resize(n,0);
  //cat.size_flags.resize(n,0);
  //cat.star_flag.resize(n,0);

  // Not used
  cat.noise.resize(n,0);
}

void ReadSXCat(ConfigFile& params, SXCAT_STRUCT& cat)
{

  int cat_hdu = 1;
  if (params.keyExists("cat_hdu")) cat_hdu = params["cat_hdu"];
  std::string file = Name(params,"cat",true);
  int minrows = 100;
  if (params.keyExists("sx_minrows")) minrows = params["sx_minrows"];

  std::cout<< "Reading cat from file: " << file << std::endl;

  FitsFile fits(file);

  //int hdutype;

  dbg<<"Moving to HDU #"<<cat_hdu<<std::endl;
  fits.GotoHDU(cat_hdu);

  // These are not very portable.
  //int nfound;
  long nrows=0;
  fits.ReadKey("NAXIS2", TLONG, (char *)&nrows);

  std::cout<<"  nrows = "<<nrows<<std::endl;
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

  // use get() since it checks for existence
  std::string id_name=params.get("sx_id_name");
  std::string x_name=params.get("sx_x_name");
  std::string y_name=params.get("sx_y_name");
  std::string local_sky_name=params.get("sx_local_sky_name");
  std::string mag_name=params.get("sx_mag_name");
  std::string mag_err_name=params.get("sx_mag_err_name");
  std::string flags_name=params.get("sx_flags_name");

  std::cout<<"Reading columns"<<std::endl;
  std::cout<<"  "<<id_name<<std::endl;
  fits.ReadScalarCol((char *)id_name.c_str(),TLONG,(char *)&cat.id[0], nrows);
  std::cout<<"  "<<x_name<<std::endl;
  fits.ReadScalarCol((char *)x_name.c_str(),TFLOAT,(char *)&cat.x[0], nrows);
  std::cout<<"  "<<y_name<<std::endl;
  fits.ReadScalarCol((char *)y_name.c_str(),TFLOAT,(char *)&cat.y[0], nrows);
  std::cout<<"  "<<local_sky_name<<std::endl;
  fits.ReadScalarCol((char *)local_sky_name.c_str(),TDOUBLE,(char *)&cat.local_sky[0], nrows);
  std::cout<<"  "<<mag_name<<std::endl;
  fits.ReadScalarCol((char *)mag_name.c_str(),TFLOAT,(char *)&cat.mag[0], nrows);
  std::cout<<"  "<<mag_err_name<<std::endl;
  fits.ReadScalarCol((char *)mag_err_name.c_str(),TFLOAT,(char *)&cat.mag_err[0], nrows);
  std::cout<<"  "<<flags_name<<std::endl;
  fits.ReadScalarCol((char *)flags_name.c_str(),TSHORT,(char *)&cat.flags[0], nrows);

  for (long i=0; i< nrows; i++) {
    cat.pos[i] = Position(cat.x[i], cat.y[i]);
  }

  fits.Close();
}







// This writes to hdu=2
void WriteFindStarsCat(ConfigFile& params, FINDSTARS_STRUCT& cat)
{

  std::string file = Name(params,"allcat");
  FitsFile fits(file.c_str(), READWRITE, true);
  fitsfile* fptr = fits.get_fptr();

  std::string id_name=params.get("fs_id_name");
  std::string sigma0_name=params.get("fs_sigma0_name");
  std::string size_flags_name=params.get("fs_size_flags_name");
  std::string star_flag_name=params.get("fs_star_flag_name");

  int nfields=4;
  char *table_names[] =  
      {(char*)id_name.c_str(),
	(char*)sigma0_name.c_str(),
	(char*)size_flags_name.c_str(),
	(char*)star_flag_name.c_str()};
  char *table_types[] =  
      {(char*)"1j",
	(char*)"1d",
	(char*)"1j",
	(char*)"1j"};
  char *table_units[] =  
      {(char*)"None",
	(char*)"pixels",
	(char*)"None",
	(char*)"None"};


  // Create a binary table
  int fits_status=0;
  int tbl_type = BINARY_TBL;
  fits_create_tbl(fptr, tbl_type, cat.id.size(), nfields, 
      table_names, table_types, table_units, NULL, &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error creating FindStars FITS table: ";
    throw FitsException(serr);
  }


  // need a more general way to write columns
  fits_status=0;
  fits_write_col(
      fptr, TLONG, 1, 1, 1, cat.id.size(), &cat.id[0], &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing 'id' to fits column";
    throw FitsException(serr);
  }

  fits_status=0;
  fits_write_col(
      fptr, TDOUBLE, 2, 1, 1, cat.sigma0.size(), &cat.sigma0[0], 
      &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing 'sigma' to fits column: ";
    throw FitsException(serr);
  }

  fits_status=0;
  fits_write_col(
      fptr, TINT, 3, 1, 1, cat.id.size(), &cat.size_flags[0], 
      &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing 'size_flags' to fits column";
    throw FitsException(serr);
  }

  fits_status=0;
  fits_write_col(
      fptr, TINT, 4, 1, 1, cat.id.size(), &cat.star_flag[0], 
      &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing 'star_flag' to fits column";
    throw FitsException(serr);
  }




  fits.Close();

}


void ReadFindStarsCat(ConfigFile& params, FINDSTARS_STRUCT& cat)
{

  int cat_hdu = 2;
  if (params.keyExists("fs_cat_hdu")) cat_hdu = params["fs_cat_hdu"];

  std::string file = Name(params,"allcat",true);
  std::cout<< "Reading FindStars cat from file: " << file << std::endl;

  FitsFile fits(file);

  //int hdutype;

  dbg<<"Moving to HDU #"<<cat_hdu<<std::endl;
  fits.GotoHDU(cat_hdu);

  // These are not very portable.
  //int nfound;
  long nrows=0;
  fits.ReadKey("NAXIS2", TLONG, (char *)&nrows);

  std::cout<<"  nrows = "<<nrows<<std::endl;
  if (nrows <= 0) {
    std::cerr<<"nrows must be >= 0"<<std::endl;
    throw FAILURE_FORMAT_ERROR;
  }

  // Allocate memory for the columns we will read
  ResizeFindStarsCat(cat, nrows);

  std::string id_name=params.get("fs_id_name");
  std::string sigma0_name=params.get("fs_sigma0_name");
  std::string size_flags_name=params.get("fs_size_flags_name");
  std::string star_flag_name=params.get("fs_star_flag_name");

  std::cout<<"Reading columns"<<std::endl;
  std::cout<<"  "<<id_name<<std::endl;
  fits.ReadScalarCol((char *)id_name.c_str(),TLONG,(char *)&cat.id[0], nrows);

  fits.ReadScalarCol((char *)sigma0_name.c_str(),TDOUBLE,(char *)&cat.sigma0[0], nrows);

  fits.ReadScalarCol((char *)size_flags_name.c_str(),TINT,(char *)&cat.size_flags[0], nrows);
  fits.ReadScalarCol((char *)star_flag_name.c_str(),TINT,(char *)&cat.star_flag[0], nrows);

  fits.Close();
}

void ResizeFindStarsCat(FINDSTARS_STRUCT& cat, size_t n)
{
  cat.id.resize(n,0);
  cat.sigma0.resize(n,0);
  cat.size_flags.resize(n,0);
  cat.star_flag.resize(n,0);

  cat.local_sky.resize(n,0);
  cat.pos.resize(n);
  cat.noise.resize(n,0);

}



void ResizePSFCat(PSF_STRUCT& cat, size_t n, int psf_order, double sigma)
{

  cat.id.resize(n,0);
  cat.psf_flags.resize(n,0);
  cat.nu.resize(n,0);

  cat.psf_order.clear();
  cat.psf_order.resize(n,psf_order);

  cat.sigma_p.resize(n,sigma);
  cat.psf.resize(n, BVec(psf_order,sigma));
}


// This writes to hdu=2
void WritePSFCat(ConfigFile& params, PSF_STRUCT& cat)
{

  int colnum;
  LONGLONG firstrow;
  LONGLONG firstel;
  LONGLONG nel;
  std::stringstream err;


  std::string file = Name(params,"psf");
  FitsFile fits(file.c_str(), READWRITE, true);
  fitsfile* fptr = fits.get_fptr();

  int ncoeff = cat.psf[0].size();

  std::stringstream coeff_form;
  coeff_form << ncoeff << "d";

  std::string id_name=params.get("psf_id_name");
  std::string psf_flags_name=params.get("psf_flags_name");
  std::string nu_name=params.get("psf_nu_name");
  std::string psf_order_name=params.get("psf_order_name");
  std::string sigma_p_name=params.get("psf_sigma_p_name");
  std::string coeffs_name=params.get("psf_coeffs_name");



  int nfields=6;
  char *table_names[] =  
      {(char*)id_name.c_str(),
	(char*)psf_flags_name.c_str(),
	(char*)nu_name.c_str(),
	(char*)psf_order_name.c_str(),
	(char*)sigma_p_name.c_str(),
	(char*)coeffs_name.c_str()};
  char *table_types[] =  
      {(char*)"1j",
	(char*)"1j",
	(char*)"1d",
	(char*)"1j",
	(char*)"1d",
	(char*)coeff_form.str().c_str()};
  char *table_units[] =  
      {(char*)"None",
	(char*)"None",
	(char*)"None",
	(char*)"None",
	(char*)"arcsec",
	(char*)"None"};


  // Create a binary table
  int fits_status=0;
  int tbl_type = BINARY_TBL;
  fits_create_tbl(fptr, tbl_type, cat.id.size(), nfields, 
      table_names, table_types, table_units, NULL, &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error creating PSF FITS table: ";
    throw FitsException(serr);
  }




  // Write some things that are the same for the whole field into a
  // header keyword.  These are also written as columns right now becaues
  // we want to allow for them to be variable.
  fits_status=0;
  fits_write_key(fptr, TINT, "PSFORDER", &cat.psf_order[0], "Order of shapelets expansion for PSF stars", &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing keyword PSFORDER";
    throw FitsException(serr);
  }

  fits_status=0;
  fits_write_key(fptr, TINT, "NCOEFFP", &ncoeff, "Number of coeffs (psforder+1)*(psforder+2)/2", &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing keyword NCOEFFP";
    throw FitsException(serr);
  }

  fits_status=0;
  fits_write_key(fptr, TDOUBLE, "SIGMA_P", &cat.sigma_p[0], "Scale used in shapelets expansion", &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing keyword SIGMA_P";
    throw FitsException(serr);
  }





  // need a more general way to write columns
  fits_status=0;
  colnum = 1;  
  firstrow = 1;
  firstel = 1;
  nel = cat.id.size();
  fits_write_col(
      fptr, TLONG, 
      colnum, firstrow, firstel, nel, 
      &cat.id[0], 
      &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing 'id' to fits column";
    throw FitsException(serr);
  }


  fits_status=0;
  colnum = 2;  
  firstrow = 1;
  firstel = 1;
  nel = cat.psf_flags.size();
  fits_write_col(
      fptr, TINT, 
      colnum, firstrow, firstel, nel, 
      &cat.psf_flags[0], 
      &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing 'psf_flags' to fits column";
    throw FitsException(serr);
  }


  fits_status=0;
  colnum = 3;  
  firstrow = 1;
  firstel = 1;
  nel = cat.nu.size();
  fits_write_col(
      fptr, TDOUBLE, 
      colnum, firstrow, firstel, nel, 
      &cat.nu[0], 
      &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing 'nu' to fits column: ";
    throw FitsException(serr);
  }



  fits_status=0;
  colnum = 4;  
  firstrow = 1;
  firstel = 1;
  nel = cat.psf_order.size();
  fits_write_col(
      fptr, TINT, 
      colnum, firstrow, firstel, nel, 
      &cat.psf_order[0], 
      &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing 'psf_order' to fits column";
    throw FitsException(serr);
  }




  fits_status=0;
  colnum = 5;  
  firstrow = 1;
  firstel = 1;
  nel = cat.sigma_p.size();
  fits_write_col(
      fptr, TDOUBLE, 
      colnum, firstrow, firstel, nel, 
      &cat.sigma_p[0], 
      &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing 'sigma_p' to fits column: ";
    throw FitsException(serr);
  }

  // Now we have to loop through each psf decomposition and write
  // separately.  This is pretty dumb.
  colnum = 6;  

  for (size_t i=0; i<cat.id.size(); i++) {
    fits_status=0;
    LONGLONG row = i+1;
    firstel = 1;
    //nel=1;
    nel=cat.psf[0].size();
    //nel = cat.sigma_p.size();
    fits_write_col(
	fptr, TDOUBLE, 
	colnum, row, firstel, nel, 
	&cat.psf[i][0], 
	&fits_status);
    if (!fits_status==0) {
      fits_report_error(stderr, fits_status); 
      std::stringstream serr;
      serr<<"Error writing row "<<row<<" for column 'coeffs'";
      throw FitsException(serr.str());
    }
  }






  fits.Close();

}

void ReadPSFCat(ConfigFile& params, PSF_STRUCT& cat)
{

  int cat_hdu = 2;
  if (params.keyExists("psf_cat_hdu")) cat_hdu = params["psf_cat_hdu"];

  std::string file = Name(params,"psf",true);
  std::cout<< "Reading PSF cat from file: " << file << std::endl;

  FitsFile fits(file);

  //int hdutype;

  dbg<<"Moving to HDU #"<<cat_hdu<<std::endl;
  fits.GotoHDU(cat_hdu);

  // These are not very portable.
  //int nfound;
  long nrows=0;
  fits.ReadKey("NAXIS2", TLONG, (char *)&nrows);

  int psforder=0;
  fits.ReadKey("PSFORDER", TINT, (char*)&psforder);

  std::cout<<"  nrows = "<<nrows<<std::endl;
  if (nrows <= 0) {
    std::cerr<<"nrows must be > 0"<<std::endl;
    throw FAILURE_FORMAT_ERROR;
  }

  // Allocate memory for the columns we will read
  ResizePSFCat(cat, nrows, psforder);

  std::string id_name=params.get("psf_id_name");
  std::string psf_flags_name=params.get("psf_flags_name");
  std::string nu_name=params.get("psf_nu_name");
  std::string psf_order_name=params.get("psf_order_name");
  std::string sigma_p_name=params.get("psf_sigma_p_name");
  std::string coeffs_name=params.get("psf_coeffs_name");

  std::cout<<"Reading columns"<<std::endl;
  std::cout<<"  "<<id_name<<std::endl;
  fits.ReadScalarCol((char *)id_name.c_str(),TLONG,(char *)&cat.id[0], nrows);
  fits.ReadScalarCol((char *)psf_flags_name.c_str(),TINT,(char *)&cat.psf_flags[0], nrows);
  fits.ReadScalarCol((char *)nu_name.c_str(),TDOUBLE,(char *)&cat.nu[0], nrows);
  fits.ReadScalarCol((char *)psf_order_name.c_str(),TINT,(char *)&cat.psf_order[0], nrows);
  fits.ReadScalarCol((char *)sigma_p_name.c_str(),TDOUBLE,(char *)&cat.sigma_p[0], nrows);

  // gotta loop for this one
  int ncoeff=(psforder+1)*(psforder+2)/2;
  for (size_t i=0; i<cat.id.size(); i++) { 

    size_t row=i+1;
    fits.ReadCell(
	(char*)coeffs_name.c_str(),
	TDOUBLE,
	(char*)&cat.psf[i][0],
	row,
	ncoeff);
  }

  fits.Close();
}




//void WriteFittedPSF(ConfigFile& params, FittedPSF& fpsf);

// This writes to hdu=2
void WriteFittedPSF(ConfigFile& params, FittedPSF& fpsf)
{

  int colnum;
  LONGLONG firstrow;
  LONGLONG firstel;
  LONGLONG nel;
  std::stringstream err;


  std::string file = Name(params,"fitpsf");
  //std::string file = "test_fittedpsf.fits";
  FitsFile fits(file.c_str(), READWRITE, true);
  fitsfile* fptr = fits.get_fptr();


  double sigma = fpsf.GetSigma();

  // Note the actual coeffs may be less than this the way Mike does things
  // but we will fill the extra with zeros
  int psforder = fpsf.GetOrder();
  int fitorder = fpsf.GetFitOrder();

  int npca=fpsf.GetNpca();

  int n_shapelet_coeff = (psforder+1)*(psforder+2)/2;
  int n_rot_matrix_max = n_shapelet_coeff*n_shapelet_coeff;
  int n_rot_matrix_actual = npca*n_shapelet_coeff;

  int n_fit_coeff = (fitorder+1)*(fitorder+2)/2;
  int n_interp_matrix_max = n_shapelet_coeff*n_fit_coeff;
  int n_interp_matrix_actual = npca*n_fit_coeff;

  std::stringstream ave_psf_form;
  ave_psf_form << n_shapelet_coeff << "d";

  std::stringstream rot_matrix_form;
  rot_matrix_form << n_rot_matrix_max << "d";

  std::stringstream interp_matrix_form;
  interp_matrix_form << n_interp_matrix_max << "d";

  std::string psf_order_name=params.get("fitpsf_psf_order_name");
  std::string sigma_name=params.get("fitpsf_sigma_name");
  std::string fit_order_name=params.get("fitpsf_fit_order_name");
  std::string npca_name=params.get("fitpsf_npca_name");

  std::string xmin_name=params.get("fitpsf_xmin_name");
  std::string xmax_name=params.get("fitpsf_xmax_name");
  std::string ymin_name=params.get("fitpsf_ymin_name");
  std::string ymax_name=params.get("fitpsf_ymax_name");

  std::string ave_psf_name=params.get("fitpsf_ave_psf_name");
  std::string rot_matrix_name=params.get("fitpsf_rot_matrix_name");
  std::string interp_matrix_name=params.get("fitpsf_interp_matrix_name");

  int nfields=11;
  char *table_names[] =  
      {(char*)psf_order_name.c_str(),
	(char*)sigma_name.c_str(),
	(char*)fit_order_name.c_str(),
	(char*)npca_name.c_str(),
	(char*)xmin_name.c_str(),
	(char*)xmax_name.c_str(),
	(char*)ymin_name.c_str(),
	(char*)ymax_name.c_str(),
	(char*)ave_psf_name.c_str(),
	(char*)rot_matrix_name.c_str(),
	(char*)interp_matrix_name.c_str()};
  char *table_types[] =  
      {(char*)"1j",   // psf_order
	(char*)"1d",  // sigma
	(char*)"1j",  // fit_order
	(char*)"1j",  // npca
	(char*)"1e",  // xmin  
	(char*)"1e",  // xmax
	(char*)"1e",  // ymin
	(char*)"1e",  // ymax
	(char*)ave_psf_form.str().c_str(),  // ave_psf
	(char*)rot_matrix_form.str().c_str(),  // rot_matrix
	(char*)interp_matrix_form.str().c_str()};  // interp_matrix
  char *table_units[] =  
      {(char*)"None",    // psf_order
	(char*)"arcsec", // sigma
	(char*)"None",   // fit_order
	(char*)"None",   // npca
	(char*)"pixels",   // xmin
	(char*)"pixels",   // xmax
	(char*)"pixels",   // ymin
	(char*)"pixels",   // ymax
	(char*)"None",   // ave_psf
	(char*)"None", // rot_matrix
	(char*)"None"};  // interp_matrix


  // Create a binary table
  int fits_status=0;
  int tbl_type = BINARY_TBL;
  int nrows=1;
  fits_create_tbl(fptr, tbl_type, nrows, nfields, 
      table_names, table_types, table_units, NULL, &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error creating FittedPSF FITS table";
    throw FitsException(serr);
  }

  // dimensions information
  std::stringstream tdim10;
  tdim10<<"("<<n_shapelet_coeff<<","<<n_shapelet_coeff<<")";
  fits_status=0;
  fits_write_key(fptr, TSTRING, "TDIM10", (void*)&tdim10.str()[0], "dimensions of rot_matrix", &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing keyword TDIM10";
    throw FitsException(serr);
  }

  std::stringstream tdim11;
  //tdim11<<"("<<n_fit_coeff<<","<<npca<<")";
  tdim11<<"("<<n_fit_coeff<<","<<n_shapelet_coeff<<")";
  fits_status=0;
  fits_write_key(fptr, TSTRING, "TDIM11", (void*)&tdim11.str()[0], "dimensions of interp_matrix", &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing keyword TDIM11";
    throw FitsException(serr);
  }




  colnum = 1;  
  firstrow = 1;
  firstel = 1;
  nel = 1;
  fits.WriteColumn(TINT, colnum, firstrow, firstel, nel, &psforder);


  colnum = 2;  
  firstrow = 1;
  firstel = 1;
  nel = 1;
  fits.WriteColumn(TDOUBLE, colnum, firstrow, firstel, nel, &sigma);

  colnum = 3;  
  firstrow = 1;
  firstel = 1;
  nel = 1;
  fits.WriteColumn(TINT, colnum, firstrow, firstel, nel, &fitorder);

  colnum = 4;  
  firstrow = 1;
  firstel = 1;
  nel = 1;
  fits.WriteColumn(TINT, colnum, firstrow, firstel, nel, &npca);


  float xmin = fpsf.GetXMin();
  float xmax = fpsf.GetXMax();
  float ymin = fpsf.GetYMin();
  float ymax = fpsf.GetYMax();

  colnum = 5;  
  firstrow = 1;
  firstel = 1;
  nel = 1;
  fits.WriteColumn(TFLOAT, colnum, firstrow, firstel, nel, &xmin);

  colnum = 6;  
  firstrow = 1;
  firstel = 1;
  nel = 1;
  fits.WriteColumn(TFLOAT, colnum, firstrow, firstel, nel, &xmax);

  colnum = 7;  
  firstrow = 1;
  firstel = 1;
  nel = 1;
  fits.WriteColumn(TFLOAT, colnum, firstrow, firstel, nel, &ymin);

  colnum = 8;  
  firstrow = 1;
  firstel = 1;
  nel = 1;
  fits.WriteColumn(TFLOAT, colnum, firstrow, firstel, nel, &ymax);



  double* avg_psf = fpsf.GetAvePSFPtr();
  colnum = 9;  
  firstrow = 1;
  firstel = 1;
  nel = n_shapelet_coeff;
  fits.WriteColumn(TDOUBLE, colnum, firstrow, firstel, nel, avg_psf);


  double* rot_matrix = fpsf.GetRotMatrixPtr();

  //for (int i=0; i<n_rot_matrix_actual; i++) {
  //for (int i=0; i<npca; i++) {
  for (int i=0; i<n_shapelet_coeff; i++) {
    for (int j=0; j<n_shapelet_coeff; j++) {
      if (i >= npca) {
	*(rot_matrix +j*n_shapelet_coeff + i) = 0;
      }
     // std::cout
     //	<<"rot_matrix["<<i<<"]["<<j<<"] = "
     //	<<*(rot_matrix +j*n_shapelet_coeff + i)<<std::endl;
    }
  }

  colnum = 10;  
  firstrow = 1;
  firstel = 1;
  nel = n_rot_matrix_actual;
  nel = n_rot_matrix_max;
  fits.WriteColumn(TDOUBLE, colnum, firstrow, firstel, nel, rot_matrix);


  double* interp_matrix = fpsf.GetInterpMatrixPtr();


  colnum = 11;  
  firstrow = 1;
  firstel = 1;
  nel = n_interp_matrix_actual;
  fits.WriteColumn(TDOUBLE, colnum, firstrow, firstel, nel, interp_matrix);

  fits.Close();

}

void ReadFittedPSF(ConfigFile& params, FittedPSF& fpsf)
{

  int cat_hdu = 2;
  if (params.keyExists("fitpsf_cat_hdu")) cat_hdu = params["fitpsf_cat_hdu"];

  std::string file = Name(params,"fitpsf",true);
  std::cout<< "Reading FittedPSF from file: " << file << std::endl;

  FitsFile fits(file);

  //int hdutype;

  dbg<<"Moving to HDU #"<<cat_hdu<<std::endl;
  fits.GotoHDU(cat_hdu);

  // Allocate memory for the columns we will read
  //ResizePSFCat(cat, nrows, psforder);

  std::string psf_order_name=params.get("fitpsf_psf_order_name");
  std::string sigma_name=params.get("fitpsf_sigma_name");
  std::string fit_order_name=params.get("fitpsf_fit_order_name");
  std::string npca_name=params.get("fitpsf_npca_name");

  std::string xmin_name=params.get("fitpsf_xmin_name");
  std::string xmax_name=params.get("fitpsf_xmax_name");
  std::string ymin_name=params.get("fitpsf_ymin_name");
  std::string ymax_name=params.get("fitpsf_ymax_name");

  std::string ave_psf_name=params.get("fitpsf_ave_psf_name");
  std::string rot_matrix_name=params.get("fitpsf_rot_matrix_name");
  std::string interp_matrix_name=params.get("fitpsf_interp_matrix_name");

  int nrows=1;

  int psf_order=0;
  double sigma;
  int fit_order;
  int npca;

  fits.ReadScalarCol(
      (char *)psf_order_name.c_str(),
      TINT,  (char *)&psf_order, nrows);
  fits.ReadScalarCol(
      (char *)sigma_name.c_str(),
      TDOUBLE,  (char *)&sigma, nrows);
  fits.ReadScalarCol(
      (char *)fit_order_name.c_str(),
      TINT,  (char *)&fit_order, nrows);
  fits.ReadScalarCol(
      (char *)npca_name.c_str(),
      TINT,  (char *)&npca, nrows);

  fpsf.Reset(psf_order, sigma, fit_order, npca);

  float xmin;
  float xmax;
  float ymin;
  float ymax;

  fits.ReadScalarCol(
      (char *)xmin_name.c_str(),
      TFLOAT,  (char *)&xmin, nrows);
  fits.ReadScalarCol(
      (char *)xmax_name.c_str(),
      TFLOAT,  (char *)&xmax, nrows);
  fits.ReadScalarCol(
      (char *)ymin_name.c_str(),
      TFLOAT,  (char *)&ymin, nrows);
  fits.ReadScalarCol(
      (char *)ymax_name.c_str(),
      TFLOAT,  (char *)&ymax, nrows);

  fpsf.SetBounds(xmin, xmax, ymin, ymax);


  double* avg_psf = fpsf.GetAvePSFPtr();
  int n_shapelet_coeff = (psf_order+1)*(psf_order+2)/2;
  fits.ReadCell(
      (char*)ave_psf_name.c_str(),
      TDOUBLE, (char*)avg_psf, 
      1, n_shapelet_coeff);

  double* rot_matrix = fpsf.GetRotMatrixPtr();
  int n_rot_matrix_max = n_shapelet_coeff*n_shapelet_coeff;
  fits.ReadCell(
      (char*)rot_matrix_name.c_str(),
      TDOUBLE, (char*)rot_matrix, 
      1, n_rot_matrix_max);

  double* interp_matrix = fpsf.GetInterpMatrixPtr();
  int n_fit_coeff = (fit_order+1)*(fit_order+2)/2;
  int n_interp_matrix_max = n_shapelet_coeff*n_fit_coeff;
  fits.ReadCell(
      (char*)interp_matrix_name.c_str(),
      TDOUBLE, (char*)interp_matrix, 
      1, n_interp_matrix_max);



  fits.Close();
}




void ResizeShearCat(SHEAR_STRUCT& cat, size_t n, int gal_order, double sigma)
{

  cat.id.resize(n,0);
  cat.shear_flags.resize(n,0);

  cat.shear1.resize(n,0);
  cat.shear2.resize(n,0);


  cat.shear_cov00.resize(n,0);
  cat.shear_cov01.resize(n,0);
  cat.shear_cov11.resize(n,0);

  cat.gal_order.clear();
  cat.gal_order.resize(n,gal_order);

  cat.shapelets_prepsf.resize(n, BVec(gal_order,sigma));

  // These are extra in order to get around Joe's stupidity
#ifdef SHEXTRA_PARS
  cat.sigma0.resize(n,0);
  cat.size_flags.resize(n,0);
  cat.star_flag.resize(n,0);
#endif

}


// This writes to hdu=2
void WriteShearCat(ConfigFile& params, SHEAR_STRUCT& cat)
{

  int colnum;
  LONGLONG firstrow;
  LONGLONG firstel;
  LONGLONG nel;



  std::string file = Name(params,"shear");
  FitsFile fits(file.c_str(), READWRITE, true);
  fitsfile* fptr = fits.get_fptr();


  std::string id_name=params.get("shear_id_name");

  std::string size_flags_name=params.get("shear_size_flags_name");
  std::string star_flag_name=params.get("shear_star_flag_name");
  std::string sigma0_name=params.get("shear_sigma0_name");

  std::string shear_flags_name=params.get("shear_flags_name");

  std::string shear1_name=params.get("shear_shear1_name");
  std::string shear2_name=params.get("shear_shear2_name");

  std::string cov00_name=params.get("shear_cov00_name");
  std::string cov01_name=params.get("shear_cov01_name");
  std::string cov11_name=params.get("shear_cov11_name");

  std::string gal_order_name=params.get("shear_gal_order_name");
  std::string coeffs_name=params.get("shear_coeffs_name");

  int ncoeff = cat.shapelets_prepsf[0].size();
  std::stringstream coeff_form;
  coeff_form << ncoeff << "d";


#ifdef SHEXTRA_PARS
  int nfields=12;
#else
  int nfields=9;
#endif
  char *table_names[] =  
      {(char*)id_name.c_str(),

#ifdef SHEXTRA_PARS
	(char*)size_flags_name.c_str(),
	(char*)star_flag_name.c_str(),
	(char*)sigma0_name.c_str(),
#endif
	(char*)shear_flags_name.c_str(),
	(char*)shear1_name.c_str(),
	(char*)shear2_name.c_str(),
	(char*)cov00_name.c_str(),
	(char*)cov01_name.c_str(),
	(char*)cov11_name.c_str(),
	(char*)gal_order_name.c_str(),
	(char*)coeffs_name.c_str()};
  char *table_types[] =  
      {(char*)"1j",        // id

#ifdef SHEXTRA_PARS
	(char*)"1j",        // size_flags
	(char*)"1j",        // star_flag
	(char*)"1d",        // sigma0
#endif

	(char*)"1j",        // shear_flags
	(char*)"1d",        // shear1
	(char*)"1d",        // shear1
	(char*)"1d",        // shear_cov00
	(char*)"1d",        // shear_cov01
	(char*)"1d",        // shear_cov11
	(char*)"1i",        // gal_order
	(char*)coeff_form.str().c_str()};

  char *table_units[] =  
      {(char*)"None",        // id

#ifdef SHEXTRA_PARS
	(char*)"None",        // size_flags
	(char*)"None",        // star_flag
	(char*)"pixels",        // sigma0
#endif

	(char*)"None",        // shear_flags
	(char*)"None",        // shear1
	(char*)"None",        // shear1
	(char*)"None",        // shear_cov00
	(char*)"None",        // shear_cov01
	(char*)"None",        // shear_cov11
	(char*)"None",        // gal_order
	(char*)"None"};       // coeffs


  // Create a binary table
  int fits_status=0;
  int tbl_type = BINARY_TBL;
  fits_create_tbl(fptr, tbl_type, cat.id.size(), nfields, 
      table_names, table_types, table_units, NULL, &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error creating FindStars FITS table: ";
    throw FitsException(serr);
  }

  fits_status=0;
  fits_write_key(fptr, TINT, "GALORDER", &cat.gal_order[0], "Order of pre-psf shapelets expansion", &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing keyword GALORDER";
    throw FitsException(serr);
  }

  fits_status=0;
  fits_write_key(fptr, TINT, "NCOEFFG", &ncoeff, "Number of coeffs (galorder+1)*(galorder+2)/2", &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing keyword NCOEFFP";
    throw FitsException(serr);
  }



  colnum = 1;  
  firstrow = 1;
  firstel = 1;
  nel = cat.id.size();
  fits.WriteColumn(TLONG, colnum, firstrow, firstel, nel, &cat.id[0]);

#ifdef SHEXTRA_PARS
  colnum++;  
  firstrow = 1;
  firstel = 1;
  nel = cat.size_flags.size();
  fits.WriteColumn(TINT, colnum, firstrow, firstel, nel, &cat.size_flags[0]);

  colnum++;  
  firstrow = 1;
  firstel = 1;
  nel = cat.star_flag.size();
  fits.WriteColumn(TINT, colnum, firstrow, firstel, nel, &cat.star_flag[0]);

  colnum++;  
  firstrow = 1;
  firstel = 1;
  nel = cat.sigma0.size();
  fits.WriteColumn(TDOUBLE, colnum, firstrow, firstel, nel, &cat.sigma0[0]);
#endif

  colnum++;  
  firstrow = 1;
  firstel = 1;
  nel = cat.shear_flags.size();
  fits.WriteColumn(TINT, colnum, firstrow, firstel, nel, &cat.shear_flags[0]);


  colnum++;  
  firstrow = 1;
  firstel = 1;
  nel = cat.shear1.size();
  fits.WriteColumn(TDOUBLE, colnum, firstrow, firstel, nel, &cat.shear1[0]);

  colnum++;  
  firstrow = 1;
  firstel = 1;
  nel = cat.shear2.size();
  fits.WriteColumn(TDOUBLE, colnum, firstrow, firstel, nel, &cat.shear2[0]);


  colnum++;  
  firstrow = 1;
  firstel = 1;
  nel = cat.shear_cov00.size();
  fits.WriteColumn(TDOUBLE, colnum, firstrow, firstel, nel, &cat.shear_cov00[0]);

  colnum++;  
  firstrow = 1;
  firstel = 1;
  nel = cat.shear_cov01.size();
  fits.WriteColumn(TDOUBLE, colnum, firstrow, firstel, nel, &cat.shear_cov01[0]);

  colnum++;  
  firstrow = 1;
  firstel = 1;
  nel = cat.shear_cov11.size();
  fits.WriteColumn(TDOUBLE, colnum, firstrow, firstel, nel, &cat.shear_cov11[0]);

  colnum++;  
  firstrow = 1;
  firstel = 1;
  nel = cat.gal_order.size();
  fits.WriteColumn(TINT, colnum, firstrow, firstel, nel, &cat.gal_order[0]);



  // Now we have to loop through each decomposition and write
  // separately.  This is pretty dumb.
  colnum++;  

  for (size_t i=0; i<cat.id.size(); i++) {
    fits_status=0;
    LONGLONG row = i+1;
    firstel = 1;
    //nel=1;
    nel=cat.shapelets_prepsf[0].size();
    fits_write_col(
	fptr, TDOUBLE, 
	colnum, row, firstel, nel, 
	&cat.shapelets_prepsf[i][0], 
	&fits_status);
    if (!fits_status==0) {
      fits_report_error(stderr, fits_status); 
      std::stringstream serr;
      serr<<"Error writing row "<<row<<" for column 'shapelets_prepsf'";
      throw FitsException(serr.str());
    }
  }






  fits.Close();

}


void ReadShearCat(ConfigFile& params, SHEAR_STRUCT& cat)
{

  int cat_hdu = 2;
  if (params.keyExists("shear_cat_hdu")) cat_hdu = params["shear_cat_hdu"];

  std::string file = Name(params,"shear",true);
  std::cout<< "Reading Shear cat from file: " << file << std::endl;

  FitsFile fits(file);

  //int hdutype;

  dbg<<"Moving to HDU #"<<cat_hdu<<std::endl;
  fits.GotoHDU(cat_hdu);

  // These are not very portable.
  //int nfound;
  long nrows=0;
  fits.ReadKey("NAXIS2", TLONG, (char *)&nrows);

  std::cout<<"  nrows = "<<nrows<<std::endl;
  if (nrows <= 0) {
    std::cerr<<"nrows must be >= 0"<<std::endl;
    throw FAILURE_FORMAT_ERROR;
  }

  int gal_order=0;
  fits.ReadKey("GALORDER", TINT, (char*)&gal_order);

  // Allocate memory for the columns we will read
  ResizeShearCat(cat, nrows, gal_order);

  std::string id_name=params.get("shear_id_name");

#ifdef SHEXTRA_PARS
  std::string size_flags_name=params.get("shear_size_flags_name");
  std::string star_flag_name=params.get("shear_star_flag_name");
  std::string sigma0_name=params.get("shear_sigma0_name");
#endif

  std::string shear_flags_name=params.get("shear_flags_name");

  std::string shear1_name=params.get("shear_shear1_name");
  std::string shear2_name=params.get("shear_shear2_name");

  std::string cov00_name=params.get("shear_cov00_name");
  std::string cov01_name=params.get("shear_cov01_name");
  std::string cov11_name=params.get("shear_cov11_name");

  std::string gal_order_name=params.get("shear_gal_order_name");
  std::string coeffs_name=params.get("shear_coeffs_name");

  /*
  std::cerr<<
    "'"<<id_name<<"'"
    "'"<<shear_flags_name<<"'"
    "'"<<shear1_name<<"'"
    "'"<<shear2_name<<"'"
    "'"<<cov00_name<<"'"
    "'"<<cov01_name<<"'"
    "'"<<cov11_name<<"'"<<std::endl;
    */

  std::cout<<"Reading columns"<<std::endl;
  fits.ReadScalarCol((char *)id_name.c_str(),TLONG,(char *)&cat.id[0], nrows);


#ifdef SHEXTRA_PARS
  fits.ReadScalarCol((char *)size_flags_name.c_str(),TINT,(char *)&cat.size_flags[0], nrows);
  fits.ReadScalarCol((char *)star_flag_name.c_str(),TINT,(char *)&cat.star_flag[0], nrows);
  fits.ReadScalarCol((char *)sigma0_name.c_str(),TDOUBLE,(char *)&cat.sigma0[0], nrows);
#endif

  fits.ReadScalarCol((char *)shear_flags_name.c_str(),TINT,(char *)&cat.shear_flags[0], nrows);

  fits.ReadScalarCol((char *)shear1_name.c_str(),TDOUBLE,(char *)&cat.shear1[0], nrows);
  fits.ReadScalarCol((char *)shear2_name.c_str(),TDOUBLE,(char *)&cat.shear2[0], nrows);

  fits.ReadScalarCol((char *)cov00_name.c_str(),TDOUBLE,(char *)&cat.shear_cov00[0], nrows);
  fits.ReadScalarCol((char *)cov01_name.c_str(),TDOUBLE,(char *)&cat.shear_cov01[0], nrows);
  fits.ReadScalarCol((char *)cov11_name.c_str(),TDOUBLE,(char *)&cat.shear_cov11[0], nrows);

  fits.ReadScalarCol((char *)gal_order_name.c_str(),TINT,(char *)&cat.gal_order[0], nrows);

  // gotta loop for this one
  int ncoeff=(gal_order+1)*(gal_order+2)/2;
  for (size_t i=0; i<cat.id.size(); i++) { 

    size_t row=i+1;
    fits.ReadCell(
	(char*)coeffs_name.c_str(),
	TDOUBLE,
	(char*)&cat.shapelets_prepsf[i][0],
	row,
	ncoeff);
  }



  fits.Close();
}



