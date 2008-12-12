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







// This writes to hdu=2
void WriteFindStarsCat(ConfigFile& params, FINDSTARS_STRUCT& cat)
{

  std::string file = Name(params,"allcat");
  FitsFile fits(file.c_str(), READWRITE, true);
  fitsfile* fptr = fits.get_fptr();

  int nfields=4;
  char *table_names[] =  
      {(char*)"id",
	(char*)"sigma0",
	(char*)"size_flags",
	(char*)"star_flag"};
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
  std::cerr<< "Reading FindStars cat from file: " << file << std::endl;

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

  // Allocate memory for the columns we will read
  ResizeFindStarsCat(cat, nrows);

  std::string id_name=params["fs_id_name"];
  std::string sigma0_name=params["fs_sigma0_name"];
  std::string size_flags_name=params["fs_size_flags_name"];
  std::string star_flag_name=params["fs_star_flag_name"];

  std::cerr<<"Reading columns"<<std::endl;
  std::cerr<<"  "<<id_name<<std::endl;
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

  int nfields=6;
  char *table_names[] =  
      {(char*)"id",
	(char*)"psf_flags",
	(char*)"nu",
	(char*)"psf_order",
	(char*)"sigma_p",
	(char*)"coeffs"};
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
    std::string serr="Error writing keyword PSF_ORDER";
    throw FitsException(serr);
  }

  fits_status=0;
  fits_write_key(fptr, TINT, "NCOEFF", &ncoeff, "Number of coeffs (psforder+1)*(psforder+2)/2", &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error writing keyword NCOEFF";
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
  std::cerr<< "Reading PSF cat from file: " << file << std::endl;

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

  std::cerr<<"  nrows = "<<nrows<<std::endl;
  if (nrows <= 0) {
    std::cerr<<"nrows must be > 0"<<std::endl;
    throw FAILURE_FORMAT_ERROR;
  }

  // Allocate memory for the columns we will read
  ResizePSFCat(cat, nrows, psforder);

  std::string id_name=params["psf_id_name"];
  std::string psf_flags_name=params["psf_flags_name"];
  std::string nu_name=params["psf_nu_name"];
  std::string psf_order_name=params["psf_order_name"];
  std::string sigma_p_name=params["psf_sigma_p_name"];
  std::string coeffs_name=params["psf_coeffs_name"];

  std::cerr<<"Reading columns"<<std::endl;
  std::cerr<<"  "<<id_name<<std::endl;
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


  int nfields=11;
  char *table_names[] =  
      {(char*)"psf_order",
	(char*)"sigma",
	(char*)"fit_order",
	(char*)"npca",
	(char*)"xmin",
	(char*)"xmax",
	(char*)"ymin",
	(char*)"ymax",
	(char*)"ave_psf",
	(char*)"rot_matrix",
	(char*)"interp_matrix"};
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
  std::cerr<< "Reading FittedPSF from file: " << file << std::endl;

  FitsFile fits(file);

  //int hdutype;

  dbg<<"Moving to HDU #"<<cat_hdu<<std::endl;
  fits.GotoHDU(cat_hdu);

  // Allocate memory for the columns we will read
  //ResizePSFCat(cat, nrows, psforder);

  std::string psf_order_name=params["fitpsf_psf_order_name"];
  std::string sigma_name=params["fitpsf_sigma_name"];
  std::string fit_order_name=params["fitpsf_fit_order_name"];
  std::string npca_name=params["fitpsf_npca_name"];

  std::string xmin_name=params["fitpsf_xmin_name"];
  std::string xmax_name=params["fitpsf_xmax_name"];
  std::string ymin_name=params["fitpsf_ymin_name"];
  std::string ymax_name=params["fitpsf_ymax_name"];

  std::string ave_psf_name=params["fitpsf_ave_psf_name"];
  std::string rot_matrix_order_name=params["fitpsf_rot_matrix_name"];
  std::string interp_matrix_name=params["fitpsf_interp_matrix_name"];

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


/*
  double* avg_psf = fpsf.GetAvePSFPtr();
  int n_shapelet_coeff = (psf_order+1)*(psf_order+2)/2;
  int n_rot_matrix_max = n_shapelet_coeff*n_shapelet_coeff;
  */


  fits.Close();
}


