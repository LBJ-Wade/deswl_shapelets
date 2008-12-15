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
  cat.ra.resize(n,0);
  cat.dec.resize(n,0);

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
  if (params.keyExists("cat_minrows")) minrows = params["cat_minrows"];

  dbg<< "Reading cat from file: " << file << std::endl;

  FitsFile fits(file);

  //int hdutype;

  dbg<<"Moving to HDU #"<<cat_hdu<<std::endl;
  fits.GotoHDU(cat_hdu);

  // These are not very portable.
  //int nfound;
  long nrows=0;
  fits.ReadKey("NAXIS2", TLONG, (char *)&nrows);

  dbg<<"  nrows = "<<nrows<<std::endl;
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
  std::string id_name=params.get("cat_id_name");
  std::string x_name=params.get("cat_x_name");
  std::string y_name=params.get("cat_y_name");
  std::string local_sky_name=params.get("cat_local_sky_name");
  std::string mag_name=params.get("cat_mag_name");
  std::string mag_err_name=params.get("cat_mag_err_name");
  std::string flags_name=params.get("cat_flags_name");
  std::string ra_name=params["cat_ra_name"];
  std::string dec_name=params["cat_dec_name"];

  dbg<<"Reading columns"<<std::endl;
  dbg<<"  "<<id_name<<std::endl;
  fits.ReadScalarCol((char *)id_name.c_str(),TLONG,(char *)&cat.id[0], nrows);
  dbg<<"  "<<x_name<<std::endl;
  fits.ReadScalarCol((char *)x_name.c_str(),TFLOAT,(char *)&cat.x[0], nrows);
  dbg<<"  "<<y_name<<std::endl;
  fits.ReadScalarCol((char *)y_name.c_str(),TFLOAT,(char *)&cat.y[0], nrows);
  dbg<<"  "<<local_sky_name<<std::endl;
  fits.ReadScalarCol((char *)local_sky_name.c_str(),TDOUBLE,(char *)&cat.local_sky[0], nrows);
  dbg<<"  "<<mag_name<<std::endl;
  fits.ReadScalarCol((char *)mag_name.c_str(),TFLOAT,(char *)&cat.mag[0], nrows);
  dbg<<"  "<<mag_err_name<<std::endl;
  fits.ReadScalarCol((char *)mag_err_name.c_str(),TFLOAT,(char *)&cat.mag_err[0], nrows);
  dbg<<"  "<<flags_name<<std::endl;
  fits.ReadScalarCol((char *)flags_name.c_str(),TSHORT,(char *)&cat.flags[0], nrows);
  dbg<<"  "<<ra_name<<std::endl;
  fits.ReadScalarCol((char *)ra_name.c_str(),TFLOAT,(char *)&cat.ra[0], nrows);
  dbg<<"  "<<dec_name<<std::endl;
  fits.ReadScalarCol((char *)dec_name.c_str(),TFLOAT,(char *)&cat.dec[0], nrows);

  for (long i=0; i< nrows; i++) {
    cat.pos[i] = Position(cat.x[i], cat.y[i]);
  }

  fits.Close();
}


// This writes to hdu=2
void WriteFindStarsCat(ConfigFile& params, FINDSTARS_STRUCT& cat)
{

  int colnum;
  LONGLONG firstrow;
  LONGLONG firstel;
  LONGLONG nel;

  std::string file = Name(params,"stars");
  FitsFile fits(file.c_str(), READWRITE, true);
  fitsfile* fptr = fits.get_fptr();

  std::string id_name=params.get("stars_id_name");
  std::string sigma0_name=params.get("stars_sigma0_name");
  std::string size_flags_name=params.get("stars_size_flags_name");
  std::string star_flag_name=params.get("stars_star_flag_name");

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


  colnum = 1;  
  firstrow = 1;
  firstel = 1;
  nel = cat.id.size();
  fits.WriteColumn(TLONG, colnum, firstrow, firstel, nel, &cat.id[0]);

  colnum = 2;  
  fits.WriteColumn(TDOUBLE, colnum, firstrow, firstel, nel, &cat.sigma0[0]);
  colnum = 3;  
  fits.WriteColumn(TINT, colnum, firstrow, firstel, nel, &cat.size_flags[0]);
  colnum = 4;  
  fits.WriteColumn(TINT, colnum, firstrow, firstel, nel, &cat.star_flag[0]);

  fits.Close();

}


void ReadFindStarsCat(ConfigFile& params, FINDSTARS_STRUCT& cat)
{

  int cat_hdu = 2;
  if (params.keyExists("stars_cat_hdu")) cat_hdu = params["stars_cat_hdu"];

  std::string file = Name(params,"stars",false,true);
  // false,true = input_prefix=false, mustexist=true.
  // It is an input here, but it is in the output_prefix directory.
  dbg<< "Reading FindStars cat from file: " << file << std::endl;

  FitsFile fits(file);

  //int hdutype;

  dbg<<"Moving to HDU #"<<cat_hdu<<std::endl;
  fits.GotoHDU(cat_hdu);

  // These are not very portable.
  //int nfound;
  long nrows=0;
  fits.ReadKey("NAXIS2", TLONG, (char *)&nrows);

  dbg<<"  nrows = "<<nrows<<std::endl;
  if (nrows <= 0) {
    std::cerr<<"nrows must be >= 0"<<std::endl;
    throw FAILURE_FORMAT_ERROR;
  }

  // Allocate memory for the columns we will read
  ResizeFindStarsCat(cat, nrows);

  std::string id_name=params.get("stars_id_name");
  std::string sigma0_name=params.get("stars_sigma0_name");
  std::string size_flags_name=params.get("stars_size_flags_name");
  std::string star_flag_name=params.get("stars_star_flag_name");

  dbg<<"Reading columns"<<std::endl;
  dbg<<"  "<<id_name<<std::endl;
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


  fits.WriteKey("PSFORDER", TINT, cat.psf_order[0], 
      "Order of shapelets expansion for PSF stars");
  fits.WriteKey("NCOEFFP", TINT, ncoeff, 
      "Number of coeffs (psforder+1)*(psforder+2)/2");
  fits.WriteKey("SIGMA_P", TDOUBLE, cat.sigma_p[0], 
      "Scale used in shapelets expansion");


  fits_status=0;
  colnum = 1;  
  firstrow = 1;
  firstel = 1;
  nel = cat.id.size();
  fits.WriteColumn(TLONG, colnum, firstrow, firstel, nel, &cat.id[0]);

  
  fits_status=0;
  colnum = 2;  
  firstrow = 1;
  firstel = 1;
  nel = cat.psf_flags.size();
  fits.WriteColumn(TINT, colnum, firstrow, firstel, nel, &cat.psf_flags[0]);


  fits_status=0;
  colnum = 3;  
  firstrow = 1;
  firstel = 1;
  nel = cat.nu.size();
  fits.WriteColumn(TDOUBLE, colnum, firstrow, firstel, nel, &cat.nu[0]);


  fits_status=0;
  colnum = 4;  
  firstrow = 1;
  firstel = 1;
  nel = cat.psf_order.size();
  fits.WriteColumn(TINT, colnum, firstrow, firstel, nel, &cat.psf_order[0]);

  


  fits_status=0;
  colnum = 5;  
  firstrow = 1;
  firstel = 1;
  nel = cat.sigma_p.size();
  fits.WriteColumn(TDOUBLE, colnum, firstrow, firstel, nel, &cat.sigma_p[0]);
  
  // Now we have to loop through each psf decomposition and write
  // separately.  This is pretty dumb.
  colnum = 6;  

  for (size_t i=0; i<cat.id.size(); i++) {
    fits_status=0;
    LONGLONG row = i+1;
    firstel = 1;
    nel=cat.psf[0].size();
    fits.WriteColumn(TDOUBLE, colnum, row, firstel, nel, &cat.psf[i][0]);
  }

  fits.Close();

}

void ReadPSFCat(ConfigFile& params, PSF_STRUCT& cat)
{

  int cat_hdu = 2;
  if (params.keyExists("psf_cat_hdu")) cat_hdu = params["psf_cat_hdu"];

  std::string file = Name(params,"psf",false,true);
  dbg<< "Reading PSF cat from file: " << file << std::endl;

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

  dbg<<"  nrows = "<<nrows<<std::endl;
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

  dbg<<"Reading columns"<<std::endl;
  dbg<<"  "<<id_name<<std::endl;
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


  fits.WriteKey("GALORDER", TINT, cat.gal_order[0], 
      "Order of pre-psf shapelets expansions");
  fits.WriteKey("NCOEFFG", TINT, ncoeff, 
      "Number of coeffs (galorder+1)*(galorder+2)/2");

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

  std::string file = Name(params,"shear",false,true);
  dbg<< "Reading Shear cat from file: " << file << std::endl;

  FitsFile fits(file);

  //int hdutype;

  dbg<<"Moving to HDU #"<<cat_hdu<<std::endl;
  fits.GotoHDU(cat_hdu);

  // These are not very portable.
  //int nfound;
  long nrows=0;
  fits.ReadKey("NAXIS2", TLONG, (char *)&nrows);

  dbg<<"  nrows = "<<nrows<<std::endl;
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

  dbg<<"Reading columns"<<std::endl;
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


