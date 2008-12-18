
#include "FittedPSF.h"
#include "TMV_Diag.h"
#include "dbg.h"
#include "Function2D.h"
#include "Legendre2D.h"
#include "Name.h"
#include <fstream>

static tmv::Vector<double> DefinePXY(size_t order, double x,
    double xmin, double xmax)
{
  tmv::Vector<double> temp(order+1);
  double newx = (2.*x-xmin-xmax)/(xmax-xmin);
  temp[0] = 1.;
  if(order>0) temp[1] = newx;
  for(size_t i=2;i<=order;i++)
    temp[i] = ((2.*i-1.)*newx*temp[i-1] - (i-1.)*temp[i-2])/i;
  return temp;
}

static void SetPRow(size_t fitorder, Position pos, const Bounds& bounds, 
    const tmv::VectorView<double>& Prow)
{
  Assert(Prow.size() == (fitorder+1)*(fitorder+2)/2);
  tmv::Vector<double> px = DefinePXY(fitorder,pos.GetX(),
      bounds.GetXMin(),bounds.GetXMax());
  tmv::Vector<double> py = DefinePXY(fitorder,pos.GetY(),
      bounds.GetYMin(),bounds.GetYMax());
  size_t pq = 0;
  for(size_t pplusq=0;pplusq<=fitorder;pplusq++)
    for(size_t p=pplusq,q=pplusq-p;q<=pplusq;p--,q++) {
      Assert(pq < Prow.size());
      Prow(pq) = px[p]*py[q];
      pq++;
    }
  Assert(pq == Prow.size());
}

FittedPSF::FittedPSF(const std::vector<BVec>& psf,
    const std::vector<long>& flagvec,
    const std::vector<Position>& pos,
    const std::vector<double>& nu,
    double _sigma, const ConfigFile& params) :
  psforder(params["psf_order"]), sigma(_sigma),
  fitorder(params["fitpsf_order"]), fitsize((fitorder+1)*(fitorder+2)/2)
  //fitorder(2), fitsize((fitorder+1)*(fitorder+2)/2)
{
  xdbg<<"FittedPSF constructor\n";
  // Do a polynomial fit of the psf shapelet vectors
  //
  // First, rotate the vectors into their eigen directions.
  // The matrix V is stored to let us get back to the original basis.
  avepsf.reset(new BVec(psforder,sigma));
  avepsf->Zero();
  size_t psfsize = avepsf->size();
  size_t ngoodpsf=0;
  for(size_t n=0;n<psf.size();n++) if (!flagvec[n]) {
    xxdbg<<"n = "<<n<<", psf[n] = "<<psf[n]<<std::endl;
    *avepsf += psf[n];
    ngoodpsf++;
  }
  *avepsf /= double(ngoodpsf);
  xdbg<<"ngoodpsf = "<<ngoodpsf<<std::endl;
  xdbg<<"average psf = "<<*avepsf<<std::endl;
  tmv::Matrix<double> M(ngoodpsf,psfsize);
  tmv::DiagMatrix<double> invsig(ngoodpsf);
  size_t i=0;
  for(size_t n=0;n<psf.size();n++) if (!flagvec[n]) {
    Assert(psf[n].size() == psfsize);
    Assert(i < ngoodpsf);
    M.row(i) = psf[n] - *avepsf;
    invsig(i) = nu[n];
    bounds += pos[n];
    i++;
  }
  Assert(i == ngoodpsf);
  xdbg<<"bounds = "<<bounds<<std::endl;
  xxdbg<<"M = "<<M<<std::endl;
  M = invsig * M;
  //tmv::Matrix<double> M1 = M;

  size_t K = std::min(ngoodpsf,psfsize);
  tmv::DiagMatrix<double> S(K);
  tmv::MatrixView<double> U = M.Cols(0,K);
  V.reset(new tmv::Matrix<double,tmv::RowMajor>(K,psfsize));
  if (ngoodpsf > psfsize) {
    SV_Decompose(U.View(),S.View(),V->View(),true);
  } else {
    *V = M;
    SV_Decompose(V->Transpose(),S.View(),U.Transpose());
  }
  xdbg<<"In FittedPSF: SVD S = "<<S.diag()<<std::endl;
  xxdbg<<"U = "<<U<<std::endl;
  xxdbg<<"V = "<<*V<<std::endl;
  //xdbg<<"M-USV = "<<M1-U*S*(*V)<<std::endl;
  if (params.keyExists("fitpsf_npca")) {
    npca = params["fitpsf_npca"];
    xdbg<<"npca = "<<npca<<" from parameter file\n";
  } else {
    double thresh = S(0);
    if (params.keyExists("fitpsf_pca_thresh")) 
      thresh *= double(params["fitpsf_pca_thresh"]);
    else thresh *= tmv::Epsilon<double>();
    xdbg<<"thresh = "<<thresh<<std::endl;
    for(npca=1;npca<int(M.rowsize());npca++) if (S(npca) < thresh) break;
    xdbg<<"npca = "<<npca<<std::endl;
  }
  U.Cols(0,npca) *= S.SubDiagMatrix(0,npca);
  xxdbg<<"US = "<<U.Cols(0,npca)<<std::endl;
  //dbg<<"MVt = "<<M1*V->Rows(0,npca).Adjoint()<<std::endl;
  // U S = M(orig) * Vt

  xxdbg<<"V.Rows(0,npca) = "<<V->Rows(0,npca)<<std::endl;
  tmv::Matrix<double> P(ngoodpsf,fitsize,0.);
  i=0;
  for(size_t n=0;n<psf.size();n++) if (!flagvec[n]) {
    SetPRow(fitorder,pos[n],bounds,P.row(i));
    i++;
  }
  Assert(i == ngoodpsf);
  xxdbg<<"P = "<<P<<std::endl;
  P = invsig * P;

  // TODO: need to implement outlier rejection.
  f.reset(new tmv::Matrix<double>(U.Cols(0,npca)/P));
  xxdbg<<"f = "<<*f<<std::endl;
  xxdbg<<"P*f = "<<P*(*f)<<std::endl;
  xdbg<<"Done making FittedPSF\n";
}

void FittedPSF::Write(std::ostream& os) const
{
  os << psforder <<"  "<< sigma <<"  ";
  os << fitorder <<"  "<< npca << std::endl;
  os << bounds << std::endl;
  os << *avepsf << std::endl;
  os << V->Rows(0,npca) <<std::endl;
  os << *f << std::endl;
}

void FittedPSF::Read(std::istream& is)
{
  xdbg<<"Reading FittedPSF:\n";
  is >> psforder >> sigma >> fitorder >> npca >> bounds;
  xdbg<<"psforder = "<<psforder<<", sigma_psf = "<<sigma<<std::endl;
  xdbg<<"fitorder = "<<fitorder<<", npca = "<<npca<<std::endl;
  xdbg<<"bounds = "<<bounds<<std::endl;
  fitsize = (fitorder+1)*(fitorder+2)/2;
  avepsf.reset(new BVec(psforder,sigma));
  is >> *avepsf;
  xxdbg<<"avepsf = "<<*avepsf<<std::endl;
  V.reset(new tmv::Matrix<double,tmv::RowMajor>(npca,avepsf->size()));
  is >> *V;
  xxdbg<<"V = "<<*V<<std::endl;
  f.reset(new tmv::Matrix<double>(fitsize,npca));
  is >> *f;
  xxdbg<<"f = "<<*f<<std::endl;
}


void FittedPSF::Write(const ConfigFile& params) const
{
  std::string file = Name(params,"fitpsf");
  // MJ -- check which kind of write from ext or other way.
  if (1) {
    WriteFits(file,params);
  } else {
    std::ofstream fout(file.c_str());
    Write(fout);
    fout.close();
  }
}

void FittedPSF::Read(const ConfigFile& params)
{
  std::string file = Name(params,"fitpsf",false,true);
  dbg<< "Reading FittedPSF from file: " << file << std::endl;
  // MJ -- check which kind of read from ext or other way.
  if (1) {
    int hdu = 2;
    if (params.keyExists("fitpsf_hdu")) hdu = params["fitpsf_hdu"];
    ReadFits(file,hdu,params);
  } else {
    std::ifstream fin(file.c_str());
    Assert(fin);
    Read(fin);
    fin.close();
  }
}


void FittedPSF::WriteFitsKeywords(FitsFile& fits, const ConfigFile& params) const
{
  fits.WriteParKey(params, "version", TSTRING);
  fits.WriteParKey(params, "noise_method", TSTRING);
  fits.WriteParKey(params, "dist_method", TSTRING);
  fits.WriteParKey(params, "fitpsf_order", TLONG);
  fits.WriteParKey(params, "fitpsf_pca_thresh", TDOUBLE);
} 
void FittedPSF::WriteFits(std::string file, const ConfigFile& params) const
{
  int colnum;
  LONGLONG firstrow;
  LONGLONG firstel;
  LONGLONG nel;
  std::stringstream err;

  FitsFile fits(file.c_str(), READWRITE, true);
  fitsfile* fptr = fits.get_fptr();

  // Note the actual coeffs may be less than this the way Mike does things
  // but we will fill the extra with zeros
  int n_shapelet_coeff = (psforder+1)*(psforder+2)/2;
  //int n_rot_matrix_max = n_shapelet_coeff*n_shapelet_coeff;
  int n_rot_matrix = npca*n_shapelet_coeff;

  int n_fit_coeff = (fitorder+1)*(fitorder+2)/2;
  //int n_interp_matrix_max = n_shapelet_coeff*n_fit_coeff;
  int n_interp_matrix = npca*n_fit_coeff;

  std::stringstream ave_psf_form;
  ave_psf_form << n_shapelet_coeff << "d";

  std::stringstream rot_matrix_form;
  //rot_matrix_form << n_rot_matrix_max << "d";
  rot_matrix_form << n_rot_matrix << "d";

  std::stringstream interp_matrix_form;
  //interp_matrix_form << n_interp_matrix_max << "d";
  interp_matrix_form << n_interp_matrix << "d";

#if 0
  // I think this should work, but it didn't on my AMD machine with pgCC.
  char* ave_psf_form_cstr = (char*) ave_psf_form.str().c_str();
  char* rot_matrix_form_cstr = (char*) rot_matrix_form.str().c_str();
  char* interp_matrix_form_cstr = (char*) interp_matrix_form.str().c_str();
#else
  Assert(ave_psf_form.str().size() < 9);       // Usually size = 3
  Assert(rot_matrix_form.str().size() < 9);    // Usually size = 4 or 5
  Assert(interp_matrix_form.str().size() < 9); // Usually size = 4 or 5
  char ave_psf_form_cstr[10];  
  strcpy(ave_psf_form_cstr,ave_psf_form.str().c_str());
  char rot_matrix_form_cstr[10];  
  strcpy(rot_matrix_form_cstr,rot_matrix_form.str().c_str());
  char interp_matrix_form_cstr[10];  
  strcpy(interp_matrix_form_cstr,interp_matrix_form.str().c_str());
#endif

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
      {(char*)"1j",   // psforder
	(char*)"1d",  // sigma
	(char*)"1j",  // fitorder
	(char*)"1j",  // npca
	(char*)"1e",  // xmin  
	(char*)"1e",  // xmax
	(char*)"1e",  // ymin
	(char*)"1e",  // ymax
	(char*)ave_psf_form_cstr,  // ave_psf
	(char*)rot_matrix_form_cstr,  // rot_matrix
	(char*)interp_matrix_form_cstr};  // interp_matrix
  char *table_units[] =  
      {(char*)"None",    // psforder
	(char*)"arcsec", // sigma
	(char*)"None",   // fitorder
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
  fits.WriteKey("TDIM10",TSTRING,tdim10.str().c_str(),
      "dimensions of rot_matrix");


  std::stringstream tdim11;
  tdim11<<"("<<n_fit_coeff<<","<<n_shapelet_coeff<<")";
  fits.WriteKey("TDIM11",TSTRING,tdim11.str().c_str(),
      "dimensions of interp_matrix");

  WriteFitsKeywords(fits, params);

  float xmin = bounds.GetXMin();
  float xmax = bounds.GetXMax();
  float ymin = bounds.GetYMin();
  float ymax = bounds.GetYMax();

  colnum = 1;  
  firstrow = 1;
  firstel = 1;
  nel = 1;
  fits.WriteColumn(TLONG, colnum, firstrow, firstel, nel, &psforder);


  colnum = 2;  
  firstrow = 1;
  firstel = 1;
  nel = 1;
  fits.WriteColumn(TDOUBLE, colnum, firstrow, firstel, nel, &sigma);

  colnum = 3;  
  firstrow = 1;
  firstel = 1;
  nel = 1;
  fits.WriteColumn(TLONG, colnum, firstrow, firstel, nel, &fitorder);

  colnum = 4;  
  firstrow = 1;
  firstel = 1;
  nel = 1;
  fits.WriteColumn(TLONG, colnum, firstrow, firstel, nel, &npca);


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

  colnum = 9;  
  firstrow = 1;
  firstel = 1;
  nel = n_shapelet_coeff;
  fits.WriteColumn(TDOUBLE, colnum, firstrow, firstel, nel, avepsf->cptr());

#if 0
  for (int i=0; i<n_shapelet_coeff; i++) {
    for (int j=0; j<n_shapelet_coeff; j++) {
      if (i >= npca) {
	*(rot_matrix +j*n_shapelet_coeff + i) = 0;
      }
     // dbg
     //	<<"rot_matrix["<<i<<"]["<<j<<"] = "
     //	<<*(rot_matrix +j*n_shapelet_coeff + i)<<std::endl;
    }
  }
#endif
  colnum = 10;  
  firstrow = 1;
  firstel = 1;
  nel = n_rot_matrix;
  //nel = n_rot_matrix_max;
  fits.WriteColumn(TDOUBLE, colnum, firstrow, firstel, nel, V->ptr());

  colnum = 11;  
  firstrow = 1;
  firstel = 1;
  nel = n_interp_matrix;
  fits.WriteColumn(TDOUBLE, colnum, firstrow, firstel, nel, f->ptr());

  fits.Close();

  // Test FITS I/O:
  if (0) {
    FittedPSF test;
    test.Read(params);
    // Compare  
    dbg<<"   orig psforder: "<<GetOrder()<<" new: "<<test.GetOrder()<<std::endl;
    dbg<<"   orig fitorder: "<<GetFitOrder()<<" new: "<<test.GetFitOrder()<<std::endl;
    dbg<<"   orig npca: "<<GetNpca()<<" new: "<<test.GetNpca()<<std::endl;
    dbg<<"   orig sigma: "<<GetSigma()<<" new: "<<test.GetSigma()<<std::endl;
    dbg<<"   orig xmin: "<<GetXMin()<<" new: "<<test.GetXMin()<<std::endl;
    dbg<<"   orig xmax: "<<GetXMax()<<" new: "<<test.GetXMax()<<std::endl;
    dbg<<"   orig ymin: "<<GetYMin()<<" new: "<<test.GetYMin()<<std::endl;
    dbg<<"   orig ymax: "<<GetYMax()<<" new: "<<test.GetYMax()<<std::endl;

    dbg<<"\n";
    dbg<<"orig avepsf = "<<(*avepsf)<<std::endl;
    dbg<<"     new = "<<(*test.avepsf)<<std::endl;
    dbg<<"     diff = "<<(*avepsf)-(*test.avepsf)<<std::endl;
    dbg<<"     Norm(diff) = "<<Norm((*avepsf)-(*test.avepsf))<<std::endl;
    if (Norm((*avepsf)-(*test.avepsf)) > 0.001) 
      throw std::runtime_error("Error in FittedPSF Fits I/O");

    dbg<<"\n";
    dbg<<"orig rot_matrix = "<<V->Rows(0,npca);
    dbg<<"     new =  "<<(*test.V);
    dbg<<"     diff =  "<<V->Rows(0,npca)-(*test.V)<<std::endl;
    dbg<<"     Norm(diff) =  "<<Norm(V->Rows(0,npca)-(*test.V))<<std::endl;
    if (Norm(V->Rows(0,npca)-(*test.V)) > 0.001) 
      throw std::runtime_error("Error in FittedPSF Fits I/O");

    dbg<<"\n";
    dbg<<"orig interp_matrix = "<<(*f);
    dbg<<"     new: "<<(*test.f);
    dbg<<"     diff: "<<(*f)-(*test.f);
    dbg<<"     Norm(diff): "<<Norm((*f)-(*test.f));
    if (Norm((*f)-(*test.f)) > 0.001) 
      throw std::runtime_error("Error in FittedPSF Fits I/O");

    dbg<<"Done testing fitted PSF file IO\n";
  }
}

void FittedPSF::ReadFits(std::string file, int hdu, const ConfigFile& params)
{
  FitsFile fits(file);

  dbg<<"Moving to HDU #"<<hdu<<std::endl;
  fits.GotoHDU(hdu);

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

  fits.ReadScalarCol(
      (char *)psf_order_name.c_str(),
      TLONG,  (char *)&psforder, nrows);
  xdbg<<"psforder = "<<psforder<<std::endl;
  fits.ReadScalarCol(
      (char *)sigma_name.c_str(),
      TDOUBLE,  (char *)&sigma, nrows);
  xdbg<<"sigma = "<<sigma<<std::endl;
  fits.ReadScalarCol(
      (char *)fit_order_name.c_str(),
      TLONG,  (char *)&fitorder, nrows);
  xdbg<<"fitorder = "<<fitorder<<std::endl;
  fits.ReadScalarCol(
      (char *)npca_name.c_str(),
      TLONG,  (char *)&npca, nrows);
  xdbg<<"npca = "<<npca<<std::endl;

  float xmin;
  float xmax;
  float ymin;
  float ymax;

  fits.ReadScalarCol(
      (char *)xmin_name.c_str(),
      TFLOAT,  (char *)&xmin, nrows);
  xdbg<<"xmin = "<<xmin<<std::endl;
  fits.ReadScalarCol(
      (char *)xmax_name.c_str(),
      TFLOAT,  (char *)&xmax, nrows);
  xdbg<<"xmax = "<<xmax<<std::endl;
  fits.ReadScalarCol(
      (char *)ymin_name.c_str(),
      TFLOAT,  (char *)&ymin, nrows);
  xdbg<<"ymin = "<<ymin<<std::endl;
  fits.ReadScalarCol(
      (char *)ymax_name.c_str(),
      TFLOAT,  (char *)&ymax, nrows);
  xdbg<<"ymax = "<<ymax<<std::endl;

  bounds.SetXMin(xmin);
  bounds.SetXMax(xmax);
  bounds.SetYMin(ymin);
  bounds.SetYMax(ymax);
  xdbg<<"bounds = "<<bounds<<std::endl;

  fitsize = (fitorder+1)*(fitorder+2)/2;
  xdbg<<"fitsize = "<<fitsize<<std::endl;
  avepsf.reset(new BVec(psforder,sigma));
  V.reset(new tmv::Matrix<double,tmv::RowMajor>(npca,avepsf->size()));
  f.reset(new tmv::Matrix<double>(fitsize,npca));

  int n_shapelet_coeff = (psforder+1)*(psforder+2)/2;
  fits.ReadCell(
      (char*)ave_psf_name.c_str(),
      TDOUBLE, (char*)(avepsf->ptr()),
      1, n_shapelet_coeff);
  xdbg<<"avepsf = "<<*avepsf<<std::endl;

  int n_rot_matrix = npca*n_shapelet_coeff;
  fits.ReadCell(
      (char*)rot_matrix_name.c_str(),
      TDOUBLE, (char*)(V->ptr()),
      1, n_rot_matrix);
  xdbg<<"V = "<<*V<<std::endl;

  int n_fit_coeff = (fitorder+1)*(fitorder+2)/2;
  int n_interp_matrix = npca*n_fit_coeff;
  fits.ReadCell(
      (char*)interp_matrix_name.c_str(),
      TDOUBLE, (char*)(f->ptr()),
      1, n_interp_matrix);
  xdbg<<"f = "<<*f<<std::endl;

  fits.Close();
}

void FittedPSF::InterpolateVector(
    Position pos, const tmv::VectorView<double>& b) const
{
  tmv::Vector<double> P(fitsize);
  SetPRow(fitorder,pos,bounds,P.View());
  tmv::Vector<double> b1 = P * (*f);
  b = b1 * V->Rows(0,npca);
  b += *avepsf;
}

