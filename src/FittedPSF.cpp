
#include "FittedPSF.h"
#include "FitsFile.h"
#include "TMV.h"
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

FittedPSF::FittedPSF(const PSFCatalog& psfcat, const ConfigFile& _params) :
  params(_params), psforder(params.get("psf_order")),
  fitorder(params.get("fitpsf_order")), fitsize((fitorder+1)*(fitorder+2)/2)
{
  xdbg<<"FittedPSF constructor\n";
  // Do a polynomial fit of the psf shapelet vectors

  Assert(psfcat.psf.size() > 0);
  sigma = psfcat.psf[0].GetSigma();

  // Calculate the average psf vector
  avepsf.reset(new BVec(psforder,sigma));
  avepsf->Zero();
  size_t psfsize = avepsf->size();
  size_t ngoodpsf=0;
  for(size_t n=0;n<psfcat.size();n++) if (!psfcat.flags[n]) {
    xxdbg<<"n = "<<n<<", psf[n] = "<<psfcat.psf[n]<<std::endl;
    Assert(psfcat.psf[n].GetSigma() == sigma);
    *avepsf += psfcat.psf[n];
    ngoodpsf++;
  }
  *avepsf /= double(ngoodpsf);
  xdbg<<"ngoodpsf = "<<ngoodpsf<<std::endl;
  xdbg<<"average psf = "<<*avepsf<<std::endl;

  // Rotate the vectors into their eigen directions.
  // The matrix V is stored to let us get back to the original basis.
  tmv::Matrix<double> M(ngoodpsf,psfsize);
  tmv::DiagMatrix<double> invsig(ngoodpsf);
  size_t i=0;
  for(size_t n=0;n<psfcat.size();n++) if (!psfcat.flags[n]) {
    Assert(psfcat.psf[n].size() == psfsize);
    Assert(i < ngoodpsf);
    M.row(i) = psfcat.psf[n] - *avepsf;
    invsig(i) = psfcat.nu[n];
    bounds += psfcat.pos[n];
    i++;
  }
  Assert(i == ngoodpsf);
  xdbg<<"bounds = "<<bounds<<std::endl;
  M = invsig * M;

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
  // U S = M(orig) * Vt

  tmv::Matrix<double> P(ngoodpsf,fitsize,0.);
  i=0;
  for(size_t n=0;n<psfcat.size();n++) if (!psfcat.flags[n]) {
    SetPRow(fitorder,psfcat.pos[n],bounds,P.row(i));
    i++;
  }
  Assert(i == ngoodpsf);
  P = invsig * P;

  // TODO: need to implement outlier rejection.
  f.reset(new tmv::Matrix<double>(U.Cols(0,npca)/P));
  xdbg<<"Done making FittedPSF\n";
}

FittedPSF::FittedPSF(const ConfigFile& _params) : params(_params)
{
  Read();
}

void FittedPSF::WriteAscii(std::string file) const
{
  std::ofstream fout(file.c_str());
  fout << psforder <<"  "<< sigma <<"  ";
  fout << fitorder <<"  "<< npca << std::endl;
  fout << bounds << std::endl;
  fout << *avepsf << std::endl;
  fout << V->Rows(0,npca) <<std::endl;
  fout << *f << std::endl;
}

void FittedPSF::ReadAscii(std::string file)
{
  std::ifstream fin(file.c_str());
  if (!fin) {
    throw std::runtime_error("Error opening fitpsf file");
  }

  xdbg<<"Reading FittedPSF:\n";
  fin >> psforder >> sigma >> fitorder >> npca >> bounds;
  xdbg<<"psforder = "<<psforder<<", sigma_psf = "<<sigma<<std::endl;
  xdbg<<"fitorder = "<<fitorder<<", npca = "<<npca<<std::endl;
  xdbg<<"bounds = "<<bounds<<std::endl;
  fitsize = (fitorder+1)*(fitorder+2)/2;
  avepsf.reset(new BVec(psforder,sigma));
  fin >> *avepsf;
  xxdbg<<"avepsf = "<<*avepsf<<std::endl;
  V.reset(new tmv::Matrix<double,tmv::RowMajor>(npca,avepsf->size()));
  fin >> *V;
  xxdbg<<"V = "<<*V<<std::endl;
  f.reset(new tmv::Matrix<double>(fitsize,npca));
  fin >> *f;
  xxdbg<<"f = "<<*f<<std::endl;
}

void FittedPSF::Write() const
{
  std::vector<std::string> files = MultiName(params, "fitpsf");  

  for(size_t i=0; i<files.size(); ++i) {
    const std::string& file = files[i];
    dbg<<"Writing fitted psf to file: "<<file<<std::endl;

    bool fitsio = false;
    if (params.keyExists("fitpsf_io")) {
      std::vector<std::string> ios = params["fitpsf_io"];
      Assert(ios.size() == files.size());
      fitsio = (ios[i] == "FITS");
    }
    else if (file.find("fits") != std::string::npos) 
      fitsio = true;

    if (fitsio) {
      //WriteFits(file);
      WriteFitsCCfits(file);
    } else {
      WriteAscii(file);
    }
  }
  dbg<<"Done Write FittedPSF\n";
}

void FittedPSF::Read()
{
  std::string file = Name(params,"fitpsf",false,true);
  // false,true = input_prefix=false, mustexist=true.
  // It is an input here, but it is in the output_prefix directory.
  dbg<< "Reading FittedPSF from file: " << file << std::endl;

  bool fitsio = false;
  if (params.keyExists("fitpsf_io")) 
    fitsio = (params["fitpsf_io"] == "FITS");
  else if (file.find("fits") != std::string::npos) 
    fitsio = true;

  if (fitsio) {
    ReadFits(file);
  } else {
    ReadAscii(file);
  }
  dbg<<"Done Read FittedPSF\n";
}

void FittedPSF::WriteFitsCCfits(std::string file) const
{
  dbg<<"Start WriteFits"<<std::endl;

  // ! means overwrite existing file
  CCfits::FITS fits("!"+file, CCfits::Write);

  dbg<<"Made fits"<<std::endl;

  // Note the actual coeffs may be less than this the way Mike does things
  // but we will fill the extra with zeros
  int n_shapelet_coeff = (psforder+1)*(psforder+2)/2;
  //int n_rot_matrix_max = n_shapelet_coeff*n_shapelet_coeff;
  int n_rot_matrix = npca*n_shapelet_coeff;

  int n_fit_coeff = (fitorder+1)*(fitorder+2)/2;
  //int n_interp_matrix_max = n_shapelet_coeff*n_fit_coeff;
  int n_interp_matrix = npca*n_fit_coeff;

  const int nfields=11;

  std::vector<string> colnames(nfields);
  std::vector<string> colfmts(nfields);
  std::vector<string> colunits(nfields);

  colnames[0] = params["fitpsf_psf_order_col"];
  colnames[1] = params["fitpsf_sigma_col"];
  colnames[2] = params["fitpsf_fit_order_col"];
  colnames[3] = params["fitpsf_npca_col"];

  colnames[4] = params["fitpsf_xmin_col"];
  colnames[5] = params["fitpsf_xmax_col"];
  colnames[6] = params["fitpsf_ymin_col"];
  colnames[7] = params["fitpsf_ymax_col"];

  colnames[8] = params["fitpsf_ave_psf_col"];
  colnames[9] = params["fitpsf_rot_matrix_col"];
  colnames[10] = params["fitpsf_interp_matrix_col"];

  colfmts[0] = "1J"; // psf order
  colfmts[1] = "1D"; // sigma
  colfmts[2] = "1J"; // fit order
  colfmts[3] = "1J"; // npca
  colfmts[4] = "1E"; // xmin
  colfmts[5] = "1E"; // xmax
  colfmts[6] = "1E"; // ymin
  colfmts[7] = "1E"; // ymax

  std::stringstream fmt;
  fmt << n_shapelet_coeff << "D";
  colfmts[8] = fmt.str(); // average psf shapelets
  fmt.str("");
  fmt << n_rot_matrix << "D";
  colfmts[9] = fmt.str(); // rotation matrix
  fmt.str("");
  fmt << n_interp_matrix << "D";
  colfmts[10] = fmt.str(); // interp matrix

  colunits[0] = "None";     // psf order
  colunits[1] = "arcsec";   // sigma
  colunits[2] = "None";     // fit order
  colunits[3] = "None";     // npca
  colunits[4] = "pixels";   // xmin
  colunits[5] = "pixels";   // xmax
  colunits[6] = "pixels";   // ymin
  colunits[7] = "pixels";   // ymax
  colunits[8] = "None";     // avg psf shapelets
  colunits[9] = "None";     // rot_matrix
  colunits[10] = "None";    // interp_matrix


  int nrows=1;
  dbg<<"Before Create table"<<std::endl;
  CCfits::Table* table;
  table = fits.addTable("fitpsf",nrows,colnames,colfmts,colunits);
  dbg<<"After Create table"<<std::endl;


  // dimensions information
  dbg<<"Adding dimensional info"<<std::endl;
  std::stringstream tdim10, tdim11;
  tdim10<<"("<<npca<<","<<n_shapelet_coeff<<")";
  tdim11<<"("<<n_fit_coeff<<","<<npca<<")";

  table->addKey("TDIM10",tdim10.str(),"dimensions of rot_matrix");
  table->addKey("TDIM11",tdim11.str(),"dimensions of interp_matrix");

  // write the data
  int startrow=1;
  int nel=1;

  // CCfits write() doesn't allow constant scalar arguments, so we have copy
  // out.  Might as well be vectors to simplify the signature
  
  std::vector<int> l_psforder(1,psforder);
  std::vector<double> l_sigma(1,sigma);
  std::vector<int> l_fitorder(1,fitorder);
  std::vector<int> l_npca(1,npca); 
  std::vector<double> xmin(1,bounds.GetXMin());
  std::vector<double> xmax(1,bounds.GetXMax());
  std::vector<double> ymin(1,bounds.GetYMin());
  std::vector<double> ymax(1,bounds.GetYMax());

  table->column(colnames[0]).write(l_psforder, startrow);
  table->column(colnames[1]).write(l_sigma, startrow);
  table->column(colnames[2]).write(l_fitorder, startrow);
  table->column(colnames[3]).write(l_npca, startrow);
  table->column(colnames[4]).write(xmin, startrow);
  table->column(colnames[5]).write(xmax, startrow);
  table->column(colnames[6]).write(ymin, startrow);
  table->column(colnames[7]).write(ymax, startrow);

  double* cptr;

  cptr = (double *) avepsf->cptr();
  table->column(colnames[8]).write(cptr, n_shapelet_coeff, nrows, startrow);
  cptr = (double *) V->cptr();
  table->column(colnames[9]).write(cptr, n_rot_matrix, nrows, startrow);
  cptr = (double *) f->cptr();
  table->column(colnames[10]).write(cptr, n_interp_matrix, nrows, startrow);


  std::string str;
  double dbl;
  int intgr;

  dbg<<"Before Write Par Keys"<<std::endl;
  WriteParKeyCCfits(params, table, "version", str);
  WriteParKeyCCfits(params, table, "noise_method", str);
  WriteParKeyCCfits(params, table, "dist_method", str);

  WriteParKeyCCfits(params, table, "fitpsf_order", intgr);
  WriteParKeyCCfits(params, table, "fitpsf_pca_thresh", dbl);
  dbg<<"After Write Par Keys"<<std::endl;

}



void FittedPSF::WriteFits(std::string file) const
{
  dbg<<"Start WriteFits"<<std::endl;

  FitsFile fits(file, READWRITE, true);
  dbg<<"Made fits"<<std::endl;

  // Note the actual coeffs may be less than this the way Mike does things
  // but we will fill the extra with zeros
  int n_shapelet_coeff = (psforder+1)*(psforder+2)/2;
  //int n_rot_matrix_max = n_shapelet_coeff*n_shapelet_coeff;
  int n_rot_matrix = npca*n_shapelet_coeff;

  int n_fit_coeff = (fitorder+1)*(fitorder+2)/2;
  //int n_interp_matrix_max = n_shapelet_coeff*n_fit_coeff;
  int n_interp_matrix = npca*n_fit_coeff;


  std::string psf_order_col=params.get("fitpsf_psf_order_col");
  std::string sigma_col=params.get("fitpsf_sigma_col");
  std::string fit_order_col=params.get("fitpsf_fit_order_col");
  std::string npca_col=params.get("fitpsf_npca_col");

  std::string xmin_col=params.get("fitpsf_xmin_col");
  std::string xmax_col=params.get("fitpsf_xmax_col");
  std::string ymin_col=params.get("fitpsf_ymin_col");
  std::string ymax_col=params.get("fitpsf_ymax_col");

  std::string ave_psf_col=params.get("fitpsf_ave_psf_col");
  std::string rot_matrix_col=params.get("fitpsf_rot_matrix_col");
  std::string interp_matrix_col=params.get("fitpsf_interp_matrix_col");

  dbg<<"Before table_cols"<<std::endl;
  const int nfields=11;
  std::string table_cols[nfields] = {
    psf_order_col,
    sigma_col,
    fit_order_col,
    npca_col,
    xmin_col, xmax_col, ymin_col, ymax_col,
    ave_psf_col,
    rot_matrix_col,
    interp_matrix_col
  };
  int table_nelem[nfields] = {
    1,1,1,1,1,1,1,1,
    n_shapelet_coeff, 
    n_rot_matrix,
    n_interp_matrix
  };
  Type_FITS table_types[nfields] = {
    XINT,
    XDOUBLE,
    XINT,
    XINT,
    XFLOAT, XFLOAT, XFLOAT, XFLOAT,
    XDOUBLE,
    XDOUBLE,
    XDOUBLE
  };
  std::string table_units[nfields] = {
    "None",
    "arcsec",
    "None",
    "None",
    "pixels", "pixels", "pixels", "pixels",
    "None",
    "None",
    "None"
  };

  dbg<<"Before Create table"<<std::endl;
  fits.CreateBinaryTable(1,nfields,
      table_cols,table_nelem,table_types,table_units);
  dbg<<"After Create table"<<std::endl;

  // dimensions information
  std::stringstream tdim10;
  tdim10<<"("<<npca<<","<<n_shapelet_coeff<<")";
  fits.WriteKey("TDIM10",XSTRING,tdim10.str().c_str(),
      "dimensions of rot_matrix");

  std::stringstream tdim11;
  tdim11<<"("<<n_fit_coeff<<","<<npca<<")";
  fits.WriteKey("TDIM11",XSTRING,tdim11.str().c_str(),
      "dimensions of interp_matrix");

  WriteParKey("version");
  WriteParKey("noise_method");
  WriteParKey("dist_method");
  WriteParKey("fitpsf_order");
  WriteParKey("fitpsf_pca_thresh");
  dbg<<"After Write keys"<<std::endl;

  float xmin = bounds.GetXMin();
  float xmax = bounds.GetXMax();
  float ymin = bounds.GetYMin();
  float ymax = bounds.GetYMax();

  fits.WriteColumn(XINT, 1, 1, 1, 1, &psforder);
  fits.WriteColumn(XDOUBLE, 2, 1, 1, 1, &sigma);
  fits.WriteColumn(XINT, 3, 1, 1, 1, &fitorder);
  fits.WriteColumn(XINT, 4, 1, 1, 1, &npca);
  fits.WriteColumn(XFLOAT, 5, 1, 1, 1, &xmin);
  fits.WriteColumn(XFLOAT, 6, 1, 1, 1, &xmax);
  fits.WriteColumn(XFLOAT, 7, 1, 1, 1, &ymin);
  fits.WriteColumn(XFLOAT, 8, 1, 1, 1, &ymax);
  fits.WriteColumn(XDOUBLE, 9, 1, 1, n_shapelet_coeff, avepsf->cptr());
  fits.WriteColumn(XDOUBLE, 10, 1, 1, n_rot_matrix, V->ptr());
  fits.WriteColumn(XDOUBLE, 11, 1, 1, n_interp_matrix, f->ptr());
}

void FittedPSF::ReadFits(std::string file)
{
  FitsFile fits(file);

  int hdu = params.read("fitpsf_hdu",2);
  dbg<<"Moving to HDU #"<<hdu<<std::endl;
  fits.GotoHDU(hdu);

  // Allocate memory for the columns we will read
  //ResizePSFCat(cat, nrows, psforder);

  std::string psf_order_col=params.get("fitpsf_psf_order_col");
  std::string sigma_col=params.get("fitpsf_sigma_col");
  std::string fit_order_col=params.get("fitpsf_fit_order_col");
  std::string npca_col=params.get("fitpsf_npca_col");

  std::string xmin_col=params.get("fitpsf_xmin_col");
  std::string xmax_col=params.get("fitpsf_xmax_col");
  std::string ymin_col=params.get("fitpsf_ymin_col");
  std::string ymax_col=params.get("fitpsf_ymax_col");

  std::string ave_psf_col=params.get("fitpsf_ave_psf_col");
  std::string rot_matrix_col=params.get("fitpsf_rot_matrix_col");
  std::string interp_matrix_col=params.get("fitpsf_interp_matrix_col");

  int nrows=1;

  fits.ReadScalarCol(psf_order_col, XINT,  &psforder, nrows);
  xdbg<<"psforder = "<<psforder<<std::endl;
  fits.ReadScalarCol(sigma_col, XDOUBLE, &sigma, nrows);
  xdbg<<"sigma = "<<sigma<<std::endl;
  fits.ReadScalarCol(fit_order_col, XINT, &fitorder, nrows);
  xdbg<<"fitorder = "<<fitorder<<std::endl;
  fits.ReadScalarCol(npca_col, XINT, &npca, nrows);
  xdbg<<"npca = "<<npca<<std::endl;

  float xmin;
  float xmax;
  float ymin;
  float ymax;

  fits.ReadScalarCol(xmin_col, XFLOAT, &xmin, nrows);
  xdbg<<"xmin = "<<xmin<<std::endl;
  fits.ReadScalarCol(xmax_col, XFLOAT, &xmax, nrows);
  xdbg<<"xmax = "<<xmax<<std::endl;
  fits.ReadScalarCol(ymin_col, XFLOAT, &ymin, nrows);
  xdbg<<"ymin = "<<ymin<<std::endl;
  fits.ReadScalarCol(ymax_col, XFLOAT, &ymax, nrows);
  xdbg<<"ymax = "<<ymax<<std::endl;

  bounds.SetXMin(xmin);
  bounds.SetXMax(xmax);
  bounds.SetYMin(ymin);
  bounds.SetYMax(ymax);
  xdbg<<"bounds = "<<bounds<<std::endl;

  avepsf.reset(new BVec(psforder,sigma));
  int n_shapelet_coeff = (psforder+1)*(psforder+2)/2;
  Assert(int(avepsf->size()) == n_shapelet_coeff);
  fits.ReadCell(ave_psf_col, XDOUBLE, avepsf->ptr(), 1, n_shapelet_coeff);
  xdbg<<"avepsf = "<<*avepsf<<std::endl;

  int n_rot_matrix = npca*n_shapelet_coeff;
  V.reset(new tmv::Matrix<double,tmv::RowMajor>(npca,avepsf->size()));
  Assert(int(V->LinearView().size()) == n_rot_matrix);
  fits.ReadCell(rot_matrix_col, XDOUBLE, V->ptr(), 1, n_rot_matrix);
  xdbg<<"V = "<<*V<<std::endl;

  fitsize = (fitorder+1)*(fitorder+2)/2;
  xdbg<<"fitsize = "<<fitsize<<std::endl;
  int n_interp_matrix = npca*fitsize;
  f.reset(new tmv::Matrix<double>(fitsize,npca));
  Assert(int(f->LinearView().size()) == n_interp_matrix);
  fits.ReadCell(interp_matrix_col, XDOUBLE, f->ptr(), 1, n_interp_matrix);
  xdbg<<"f = "<<*f<<std::endl;
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

