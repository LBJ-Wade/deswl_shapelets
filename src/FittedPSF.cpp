
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

FittedPSF::FittedPSF(const PSFCatalog& psfcat,
    const ConfigFile& _params, std::string _prefix) :
  params(_params), prefix(_prefix),
  psforder(params["psf_order"]),
  fitorder(params[prefix + "order"]), fitsize((fitorder+1)*(fitorder+2)/2)
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
  if (params.keyExists(prefix + "npca")) {
    npca = params[prefix + "npca"];
    xdbg<<"npca = "<<npca<<" from parameter file\n";
  } else {
    double thresh = S(0);
    if (params.keyExists(prefix + "pca_thresh")) 
      thresh *= double(params[prefix + "pca_thresh"]);
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
  for(size_t n=0;n<psfcat.size();n++) if (!psfcat.flags[n]) {
    SetPRow(fitorder,psfcat.pos[n],bounds,P.row(i));
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

FittedPSF::FittedPSF(const ConfigFile& _params, std::string _prefix) :
  params(_params), prefix(_prefix)
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
  std::string file = Name(params,"fitpsf");
  dbg<<"Writing to fitpsf file: "<<file<<std::endl;

  if (file.find("fits") != std::string::npos) {
    WriteFits(file);
  } else {
    WriteAscii(file);
  }
  dbg<<"Done Write\n";
}

void FittedPSF::Read()
{
  std::string file = Name(params,"fitpsf",false,true);
  // false,true = input_prefix=false, mustexist=true.
  // It is an input here, but it is in the output_prefix directory.
  dbg<< "Reading FittedPSF from file: " << file << std::endl;

  if (file.find("fits") != std::string::npos) {
    ReadFits(file);
  } else {
    ReadAscii(file);
  }
  dbg<<"Done Read PSFCatalog\n";
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

  std::stringstream ave_psf_form;
  ave_psf_form << n_shapelet_coeff << "d";

  std::stringstream rot_matrix_form;
  //rot_matrix_form << n_rot_matrix_max << "d";
  rot_matrix_form << n_rot_matrix << "d";

  std::stringstream interp_matrix_form;
  //interp_matrix_form << n_interp_matrix_max << "d";
  interp_matrix_form << n_interp_matrix << "d";

  std::string psf_order_col=params.get(prefix + "psf_order_col");
  std::string sigma_col=params.get(prefix + "sigma_col");
  std::string fit_order_col=params.get(prefix + "fit_order_col");
  std::string npca_col=params.get(prefix + "npca_col");

  std::string xmin_col=params.get(prefix + "xmin_col");
  std::string xmax_col=params.get(prefix + "xmax_col");
  std::string ymin_col=params.get(prefix + "ymin_col");
  std::string ymax_col=params.get(prefix + "ymax_col");

  std::string ave_psf_col=params.get(prefix + "ave_psf_col");
  std::string rot_matrix_col=params.get(prefix + "rot_matrix_col");
  std::string interp_matrix_col=params.get(prefix + "interp_matrix_col");

  dbg<<"Before table_cols"<<std::endl;
  const int nfields=11;
  std::string table_cols[nfields] = {
    psf_order_col,
    sigma_col,
    fit_order_col,
    npca_col,
    xmin_col,
    xmax_col,
    ymin_col,
    ymax_col,
    ave_psf_col,
    rot_matrix_col,
    interp_matrix_col
  };
  std::string table_types[nfields] = {
    "1j",   // psforder
    "1d",  // sigma
    "1j",  // fitorder
    "1j",  // npca
    "1e",  // xmin  
    "1e",  // xmax
    "1e",  // ymin
    "1e",  // ymax
    ave_psf_form.str(),  // ave_psf
    rot_matrix_form.str(),  // rot_matrix
    interp_matrix_form.str()
  };  // interp_matrix
  std::string table_units[nfields] = {
    "None",    // psforder
    "arcsec", // sigma
    "None",   // fitorder
    "None",   // npca
    "pixels",   // xmin
    "pixels",   // xmax
    "pixels",   // ymin
    "pixels",   // ymax
    "None",   // ave_psf
    "None", // rot_matrix
    "None"
  };  // interp_matrix

  dbg<<"Before Create table"<<std::endl;
  fits.CreateBinaryTable(1,nfields,table_cols,table_types,table_units);
  dbg<<"After Create table"<<std::endl;

  // dimensions information
  std::stringstream tdim10;
  tdim10<<"("<<n_shapelet_coeff<<","<<n_shapelet_coeff<<")";
  fits.WriteKey("TDIM10",TSTRING,tdim10.str().c_str(),
      "dimensions of rot_matrix");

  std::stringstream tdim11;
  tdim11<<"("<<n_fit_coeff<<","<<n_shapelet_coeff<<")";
  fits.WriteKey("TDIM11",TSTRING,tdim11.str().c_str(),
      "dimensions of interp_matrix");

  fits.WriteParKey(params, "version", TSTRING);
  fits.WriteParKey(params, "noise_method", TSTRING);
  fits.WriteParKey(params, "dist_method", TSTRING);
  fits.WriteParKey(params, prefix + "order", TLONG);
  fits.WriteParKey(params, prefix + "pca_thresh", TDOUBLE);
  dbg<<"After Write keys"<<std::endl;

  float xmin = bounds.GetXMin();
  float xmax = bounds.GetXMax();
  float ymin = bounds.GetYMin();
  float ymax = bounds.GetYMax();

  fits.WriteColumn(TINT, 1, 1, 1, 1, &psforder);

  fits.WriteColumn(TDOUBLE, 2, 1, 1, 1, &sigma);

  fits.WriteColumn(TINT, 3, 1, 1, 1, &fitorder);

  fits.WriteColumn(TINT, 4, 1, 1, 1, &npca);

  fits.WriteColumn(TFLOAT, 5, 1, 1, 1, &xmin);

  fits.WriteColumn(TFLOAT, 6, 1, 1, 1, &xmax);

  fits.WriteColumn(TFLOAT, 7, 1, 1, 1, &ymin);

  fits.WriteColumn(TFLOAT, 8, 1, 1, 1, &ymax);

  fits.WriteColumn(TDOUBLE, 9, 1, 1, n_shapelet_coeff, avepsf->cptr());

  fits.WriteColumn(TDOUBLE, 10, 1, 1, n_rot_matrix, V->ptr());

  fits.WriteColumn(TDOUBLE, 11, 1, 1, n_interp_matrix, f->ptr());

  // Test FITS I/O:
  if (0) {
    FittedPSF test(params,"fitpsf_");
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
    dbg<<"     Norm(diff) = "<<Norm((*avepsf)-(*test.avepsf))<<std::endl;
    if (Norm((*avepsf)-(*test.avepsf)) > 0.001) {
      dbg<<"orig avepsf = "<<(*avepsf)<<std::endl;
      dbg<<"     new = "<<(*test.avepsf)<<std::endl;
      dbg<<"     diff = "<<(*avepsf)-(*test.avepsf)<<std::endl;
      dbg<<"     Norm(diff) = "<<Norm((*avepsf)-(*test.avepsf))<<std::endl;
      throw std::runtime_error("Error in FittedPSF Fits I/O");
    }

    dbg<<"     Norm(diff) =  "<<Norm(V->Rows(0,npca)-(*test.V))<<std::endl;
    if (Norm(V->Rows(0,npca)-(*test.V)) > 0.001) {
      dbg<<"orig rot_matrix = "<<V->Rows(0,npca)<<std::endl;
      dbg<<"     new =  "<<(*test.V)<<std::endl;
      dbg<<"     diff =  "<<V->Rows(0,npca)-(*test.V)<<std::endl;
      dbg<<"     Norm(diff) =  "<<Norm(V->Rows(0,npca)-(*test.V))<<std::endl;
      throw std::runtime_error("Error in FittedPSF Fits I/O");
    }

    dbg<<"     Norm(diff): "<<Norm((*f)-(*test.f))<<std::endl;
    if (Norm((*f)-(*test.f)) > 0.001) {
      dbg<<"orig interp_matrix = "<<(*f)<<std::endl;
      dbg<<"     new: "<<(*test.f)<<std::endl;
      dbg<<"     diff: "<<(*f)-(*test.f)<<std::endl;
      dbg<<"     Norm(diff): "<<Norm((*f)-(*test.f))<<std::endl;
      throw std::runtime_error("Error in FittedPSF Fits I/O");
    }

    dbg<<"Done testing fitted PSF file IO"<<std::endl;
  }
}

void FittedPSF::ReadFits(std::string file)
{
  FitsFile fits(file);

  int hdu = 2;
  if (params.keyExists(prefix + "hdu")) hdu = params[prefix + "hdu"];
  dbg<<"Moving to HDU #"<<hdu<<std::endl;
  fits.GotoHDU(hdu);

  // Allocate memory for the columns we will read
  //ResizePSFCat(cat, nrows, psforder);

  std::string psf_order_col=params.get(prefix + "psf_order_col");
  std::string sigma_col=params.get(prefix + "sigma_col");
  std::string fit_order_col=params.get(prefix + "fit_order_col");
  std::string npca_col=params.get(prefix + "npca_col");

  std::string xmin_col=params.get(prefix + "xmin_col");
  std::string xmax_col=params.get(prefix + "xmax_col");
  std::string ymin_col=params.get(prefix + "ymin_col");
  std::string ymax_col=params.get(prefix + "ymax_col");

  std::string ave_psf_col=params.get(prefix + "ave_psf_col");
  std::string rot_matrix_col=params.get(prefix + "rot_matrix_col");
  std::string interp_matrix_col=params.get(prefix + "interp_matrix_col");

  int nrows=1;

  fits.ReadScalarCol(psf_order_col, TINT,  &psforder, nrows);
  xdbg<<"psforder = "<<psforder<<std::endl;
  fits.ReadScalarCol(sigma_col, TDOUBLE, &sigma, nrows);
  xdbg<<"sigma = "<<sigma<<std::endl;
  fits.ReadScalarCol(fit_order_col, TINT, &fitorder, nrows);
  xdbg<<"fitorder = "<<fitorder<<std::endl;
  fits.ReadScalarCol(npca_col, TINT, &npca, nrows);
  xdbg<<"npca = "<<npca<<std::endl;

  float xmin;
  float xmax;
  float ymin;
  float ymax;

  fits.ReadScalarCol(xmin_col, TFLOAT, &xmin, nrows);
  xdbg<<"xmin = "<<xmin<<std::endl;
  fits.ReadScalarCol(xmax_col, TFLOAT, &xmax, nrows);
  xdbg<<"xmax = "<<xmax<<std::endl;
  fits.ReadScalarCol(ymin_col, TFLOAT, &ymin, nrows);
  xdbg<<"ymin = "<<ymin<<std::endl;
  fits.ReadScalarCol(ymax_col, TFLOAT, &ymax, nrows);
  xdbg<<"ymax = "<<ymax<<std::endl;

  bounds.SetXMin(xmin);
  bounds.SetXMax(xmax);
  bounds.SetYMin(ymin);
  bounds.SetYMax(ymax);
  xdbg<<"bounds = "<<bounds<<std::endl;

  avepsf.reset(new BVec(psforder,sigma));
  int n_shapelet_coeff = (psforder+1)*(psforder+2)/2;
  Assert(int(avepsf->size()) == n_shapelet_coeff);
  fits.ReadCell(ave_psf_col, TDOUBLE, avepsf->ptr(), 1, n_shapelet_coeff);
  xdbg<<"avepsf = "<<*avepsf<<std::endl;

  int n_rot_matrix = npca*n_shapelet_coeff;
  V.reset(new tmv::Matrix<double,tmv::RowMajor>(npca,avepsf->size()));
  Assert(int(V->LinearView().size()) == n_rot_matrix);
  fits.ReadCell(rot_matrix_col, TDOUBLE, V->ptr(), 1, n_rot_matrix);
  xdbg<<"V = "<<*V<<std::endl;

  fitsize = (fitorder+1)*(fitorder+2)/2;
  xdbg<<"fitsize = "<<fitsize<<std::endl;
  int n_interp_matrix = npca*fitsize;
  f.reset(new tmv::Matrix<double>(fitsize,npca));
  Assert(int(f->LinearView().size()) == n_interp_matrix);
  fits.ReadCell(interp_matrix_col, TDOUBLE, f->ptr(), 1, n_interp_matrix);
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

