
#include <valarray>
#include <fstream>
#include "TMV.h"
#include <CCfits/CCfits>

#include "FittedPSF.h"
#include "dbg.h"
#include "Function2D.h"
#include "Legendre2D.h"
#include "Name.h"
#include "WlVersion.h"
#include "TMV.h"
#include "WriteParKey.h"

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
  params(_params), psforder(params.read<int>("psf_order")),
  fitorder(params.read<int>("fitpsf_order")),
  fitsize((fitorder+1)*(fitorder+2)/2)
{
  xdbg<<"FittedPSF constructor\n";
  // Do a polynomial fit of the psf shapelet vectors

  Assert(psfcat.psf.size() > 0);

  // This used to be set to the sigma of psfcat.psf[0].
  // This doesn't work if the first object failed. 
  // So now I set it to -1, and then update it when we get to the first
  // object without an error flag (usually n = 0).
  sigma = -1.;
  for(size_t n=0;n<psfcat.size();n++) if (!psfcat.flags[n]) {
    sigma = psfcat.psf[n].GetSigma();
    break;
  }

  // Calculate the average psf vector
  avepsf.reset(new BVec(psforder,sigma));
  avepsf->Zero();
  size_t psfsize = avepsf->size();
  size_t ngoodpsf=0;
  for(size_t n=0;n<psfcat.size();n++) if (!psfcat.flags[n]) {
    xxdbg<<"n = "<<n<<", psf[n] = "<<psfcat.psf[n]<<std::endl;
    xxdbg<<"sigma = "<<sigma<<std::endl;
    xxdbg<<"psfcat.psf["<<n<<"].sigma = "<<psfcat.psf[n].GetSigma()<<std::endl;
    xxdbg<<"diff = "<<std::abs(sigma-psfcat.psf[n].GetSigma())<<std::endl;
    Assert(psfcat.psf[n].GetSigma() == sigma);
    *avepsf += psfcat.psf[n];
    ngoodpsf++;
  }
  *avepsf /= double(ngoodpsf);
  xxdbg<<"ngoodpsf = "<<ngoodpsf<<std::endl;
  xxdbg<<"average psf = "<<*avepsf<<std::endl;

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

// With this one we don't have to use the whole root= thing
FittedPSF::FittedPSF(const ConfigFile& _params, std::string file) : params(_params)
{
  Read(file);
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
    throw ReadError("Error opening fitpsf file"+file);
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
  if (!fin) {
    throw ReadError("Error reading fitpsf file"+file);
  }
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

    try
    {
      if (fitsio) {
	WriteFits(file);
      } else {
	WriteAscii(file);
      }
    }
    catch (CCfits::FitsException& e)
    {
      throw WriteError("Error writing to "+file+" -- caught error\n" +
	  e.message());
    }
    catch (std::exception& e)
    {
      throw WriteError("Error writing to "+file+" -- caught error\n" +
	  e.what());
    }
    catch (...)
    {
      throw WriteError("Error writing to "+file+" -- caught unknown error");
    }
  }
  dbg<<"Done Write FittedPSF\n";
}

void FittedPSF::Read()
{
	std::string file = Name(params,"fitpsf",false,true);
	Read(file);
}
void FittedPSF::Read(std::string file)
{
  // false,true = input_prefix=false, mustexist=true.
  // It is an input here, but it is in the output_prefix directory.
  dbg<< "Reading FittedPSF from file: " << file << std::endl;

  bool fitsio = false;
  if (params.keyExists("fitpsf_io")) 
    fitsio = (params["fitpsf_io"] == "FITS");
  else if (file.find("fits") != std::string::npos) 
    fitsio = true;

  if (!FileExists(file))
  {
    throw FileNotFound(file);
  }
  try 
  {
    if (fitsio) {
      ReadFits(file);
    } else {
      ReadAscii(file);
    }
  }
  catch (CCfits::FitsException& e)
  {
    throw ReadError("Error reading from "+file+" -- caught error\n" +
	e.message());
  }
  catch (std::exception& e)
  {
    throw ReadError("Error reading from "+file+" -- caught error\n" +
	e.what());
  }
  catch (...)
  {
    throw ReadError("Error reading from "+file+" -- caught unknown error");
  }
  dbg<<"Done Read FittedPSF\n";
}

void FittedPSF::WriteFits(std::string file) const
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


  // Header keywords

  // dimensions information
  dbg<<"Adding dimensional info"<<std::endl;
  std::stringstream tdim10, tdim11;
  tdim10<<"("<<npca<<","<<n_shapelet_coeff<<")";
  tdim11<<"("<<n_fit_coeff<<","<<npca<<")";

  table->addKey("TDIM10",tdim10.str(),"dimensions of rot_matrix");
  table->addKey("TDIM11",tdim11.str(),"dimensions of interp_matrix");


  std::string tmvvers = tmv::TMV_Version();
  std::string wlvers = WlVersion();

  table->addKey("tmvvers", tmvvers, "version of TMV code");
  table->addKey("wlvers", wlvers, "version of weak lensing code");

  std::string str;
  double dbl;
  int intgr;

  dbg<<"Before Write Par Keys"<<std::endl;

  // if serun= is sent we'll put it in the header.  This allows us to 
  // associate some more, possibly complicated, metadata with this file
  if ( params.keyExists("serun") ) {
	  CCfitsWriteParKey(params, table, "serun", str);
  }

  //CCfitsWriteParKey(params, table, "version", str);
  CCfitsWriteParKey(params, table, "noise_method", str);
  CCfitsWriteParKey(params, table, "dist_method", str);

  CCfitsWriteParKey(params, table, "fitpsf_order", intgr);
  CCfitsWriteParKey(params, table, "fitpsf_pca_thresh", dbl);
  dbg<<"After Write Par Keys"<<std::endl;


  // data

  // write the data
  int startrow=1;
  //int nel=1;

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
}

void FittedPSF::ReadFits(std::string file)
{
  // must do this way because of the const thing
  int hdu=2;
  if (params.keyExists("fitpsf_hdu")) {
    hdu = params.read<int>("fitpsf_hdu");
  }

  dbg<<"Opening FITS file at hdu "<<hdu<<std::endl;
  // true means read all as part of the construction
  CCfits::FITS fits(file, CCfits::Read, hdu-1, true);

  CCfits::ExtHDU& table=fits.extension(hdu-1);

  //long nrows=table.rows();


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

  long start=1;
  long end=1;

  std::vector<int> psforderV;
  table.column(psf_order_col).read(psforderV, start, end);
  psforder=psforderV[0];
  xdbg<<"psforder = "<<psforder<<std::endl;

  std::vector<double> sigmaV;
  table.column(sigma_col).read(sigmaV, start, end);
  sigma=sigmaV[0];
  xdbg<<"sigma = "<<sigma<<std::endl;

  std::vector<int> fitorderV;
  table.column(fit_order_col).read(fitorderV, start, end);
  fitorder=fitorderV[0];
  xdbg<<"fitorder = "<<fitorder<<std::endl;

  std::vector<int> npcaV;
  table.column(npca_col).read(npcaV, start, end);
  npca=npcaV[0];
  xdbg<<"npca = "<<npca<<std::endl;

  std::vector<float> xmin;
  std::vector<float> xmax;
  std::vector<float> ymin;
  std::vector<float> ymax;

  table.column(xmin_col).read(xmin, start, end);
  xdbg<<"xmin = "<<xmin[0]<<std::endl;
  table.column(xmax_col).read(xmax, start, end);
  xdbg<<"xmax = "<<xmax[0]<<std::endl;
  table.column(ymin_col).read(ymin, start, end);
  xdbg<<"ymin = "<<ymin[0]<<std::endl;
  table.column(ymax_col).read(ymax, start, end);
  xdbg<<"ymax = "<<ymax[0]<<std::endl;

  bounds.SetXMin(xmin[0]);
  bounds.SetXMax(xmax[0]);
  bounds.SetYMin(ymin[0]);
  bounds.SetYMax(ymax[0]);
  xdbg<<"bounds = "<<bounds<<std::endl;



  // vector columns
  double* dptr=NULL;
  std::valarray<double> dvec;

  avepsf.reset(new BVec(psforder,sigma));
  int n_shapelet_coeff = (psforder+1)*(psforder+2)/2;
  Assert(int(avepsf->size()) == n_shapelet_coeff);

  dvec.resize(0);
  table.column(ave_psf_col).read(dvec, 1);
  dptr=(double*) avepsf->cptr();
  for (size_t j=0; j<dvec.size(); j++) {
    dptr[j] = dvec[j];
  }
  xxdbg<<"avepsf = "<<*avepsf<<std::endl;



  int n_rot_matrix = npca*n_shapelet_coeff;
  V.reset(new tmv::Matrix<double,tmv::RowMajor>(npca,avepsf->size()));
  Assert(int(V->LinearView().size()) == n_rot_matrix);

  dvec.resize(0);
  table.column(rot_matrix_col).read(dvec, 1);
  dptr=(double*) V->cptr();
  for (size_t j=0; j<dvec.size(); j++) {
    dptr[j] = dvec[j];
  }
  xxdbg<<"V = "<<*V<<std::endl;



  fitsize = (fitorder+1)*(fitorder+2)/2;
  xdbg<<"fitsize = "<<fitsize<<std::endl;
  int n_interp_matrix = npca*fitsize;
  f.reset(new tmv::Matrix<double>(fitsize,npca));
  Assert(int(f->LinearView().size()) == n_interp_matrix);

  dvec.resize(0);
  table.column(interp_matrix_col).read(dvec, 1);
  dptr=(double*) f->cptr();
  for (size_t j=0; j<dvec.size(); j++) {
    dptr[j] = dvec[j];
  }
  xxdbg<<"f = "<<*f<<std::endl;
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

