
#include "FittedPSF.h"
#include "TMV_Diag.h"
#include "dbg.h"
#include "Function2D.h"
#include "Legendre2D.h"
#include <fstream>

tmv::Vector<double> DefinePXY(size_t order, double x,
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

void SetPRow(size_t fitorder, Position pos, const Bounds& bounds, 
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
    const std::vector<int32>& flagvec,
    const std::vector<Position>& pos,
    const std::vector<double>& nu,
    double _sigma, ConfigFile& params) :
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
  V.reset(new tmv::Matrix<double>(K,psfsize));
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

FittedPSF::FittedPSF(std::istream& is)
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
  V.reset(new tmv::Matrix<double>(npca,avepsf->size()));
  is >> *V;
  xxdbg<<"V = "<<*V<<std::endl;
  f.reset(new tmv::Matrix<double>(fitsize,npca));
  is >> *f;
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

