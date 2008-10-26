
#include "Transformation.h"
#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Tri.h"
#include "dbg.h"
#include  <fitsio.h>
#include "Legendre2D.h"
#include "NLSolver.h"

// 
// Constructor, and basic definition routines
//

Transformation::Transformation() : 
  up(0), vp(0), dudxp(0), dudyp(0), dvdxp(0), dvdyp(0)
{}

//
// I/O
//

void Transformation::Read(std::istream& is) 
{
  is >> up >> vp;
  dudxp = up->DFDX();
  dudyp = up->DFDY();
  dvdxp = vp->DFDX();
  dvdyp = vp->DFDY();
}

void Transformation::Write(std::ostream& os) const
{
  os << *up << std::endl;
  os << *vp << std::endl;
}


// 
// Basic routines to apply the transformation
//

void Transformation::Transform(Position pxy, Position& puv) const
{
  if (up.get()) {
    double u1 = (*up)(pxy);
    double v1 = (*vp)(pxy);
    puv = Position(u1,v1);
  } else {
    puv = pxy;
  }
}

void Transformation::Distort(Position pos,
    double& ixx, double& ixy, double& iyy) const
{
  tmv::SmallMatrix<double,2,2> D;
  GetDistortion(pos,D);
  tmv::SmallMatrix<double,2,2> S;
  S(0,0) = ixx; S(0,1) = ixy;
  S(1,0) = ixy; S(1,1) = iyy;
  tmv::SmallMatrix<double,2,2> S2 = D * S * D.Transpose();
  ixx = S2(0,0);
  ixy = S2(0,1);
  iyy = S2(1,1);
}

void Transformation::GetDistortion(Position pos,
    double& dudx, double& dudy, double& dvdx, double& dvdy) const
{
  if (up.get()) {
    dudx = (*dudxp)(pos);
    dudy = (*dudyp)(pos);
    dvdx = (*dvdxp)(pos);
    dvdy = (*dvdyp)(pos);
  } else {
    dudx = 1.;
    dudy = 0.;
    dvdx = 0.;
    dvdy = 1.;
  }
}

void Transformation::GetDistortion(Position pos,
    tmv::SmallMatrix<double,2,2>& J) const
{ GetDistortion(pos,J(0,0),J(0,1),J(1,0),J(1,1)); }

void Transformation::GetDistortion(Position pos,
    tmv::Matrix<double>& J) const
{ GetDistortion(pos,J(0,0),J(0,1),J(1,0),J(1,1)); }


// 
// InverseTransform with accompanying helper class
//

class InverseSolver : public NLSolver
{
  public :

    InverseSolver(const Transformation& _t, const Position& _target) :
      t(_t), target(_target) {}

    void F(const tmv::Vector<double>& x, tmv::Vector<double>& f) const
    {
      Assert(x.size() == 2);
      Assert(f.size() == 2);
      Position pxy(x[0],x[1]);
      Position puv;
      t.Transform(pxy,puv);
      std::complex<double> diff = puv - target;
      f(0) = real(diff);
      f(1) = imag(diff);
    }

    void J(const tmv::Vector<double>& x, const tmv::Vector<double>& ,
	tmv::Matrix<double>& j) const
    {
      Assert(x.size() == 2);
      Assert(j.colsize() == 2);
      Assert(j.rowsize() == 2);
      Position pxy(x[0],x[1]);
      t.GetDistortion(pxy,j);
    }

  private :

    const Transformation& t;
    Position target;
};

bool Transformation::InverseTransform(Position puv, Position& pxy) const
{
  InverseSolver solver(*this, puv);
  tmv::Vector<double> x(2);
  x[0] = pxy.GetX();
  x[1] = pxy.GetY();
  tmv::Vector<double> f(2);

  solver.method = NLSolver::Dogleg;
  solver.ftol = 1.e-8;
  solver.min_step = 1.e-15;
  //solver.nlout = dbgout;
  //solver.verbose = true;
  bool success = solver.Solve(x,f);
  xdbg<<"In inverseTransform, final f = "<<f<<std::endl;
  pxy = Position(x[0],x[1]);
  return success;
}


//
// MakeInverseOf
//

void Transformation::MakeInverseOf(const Transformation& t2,
    const Bounds& bounds, int order)
{
  //int ngrid = 5*order;
  int ngrid = order;
  double dx = (bounds.GetXMax() - bounds.GetXMin()) / ngrid;
  double dy = (bounds.GetYMax() - bounds.GetXMin()) / ngrid;
  int ntot = (ngrid+1)*(ngrid+1);
  Position puv;
  std::vector<Position> v_pos;
  std::vector<double> v_x;
  std::vector<double> v_y;
  v_pos.reserve(ntot);
  v_x.reserve(ntot);
  v_y.reserve(ntot);
  std::vector<bool> use(ntot,true);
  Bounds newbounds;
  dbg<<"MakeInverse\n";
  dbg<<"bounds = "<<bounds<<std::endl;
  dbg<<"order = "<<order<<std::endl;
  double x = bounds.GetXMin();
  for(int i=0; i<=ngrid; i++, x+=dx) {
    double y = bounds.GetYMin();
    for(int j=0; j<=ngrid; j++, y+=dy) {
      dbg<<"i,j = "<<i<<','<<j<<std::endl;
      dbg<<"x,y = "<<x<<','<<y<<std::endl;
      Position pxy(x,y);
      t2.Transform(pxy,puv);
      dbg<<"u,v = "<<puv<<std::endl;
      v_pos.push_back(puv);
      v_x.push_back(x);
      v_y.push_back(y);
      newbounds += puv;
    }
  }
  dbg<<"newbounds = "<<newbounds<<std::endl;
  up.reset(new Legendre2D<double>(newbounds));
  up->SimpleFit(order,v_pos,v_x,use);
  vp.reset(new Legendre2D<double>(newbounds));
  vp->SimpleFit(order,v_pos,v_y,use);

  dudxp = up->DFDX();
  dudyp = up->DFDY();
  dvdxp = vp->DFDX();
  dvdyp = vp->DFDY();
}


//
// ReadWCS, with accompanying helper functions
//

enum TNXFUNC {TNX_LEGENDRE=1,TNX_CHEBYSHEV,TNX_POLYNOMIAL};
enum TNXCROSS {TNX_XNONE,TNX_XFULL,TNX_XHALF};

void ReadWCSFits(const std::string& filename, 
    std::string* lng, std::string* lat, 
    tmv::Matrix<double>* cd, tmv::Vector<double>* crpix, 
    tmv::Vector<double>* crval)
{
  Assert(crpix->size() == 2);
  Assert(crval->size() == 2);
  Assert(cd->rowsize() == 2);
  Assert(cd->colsize() == 2);

  dbg<<"starting read fits\n";

  int status=0;
  fitsfile *fitsptr;
  char temp[80];

  if (fits_open_file(&fitsptr,filename.c_str(),READONLY,&status))
    dbg<<"fits open file: "<<status<<std::endl;
  if (status) myerror("opening fits file: ",filename.c_str());

  if (fits_read_key(fitsptr,TSTRING,(char*)"CTYPE1",temp,NULL,&status))
    dbg<<"fits read ctype1: "<<status<<std::endl;
  if (std::string(temp) != "RA---TNX") {
    dbg << "ctype = "<<temp<<std::endl;
    myerror("ctype1 is not as expected");
  }

  if (fits_read_key(fitsptr,TSTRING,(char*)"CTYPE2",temp,NULL,&status))
    dbg<<"fits read ctype2: "<<status<<std::endl;
  if (std::string(temp) != "DEC--TNX") {
    dbg << "ctype = "<<temp<<std::endl;
    myerror("ctype2 is not as expected");
  }
 
  if (fits_read_key(fitsptr,TDOUBLE,(char*)"CRVAL1",&(*crval)[0],NULL,&status))
    dbg<<"fits read crval1: "<<status<<std::endl;
  if (fits_read_key(fitsptr,TDOUBLE,(char*)"CRVAL2",&(*crval)[1],NULL,&status))
    dbg<<"fits read crval2: "<<status<<std::endl;

  if (fits_read_key(fitsptr,TDOUBLE,(char*)"CRPIX1",&(*crpix)[0],NULL,&status))
    dbg<<"fits read crpixl1: "<<status<<std::endl;
  if (fits_read_key(fitsptr,TDOUBLE,(char*)"CRPIX2",&(*crpix)[1],NULL,&status))
    dbg<<"fits read crpixl2: "<<status<<std::endl;
  
  if (fits_read_key(fitsptr,TDOUBLE,(char*)"CD1_1",&(*cd)(0,0),NULL,&status))
    dbg<<"fits read cd1_1: "<<status<<std::endl;
  if (fits_read_key(fitsptr,TDOUBLE,(char*)"CD2_1",&(*cd)(1,0),NULL,&status))
    dbg<<"fits read cd2_1: "<<status<<std::endl;
  if (fits_read_key(fitsptr,TDOUBLE,(char*)"CD1_2",&(*cd)(0,1),NULL,&status))
    dbg<<"fits read cd1_2: "<<status<<std::endl;
  if (fits_read_key(fitsptr,TDOUBLE,(char*)"CD2_2",&(*cd)(1,1),NULL,&status))
    dbg<<"fits read cd2_2: "<<status<<std::endl;

  if (status) myerror("reading wcs info");

  char wat[10] = "WAT0_001";
  std::string lnglat[2];
  for(int i=0; i<=1; i++) {
    wat[3]='1'+i;
    sprintf(wat+5,"%03d",1);
    int watnum = 1;
    while(fits_read_key(fitsptr,TSTRING,wat,temp,NULL,&status) == 0) {
      lnglat[i] += temp;
      watnum++;
      sprintf(wat+5,"%03d",watnum);
    }
    if (watnum > 1) status = 0;
    else dbg<<"fits read lnglat: "<<i<<' '<<watnum<<
      ", status = "<<status<<std::endl;
  }

  *lng = lnglat[0];
  *lat = lnglat[1];

  if (status) {
    char errmsg[80];
    fits_get_errstatus(status,errmsg);
    std::cerr << "fits error status "<<status<<" = "<<errmsg<<std::endl;
    while (fits_read_errmsg(errmsg))
      std::cerr << "fits error: "<<errmsg<<std::endl;
    myerror("fits errors");
  }
}

Function2D<double>* WCSConvert(const std::string& wcsstr,
    const tmv::Matrix<double>& cd, const tmv::Vector<double>& crpix)
{
  std::istringstream wcsin(wcsstr);
  wcsin.ignore(80,'\"');

  TNXFUNC functype;
  int xorder,yorder;
  TNXCROSS crossterms;
  float x;

  wcsin >> x;
  functype = TNXFUNC(int(x+0.5));
  dbg<<"functype: "<<x<<' '<<functype<<std::endl;
  wcsin >> x;
  xorder = int(x+0.5);
  dbg<<"xorder: "<<x<<' '<<xorder<<std::endl;
  wcsin >> x;
  yorder = int(x+0.5);
  dbg<<"yorder: "<<x<<' '<<yorder<<std::endl;
  wcsin >> x;
  crossterms = TNXCROSS(int(x+0.5));
  dbg<<"crossterms: "<<x<<' '<<crossterms<<std::endl;

  double xmin,xmax,ymin,ymax;
  wcsin >> xmin >> xmax >> ymin >> ymax;
  Bounds b(xmin,xmax,ymin,ymax);
  dbg<<"bounds = "<<b<<std::endl;

  tmv::Matrix<double> a(xorder,yorder,0.);
  int xorder1 = xorder;
  int maxorder = std::max(xorder,yorder);

  for(int j=0;j<yorder;j++) {
    for(int i=0;i<xorder1;i++) {
      wcsin >> a[i][j];
    }
    if (crossterms == TNX_XNONE) xorder1 = 1;
    else if (crossterms == TNX_XHALF) 
      if (j+1+xorder > maxorder) xorder1--;
  }
  if (!wcsin) myerror("reading wat line");
  
  // This next bit is only right if we have a polynomial function
  Assert(functype == TNX_POLYNOMIAL);

  // The function here is actually a function of 
  // z' = [CD](z-z0), 
  // so we need to convert the function to that.

  //  x' = (-cd00 x0 - cd01 y0) + cd00 x + cd01 y = xc + cd00 x + cd01 y
  //  y' = (-cd10 x0 - cd11 y0) + cd10 x + cd11 y = yc + cd10 x + cd11 y

  // (x')^p (y')^q = (c00 x + c01 y + x0)^i (c10 x + c11 y + y0)^j
  // = Sum_i=0..p Sum_m=0..p-i Sum_j=0..q Sum_n=0..q-j
  //     (pCi)(p-iCm)(qCj)(q-jCn) c00^i c01^m c10^j c11^n
  //     xc^(p-i-m) yc^(q-j-n) x^(i+j) y^(m+n)

  tmv::Matrix<double> newar(xorder,yorder,0.);

  double xconst = -cd(0,0)*crpix[0] - cd(0,1)*crpix[1];
  double yconst = -cd(1,0)*crpix[0] - cd(1,1)*crpix[1];

  std::vector<double> cd00tothe(maxorder,1.);
  std::vector<double> cd01tothe(maxorder,1.);
  std::vector<double> cd10tothe(maxorder,1.);
  std::vector<double> cd11tothe(maxorder,1.);
  std::vector<double> xctothe(maxorder,1.);
  std::vector<double> yctothe(maxorder,1.);
  for(int i=1;i<maxorder;i++) {
    cd00tothe[i] = cd00tothe[i-1]*cd(0,0);
    cd01tothe[i] = cd01tothe[i-1]*cd(0,1);
    cd10tothe[i] = cd10tothe[i-1]*cd(1,0);
    cd11tothe[i] = cd11tothe[i-1]*cd(1,1);
    xctothe[i] = xctothe[i-1]*xconst;
    yctothe[i] = yctothe[i-1]*yconst;
  }

  tmv::LowerTriMatrix<double,tmv::UnitDiag> binom(maxorder+1);
  for(int n=1;n<=maxorder;++n) {
    binom(n,0) = 1.;
    for(int m=1;m<n;++m) {
      binom(n,m) = binom(n-1,m-1) + binom(n-1,m);
    }
  }

  for(int p=0;p<xorder;p++) for(int q=0;q<std::min(maxorder-p,yorder);q++) {
    // This is the term for (x')^i(y')^j 
    // = (cd00 x + cd01 y)^i (cd10 x + cd11 y)^j
    // = Sum_k=0..i (iCk) cd00^k x^k cd01^(i-k) y^(i-k)
    // * Sum_m=0..j (jCm) cd10^m x^m cd11^(j-m) y^(j-m)
    for(int i=0;i<=p;i++) for(int m=0;m<=p-i;m++)
      for(int j=0;j<=q;j++) for(int n=0;n<=q-j;n++) {
        newar[i+j][m+n] += 
          binom(p,i)*binom(p-i,m)*binom(q,j)*binom(q-j,n)*
          cd00tothe[i]*cd01tothe[m]*cd10tothe[j]*cd11tothe[n]*
          xctothe[p-i-m]*yctothe[q-j-n]*
          a[p][q];
    }
  }

  Function2D<double>* f=0;

  dbg<<"before switch\n";
  switch(functype) {
    case TNX_CHEBYSHEV:
      dbg<<"cheby\n";
      myerror("No Chebyshev functions");
      //f = new Cheby2D<double>(b,newar);
      break;
    case TNX_LEGENDRE:
      dbg<<"legendre\n";
      f = new Legendre2D<double>(b,newar);
      break;
    case TNX_POLYNOMIAL:
      dbg<<"poly\n";
      f = new Polynomial2D<double>(newar);
      break;
    default:
      myerror("Unknown TNX surface type");
  }
  dbg<<"done WCSConvert\n";
  return f;
}

void Transformation::ReadWCS(std::string fitsfile)
{
  std::string lngstr,latstr;
  tmv::Matrix<double> cd(2,2);
  tmv::Vector<double> crpix(2), crval(2);

  ReadWCSFits(fitsfile,&lngstr,&latstr,&cd,&crpix,&crval);
  dbg<<"done read fits: \n";
  dbg<<"lngstr = \n"<<lngstr<<std::endl;
  dbg<<"latstr = \n"<<latstr<<std::endl;
  dbg<<"cd = \n"<<cd<<std::endl;
  dbg<<"crpix = \n"<<crpix<<std::endl;
  dbg<<"crval = \n"<<crval<<std::endl;
  
  up.reset(WCSConvert(lngstr,cd,crpix));
  vp.reset(WCSConvert(latstr,cd,crpix));
  
  // These are just the distortion terms so far. 
  // There is also a rotation and vector shift.
  // z' = [CD](z-z0) + dist
  // Linear adjustments are:
  //   X:  (-cd[0,0]x0 - cd[0,1]y0) + cd[0,0]x + cd[0,1]y
  //   Y:  (-cd[1,0]x0 - cd[1,1]y0) + cd[1,0]x + cd[1,1]y

  // First though, the cd terms give units of degrees.
  // We want units of arcsec, so multiply each by 3600.

  cd(0,0) *= 3600.; cd(0,1) *= 3600.;
  cd(1,0) *= 3600.; cd(1,1) *= 3600.;

  double xconst = -cd(0,0)*crpix[0] - cd(0,1)*crpix[1];
  double yconst = -cd(1,0)*crpix[0] - cd(1,1)*crpix[1];

  up->AddLinear(xconst,cd(0,0),cd(0,1));
  vp->AddLinear(yconst,cd(1,0),cd(1,1));

  dudxp = up->DFDX();
  dudyp = up->DFDY();
  dvdxp = vp->DFDX();
  dvdyp = vp->DFDY();
}


