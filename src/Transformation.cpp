
#include "Transformation.h"
#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Tri.h"
#include "dbg.h"
#include  <fitsio.h>
#include "Legendre2D.h"
#include "NLSolver.h"
#include "Name.h"
#include <stdexcept>

const double PI = 4.*atan(1.);

// 
// Constructor, and basic definition routines
//

Transformation::Transformation() : 
  up(0), vp(0), dudxp(0), dudyp(0), dvdxp(0), dvdyp(0)
{}

Transformation::Transformation(const ConfigFile& params) : 
  up(0), vp(0), dudxp(0), dudyp(0), dvdxp(0), dvdyp(0)
{
  Assert(params.keyExists("dist_method"));

  if (params["dist_method"] == "SCALE") {
    Assert(params.keyExists("pixel_scale"));
    double pixel_scale = params["pixel_scale"];
    SetToScale(pixel_scale);
  } else if (params["dist_method"] == "JACOBIAN") {
    Assert(params.keyExists("dudx"));
    Assert(params.keyExists("dudy"));
    Assert(params.keyExists("dvdx"));
    Assert(params.keyExists("dvdy"));
    double dudx = params["dudx"];
    double dudy = params["dudy"];
    double dvdx = params["dvdx"];
    double dvdy = params["dvdy"];
    SetToJacobian(dudx,dudy,dvdx,dvdy);
  } else if (params["dist_method"] == "FUNC2D") {
    std::string distfile = Name(params,"dist",true);
    std::ifstream distin(distfile.c_str());
    Assert(distin);
    ReadFunc2D(distin);
    xdbg<<"Done read distortion "<<distfile<<std::endl;
  } else if (params["dist_method"] == "WCS") {
    std::string distfile = Name(params,"dist",true);
    int hdu = 1;
    if (params.keyExists("dist_hdu")) hdu = params["dist_hdu"];
    ReadWCS(distfile,hdu);
    xdbg<<"Done read WCS distortion "<<distfile<<std::endl;
  } else {
    dbg<<"Unknown transformation method\n";
    throw std::runtime_error("Unknown distortion method");
  }
}

//
// I/O
//

void Transformation::SetToScale(double pixel_scale)
{
  tmv::Matrix<double> m(2,2,0.); 
  m(1,0) = pixel_scale; m(0,1) = 0.;
  up.reset(new Polynomial2D<double>(m));
  m(1,0) = 0.; m(0,1) = pixel_scale;
  vp.reset(new Polynomial2D<double>(m));
  dudxp.reset(new Constant2D<double>(pixel_scale));
  dudyp.reset(new Constant2D<double>(0.));
  dvdxp.reset(new Constant2D<double>(0.));
  dvdyp.reset(new Constant2D<double>(pixel_scale));
}

void Transformation::SetToJacobian(
    double dudx, double dudy, double dvdx, double dvdy) 
{
  tmv::Matrix<double> m(2,2,0.); 
  m(1,0) = dudx; m(0,1) = dudy;
  up.reset(new Polynomial2D<double>(m));
  m(1,0) = dvdx; m(0,1) = dvdy;
  vp.reset(new Polynomial2D<double>(m));
  dudxp.reset(new Constant2D<double>(dudx));
  dudyp.reset(new Constant2D<double>(dudy));
  dvdxp.reset(new Constant2D<double>(dvdx));
  dvdyp.reset(new Constant2D<double>(dvdy));
}

void Transformation::ReadFunc2D(std::istream& is) 
{
  is >> up >> vp;
  dudxp = up->DFDX();
  dudyp = up->DFDY();
  dvdxp = vp->DFDX();
  dvdyp = vp->DFDY();
}

void Transformation::WriteFunc2D(std::ostream& os) const
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
  tmv::SmallMatrix<double,2,2> S2 = D * S;
  S2 *= D.Transpose();
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
  xdbg<<"MakeInverse\n";
  xdbg<<"bounds = "<<bounds<<std::endl;
  xdbg<<"order = "<<order<<std::endl;
  double x = bounds.GetXMin();
  for(int i=0; i<=ngrid; i++, x+=dx) {
    double y = bounds.GetYMin();
    for(int j=0; j<=ngrid; j++, y+=dy) {
      xdbg<<"i,j = "<<i<<','<<j<<std::endl;
      xdbg<<"x,y = "<<x<<','<<y<<std::endl;
      Position pxy(x,y);
      t2.Transform(pxy,puv);
      xdbg<<"u,v = "<<puv<<std::endl;
      v_pos.push_back(puv);
      v_x.push_back(x);
      v_y.push_back(y);
      newbounds += puv;
    }
  }
  xdbg<<"newbounds = "<<newbounds<<std::endl;
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

enum WCSType { TNX , TAN };
enum TNXFUNC { TNX_LEGENDRE=1 , TNX_CHEBYSHEV , TNX_POLYNOMIAL};
enum TNXCROSS { TNX_XNONE , TNX_XFULL , TNX_XHALF};

static void ReadWCSFits(const std::string& filename, int hdu,
    WCSType& wcstype,
    tmv::Matrix<double>& cd, tmv::Vector<double>& crpix, 
    tmv::Vector<double>& crval)
{
  Assert(crpix.size() == 2);
  Assert(crval.size() == 2);
  Assert(cd.rowsize() == 2);
  Assert(cd.colsize() == 2);

  xdbg<<"starting read fits\n";

  int status=0;
  fitsfile *fitsptr;
  char temp[80];

  if (fits_open_file(&fitsptr,filename.c_str(),READONLY,&status))
    dbg<<"fits open file: "<<status<<std::endl;
  if (status) throw std::runtime_error(
      std::string("opening fits file: ") + filename);

  if (fits_movabs_hdu(fitsptr,hdu,0,&status))
    dbg<<"fits moveabs hdu: "<<status<<std::endl;
  if (status) throw std::runtime_error("Changing hdu");
  xdbg<<"Moved to hdu "<<hdu<<std::endl;

  if (fits_read_key(fitsptr,TSTRING,(char*)"CTYPE1",temp,NULL,&status))
    dbg<<"fits read ctype1: "<<status<<std::endl;
  if (std::string(temp) == "RA---TAN") {
    wcstype = TAN;
  } else if (std::string(temp) == "RA---TNX") {
    wcstype = TNX;
  } else {
    dbg << "ctype = "<<temp<<std::endl;
    throw std::runtime_error("ctype1 is not as expected");
  }

  if (fits_read_key(fitsptr,TSTRING,(char*)"CTYPE2",temp,NULL,&status))
    dbg<<"fits read ctype2: "<<status<<std::endl;
  if ((wcstype == TAN && std::string(temp) != "DEC--TAN") ||
      (wcstype == TNX && std::string(temp) != "DEC--TNX")) {
    dbg << "ctype = "<<temp<<std::endl;
    throw std::runtime_error("ctype2 is not as expected");
  }
 
  if (fits_read_key(fitsptr,TDOUBLE,(char*)"CRVAL1",&crval[0],NULL,&status))
    dbg<<"fits read crval1: "<<status<<std::endl;
  if (fits_read_key(fitsptr,TDOUBLE,(char*)"CRVAL2",&crval[1],NULL,&status))
    dbg<<"fits read crval2: "<<status<<std::endl;

  if (fits_read_key(fitsptr,TDOUBLE,(char*)"CRPIX1",&crpix[0],NULL,&status))
    dbg<<"fits read crpixl1: "<<status<<std::endl;
  if (fits_read_key(fitsptr,TDOUBLE,(char*)"CRPIX2",&crpix[1],NULL,&status))
    dbg<<"fits read crpixl2: "<<status<<std::endl;
  
  if (fits_read_key(fitsptr,TDOUBLE,(char*)"CD1_1",&cd(0,0),NULL,&status))
    dbg<<"fits read cd1_1: "<<status<<std::endl;
  if (fits_read_key(fitsptr,TDOUBLE,(char*)"CD2_1",&cd(1,0),NULL,&status))
    dbg<<"fits read cd2_1: "<<status<<std::endl;
  if (fits_read_key(fitsptr,TDOUBLE,(char*)"CD1_2",&cd(0,1),NULL,&status))
    dbg<<"fits read cd1_2: "<<status<<std::endl;
  if (fits_read_key(fitsptr,TDOUBLE,(char*)"CD2_2",&cd(1,1),NULL,&status))
    dbg<<"fits read cd2_2: "<<status<<std::endl;

  if (status) {
    char errmsg[80];
    fits_get_errstatus(status,errmsg);
    std::cerr << "fits error status "<<status<<" = "<<errmsg<<std::endl;
    while (fits_read_errmsg(errmsg))
      std::cerr << "fits error: "<<errmsg<<std::endl;
    throw std::runtime_error("fits errors");
  }
}

static void ReCenterDistortion(tmv::Matrix<double>& a, 
    const tmv::Matrix<double>& cd, const tmv::Vector<double>& crpix)
{
  // The function here is actually a function of 
  // z' = [CD](z-z0), 
  // so we need to convert the function to that.

  //  x' = (-cd00 x0 - cd01 y0) + cd00 x + cd01 y = xc + cd00 x + cd01 y
  //  y' = (-cd10 x0 - cd11 y0) + cd10 x + cd11 y = yc + cd10 x + cd11 y

  // (x')^p (y')^q = (c00 x + c01 y + x0)^i (c10 x + c11 y + y0)^j
  // = Sum_i=0..p Sum_m=0..p-i Sum_j=0..q Sum_n=0..q-j
  //     (pCi)(p-iCm)(qCj)(q-jCn) c00^i c01^m c10^j c11^n
  //     xc^(p-i-m) yc^(q-j-n) x^(i+j) y^(m+n)

  xdbg<<"Start recenter distortion:\n";
  xdbg<<"a = "<<a<<std::endl;
  xdbg<<"cd = "<<cd<<std::endl;
  xdbg<<"crpix = "<<crpix<<std::endl;

  int xorder = a.colsize();
  int yorder = a.rowsize();
  int maxorder = std::max(xorder,yorder);

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

  tmv::LowerTriMatrix<double,tmv::NonUnitDiag> binom(maxorder+1);
  binom(0,0) = 1.;
  for(int n=1;n<=maxorder;++n) {
    binom(n,0) = 1.; binom(n,n) = 1.;
    for(int m=1;m<n;++m) {
      binom(n,m) = binom(n-1,m-1) + binom(n-1,m);
    }
  }
  xdbg<<"binom = "<<binom<<std::endl;

  for(int p=0;p<xorder;p++) for(int q=0;q<std::min(maxorder-p,yorder);q++) {
    // (x')^p (y')^q = (c00 x + c01 y + x0)^i (c10 x + c11 y + y0)^j
    // = Sum_i=0..p Sum_m=0..p-i Sum_j=0..q Sum_n=0..q-j
    //     (pCi)(p-iCm)(qCj)(q-jCn) c00^i c01^m c10^j c11^n
    //     xc^(p-i-m) yc^(q-j-n) x^(i+j) y^(m+n)
    for(int i=0;i<=p;i++) for(int m=0;m<=p-i;m++)
      for(int j=0;j<=q;j++) for(int n=0;n<=q-j;n++) {
	newar(i+j,m+n) += 
	  binom(p,i)*binom(p-i,m)*binom(q,j)*binom(q-j,n)*
	  cd00tothe[i]*cd01tothe[m]*cd10tothe[j]*cd11tothe[n]*
	  xctothe[p-i-m]*yctothe[q-j-n]*
	  a(p,q);
      }
  }
  xdbg<<"new ar = "<<newar<<std::endl;
  a = newar;
}

static void ReadTANFits(const std::string& filename, int hdu,
    std::vector<tmv::Matrix<double> >& pv)
{
  int status=0;
  fitsfile *fitsptr;

  if (fits_open_file(&fitsptr,filename.c_str(),READONLY,&status))
    dbg<<"fits open file: "<<status<<std::endl;
  if (status) throw std::runtime_error(
      std::string("opening fits file: ") + filename);

  if (fits_movabs_hdu(fitsptr,hdu,0,&status))
    dbg<<"fits moveabs hdu: "<<status<<std::endl;
  if (status) throw std::runtime_error("Changing hdu");
  xdbg<<"Moved to hdu "<<hdu<<std::endl;

  char pvstr[10] = "PV1_1";
  for(int pvnum=0; pvnum<=1; pvnum++) {
    pv[pvnum].Zero();

    //pvstr[2]='1'+pvnum;
    pvstr[2]='2'-pvnum;
    // MJ -- This is backwards from what Erin told me for the 
    //       definitions of the PV1 and PV2 parameters.
    //       The current code makes the WCS output match the alpha
    //       and delta values given in the sextractor catalog.
    //       So there are three possibilities:
    //       1) Erin just mixed up the definitions, although it seems
    //          odd that 1 would be dec and 2 would be ra when the 
    //          opposite is true for every other parameter.
    //       2) Sextractor has a bug where it calculates the wrong ra/dec.
    //          Also odd, since I would have thought other people would be
    //          testing this.
    //       3) I have a bug somewhere else that switches the two matrices.
    //          This is the most plausible explanation, but I can't find it.

    int k = 0;
    pvstr[5] = '\0';
    for(int n=0;n<4;++n) {
      for(int i=n,j=0;j<=n;--i,++j,++k) {
	if (k == 3) ++k; // I don't know why they skip 3.
	if (k < 10) { pvstr[4] = '0'+k; pvstr[5] = '\0'; }
	else { pvstr[4] = '0'+(k/10); pvstr[5] = '0'+(k%10); pvstr[6] = '\0'; }
	xdbg<<"Try reading pv key |"<<pvstr<<"|\n";
	float temp;
	if (fits_read_key(fitsptr,TFLOAT,pvstr,&temp,NULL,&status))
	  dbg<<"Problem reading key: "<<pvstr<<" for n,i,j,k = "<<n<<','<<i<<','<<j<<','<<k<<std::endl;
	if (status) break;
	pv[pvnum](i,j) = temp;
      }
      if (status) break;
    }
    if (status) break;
  }
  pv[1].TransposeSelf(); // The x,y terms are opposite for pv[1]

  if (fits_close_file(fitsptr,&status))
    dbg<<"fits close file: "<<status<<std::endl;

  if (status) {
    char errmsg[80];
    fits_get_errstatus(status,errmsg);
    std::cerr << "fits error status "<<status<<" = "<<errmsg<<std::endl;
    while (fits_read_errmsg(errmsg))
      std::cerr << "fits error: "<<errmsg<<std::endl;
    throw std::runtime_error("fits errors");
  }
}

static void ReadTNXFits(const std::string& filename, int hdu,
    std::string& lng, std::string& lat)
{
  int status=0;
  fitsfile *fitsptr;
  char temp[80];

  if (fits_open_file(&fitsptr,filename.c_str(),READONLY,&status))
    dbg<<"fits open file: "<<status<<std::endl;
  if (status) throw std::runtime_error(
      std::string("opening fits file: ") + filename);

  if (fits_movabs_hdu(fitsptr,hdu,0,&status))
    dbg<<"fits moveabs hdu: "<<status<<std::endl;
  if (status) throw std::runtime_error("Changing hdu");
  xdbg<<"Moved to hdu "<<hdu<<std::endl;

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

  lng = lnglat[0];
  lat = lnglat[1];

  if (fits_close_file(fitsptr,&status))
    dbg<<"fits close file: "<<status<<std::endl;

  if (status) {
    char errmsg[80];
    fits_get_errstatus(status,errmsg);
    std::cerr << "fits error status "<<status<<" = "<<errmsg<<std::endl;
    while (fits_read_errmsg(errmsg))
      std::cerr << "fits error: "<<errmsg<<std::endl;
    throw std::runtime_error("fits errors");
  }
}

static Function2D<double>* TNXConvert(const std::string& wcsstr,
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
  xdbg<<"functype: "<<x<<' '<<functype<<std::endl;
  wcsin >> x;
  xorder = int(x+0.5);
  xdbg<<"xorder: "<<x<<' '<<xorder<<std::endl;
  wcsin >> x;
  yorder = int(x+0.5);
  xdbg<<"yorder: "<<x<<' '<<yorder<<std::endl;
  wcsin >> x;
  crossterms = TNXCROSS(int(x+0.5));
  xdbg<<"crossterms: "<<x<<' '<<crossterms<<std::endl;

  double xmin,xmax,ymin,ymax;
  wcsin >> xmin >> xmax >> ymin >> ymax;
  Bounds b(xmin,xmax,ymin,ymax);
  xdbg<<"bounds = "<<b<<std::endl;

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
  if (!wcsin) throw std::runtime_error("reading wat line");
  
  // This next bit is only right if we have a polynomial function
  Assert(functype == TNX_POLYNOMIAL);

  ReCenterDistortion(a,cd,crpix);

  Function2D<double>* f=0;

  xdbg<<"before switch\n";
  switch(functype) {
    case TNX_CHEBYSHEV:
      xdbg<<"cheby\n";
      throw std::runtime_error("No Chebyshev functions");
      //f = new Cheby2D<double>(b,a);
      //break;
    case TNX_LEGENDRE:
      xdbg<<"legendre\n";
      f = new Legendre2D<double>(b,a);
      break;
    case TNX_POLYNOMIAL:
      xdbg<<"poly\n";
      f = new Polynomial2D<double>(a);
      break;
    default:
      throw std::runtime_error("Unknown TNX surface type");
  }
  xdbg<<"done WCSConvert\n";

  return f;
}

void Transformation::ReadWCS(std::string fitsfile, int hdu)
{
  tmv::Matrix<double> cd(2,2);
  tmv::Vector<double> crpix(2), crval(2);
  WCSType wcstype;

  ReadWCSFits(fitsfile,hdu,wcstype,cd,crpix,crval);

  if (wcstype == TAN) {
    std::vector<tmv::Matrix<double> > pv(2,tmv::Matrix<double>(4,4,0.));
    ReadTANFits(fitsfile,hdu,pv);

    if (XDEBUG) {
      // Do the direct calculation here:
      double x = 951.104;
      double y = 39.3815; 
      double ra = 334.154;
      double dec = -10.2984;

      xdbg<<"x,y = "<<x<<','<<y<<std::endl;
      tmv::Vector<double> xy(2); xy(0) = x; xy(1) = y;
      xy = cd * (xy - crpix);
      double x2 = xy(0);
      double y2 = xy(1);
      xdbg<<"x2,y2 = "<<x2<<','<<y2<<std::endl;
      tmv::Vector<double> xv(4);
      tmv::Vector<double> yv(4);
      xv(0) = 1.; xv(1) = x2; xv(2) = xv(1)*x2; xv(3) = xv(2)*x2;
      yv(0) = 1.; yv(1) = y2; yv(2) = yv(1)*y2; yv(3) = yv(2)*y2;
      double x3 = xv * pv[0] * yv;
      double y3 = xv * pv[1] * yv;
      xdbg<<"x3,y3 = "<<x3<<','<<y3<<std::endl;
      x3 *= PI/180.;
      y3 *= PI/180.;
      double cenr = crval[0] * PI/180.;
      double cend = crval[1] * PI/180.;
      double sd = cos(y3)*cos(x3)*sin(cend) + sin(y3)*cos(cend);
      double cdcr = cos(y3)*cos(x3)*cos(cend) - sin(y3)*sin(cend);
      double cdsr = cos(y3)*sin(x3);
      double d = asin(sd);
      double sr = cdsr/cos(d);
      double cr = cdcr/cos(d);
      double r = atan2(sr,cr);
      r += cenr;
      d *= 180./PI;
      r *= 180./PI;
      xdbg<<"r,d = "<<r<<','<<d<<std::endl;
      xdbg<<"Should be "<<ra<<','<<dec<<std::endl;
    }

    ReCenterDistortion(pv[0],cd,crpix);
    ReCenterDistortion(pv[1],cd,crpix);
    up.reset(new Polynomial2D<double>(pv[0]));
    vp.reset(new Polynomial2D<double>(pv[1]));
  }
  else if (wcstype == TNX) {
    std::string lngstr,latstr;
    ReadTNXFits(fitsfile,hdu,lngstr,latstr);
    xdbg<<"done read fits: \n";
    xdbg<<"lngstr = \n"<<lngstr<<std::endl;
    xdbg<<"latstr = \n"<<latstr<<std::endl;
    xdbg<<"cd = \n"<<cd<<std::endl;
    xdbg<<"crpix = \n"<<crpix<<std::endl;
    xdbg<<"crval = \n"<<crval<<std::endl;

    up.reset(TNXConvert(lngstr,cd,crpix));
    vp.reset(TNXConvert(latstr,cd,crpix));
  }
  else {
    throw std::runtime_error("Unknown wcs type");
  }
  
  // These are just the distortion terms so far. 
  // Next we have to convert from the tangent plane back to 
  // spherical coordinates.  
  // The zeroth order approximation to this step is to simply recenter
  // the result by adding a constant offset.

  // The exact formulae are:
  //
  //             sin(dec) = cos(v) cos(u) sin(dec0) + sin(v) cos(dec0)
  // cos(dec) cos(ra-ra0) = cos(v) cos(u) cos(dec0) - sin(v) sin(dec0)
  // cos(dec) sin(ra-ra0) = cos(v) sin(u)
#if 0
  double sd = cos(y)*cos(x)*sin(cend) + sin(y)*cos(cend);
  double cdcr = cos(y)*cos(x)*cos(cend) - sin(y)*sin(cend);
  double cdsr = cos(y)*sin(x);
  double d = asin(sd);
  double sr = cdsr/cos(d);
  double cr = cdcr/cos(d);
  double r = atan2(sr,cr);
  r += cenr;
#endif
  // But it is probably sufficient to do a Taylor expansion of these
  // equations around the tangent point.  
  //
  // Expanding to 2nd order, and letting dec = dec0 + d, ra = ra0 + r:
  //
  // sin(dec) = sin(dec0 + d) = sin(dec0) + d cos(dec0)
  // d cos(dec0) = (-1/2 v^2 - 1/2 u^2) sin(dec0) + v cos(dec0)
  // d = v - 1/2 (u^2+v^2) tan(dec0)
  //
  // (cos(dec0) - d sin(dec0)) r = (1-1/2 v^2) u
  // r = u/cos(dec0) - 1/2 v^2 u/cos(dec0) + d u tan(dec0)
  // r = u/cos(dec0) + u v tan(dec0)
  // 
  // We should eventually include the second order terms.  But for now,
  // We just do the first order terms which are pretty simple to apply:
  (*up) *= 1./cos(crval(1) * PI/180.);
  up->AddLinear(crval(0),0.,0.);
  vp->AddLinear(crval(1),0.,0.);

  // The u and v functions give the result in degrees. 
  // We need it to be arcsec.  So convert it here:
  (*up) *= 3600.;
  (*vp) *= 3600.;

  dudxp = up->DFDX();
  dudyp = up->DFDY();
  dvdxp = vp->DFDX();
  dvdyp = vp->DFDY();
}


