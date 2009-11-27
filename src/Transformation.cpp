
#include <valarray>
#include "TMV.h"
#include <fitsio.h>
#include "Transformation.h"
#include "dbg.h"
#include "Legendre2D.h"
#include "NLSolver.h"
#include "Name.h"
#include "Params.h"
#include "TMV.h"
#include "TMV_Small.h"

const double PI = 4.*atan(1.);

// 
// Constructor, and basic definition routines
//

Transformation::Transformation() : 
  _isRaDec(false), up(0), vp(0), dudxp(0), dudyp(0), dvdxp(0), dvdyp(0)
{}

Transformation::Transformation(const ConfigFile& params) : 
  _isRaDec(false), up(0), vp(0), dudxp(0), dudyp(0), dvdxp(0), dvdyp(0)
{
  Assert(params.keyExists("dist_method"));

  std::string dist_method = params.get("dist_method");

  if (dist_method == "SCALE") {
    double pixel_scale = params.read<double>("pixel_scale");
    SetToScale(pixel_scale);
  } else if (dist_method == "JACOBIAN") {
    double dudx = params.read<double>("dudx");
    double dudy = params.read<double>("dudy");
    double dvdx = params.read<double>("dvdx");
    double dvdy = params.read<double>("dvdy");
    SetToJacobian(dudx,dudy,dvdx,dvdy);
  } else if (dist_method == "FUNC2D") {
    std::string dist_file = makeName(params,"dist",true,false);
    std::ifstream distin(dist_file.c_str());
    Assert(distin);
    ReadFunc2D(distin);
    xdbg<<"Done read distortion "<<dist_file<<std::endl;
  } else if (dist_method == "WCS") {
    std::string dist_file = makeName(params,"dist",true,false);
    int hdu = params.read("dist_hdu",1);
    ReadWCS(dist_file,hdu);
    xdbg<<"Done read WCS distortion "<<dist_file<<std::endl;
  } else {
    dbg<<"Unknown transformation method: "<<dist_method<<std::endl;
    throw ParameterException(
        "Unknown distortion method: "+dist_method);
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
  _isRaDec = false;
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
  _isRaDec = false;
}

void Transformation::ReadFunc2D(std::istream& is) 
{
  try {
    is >> up >> vp;
  }
  catch (std::exception& e)
  {
    dbg<<"In ReadFunc2D: caught error "<<e.what()<<std::endl;
    throw ReadException(
        "Error reading Func2D format transformation.\n"
        "Caught error: " + std::string(e.what()));
  }
  dudxp = up->dFdX();
  dudyp = up->dFdY();
  dvdxp = vp->dFdX();
  dvdyp = vp->dFdY();
  _isRaDec = false;
}

void Transformation::WriteFunc2D(std::ostream& os) const
{
  try {
    os << *up << std::endl;
    os << *vp << std::endl;
  }
  catch (std::exception& e)
  {
    dbg<<"In WriteFunc2D: caught error "<<e.what()<<std::endl;
    throw WriteException(
        "Error writing Func2D format transformation.\n"
        "Caught error: " + std::string(e.what()));
  }
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
  const double ARCSEC_PER_RAD = 206264.806247;

  if (up.get()) {
    dudx = (*dudxp)(pos);
    dudy = (*dudyp)(pos);
    dvdx = (*dvdxp)(pos);
    dvdy = (*dvdyp)(pos);

    if (_isRaDec) {
      // Then u is RA, v is Dec.
      // Need to scale du (RA) by cos(Dec).
      double dec = (*vp)(pos); // this is in arcsec.
      dec /= ARCSEC_PER_RAD; // now it is in radians.
      double cosdec = cos(dec);
      dudx *= cosdec;
      dudy *= cosdec;
    }
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
  x[0] = pxy.getX();
  x[1] = pxy.getY();
  tmv::Vector<double> f(2);

  solver.method = NLSolver::Dogleg;
  solver.ftol = 1.e-8;
  // E.S.S.
  solver.gtol = 0.;
  solver.min_step = 1.e-15;
  //solver.nlout = dbgout;
  //solver.verbose = true;
  bool success = solver.Solve(x,f);
  xdbg<<"In inverseTransform, final f = "<<f<<std::endl;
  if (success) pxy = Position(x[0],x[1]);
  return success;
}


//
// MakeInverseOf
//

Bounds Transformation::MakeInverseOf(const Transformation& t2,
    const Bounds& bounds, int order)
{
  //int ngrid = 5*order;
  int ngrid = order;
  double dx = (bounds.getXMax() - bounds.getXMin()) / ngrid;
  double dy = (bounds.getYMax() - bounds.getXMin()) / ngrid;
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
  double x = bounds.getXMin();
  for(int i=0; i<=ngrid; i++, x+=dx) {
    double y = bounds.getYMin();
    for(int j=0; j<=ngrid; j++, y+=dy) {
      Position pxy(x,y);
      t2.Transform(pxy,puv);
      v_pos.push_back(puv);
      v_x.push_back(x);
      v_y.push_back(y);
      newbounds += puv;
    }
  }
  xdbg<<"newbounds = "<<newbounds<<std::endl;
  newbounds.addBorder(10.);
  xdbg<<"newbounds (add border of 10 pixels) => "<<newbounds<<std::endl;
  up.reset(new Legendre2D<double>(newbounds));
  up->simpleFit(order,v_pos,v_x,use);
  vp.reset(new Legendre2D<double>(newbounds));
  vp->simpleFit(order,v_pos,v_y,use);

  dudxp = up->dFdX();
  dudyp = up->dFdY();
  dvdxp = vp->dFdX();
  dvdyp = vp->dFdY();
  _isRaDec = false;

  return newbounds;
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

  xdbg<<"starting read wcs from fits "<<filename<<"\n";

  // TODO: Use CCFits
  int status=0;
  fitsfile *fitsptr;
  char temp[80];

  dbg<<"Reading WCS from file: "<<filename<<std::endl;
  if (fits_open_file(&fitsptr,filename.c_str(),READONLY,&status))
    dbg<<"fits open file: "<<status<<std::endl;
  if (status) throw ReadException(
      "opening fits file for transformation: "+filename);

  if (fits_movabs_hdu(fitsptr,hdu,0,&status))
    dbg<<"fits moveabs hdu: "<<status<<std::endl;
  if (status) throw ReadException(
      "Error changing hdu in transformation fits file: "+filename);
  xdbg<<"Moved to hdu "<<hdu<<std::endl;

  if (fits_read_key(fitsptr,TSTRING,(char*)"CTYPE1",temp,NULL,&status))
    dbg<<"fits read ctype1: "<<status<<std::endl;
  if (std::string(temp) == "RA---TAN") {
    wcstype = TAN;
  } else if (std::string(temp) == "RA---TNX") {
    wcstype = TNX;
  } else {
    dbg << "ctype = "<<temp<<std::endl;
    throw ReadException(
        std::string("Error ctype1 (")+temp+") is not as expected");
  }

  if (fits_read_key(fitsptr,TSTRING,(char*)"CTYPE2",temp,NULL,&status))
    dbg<<"fits read ctype2: "<<status<<std::endl;
  if ((wcstype == TAN && std::string(temp) != "DEC--TAN") ||
      (wcstype == TNX && std::string(temp) != "DEC--TNX")) {
    dbg << "ctype = "<<temp<<std::endl;
    throw ReadException(
        std::string("Error ctype2 (")+temp+") is not as expected");
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

  if (fits_close_file(fitsptr,&status))
    dbg<<"fits close file: "<<status<<std::endl;

  if (status) {
    char errmsg[80];
    fits_get_errstatus(status,errmsg);
    std::cerr << "fits error status "<<status<<" = "<<errmsg<<std::endl;
    while (fits_read_errmsg(errmsg))
      std::cerr << "fits error: "<<errmsg<<std::endl;
    throw ReadException(
        "Fits errors encountered reading WCS transformation");
  }
  xdbg<<"Done reading wcs from fits."<<std::endl;
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
  a = newar;
}

static void ReadTANFits(const std::string& filename, int hdu,
    std::vector<tmv::Matrix<double> >& pv)
{
  xdbg<<"Start read TAN specific wcs parameters."<<std::endl;
  int status=0;
  fitsfile *fitsptr;

  if (fits_open_file(&fitsptr,filename.c_str(),READONLY,&status))
    dbg<<"fits open file: "<<status<<std::endl;
  if (status)
  {
    char errmsg[80];
    fits_get_errstatus(status,errmsg);
    std::cerr << "fits error status "<<status<<" = "<<errmsg<<std::endl;
    while (fits_read_errmsg(errmsg))
      std::cerr << "fits error: "<<errmsg<<std::endl;
    throw ReadException(
        "Opening fits file "+filename+" for TAN Transformation");
  }

  if (fits_movabs_hdu(fitsptr,hdu,0,&status))
    dbg<<"fits moveabs hdu: "<<status<<std::endl;
  if (status)
  {
    char errmsg[80];
    fits_get_errstatus(status,errmsg);
    std::cerr << "fits error status "<<status<<" = "<<errmsg<<std::endl;
    while (fits_read_errmsg(errmsg))
      std::cerr << "fits error: "<<errmsg<<std::endl;
    throw ReadException(
        "Changing hdu in fits file "+filename+" for Transformation read");
  }
  xdbg<<"Moved to hdu "<<hdu<<std::endl;

  char pvstr[10] = "PV1_1";
  for(int pvnum=0; pvnum<=1; pvnum++) {
    Assert(pvnum < int(pv.size()));
    Assert(pv[pvnum].colsize() >= 4);
    Assert(pv[pvnum].rowsize() >= 4);
    pv[pvnum].Zero();

    pvstr[2]='1'+pvnum;
    pvstr[5] = '\0';
    int k = 0;
    for(int n=0;n<4;++n) {
      for(int i=n,j=0;j<=n;--i,++j,++k) {
	if (k == 3) ++k; // I don't know why they skip 3.
	if (k < 10) { pvstr[4] = '0'+k; pvstr[5] = '\0'; }
	else { pvstr[4] = '0'+(k/10); pvstr[5] = '0'+(k%10); pvstr[6] = '\0'; }
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
    throw ReadException(
        "Fits errors encountered reading WCS transformation");
  }
  xdbg<<"Done read TAN specific wcs parameters."<<std::endl;
}

static void ReadTNXFits(const std::string& filename, int hdu,
    std::string& lng, std::string& lat)
{
  int status=0;
  fitsfile *fitsptr;
  char temp[80];

  if (fits_open_file(&fitsptr,filename.c_str(),READONLY,&status))
    dbg<<"fits open file: "<<status<<std::endl;
  if (status) 
    throw ReadException(
        "Opening fits file "+filename+" for TNX Transformation");

  if (fits_movabs_hdu(fitsptr,hdu,0,&status))
    dbg<<"fits moveabs hdu: "<<status<<std::endl;
  if (status) 
      throw ReadException(
          "Changing hdu in fits file "+filename+" for Transformation read");
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
    throw ReadException(
        "Fits errors encountered reading WCS transformation");
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
  wcsin >> x;
  xorder = int(x+0.5);
  wcsin >> x;
  yorder = int(x+0.5);
  wcsin >> x;
  crossterms = TNXCROSS(int(x+0.5));

  double xmin,xmax,ymin,ymax;
  wcsin >> xmin >> xmax >> ymin >> ymax;
  Bounds b(xmin,xmax,ymin,ymax);

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
  if (!wcsin) 
    throw ReadException(
        "Error reading wat line for WCS transformation");
  
  // This next bit is only right if we have a polynomial function
  Assert(functype == TNX_POLYNOMIAL);

  ReCenterDistortion(a,cd,crpix);

  Function2D<double>* f=0;

  switch(functype) {
    case TNX_CHEBYSHEV:
      throw ReadException("Error TNX_CHEBYSHEV not implemented yet.");
    case TNX_LEGENDRE:
      f = new Legendre2D<double>(b,a);
      break;
    case TNX_POLYNOMIAL:
      f = new Polynomial2D<double>(a);
      break;
    default:
      throw ReadException( "Unknown TNX surface type");
  }

  return f;
}

void Transformation::ReadWCS(std::string fitsfile, int hdu)
{
  tmv::Matrix<double> cd(2,2);
  tmv::Vector<double> crpix(2), crval(2);
  WCSType wcstype;

  ReadWCSFits(fitsfile,hdu,wcstype,cd,crpix,crval);
  xdbg<<"wcstype = "<<wcstype<<std::endl;
  xdbg<<"cd = "<<cd<<std::endl;
  xdbg<<"crpix = "<<crpix<<std::endl;
  xdbg<<"crval = "<<crval<<std::endl;

  if (wcstype == TAN) {
    std::vector<tmv::Matrix<double> > pv(2,tmv::Matrix<double>(4,4,0.));
    ReadTANFits(fitsfile,hdu,pv);

#if 0
    if (XDEBUG) {
      // Do the direct calculation here:
      double x = 457.5703735;
      double y = 60.30065536;
      double ra = 334.1531982;
      double dec = -10.3365593;

      tmv::Vector<double> xy(2); xy(0) = x; xy(1) = y;
      xdbg.precision(10);
      xdbg<<"xy0 = "<<xy<<std::endl;
      xdbg<<"xy1 = "<<xy-crpix<<std::endl;
      xy = cd * (xy - crpix);
      xdbg<<"xy2 = "<<xy<<std::endl;
      double x2 = xy(0);
      double y2 = xy(1);
      xdbg<<"xy2 = "<<x2<<','<<y2<<std::endl;
      tmv::Vector<double> xv(4);
      tmv::Vector<double> yv(4);
      xv(0) = 1.; xv(1) = x2; xv(2) = xv(1)*x2; xv(3) = xv(2)*x2;
      yv(0) = 1.; yv(1) = y2; yv(2) = yv(1)*y2; yv(3) = yv(2)*y2;
      double x3 = xv * pv[0] * yv;
      double y3 = xv * pv[1] * yv;
      xdbg<<"xy3 = "<<x3<<','<<y3<<std::endl;
      x3 *= PI/180.;
      y3 *= PI/180.;
      double ra0 = crval[0] * PI/180.;
      double dec0 = crval[1] * PI/180.;
      double b0 = -cos(ra0)*sin(dec0)*y3 - sin(ra0)*x3 + cos(ra0)*cos(dec0);
      double b1 = -sin(ra0)*sin(dec0)*y3 + cos(ra0)*x3 + sin(ra0)*cos(dec0);
      double b2 = cos(dec0)*y3 + sin(dec0);
      double d = asin(b2/sqrt(1+x3*x3+y3*y3));
      double r = atan2(b1,b0);
      d *= 180./PI;
      r *= 180./PI;
      if (r < 0.) r += 360.;
      xdbg<<"r,d = "<<r<<','<<d<<std::endl;
      xdbg<<"Should be "<<ra<<','<<dec<<std::endl;
      xdbg<<"diff = "<<(r-ra)*3600.<<','<<(d-dec)*3600.<<" arcsec\n";
    }
#endif

    ReCenterDistortion(pv[0],cd,crpix);
    ReCenterDistortion(pv[1],cd,crpix);
    up.reset(new Polynomial2D<double>(pv[0]));
    vp.reset(new Polynomial2D<double>(pv[1]));
  }
  else if (wcstype == TNX) {
    std::string lngstr,latstr;
    ReadTNXFits(fitsfile,hdu,lngstr,latstr);
    xdbg<<"done read fits: \n";

    up.reset(TNXConvert(lngstr,cd,crpix));
    vp.reset(TNXConvert(latstr,cd,crpix));
  }
  else {
    throw ReadException( "Unknown WCS type");
  }
  
  // These are just the distortion terms so far. 
  // Next we have to convert from the tangent plane back to 
  // spherical coordinates.  
  // The zeroth order approximation to this step is to simply recenter
  // the result by adding a constant offset.

  // The exact formulae are:
  //
  // tan(ra) = (-sin(ra0) sin(dec0) v + cos(ra0) u + sin(ra0) cos(dec0)) /
  //             (-cos(ra0) sin(dec0) v - sin(ra0) u + cos(ra0) cos(dec0))
  // sin(dec) = (cos(dec0) v + sin(dec0)) / (sqrt(1+u^2+v^2))
  //
  // But it is sufficient to do a Taylor expansion of these
  // equations around the tangent point.  
  //
  // Expanding to 4th order, and letting dec = dec0 + d, ra = ra0 + r:
  //
  // r cos(dec0) ~= u + u v tan(dec0) - 1/3 u^3 sec^2(dec0) + u v^2 tan^2(dec0)
  //                - u^3 v tan(dec0) sec^2(dec0) + u v^3 tan^3(dec0)
  // d ~= v - 1/2 u^2 tan(dec0) - 1/3 v^3 - 1/2 u^2 v sec^2(dec0)
  //      - 1/2 u^2 v^2 tan^3(dec0) + 1/8 u^4 tan(dec0) (tan^2(dec0) + 3)
  //
  // (I think 3rd order is sufficient, but this takes essentially no time
  //  so might as well go slightly overkill with 4th order.)

  Polynomial2D<double> uv(6,6);
  Polynomial2D<double> u2(6,6);
  Polynomial2D<double> v2(6,6);
  Polynomial2D<double> u3(9,9);
  Polynomial2D<double> u2v(9,9);
  Polynomial2D<double> uv2(9,9);
  Polynomial2D<double> v3(9,9);
  Polynomial2D<double> u4(12,12);
  Polynomial2D<double> u3v(12,12);
  Polynomial2D<double> u2v2(12,12);
  Polynomial2D<double> uv3(12,12);

  const Polynomial2D<double>& ux = static_cast<const Polynomial2D<double>&>(*up);
  const Polynomial2D<double>& vx = static_cast<const Polynomial2D<double>&>(*vp);

  u2.makeProductOf(ux,ux);
  uv.makeProductOf(ux,vx);
  v2.makeProductOf(vx,vx);
  u3.makeProductOf(ux,u2);
  u2v.makeProductOf(ux,uv);
  uv2.makeProductOf(ux,v2);
  v3.makeProductOf(vx,v2);
  u4.makeProductOf(ux,u3);
  u3v.makeProductOf(ux,u2v);
  u2v2.makeProductOf(ux,uv2);
  uv3.makeProductOf(ux,v3);

  // Above formula assume u,v in radians, but currently u,v are in degrees.
  u2 *= PI/180.;
  uv *= PI/180.;
  u3 *= std::pow(PI/180.,2);
  u2v *= std::pow(PI/180.,2);
  uv2 *= std::pow(PI/180.,2);
  v3 *= std::pow(PI/180.,2);
  u4 *= std::pow(PI/180.,3);
  u3v *= std::pow(PI/180.,3);
  u2v2 *= std::pow(PI/180.,3);
  uv3 *= std::pow(PI/180.,3);

  // r cos(dec0) ~= u + u v tan(dec0) - 1/3 u^3 sec^2(dec0) + u v^2 tan^2(dec0)
  //                - u^3 v tan(dec0) sec^2(dec0) + u v^3 tan^3(dec0)
  double dec0 = crval[1] * PI/180.;
  double t = tan(dec0);
  double s = 1./cos(dec0);
  uv *= t;
  u3 *= -(1./3.) * s*s;
  uv2 *= t*t;
  u3v *= -t*s*s;
  uv3 *= t*t*t;

  (*up) += uv;
  (*up) += u3;
  (*up) += uv2;
  (*up) += u3v;
  (*up) += uv3;
  (*up) *= s;

  // d ~= v - 1/2 u^2 tan(dec0) - 1/3 v^3 - 1/2 u^2 v sec^2(dec0)
  //      - 1/2 u^2 v^2 tan^3(dec0) + 1/8 u^4 tan(dec0) (tan^2(dec0) + 3)
  u2 *= -0.5*t;
  v3 *= -(1./3.);
  u2v *= -0.5 * s*s;
  u2v2 *= -0.5*t*t*t;
  u4 *= (1./8.)*t*(t*t+3.);

  (*vp) += u2;
  (*vp) += v3;
  (*vp) += u2v;
  (*vp) += u2v2;
  (*vp) += u4;

  up->addLinear(crval(0),0.,0.);
  vp->addLinear(crval(1),0.,0.);

  // The u and v functions give the result in degrees. 
  // We need it to be arcsec.  So convert it here:
  (*up) *= 3600.;
  (*vp) *= 3600.;

  dudxp = up->dFdX();
  dudyp = up->dFdY();
  dvdxp = vp->dFdX();
  dvdyp = vp->dFdY();

  // The transformed variables u,v are RA and Dec.
  // This has implications for how we want to use dudx, etc. so keep
  // track of this fact in the variable _isRaDec.
  _isRaDec = true;
}
