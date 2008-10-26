#include "Image.h"
#include "Function2D.h"
#include <fitsio.h>

using std::complex;
using std::vector;
using std::endl;
using tmv::Matrix;
using tmv::MatrixView;
using tmv::Vector;
using tmv::VectorView;
using tmv::ColMajor;
using tmv::VIt;
using tmv::CVIt;
using tmv::Unit;
using tmv::NonConj;

template <class T> inline T SQR(const T& x) { return x*x; }
template <class T> inline void SWAP(T& a, T& b) { T temp = a; a = b; b = temp; }

template <class T> Image<T>::Image(const std::string& filename, int hdu) 
{
  // TODO: figure out how cfitsio deals with different hdu.
  // I'm sure it's easy - I just haven't looked yet.
  if (hdu != 1) std::cout<<"hdu != 1 not implemented yet.\n";
  Assert(hdu == 1);
  xxdbg<<"Start read fitsimage"<<endl;
  fitsfile *fptr;
  int fitserr=0;

  fits_open_file(&fptr,filename.c_str(),READONLY,&fitserr);
  xxdbg<<"Done open"<<endl;
  Assert(fitserr==0);
  int bitpix, naxes;
  long sizes[2];
  fits_get_img_param(fptr, int(2), &bitpix, &naxes, sizes, &fitserr);
  xxdbg<<"done getimgparam"<<endl;
  Assert(fitserr==0);
  Assert(bitpix == FLOAT_IMG);
  Assert(naxes == 2);
  xxdbg<<"sizes = "<<sizes[0]<<"  "<<sizes[1]<<endl;

  xmin = 0;
  xmax = sizes[0];
  ymin = 0;
  ymax = sizes[1];
  sourcem.reset(new Matrix<T,ColMajor>(xmax,ymax));
  xxdbg<<"done make matrix of image"<<endl;
 
  long fpixel[2] = {1,1};
  int anynul;
  xxdbg<<"Before read_pix\n";
  fits_read_pix(fptr,TDOUBLE,fpixel,long(xmax*ymax),0,sourcem->ptr(),&anynul,
      &fitserr);
  xxdbg<<"done readpix  "<<fitserr<<endl;
  Assert(fitserr==0);

  itsm.reset(new MatrixView<T>(sourcem->View()));
  xxdbg<<"Done make matrixview"<<endl;

  fits_close_file(fptr, &fitserr);
  Assert(fitserr==0);
}

template <class T> void Image<T>::Flush(const std::string& filename)
{
  fitsfile *fptr;
  int fitserr=0;
  fits_open_file(&fptr,filename.c_str(),READWRITE,&fitserr);
  Assert(fitserr==0);

  long fpixel[2] = {1,1};
  fits_write_pix(fptr,TDOUBLE,fpixel,long(xmax*ymax),itsm->ptr(),&fitserr);
  Assert(fitserr==0);

  fits_close_file(fptr, &fitserr);
  Assert(fitserr==0);
}

template <class T> vector<vector<Image<T>*> > Image<T>::Divide(
    size_t nx, size_t ny) const
{
  vector<size_t> x(nx+1);
  vector<size_t> y(ny+1);
  x[0] = xmin;  x[nx] = xmax;
  y[0] = ymin;  y[ny] = ymax;
  size_t xstep = (xmax-xmin)/nx;
  size_t ystep = (ymax-ymin)/ny;
  for(size_t i=1;i<nx;i++) x[i] = x[i-1]+xstep;
  for(size_t j=1;j<ny;j++) y[j] = y[j-1]+ystep;
  vector<vector<Image*> > blockimages(nx,vector<Image*>(ny));
  for(size_t i=0;i<nx;i++) for(size_t j=0;j<ny;j++) 
    blockimages[i][j] = new Image(itsm->SubMatrix(x[i],x[i+1],y[j],y[j+1]),
	  x[i],x[i+1],y[j],y[j+1]);
  return blockimages;
}

template <class T> double Image<T>::Interpolate(double x, double y) const
{
  Assert(x>=double(xmin) && x<double(xmax));
  Assert(y>=double(ymin) && y<double(ymax));
  int i = int(floor(x-0.5));
  double dx = x - (i+0.5);
  int j = int(floor(y-0.5));
  double dy = y - (j+0.5);
  Assert(dx >= 0. && dx < 1.);
  Assert(dy >= 0. && dy < 1.);

/* 
   2    3
     x      the point (x,y) is within square of points 0,1,2,3 as shown
    
   0    1 

  Since the points are really the values at the center of the pixel, it
  is possible for x to fall closer to the edge of the chip than any points.
  In this case, the function does an extrapolation given the nearest
  square of pixels.
*/
  if (i==int(xmax-1)) {i--; dx += 1.;}
  if (j==int(ymax-1)) {j++; dy += 1.;}
  if (i==-1) {i++; dx -= 1.;}
  if (j==-1) {j++; dy -= 1.;}
  Assert(i>=0 && j>=0 && i+1<int(itsm->colsize()) && j+1<=int(itsm->rowsize()));

  double f0 = (*itsm)(i,j);
  double f1 = (*itsm)(i+1,j);
  double f2 = (*itsm)(i+1,j);
  double f3 = (*itsm)(i+1,j+1);
  double dfdx = f1-f0;
  double dfdy = f2-f0;
  double d2fdxdy = f3+f0-f1-f2;
  return f0 + dfdx*dx + dfdy*dy + d2fdxdy*dx*dy;
}
  
template <class T> double Image<T>::QuadInterpolate(double x, double y) const
{
  //static size_t count=0;
  //++count;
  Assert(x>=xmin && x< xmax);
  Assert(y>=ymin && y< ymax);
  size_t i = size_t (floor(x));
  double dx = x - (i+0.5);
  size_t j = size_t (floor(y));
  double dy = y - (j+0.5);
  Assert(i<itsm->colsize());
  Assert(j<itsm->rowsize());
  Assert (fabs(dx) <= 0.5);
  Assert (fabs(dy) <= 0.5);

/*
  7   4   8

       x          (x,y) is closer to point 0 than any other.
  1   0   2 


  5   3   6

  If any points are off the edge, we set them to the value they would
  have if the second derivative were 0 there.
*/
  double f0 = (*itsm)(i,j);
  double f1 = (i > 0) ? (*itsm)(i-1,j) : 0.;
  double f2 = (i < itsm->colsize()-1) ? (*itsm)(i+1,j) : 0.;
  double f3 = (j > 0) ? (*itsm)(i,j-1) : 0.;
  double f4 = (j < itsm->rowsize()-1) ? (*itsm)(i,j+1) : 0.;
  double f5 = (i > 0 && j > 0) ? (*itsm)(i-1,j-1) : 0.;
  double f6 = (i < itsm->colsize()-1 && j > 0) ? (*itsm)(i+1,j-1) : 0.;
  double f7 = (i > 0 && j < itsm->rowsize()-1) ? (*itsm)(i-1,j+1) : 0.;
  double f8 = (i < itsm->colsize()-1 && j < itsm->rowsize()-1) ?
    (*itsm)(i+1,j+1) : 0.;
  if (i == 0) {
    f1 = 2*f0 - f2;
    f5 = 2*f3 - f6;
    f7 = 2*f4 - f8;
  }
  if (i == itsm->colsize()-1) {
    f2 = 2*f0 - f1;
    f6 = 2*f3 - f5;
    f8 = 2*f4 - f7;
  }
  if (j == 0) {
    f3 = 2*f0 - f4;
    f5 = 2*f1 - f7;
    f6 = 2*f2 - f8;
  }
  if (j == itsm->rowsize()-1) {
    f4 = 2*f0 - f3;
    f7 = 2*f1 - f5;
    f8 = 2*f2 - f6;
  }
  double dfdx = (f2-f1)/2.;
  double dfdy = (f4-f3)/2.;
  double d2fdx2 = (f1+f2-2.*f0);
  double d2fdy2 = (f3+f4-2.*f0);
  double d2fdxdy = (f5+f8-f7-f6)/4.;
  double temp = f0 + dfdx*dx + dfdy*dy + 0.5*d2fdx2*dx*dx + 0.5*d2fdy2*dy*dy +
	d2fdxdy*dx*dy;
  return temp;
}

template <class T> inline void MultF(T f, double& x, double& y)
  // (x+iy) *= f  
  // This one is for real f
{ y *= f; x *= f; }

template <class T> inline void MultF2(T f, double& x, double& y)
  // x+iy = f*x
  // This one is for real f
{ x *= f; }

template <class T> inline void MultF(complex<T> f, double& x, double& y)
  // (x+iy) *= f
{
  double origy = y;
  y = y*real(f)+x*imag(f);
  x = x*real(f)-origy*imag(f);
}

template <class T> inline void MultF2(complex<T> f, double& x, double& y)
  // x+iy = f*x
{
  y = x*imag(f);
  x *= real(f);
}

template class Image<double>;
