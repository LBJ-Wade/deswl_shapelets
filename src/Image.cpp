
#include "Image.h"
#include "Function2D.h"
#include "ConfigFile.h"
#include "Name.h"

#include <fitsio.h>

template <class T> inline T SQR(const T& x) { return x*x; }
template <class T> inline void SWAP(T& a, T& b) { T temp = a; a = b; b = temp; }

template <class T> inline int BitPix() { return Assert(false), 0; }
template <> inline int BitPix<double>() { return DOUBLE_IMG; }
template <> inline int BitPix<float>() { return FLOAT_IMG; }

template <class T> inline int DataType() { return Assert(false), 0; }
template <> inline int DataType<double>() { return TDOUBLE; }
template <> inline int DataType<float>() { return TFLOAT; }

template <class T>
Image<T>::Image(const ConfigFile& params, std::auto_ptr<Image<T> >& weight_im)
{
  filename = Name(params,"image",true);
  hdu = params.read("image_hdu",1);
  ReadFits(filename,hdu);
  xdbg<<"Opened image "<<Name(params,"image",true)<<std::endl;

  // Load weight image (if necessary)
  if (params["noise_method"] == "WEIGHT_IMAGE") {
    int weight_hdu = params.read("weight_hdu",1);
    weight_im.reset(
	new Image<T>(Name(params,"weight",true),weight_hdu));
    dbg<<"Opened weight image.\n";
  }
}

template <class T>
Image<T>::Image(const ConfigFile& params)
{
  filename = Name(params,"image",true);
  hdu = params.read("image_hdu",1);
  ReadFits(Name(params,"image",true),hdu);
  xdbg<<"Opened image "<<Name(params,"image",true)<<std::endl;
}

template <class T> 
Image<T>::Image(std::string _filename, int _hdu) :
  filename(_filename), hdu(_hdu)
{
  ReadFits(filename,hdu);
}

template <class T> 
void Image<T>::ReadFits(std::string filename, int hdu) 
{
  xxdbg<<"Start read fitsimage"<<std::endl;
  fitsfile *fptr;
  int fitserr=0;

  fits_open_file(&fptr,filename.c_str(),READONLY,&fitserr);
  xxdbg<<"Done open"<<std::endl;
  if (fitserr != 0) fits_report_error(stderr,fitserr);
  Assert(fitserr==0);

  fits_movabs_hdu(fptr,hdu,0,&fitserr);
  xdbg<<"Moved to hdu "<<hdu<<std::endl;
  if (fitserr != 0) fits_report_error(stderr,fitserr);
  Assert(fitserr==0);

  int bitpix, naxes;
  long sizes[2];
  fits_get_img_param(fptr, int(2), &bitpix, &naxes, sizes, &fitserr);
  xxdbg<<"done getimgparam"<<std::endl;
  xxdbg<<"naxes = "<<naxes<<std::endl;
  xxdbg<<"bitpix = "<<bitpix<<std::endl;
  xxdbg<<"FLOAT_IMG = "<<FLOAT_IMG<<std::endl;
  if (fitserr != 0) fits_report_error(stderr,fitserr);
  Assert(fitserr==0);
  //Assert(bitpix == FLOAT_IMG);
  Assert(naxes == 2);
  xxdbg<<"sizes = "<<sizes[0]<<"  "<<sizes[1]<<std::endl;

  xmin = 0;
  xmax = sizes[0];
  ymin = 0;
  ymax = sizes[1];
  sourcem.reset(new tmv::Matrix<T,tmv::ColMajor>(xmax,ymax));
  xxdbg<<"done make matrix of image"<<std::endl;
 
  long fpixel[2] = {1,1};
  int anynul;
  xdbg<<"Before read_pix\n";
  fits_read_pix(fptr,DataType<T>(),fpixel,long(xmax*ymax),0,
      sourcem->ptr(),&anynul,&fitserr);
  xxdbg<<"done readpix  "<<fitserr<<std::endl;
  xxdbg<<"anynul = "<<anynul<<std::endl;
  xxdbg<<"fitserr = "<<fitserr<<std::endl;
  if (fitserr != 0) fits_report_error(stderr,fitserr);
  Assert(fitserr==0);

  itsm.reset(new tmv::MatrixView<T>(sourcem->View()));
  xxdbg<<"Done make matrixview"<<std::endl;

  fits_close_file(fptr, &fitserr);
  if (fitserr != 0) fits_report_error(stderr,fitserr);
  Assert(fitserr==0);
}

template <class T> 
void Image<T>::Flush(std::string _filename, int _hdu) const
{
  filename = _filename;
  hdu = _hdu;
  Flush();
}

template <class T> 
void Image<T>::Flush() const
{
  Assert(filename != "");

  fitsfile *fptr;
  int fitserr=0;
  fits_open_file(&fptr,filename.c_str(),READWRITE,&fitserr);
  if (fitserr != 0) fits_report_error(stderr,fitserr);
  Assert(fitserr==0);

  if (hdu != 1) {
    int hdutype;
    fits_movabs_hdu(fptr,hdu,&hdutype,&fitserr);
    if (fitserr != 0) fits_report_error(stderr,fitserr);
    Assert(fitserr==0);
  }

  long fpixel[2] = {1,1};
  fits_write_pix(fptr,DataType<T>(),fpixel,long(xmax*ymax),
      itsm->ptr(),&fitserr);
  if (fitserr != 0) fits_report_error(stderr,fitserr);
  Assert(fitserr==0);

  fits_close_file(fptr, &fitserr);
  if (fitserr != 0) fits_report_error(stderr,fitserr);
  Assert(fitserr==0);
}

template <class T> 
void Image<T>::Write(std::string _filename) const
{
  filename = _filename;
  hdu = 1;

  fitsfile *fptr;
  int fitserr=0;
  fits_create_file(&fptr,("!"+filename).c_str(),&fitserr);
  if (fitserr != 0) fits_report_error(stderr,fitserr);
  Assert(fitserr==0);

  int bitpix = BitPix<T>();
  int naxes = 2;
  long sizes[2] = { itsm->colsize(), itsm->rowsize() };
  fits_create_img(fptr, bitpix, naxes, sizes, &fitserr);
  if (fitserr != 0) fits_report_error(stderr,fitserr);
  Assert(fitserr==0);

  long fpixel[2] = {1,1};
  fits_write_pix(fptr,DataType<T>(),fpixel,long(xmax*ymax),
      itsm->ptr(),&fitserr);
  if (fitserr != 0) fits_report_error(stderr,fitserr);
  Assert(fitserr==0);

  fits_close_file(fptr, &fitserr);
  if (fitserr != 0) fits_report_error(stderr,fitserr);
  Assert(fitserr==0);
}

template <class T> 
std::vector<std::vector<Image<T>*> > Image<T>::Divide(
    size_t nx, size_t ny) const
{
  std::vector<size_t> x(nx+1);
  std::vector<size_t> y(ny+1);
  x[0] = xmin;  x[nx] = xmax;
  y[0] = ymin;  y[ny] = ymax;
  size_t xstep = (xmax-xmin)/nx;
  size_t ystep = (ymax-ymin)/ny;
  for(size_t i=1;i<nx;i++) x[i] = x[i-1]+xstep;
  for(size_t j=1;j<ny;j++) y[j] = y[j-1]+ystep;
  std::vector<std::vector<Image*> > blockimages(nx,std::vector<Image*>(ny));
  for(size_t i=0;i<nx;i++) for(size_t j=0;j<ny;j++) 
    blockimages[i][j] = new Image(itsm->SubMatrix(x[i],x[i+1],y[j],y[j+1]),
	  x[i],x[i+1],y[j],y[j+1]);
  return blockimages;
}

template <class T> 
T Image<T>::Interpolate(double x, double y) const
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

  T f0 = (*itsm)(i,j);
  T f1 = (*itsm)(i+1,j);
  T f2 = (*itsm)(i+1,j);
  T f3 = (*itsm)(i+1,j+1);
  T dfdx = f1-f0;
  T dfdy = f2-f0;
  T d2fdxdy = f3+f0-f1-f2;
  return f0 + dfdx*dx + dfdy*dy + d2fdxdy*dx*dy;
}
  
template <class T> 
T Image<T>::QuadInterpolate(double x, double y) const
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
  Assert (std::abs(dx) <= 0.5);
  Assert (std::abs(dy) <= 0.5);

/*
  7   4   8

       x          (x,y) is closer to point 0 than any other.
  1   0   2 


  5   3   6

  If any points are off the edge, we set them to the value they would
  have if the second derivative were 0 there.
*/
  T f0 = (*itsm)(i,j);
  T f1 = (i > 0) ? (*itsm)(i-1,j) : 0.;
  T f2 = (i < itsm->colsize()-1) ? (*itsm)(i+1,j) : 0.;
  T f3 = (j > 0) ? (*itsm)(i,j-1) : 0.;
  T f4 = (j < itsm->rowsize()-1) ? (*itsm)(i,j+1) : 0.;
  T f5 = (i > 0 && j > 0) ? (*itsm)(i-1,j-1) : 0.;
  T f6 = (i < itsm->colsize()-1 && j > 0) ? (*itsm)(i+1,j-1) : 0.;
  T f7 = (i > 0 && j < itsm->rowsize()-1) ? (*itsm)(i-1,j+1) : 0.;
  T f8 = (i < itsm->colsize()-1 && j < itsm->rowsize()-1) ?
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
  T dfdx = (f2-f1)/2.;
  T dfdy = (f4-f3)/2.;
  T d2fdx2 = (f1+f2-2.*f0);
  T d2fdy2 = (f3+f4-2.*f0);
  T d2fdxdy = (f5+f8-f7-f6)/4.;
  T temp = f0 + dfdx*dx + dfdy*dy + 0.5*d2fdx2*dx*dx + 0.5*d2fdy2*dy*dy +
	d2fdxdy*dx*dy;
  return temp;
}

template <class T> 
T Image<T>::Median() const
{
  std::vector<T> pixels;
  pixels.reserve(itsm->colsize()*itsm->rowsize());
  for(size_t i=0;i<itsm->colsize();i++)
    for(size_t j=0;j<itsm->rowsize();j++)
      pixels.push_back((*itsm)(i,j));
  sort(pixels.begin(),pixels.end());
  return pixels[pixels.size()/2];
}

template class Image<double>;
template class Image<float>;
