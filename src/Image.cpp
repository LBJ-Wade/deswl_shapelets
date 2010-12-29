
#include <sstream>
#include <valarray>
#include <CCfits/CCfits>

#include "Image.h"
#include "Function2D.h"
#include "ConfigFile.h"
#include "Name.h"

#include <fitsio.h>

template <typename T>
Image<T>::~Image()  
{}

template <typename T>
Image<T>::Image(
    const ConfigFile& params, std::auto_ptr<Image<T> >& weightIm)
{
    _fileName = makeName(params,"image",true,true);
    _hdu = getHdu(params,"image",_fileName,1);
    readFits();
    xdbg<<"Opened image "<<_fileName<<" at hdu "<<_hdu<<std::endl;

    // Load weight image (if necessary)
    if (params["noise_method"] == "WEIGHT_IMAGE") {
        std::string weightName = makeName(params,"weight",true,true);
        int weightHdu = getHdu(params,"weight",weightName,1);
        try {
            weightIm.reset(new Image<T>(weightName,weightHdu));
        } catch (ReadException& e) {
            xdbg<<"Caught ReadException: \n"<<e.what()<<std::endl;
            throw ReadException(
                "Error reading weight image for " + _fileName + "\n" +
                e.what());
        }
        dbg<<"Opened weight image.\n";

        // Make sure any bad pixels are marked with 0 variance.
        if (params.keyExists("badpix_file") || params.keyExists("badpix_ext")) {
            std::string badpixName = makeName(params,"badpix",true,true);
            int badpixHdu = getHdu(params,"badpix",badpixName,1);
            dbg<<"badpix name = "<<badpixName<<std::endl;
            dbg<<"hdu = "<<badpixHdu<<std::endl;
            try {
                weightIm.reset(new Image<T>(weightName,weightHdu));
                Image<double> badpixIm(badpixName,badpixHdu);
                dbg<<"Opened badpix image.\n";

                for(int i=0;i<=weightIm->getMaxI();++i) {
                    for(int j=0;j<=weightIm->getMaxJ();++j) {
                        if (badpixIm(i,j) > 0.0) (*weightIm)(i,j) = 0.0;
                    }
                }
            } catch (ReadException& e) {
                xdbg<<"Caught ReadException: \n"<<e.what()<<std::endl;
                throw ReadException(
                    "Error reading badpix image for " + _fileName + "\n" +
                    e.what());
            }
        }
    }
    _loaded=true;
}

template <typename T>
Image<T>::Image(const ConfigFile& params)
{
    _fileName = makeName(params,"image",true,true);
    _hdu = getHdu(params,"image",_fileName,1);
    readFits();
    xdbg<<"Opened image "<<_fileName<<" at hdu "<<_hdu<<std::endl;
    _loaded=true;
}

template <typename T> 
Image<T>::Image(std::string fileName, int hdu) :
    _fileName(fileName), _hdu(hdu)
{
    readFits();
    _loaded=true;
}

template <typename T> 
void Image<T>::load(std::string fileName, int hdu) 
{
    _fileName = fileName;
    _hdu = hdu;
    readFits();
    _loaded=true;
}

template <typename T> 
inline int getBitPix() { return 0; }

template <> 
inline int getBitPix<double>() { return DOUBLE_IMG; }

template <> 
inline int getBitPix<float>() { return FLOAT_IMG; }

template <typename T> 
inline int getDataType() { return 0; }

template <> 
inline int getDataType<double>() { return TDOUBLE; }

template <> 
inline int getDataType<float>() { return TFLOAT; }

template <typename T> 
void Image<T>::readFits()
{
    std::stringstream err_msg;

    xdbg<<"Start Image::readFits"<<std::endl;
    xdbg<<"filename = "<<_fileName<<std::endl;
    xdbg<<"hdu = "<<_hdu<<std::endl;

    if (!doesFileExist(_fileName)) {
        throw FileNotFoundException(_fileName);
    }

#if 1
    // New CCFits implementation
    try {
        CCfits::FITS fits(_fileName, CCfits::Read);
        xdbg<<"Made fits object\n";
        std::valarray<double> data;
        if (_hdu == 1) {
            CCfits::PHDU& image=fits.pHDU();
            xdbg<<"Got primary hdu object"<<std::endl;
            image.readAllKeys();
            xdbg<<"read all keys:"<<std::endl;
            xdbg<<image<<std::endl;
            xdbg<<"axes = "<<image.axes()<<std::endl;
            if (image.axes() != 2)
                throw std::runtime_error("Number of axes != 2");
            image.read(data);
            xdbg<<"read data into valarray\n";
            _xMax = image.axis(0);
            _yMax = image.axis(1);
        } else {
            fits.read(_hdu-1);
            xdbg<<"read extension"<<std::endl;
            CCfits::ExtHDU& image=fits.extension(_hdu-1);
            xdbg<<"Got extension hdu object"<<std::endl;
            image.readAllKeys();
            xdbg<<"read all keys:"<<std::endl;
            xdbg<<image<<std::endl;
            xdbg<<"axes = "<<image.axes()<<std::endl;
            if (image.axes() != 2)
                throw std::runtime_error("Number of axes != 2");
            image.read(data);
            xdbg<<"read data into valarray\n";
            _xMax = image.axis(0);
            _yMax = image.axis(1);
        }
        _xMin = 0;
        _yMin = 0;
        xdbg<<"size = "<<_xMax<<" , "<<_yMax<<std::endl;
        _source.reset(new TMatrix(T)(_xMax,_yMax));
        xdbg<<"done make matrix of image"<<std::endl;
        xdbg<<"data.size = "<<data.size()<<" =? "<<_xMax*_yMax<<std::endl;
        Assert(int(data.size()) == _xMax*_yMax);
        std::copy(&data[0],&data[data.size()-1],TMV_ptr(*_source));
        xdbg<<"done copy data to _source\n";
        _m.reset(new TMatrixView(T)(TMV_view(*_source)));
        xdbg<<"Done make matrixview"<<std::endl;
    } catch (CCfits::FitsException& e) {
        xdbg<<"Caught FitsException: "<<e.message()<<std::endl;
        throw ReadException(
            "Error reading from " + _fileName + 
            " hdu " + ConvertibleString(_hdu) +
            " -- caught error\n" + e.message());
    } catch (std::exception& e) {
        xdbg<<"Caught std::exception: "<<e.what()<<std::endl;
        throw ReadException(
            "Error reading from " + _fileName + 
            " hdu " + ConvertibleString(_hdu) +
            " -- caught error\n" + e.what());
    } catch (...) {
        xdbg<<"Caught other exception: "<<std::endl;
        throw ReadException(
            "Error reading from " + _fileName + 
            " hdu " + ConvertibleString(_hdu) +
            " -- caught unknown error\n");
    }
#else
    // Old cfitsio implementation.
    // Remove once above is confirmed to be working correctly.
    fitsfile *fPtr;
    int fitsErr=0;

    fits_open_file(&fPtr,_fileName.c_str(),READONLY,&fitsErr);
    xdbg<<"Done open"<<std::endl;
    if (fitsErr != 0) fits_report_error(stderr,fitsErr);
    Assert(fitsErr==0);

    fits_movabs_hdu(fPtr,_hdu,0,&fitsErr);
    xdbg<<"Moved to hdu "<<_hdu<<std::endl;
    if (fitsErr != 0) fits_report_error(stderr,fitsErr);
    Assert(fitsErr==0);

    int bitPix, nAxes;
    long sizes[2];
    fits_get_img_param(fPtr, int(2), &bitPix, &nAxes, sizes, &fitsErr);
    xdbg<<"done getimgparam"<<std::endl;
    xdbg<<"naxes = "<<nAxes<<std::endl;
    xdbg<<"bitpix = "<<bitPix<<std::endl;
    xdbg<<"FLOAT_IMG = "<<FLOAT_IMG<<std::endl;
    if (fitsErr != 0) fits_report_error(stderr,fitsErr);
    Assert(fitsErr==0);
    Assert(nAxes == 2);
    xdbg<<"sizes = "<<sizes[0]<<"  "<<sizes[1]<<std::endl;

    _xMin = 0;
    _xMax = sizes[0];
    _yMin = 0;
    _yMax = sizes[1];
    _source.reset(new TMatrix(T)(_xMax,_yMax));
    xdbg<<"done make matrix of image"<<std::endl;

    long fPixel[2] = {1,1};
    int anynul;
    xdbg<<"Before read_pix\n";
    Assert(getDataType<T>());
    fits_read_pix(fPtr,getDataType<T>(),fPixel,long(_xMax*_yMax),0,
                  TMV_ptr(*_source),&anynul,&fitsErr);
    xdbg<<"done readpix  "<<fitsErr<<std::endl;
    xdbg<<"anynul = "<<anynul<<std::endl;
    xdbg<<"fitserr = "<<fitsErr<<std::endl;
    if (fitsErr != 0) fits_report_error(stderr,fitsErr);
    Assert(fitsErr==0);

    _m.reset(new TMatrixView(T)(TMV_view(*_source)));
    xdbg<<"Done make matrixview"<<std::endl;

    fits_close_file(fPtr, &fitsErr);
    if (fitsErr != 0) fits_report_error(stderr,fitsErr);
    Assert(fitsErr==0);
#endif

    xdbg<<"Done Image ReadFits"<<std::endl;
}

// Load Subimage
// These next three are basically the same as the above constructors,
// and readFits function, but with a range of pixels to read in.
// I could probably combine these pretty easily, since most of the
// code is identical, but for now there is some significant code 
// repetition here.
template <typename T>
Image<T>::Image(
    const ConfigFile& params, std::auto_ptr<Image<T> >& weightIm,
    int x1, int x2, int y1, int y2)
{
    _fileName = makeName(params,"image",true,true);
    _hdu = getHdu(params,"image",_fileName,1);
    readFits(x1,x2,y1,y2);
    xdbg<<"Opened image "<<_fileName<<std::endl;

    // Load weight image (if necessary)
    if (params["noise_method"] == "WEIGHT_IMAGE") {
        std::string weightName = makeName(params,"weight",true,true);
        int weightHdu = getHdu(params,"weight",weightName,1);
        weightIm.reset(new Image<T>(weightName,weightHdu,x1,x2,y1,y2));
        dbg<<"Opened weight image.\n";

        // Make sure any bad pixels are marked with 0 variance.
        if (params.keyExists("badpix_file") || params.keyExists("badpix_ext")) {
            std::string badpixName = makeName(params,"badpix",true,true);
            int badpixHdu = getHdu(params,"badpix",badpixName,1);
            dbg<<"badpix name = "<<badpixName<<std::endl;
            dbg<<"hdu = "<<badpixHdu<<std::endl;
            Image<double> badpixIm(badpixName,badpixHdu,x1,x2,y1,y2);
            dbg<<"Opened badpix image.\n";

            for(int i=0;i<=weightIm->getMaxI();++i) {
                for(int j=0;j<=weightIm->getMaxJ();++j) {
                    if (badpixIm(i,j) > 0.0) (*weightIm)(i,j) = 0.0;
                }
            }
        }
    }
}

template <typename T>
Image<T>::Image(const ConfigFile& params, std::auto_ptr<Image<T> >& weightIm,
                const Bounds& bounds)
{
    int x1 = int(floor(bounds.getXMin()));
    int x2 = int(ceil(bounds.getXMax()));
    int y1 = int(floor(bounds.getYMin()));
    int y2 = int(ceil(bounds.getYMax()));
    _fileName = makeName(params,"image",true,true);
    _hdu = getHdu(params,"image",_fileName,1);
    readFits(x1,x2,y1,y2);
    xdbg<<"Opened image "<<_fileName<<std::endl;

    // Load weight image (if necessary)
    if (params["noise_method"] == "WEIGHT_IMAGE") {
        std::string weightName = makeName(params,"weight",true,true);
        int weightHdu = getHdu(params,"weight",weightName,1);
        weightIm.reset(new Image<T>(weightName,weightHdu,bounds));
        dbg<<"Opened weight image.\n";

        // Make sure any bad pixels are marked with 0 variance.
        if (params.keyExists("badpix_file") || params.keyExists("badpix_ext")) {
            std::string badpixName = makeName(params,"badpix",true,true);
            int badpixHdu = getHdu(params,"badpix",badpixName,1);
            dbg<<"badpix name = "<<badpixName<<std::endl;
            dbg<<"hdu = "<<badpixHdu<<std::endl;
            Image<double> badpixIm(badpixName,badpixHdu,bounds);
            dbg<<"Opened badpix image.\n";

            for(int i=0;i<=weightIm->getMaxI();++i) {
                for(int j=0;j<=weightIm->getMaxJ();++j) {
                    if (badpixIm(i,j) > 0.0) (*weightIm)(i,j) = 0.0;
                }
            }
        }
    }
}

template <typename T>
Image<T>::Image(const ConfigFile& params, int x1, int x2, int y1, int y2)
{
    _fileName = makeName(params,"image",true,true);
    _hdu = getHdu(params,"image",_fileName,1);
    readFits(x1,x2,y1,y2);
    xdbg<<"Opened image "<<_fileName<<std::endl;
}

template <typename T>
Image<T>::Image(const ConfigFile& params, const Bounds& bounds)
{
    int x1 = int(floor(bounds.getXMin()));
    int x2 = int(ceil(bounds.getXMax()));
    int y1 = int(floor(bounds.getYMin()));
    int y2 = int(ceil(bounds.getYMax()));
    _fileName = makeName(params,"image",true,true);
    _hdu = getHdu(params,"image",_fileName,1);
    readFits(x1,x2,y1,y2);
    xdbg<<"Opened image "<<_fileName<<std::endl;
}

template <typename T> 
Image<T>::Image(
    std::string fileName, int hdu,
    int x1, int x2, int y1, int y2) :
    _fileName(fileName), _hdu(hdu)
{
    readFits(x1,x2,y1,y2);
}

template <typename T> 
Image<T>::Image(std::string fileName, int hdu, const Bounds& bounds) :
    _fileName(fileName), _hdu(hdu)
{
    int x1 = int(floor(bounds.getXMin()));
    int x2 = int(ceil(bounds.getXMax()));
    int y1 = int(floor(bounds.getYMin()));
    int y2 = int(ceil(bounds.getYMax()));
    readFits(x1,x2,y1,y2);
}

template <typename T> 
void Image<T>::readFits(int x1, int x2, int y1, int y2) 
{
    xdbg<<"Start read fitsimage"<<std::endl;
    xdbg<<"filename = "<<_fileName<<std::endl;
    xdbg<<"hdu = "<<_hdu<<std::endl;
    // TODO: Use CCFits
    fitsfile *fPtr;
    int fitsErr=0;

    fits_open_file(&fPtr,_fileName.c_str(),READONLY,&fitsErr);
    xdbg<<"Done open"<<std::endl;
    if (fitsErr != 0) fits_report_error(stderr,fitsErr);
    Assert(fitsErr==0);

    fits_movabs_hdu(fPtr,_hdu,0,&fitsErr);
    xdbg<<"Moved to hdu "<<_hdu<<std::endl;
    if (fitsErr != 0) fits_report_error(stderr,fitsErr);
    Assert(fitsErr==0);

    int bitPix, nAxes;
    long sizes[2];
    fits_get_img_param(fPtr, int(2), &bitPix, &nAxes, sizes, &fitsErr);
    xdbg<<"done getimgparam"<<std::endl;
    xdbg<<"naxes = "<<nAxes<<std::endl;
    xdbg<<"bitpix = "<<bitPix<<std::endl;
    xdbg<<"FLOAT_IMG = "<<FLOAT_IMG<<std::endl;
    if (fitsErr != 0) fits_report_error(stderr,fitsErr);
    Assert(fitsErr==0);
    Assert(nAxes == 2);
    xdbg<<"sizes = "<<sizes[0]<<"  "<<sizes[1]<<std::endl;

    _xMin = 0;
    _xMax = sizes[0];
    _yMin = 0;
    _yMax = sizes[1];

    if (_xMin < x1) _xMin = x1;
    if (_xMax > x2) _xMax = x2; if (_xMax < _xMin+1) _xMax = _xMin+1;
    if (_yMin < y1) _yMin = y1;
    if (_yMax > y2) _yMax = y2; if (_yMax < _yMin+1) _yMax = _yMin+1;
    _source.reset(new TMatrix(T)(_xMax-_xMin,_yMax-_yMin));
    xdbg<<"done make matrix of image"<<std::endl;

    long fPixel[2] = {_xMin+1,_yMin+1};
    long lPixel[2] = {_xMax,_yMax};
    long inc[2] = {1,1};
    int anynul;
    xdbg<<"Before read_subset\n";
    Assert(getDataType<T>());
    fits_read_subset(fPtr,getDataType<T>(),fPixel,lPixel,inc,
                     0,TMV_ptr(*_source),&anynul,&fitsErr);
    xdbg<<"done read_subset "<<fitsErr<<std::endl;
    xdbg<<"anynul = "<<anynul<<std::endl;
    xdbg<<"fitserr = "<<fitsErr<<std::endl;
    if (fitsErr != 0) fits_report_error(stderr,fitsErr);
    Assert(fitsErr==0);

    _m.reset(new TMatrixView(T)(TMV_view(*_source)));
    xdbg<<"Done make matrixview"<<std::endl;

    fits_close_file(fPtr, &fitsErr);
    if (fitsErr != 0) fits_report_error(stderr,fitsErr);
    Assert(fitsErr==0);
    xdbg<<"Leaving Image ReadFits"<<std::endl;
}

template <typename T> 
void Image<T>::flush(std::string fileName, int hdu) const
{
    _fileName = fileName;
    _hdu = hdu;
    flush();
}

template <typename T> 
void Image<T>::flush() const
{
    Assert(_fileName != "");

    fitsfile *fPtr;
    int fitsErr=0;
    fits_open_file(&fPtr,_fileName.c_str(),READWRITE,&fitsErr);
    if (fitsErr != 0) fits_report_error(stderr,fitsErr);
    Assert(fitsErr==0);

    if (_hdu != 1) {
        int hduType;
        fits_movabs_hdu(fPtr,_hdu,&hduType,&fitsErr);
        if (fitsErr != 0) fits_report_error(stderr,fitsErr);
        Assert(fitsErr==0);
    }

    long fPixel[2] = {1,1};
    Assert(getDataType<T>());
    fits_write_pix(fPtr,getDataType<T>(),fPixel,long(_xMax*_yMax),
                   TMV_ptr(*_m),&fitsErr);
    if (fitsErr != 0) fits_report_error(stderr,fitsErr);
    Assert(fitsErr==0);

    fits_close_file(fPtr, &fitsErr);
    if (fitsErr != 0) fits_report_error(stderr,fitsErr);
    Assert(fitsErr==0);
}

template <typename T> 
void Image<T>::write(std::string fileName) const
{
    _fileName = fileName;
    _hdu = 1;

    fitsfile *fPtr;
    int fitsErr=0;
    fits_create_file(&fPtr,("!"+_fileName).c_str(),&fitsErr);
    if (fitsErr != 0) fits_report_error(stderr,fitsErr);
    Assert(fitsErr==0);

    Assert(getBitPix<T>());
    int bitPix = getBitPix<T>();
    int nAxes = 2;
    long sizes[2] = { _m->TMV_colsize(), _m->TMV_rowsize() };
    fits_create_img(fPtr, bitPix, nAxes, sizes, &fitsErr);
    if (fitsErr != 0) fits_report_error(stderr,fitsErr);
    Assert(fitsErr==0);

    long fPixel[2] = {1,1};
    Assert(getDataType<T>());
    fits_write_pix(fPtr,getDataType<T>(),fPixel,long(_xMax*_yMax),
                   TMV_ptr(*_m),&fitsErr);
    if (fitsErr != 0) fits_report_error(stderr,fitsErr);
    Assert(fitsErr==0);

    fits_close_file(fPtr, &fitsErr);
    if (fitsErr != 0) fits_report_error(stderr,fitsErr);
    Assert(fitsErr==0);
}

template <typename T> 
std::vector<Image<T>*> Image<T>::divide(int nX, int nY) const
{
    std::vector<int> x(nX+1);
    std::vector<int> y(nY+1);
    x[0] = _xMin;  x[nX] = _xMax;
    y[0] = _yMin;  y[nY] = _yMax;
    int xstep = (_xMax-_xMin)/nX;
    int ystep = (_yMax-_yMin)/nY;
    for(int i=1;i<nX;++i) x[i] = x[i-1]+xstep;
    for(int j=1;j<nY;++j) y[j] = y[j-1]+ystep;
    std::vector<Image*> blockImages;
    blockImages.reserve(nX*nY);
    for(int i=0;i<nX;++i) for(int j=0;j<nY;++j) {
#ifdef USE_TMV
        blockImages.push_back(
            new Image(_m->TMV_subMatrix(x[i],x[i+1],y[j],y[j+1]),
                      x[i],x[i+1],y[j],y[j+1]));
#else
        Assert(_m->TMV_colsize() == _source->TMV_colsize());
        Assert(_m->TMV_rowsize() == _source->TMV_rowsize());
        TMatrixView(T) subm = _source->TMV_subMatrix(x[i],x[i+1],y[j],y[j+1]);
        blockImages.push_back(new Image(subm,x[i],x[i+1],y[j],y[j+1]));
#endif
    }
    return blockImages;
}

template <typename T> 
T Image<T>::interpolate(double x, double y) const
{
    Assert(x>=double(_xMin) && x<double(_xMax));
    Assert(y>=double(_yMin) && y<double(_yMax));
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
    if (i==int(_xMax-1)) {--i; dx += 1.;}
    if (j==int(_yMax-1)) {--j; dy += 1.;}
    if (i==-1) {++i; dx -= 1.;}
    if (j==-1) {++j; dy -= 1.;}
    Assert(i>=0 && j>=0 && i+1<int(_m->TMV_colsize()) && 
           j+1<=int(_m->TMV_rowsize()));

    T f0 = (*_m)(i,j);
    T f1 = (*_m)(i+1,j);
    T f2 = (*_m)(i+1,j);
    T f3 = (*_m)(i+1,j+1);
    T dfdx = f1-f0;
    T dfdy = f2-f0;
    T d2fdxdy = f3+f0-f1-f2;
    return f0 + dfdx*dx + dfdy*dy + d2fdxdy*dx*dy;
}

template <typename T> 
T Image<T>::quadInterpolate(double x, double y) const
{
    Assert(x>=_xMin && x< _xMax);
    Assert(y>=_yMin && y< _yMax);
    int i = int (floor(x));
    double dx = x - (i+0.5);
    int j = int (floor(y));
    double dy = y - (j+0.5);
    Assert(i<int(_m->TMV_colsize()));
    Assert(j<int(_m->TMV_rowsize()));
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
    T f0 = (*_m)(i,j);
    T f1 = (i > 0) ? (*_m)(i-1,j) : 0.;
    T f2 = (i < int(_m->TMV_colsize())-1) ? (*_m)(i+1,j) : 0.;
    T f3 = (j > 0) ? (*_m)(i,j-1) : 0.;
    T f4 = (j < int(_m->TMV_rowsize())-1) ? (*_m)(i,j+1) : 0.;
    T f5 = (i > 0 && j > 0) ? (*_m)(i-1,j-1) : 0.;
    T f6 = (i < int(_m->TMV_colsize())-1 && j > 0) ? (*_m)(i+1,j-1) : 0.;
    T f7 = (i > 0 && j < int(_m->TMV_rowsize())-1) ? (*_m)(i-1,j+1) : 0.;
    T f8 = (i < int(_m->TMV_colsize())-1 && j < int(_m->TMV_rowsize())-1) ?
        (*_m)(i+1,j+1) : 0.;
    if (i == 0) {
        f1 = 2*f0 - f2;
        f5 = 2*f3 - f6;
        f7 = 2*f4 - f8;
    }
    if (i == int(_m->TMV_colsize())-1) {
        f2 = 2*f0 - f1;
        f6 = 2*f3 - f5;
        f8 = 2*f4 - f7;
    }
    if (j == 0) {
        f3 = 2*f0 - f4;
        f5 = 2*f1 - f7;
        f6 = 2*f2 - f8;
    }
    if (j == int(_m->TMV_rowsize())-1) {
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

template <typename T> 
T Image<T>::median() const
{
    std::vector<T> pixels;
    const int n1 = _m->TMV_colsize();
    const int n2 = _m->TMV_rowsize();
    pixels.reserve(n1*n2);
    for(int i=0;i<n1;++i) for(int j=0;j<n2;++j) {
        pixels.push_back((*_m)(i,j));
    }
    sort(pixels.begin(),pixels.end());
    return pixels[pixels.size()/2];
}

template class Image<double>;
//template class Image<float>;
