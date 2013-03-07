#ifndef ImageH
#define ImageH

#include <memory>
#include "MyMatrix.h"
#include "dbg.h"
#include "Bounds.h"
#include "ConfigFile.h"

#include <CCfits/CCfits>
#include <fitsio.h>


template <typename T> 
class Image 
{

public:

    Image() {};

    // Read new image from file
    Image(std::string fits_file, int hdu=1); 
    void load(std::string fits_file, int hdu=1);

    Image(std::string fits_file, int hdu,
          int x1, int x2, int y1, int y2);

    Image(std::string fits_file, int hdu, const Bounds& b);

    // Create blank image
    Image(int x_size, int y_size) : 
        _filename(""), _hdu(0),
        _xmin(0),_xmax(x_size),_ymin(0),_ymax(y_size),
        _source(new TMatrix(T)(x_size,y_size)),
        _m(new TMatrixView(T)(TMV_view(*_source))) 
    { _source->setZero(); }

    // New copy of image
    Image(const Image& rhs) : 
        _filename(""), _hdu(0),
        _xmin(rhs._xmin), _xmax(rhs._xmax), _ymin(rhs._ymin), _ymax(rhs._ymax),
        _source(new TMatrix(T)(*rhs._m)),
        _m(new TMatrixView(T)(TMV_view(*_source))) {}

    // subImage (with new storage)
    Image(const Image& rhs, int x1, int x2, int y1, int y2) :
        _filename(""), _hdu(0),
        _xmin(x1), _xmax(x2), _ymin(y1), _ymax(y2),
        _source(new TMatrix(T)(rhs._m->TMV_subMatrix(x1,x2,y1,y2))),
        _m(new TMatrixView(T)(TMV_view(*_source))) {}
    Image(const Image& rhs, const Bounds& b) :
        _filename(""), _hdu(0),
        _xmin(int(floor(b.getXMin()))), _xmax(int(ceil(b.getXMax()))), 
        _ymin(int(floor(b.getYMin()))), _ymax(int(ceil(b.getYMax()))),
        _source(new TMatrix(T)(rhs._m->TMV_subMatrix(_xmin,_xmax,_ymin,_ymax))),
        _m(new TMatrixView(T)(TMV_view(*_source))) {}

    // Read image given configuration parameters
    Image(const ConfigFile& params);
    Image(const ConfigFile& params, std::auto_ptr<Image<T> >& weight_im);

    // Read partial image 
    Image(const ConfigFile& params, int x1, int x2, int y1, int y2);
    Image(const ConfigFile& params, const Bounds& b);
    Image(const ConfigFile& params, std::auto_ptr<Image<T> >& weight_im,
          int x1, int x2, int y1, int y2);
    Image(const ConfigFile& params, std::auto_ptr<Image<T> >& weight_im,
          const Bounds& b);

    ~Image();

    // Copy image
    void operator=(const Image& rhs) { *_m = *rhs._m; }

    // Copy rhs to a subimage of this
    void copy(const Image& rhs, int x1, int x2, int y1, int y2)
    { _m->TMV_subMatrix(x1,x2,y1,y2) = *rhs._m; }

    // Write back to existing file
    void flush() const;
    void flush(std::string fits_file, int hdu=1) const; 

    // Write to new file.  Append if create is false
    void write(std::string fits_file, bool create=true) const; 
    // Write to existing fits object.  *much* faster than re-opening
    // each time
    void write(fitsfile *fits) const; 

    TConstMatrixView(T) getM() const { return *_m; }
    TMV_const TMatrixView(T)& getM() { return *_m; }
    int getXMin() const { return _xmin; }
    int getXMax() const { return _xmax; }
    int getYMin() const { return _ymin; }
    int getYMax() const { return _ymax; }
    int getMaxI() const { return _m->TMV_colsize()-1; }
    int getMaxJ() const { return _m->TMV_rowsize()-1; }

    // Access elements
    T& operator()(int i,int j) { return (*_m)(i,j); }
    T operator()(int i,int j) const { return (*_m)(i,j); }

    // Add two images
    void operator+=(const Image& rhs) { *_m += *rhs._m; }

    // Erase all values
    void clear() { _m->setZero(); }

    bool isSquare() { return _m->TMV_colsize() == _m->TMV_rowsize(); }

    Bounds getBounds() const 
    { return Bounds(_xmin,_xmax,_ymin,_ymax); }

    // subImage refers to the same storage as this.
    Image subImage(int x1, int x2, int y1, int y2)
    {
#ifdef USE_TMV
        return Image(_m->TMV_subMatrix(x1,x2,y1,y2),x1,x2,y1,y2); 
#else
        // I only ever use this function when *this is a full Image.
        // But in fact, the functionality is possible for SubImages
        // if we are using TMV rather than Eigen
        Assert(_m->TMV_colsize() == _source->TMV_colsize());
        Assert(_m->TMV_rowsize() == _source->TMV_rowsize());
        TMatrixView(T) subm = _source->TMV_subMatrix(x1,x2,y1,y2);
        return Image(subm,x1,x2,y1,y2); 
#endif
    }

    // Split into nx x ny subimages
    std::vector<Image*> divide(int nx, int ny) const; 

    // Iterpolate between integral pixel values
    T interpolate(double x, double y) const;
    T quadInterpolate(double x, double y) const;

    // Median value of all pixels
    T median() const;

    bool loaded() { return _loaded; };
private:

    mutable std::string _filename;
    mutable int _hdu;

    int _xmin,_xmax,_ymin,_ymax;
    std::auto_ptr<TMatrix(T)> _source;
    std::auto_ptr<TMatrixView(T)> _m;

    template <class M>
    Image(TMV_const M& m, int x1, int x2, int y1, int y2) :
        _xmin(x1), _xmax(x2), _ymin(y1), _ymax(y2),
        _source(0), _m(new TMatrixView(T)(m)) {}

    void readFits();
    void readFits(int x1, int x2, int y1, int y2);

    bool _loaded;
};

extern template class Image<double>;
extern template class Image<float>;

#endif

