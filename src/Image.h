#ifndef ImageH
#define ImageH

#include "TMV.h"
#include "Bounds.h"
#include "ConfigFile.h"

template <typename T> 
class Image 
{

public:

    // Read new image from file
    Image(std::string fitsFile, int hdu=1); 

    Image(std::string fitsFile, int hdu,
          int x1, int x2, int y1, int y2);

    Image(std::string fitsFile, int hdu, const Bounds& b);

    // Create blank image
    Image(int xSize, int ySize) : 
        _fileName(""), _hdu(0),
        _xMin(0),_xMax(xSize),_yMin(0),_yMax(ySize),
        _source(new tmv::Matrix<T,tmv::ColMajor>(xSize,ySize,0.)),
        _m(new tmv::MatrixView<T>(_source->view())) {}

    // New copy of image
    Image(const Image& rhs) : 
        _fileName(""), _hdu(0),
        _xMin(rhs._xMin), _xMax(rhs._xMax), _yMin(rhs._yMin), _yMax(rhs._yMax),
        _source(new tmv::Matrix<T,tmv::ColMajor>(*rhs._m)),
        _m(new tmv::MatrixView<T>(_source->view())) {}

    // subImage (with new storage)
    Image(const Image& rhs, int x1, int x2, int y1, int y2) :
        _fileName(""), _hdu(0),
        _xMin(x1), _xMax(x2), _yMin(y1), _yMax(y2),
        _source(new tmv::Matrix<T,tmv::ColMajor>(
                rhs._m->subMatrix(x1,x2,y1,y2))),
        _m(new tmv::MatrixView<T>(_source->view())) {}
    Image(const Image& rhs, const Bounds& b) :
        _fileName(""), _hdu(0),
        _xMin(int(floor(b.getXMin()))), _xMax(int(ceil(b.getXMax()))), 
        _yMin(int(floor(b.getYMin()))), _yMax(int(ceil(b.getYMax()))),
        _source(new tmv::Matrix<T,tmv::ColMajor>(
                rhs._m->subMatrix(_xMin,_xMax,_yMin,_yMax))),
        _m(new tmv::MatrixView<T>(_source->view())) {}

    // Read image given configuration parameters
    Image(const ConfigFile& params);
    Image(const ConfigFile& params, std::auto_ptr<Image<T> >& weightIm);

    // Read partial image 
    Image(const ConfigFile& params, int x1, int x2, int y1, int y2);
    Image(const ConfigFile& params, const Bounds& b);
    Image(const ConfigFile& params, std::auto_ptr<Image<T> >& weightIm,
          int x1, int x2, int y1, int y2);
    Image(const ConfigFile& params, std::auto_ptr<Image<T> >& weightIm,
          const Bounds& b);

    ~Image() {}

    // Copy image
    void operator=(const Image& rhs) { *_m = *rhs._m; }

    // Copy rhs to a subimage of this
    void copy(const Image& rhs, int x1, int x2, int y1, int y2)
    { _m->subMatrix(x1,x2,y1,y2) = *rhs._m; }

    // Write back to existing file
    void flush() const;
    void flush(std::string fitsFile, int hdu=1) const; 

    // Write to new file
    void write(std::string fitsFile) const; 

    tmv::ConstMatrixView<T> getM() const { return *_m; }
    const tmv::MatrixView<T>& getM() { return *_m; }
    int getXMin() const { return _xMin; }
    int getXMax() const { return _xMax; }
    int getYMin() const { return _yMin; }
    int getYMax() const { return _yMax; }
    int getMaxI() const { return _m->colsize()-1; }
    int getMaxJ() const { return _m->rowsize()-1; }

    // Access elements
    T& operator()(int i,int j) { return (*_m)(i,j); }
    T operator()(int i,int j) const { return (*_m)(i,j); }

    // Add two images
    void operator+=(const Image& rhs) { *_m += *rhs._m; }

    // Erase all values
    void clear() { _m->setZero(); }

    bool isSquare() { return _m->isSquare(); }

    Bounds getBounds() const 
    { return Bounds(_xMin,_xMax,_yMin,_yMax); }

    // subImage refers to the same storage as this.
    Image subImage(int x1, int x2, int y1, int y2)
    { return Image(_m->subMatrix(x1,x2,y1,y2),x1,x2,y1,y2); }

    // Split into nx x ny subimages
    std::vector<Image*> divide(int nx, int ny) const; 

    // Iterpolate between integral pixel values
    T interpolate(double x, double y) const;
    T quadInterpolate(double x, double y) const;

    // Median value of all pixels
    T median() const;

private:

    mutable std::string _fileName;
    mutable int _hdu;

    int _xMin,_xMax,_yMin,_yMax;
    std::auto_ptr<tmv::Matrix<T,tmv::ColMajor> > _source;
    std::auto_ptr<tmv::MatrixView<T> > _m;

    Image(const tmv::MatrixView<T>& m, 
          int x1, int x2, int y1, int y2) :
        _xMin(x1), _xMax(x2), _yMin(y1), _yMax(y2),
        _source(0), _m(new tmv::MatrixView<T>(m)) {}

    void readFits();
    void readFits(int x1, int x2, int y1, int y2);

};

extern template class Image<double>;
extern template class Image<float>;

#endif

