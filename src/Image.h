//---------------------------------------------------------------------------
#ifndef ImageH
#define ImageH
//---------------------------------------------------------------------------

#include "TMV.h"
#include "Bounds.h"
#include "ConfigFile.h"

template <class T> class Image {

  public:

    // Read new image from file
    Image(std::string fitsfile, int hdu=1); 

    // Create blank image
    Image(size_t xsize, size_t ysize) : 
      filename(""), hdu(0),
      xmin(0),xmax(xsize),ymin(0),ymax(ysize),
      sourcem(new tmv::Matrix<T,tmv::ColMajor>(xsize,ysize,0.)),
      itsm(new tmv::MatrixView<T>(sourcem->View())) {}

    // New copy of image
    Image(const Image& rhs) : 
      filename(""), hdu(0),
      xmin(rhs.xmin), xmax(rhs.xmax), ymin(rhs.ymin), ymax(rhs.ymax),
      sourcem(new tmv::Matrix<T,tmv::ColMajor>(*rhs.itsm)),
      itsm(new tmv::MatrixView<T>(sourcem->View())) {}

    // Subimage (with new storage)
    Image(const Image& rhs, size_t x1, size_t x2, size_t y1, size_t y2) :
      filename(""), hdu(0),
      xmin(x1), xmax(x2), ymin(y1), ymax(y2),
      sourcem(new tmv::Matrix<T,tmv::ColMajor>(
	    rhs.itsm->SubMatrix(x1,x2,y1,y2))),
      itsm(new tmv::MatrixView<T>(sourcem->View())) {}

    // Read image given configuration parameters
    Image(const ConfigFile& params);
    Image(const ConfigFile& params, std::auto_ptr<Image<T> >& weight_im);

    ~Image() {}

    // Copy image
    void operator=(const Image& rhs) { *itsm = *rhs.itsm; }

    // Copy rhs to a subimage of this
    void Copy(const Image& rhs, size_t x1, size_t x2, size_t y1, size_t y2)
    { itsm->SubMatrix(x1,x2,y1,y2) = *rhs.itsm; }

    // Write back to existing file
    void Flush() const;
    void Flush(std::string fitsfile, int hdu=1) const; 

    // Write to new file
    void Write(std::string fitsfile) const; 

    tmv::ConstMatrixView<T> GetM() const { return *itsm; }
    const tmv::MatrixView<T>& GetM() { return *itsm; }
    int GetXMin() const { return xmin; }
    int GetXMax() const { return xmax; }
    int GetYMin() const { return ymin; }
    int GetYMax() const { return ymax; }
    size_t GetMaxI() const { return itsm->colsize()-1; }
    size_t GetMaxJ() const { return itsm->rowsize()-1; }

    // Access elements
    T& operator()(size_t i,size_t j) { return (*itsm)(i,j); }
    T operator()(size_t i,size_t j) const { return (*itsm)(i,j); }

    // Add two images
    void operator+=(const Image& rhs) { *itsm += *rhs.itsm; }

    // Erase all values
    void Clear() { itsm->Zero(); }

    bool IsSquare() { return itsm->IsSquare(); }

    Bounds GetBounds() const 
    { return Bounds(xmin,xmax,ymin,ymax); }

    // SubImage refers to the same storage as this.
    Image SubImage(size_t x1, size_t x2, size_t y1, size_t y2)
    { return Image(itsm->SubMatrix(x1,x2,y1,y2),x1,x2,y1,y2); }
      
    // Split into nx x ny subimages
    std::vector<Image*> Divide(size_t nx, size_t ny) const; 

    // Iterpolate between integral pixel values
    T Interpolate(double x, double y) const;
    T QuadInterpolate(double x, double y) const;

    // Median value of all pixels
    T Median() const;

  private:

    mutable std::string filename;
    mutable int hdu;

    int xmin,xmax,ymin,ymax;
    std::auto_ptr<tmv::Matrix<T,tmv::ColMajor> > sourcem;
    std::auto_ptr<tmv::MatrixView<T> > itsm;

    Image(const tmv::MatrixView<T>& m, 
	size_t x1, size_t x2, size_t y1, size_t y2) :
      xmin(x1), xmax(x2), ymin(y1), ymax(y2),
      sourcem(0), itsm(new tmv::MatrixView<T>(m)) {}

    void ReadFits(std::string fitsfile, int hdu=1); 

};

#endif

