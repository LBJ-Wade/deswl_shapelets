#ifndef PIXEL_H
#define PIXEL_H

#include <complex>
#include <string>
#include "Image.h"
#include "Transformation.h"
#include "ConfigFile.h"

struct Pixel 
{ 
  Pixel() : z(0.), I(0.), wt(0.) {}
  Pixel(double _u, double _v, double _I, double _wt) :
    z(_u,_v), I(_I), wt(_wt) {}
  ~Pixel() {}

  std::complex<double> z;
  double I,wt;
};


#ifdef __INTEL_COMPILER
#pragma warning (disable : 1418)
#endif
#include "boost/shared_ptr.hpp"
#ifdef __INTEL_COMPILER
#pragma warning (default : 1418)
#endif

#include "pool_allocator.h"
#define PIXELLIST_BLOCK 1024*1024*100  // 100 MB per block

// Most of these methods are not (intrinsically) thread-safe, 
// since they might be using my pool allocator, so they need to 
// be wrapped in a critical block.
// Therefore, all the methods for PixelList are defined in Pixel_omp.cpp.

class PixelList
{
  public :

    PixelList();
    PixelList(const int n);
    ~PixelList();

    void UseBlockMem();

    size_t size() const;
    void reserve(const int n);
    void resize(const int n);
    void clear();
    void push_back(const Pixel& p);
    Pixel& operator[](const int i);
    const Pixel& operator[](const int i) const;

  private :

    bool use_block_mem;
    boost::shared_ptr<std::vector<Pixel> > v1;
    typedef pool_allocator<Pixel,PIXELLIST_BLOCK> pool_alloc;
    boost::shared_ptr<std::vector<Pixel,pool_alloc> > v2;
};

void GetPixList(const Image<double>& im, PixelList& pix,
    const Position cen, double sky, double noise, double gain,
    const Image<double>* wt_im, const Transformation& trans,
    double aperture, long& flag);

double GetLocalSky(const Image<float>& bkg, 
    const Position cen, const Transformation& trans,
    double aperture, long& flag);

void GetSubPixList(PixelList& pix,
    const PixelList& allpix,
    double aperture, long& flag);

#endif
