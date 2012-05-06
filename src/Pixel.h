#ifndef PIXEL_H
#define PIXEL_H

#define PIXELLIST_USE_POOL

#include <complex>
#include <string>

#ifdef __INTEL_COMPILER
#pragma warning (disable : 1418)
#endif
#include "boost/shared_ptr.hpp"
#ifdef __INTEL_COMPILER
#pragma warning (default : 1418)
#endif

#include "dbg.h"
#include "Image.h"
#include "Transformation.h"
#include "ConfigFile.h"
#include "Bounds.h"

#ifdef PIXELLIST_USE_POOL
#define PIXELLIST_BLOCK 1024*1024*100  // 100 MB per block
#include "PoolAllocator.h"
#endif

class Pixel 
{ 
public :

    Pixel() : _pos(0.), _flux(0.), _inverse_sigma(0.) {}

    Pixel(double u, double v, double flux, double inverse_sigma) :
        _pos(u,v), _flux(flux), _inverse_sigma(inverse_sigma) {}

    Pixel(std::complex<double> z, double flux, double inverse_sigma) :
        _pos(z), _flux(flux), _inverse_sigma(inverse_sigma) {}

    ~Pixel() {}

    std::complex<double> getPos() const { return _pos; }

    double getFlux() const { return _flux; }

    double getInverseSigma() const { return _inverse_sigma; }

    void setPos(const std::complex<double>& pos) { _pos = pos; }

    void setFlux(const double flux) { _flux = flux; }

    void setInverseSigma(const double inverse_sigma) 
    { _inverse_sigma = inverse_sigma; }

private :

    std::complex<double> _pos;
    double _flux;
    double _inverse_sigma;
};

// Most of these methods are not (intrinsically) thread-safe, 
// since they might be using my pool allocator, so they need to 
// be wrapped in a critical block.
// Therefore, all the methods for PixelList are defined in Pixel_omp.cpp.

class PixelList
{
public :

    PixelList();
    PixelList(const int n);

    // Note: Copy constructor and op= use shared ownership semantics.
    // This means a vector<PixelList> is efficient when it does
    // push_back, etc.
    PixelList(const PixelList& rhs);
    PixelList& operator=(const PixelList& rhs);
    ~PixelList();

    // These mimic the same functionality of a std::vector<Pixel>
    int size() const;
    void reserve(const int n);
    int capacity() const;
    void resize(const int n);
    void clear();
    void push_back(const Pixel& p);
    Pixel& operator[](const int i);
    const Pixel& operator[](const int i) const;
    void sort(const Position& cen);

    // Start not using Pool allocator.  Turn it on with this:
    void usePool();
    static void dumpPool(std::ostream& os);
    static void reclaimMemory();

private :

    bool _use_pool;
    boost::shared_ptr<std::vector<Pixel> > _v1;
#ifdef PIXELLIST_USE_POOL
    typedef PoolAllocator<Pixel,PIXELLIST_BLOCK> PoolAllocPixel;
    boost::shared_ptr<std::vector<Pixel,PoolAllocPixel> > _v2;
#else
    boost::shared_ptr<std::vector<Pixel> > _v2;
#endif

};

void GetPixList(
    const Image<double>& im, PixelList& pix,
    const Position cen, double sky, double noise,
    const Image<double>* weight_image, const Transformation& trans,
    double aperture, const ConfigFile& params, long& flag);

double GetLocalSky(
    const Image<double>& bkg, 
    const Position cen, const Transformation& trans, double aperture,
    const ConfigFile& params, long& flag);

void GetSubPixList(
    PixelList& pix, const PixelList& allpix,
    std::complex<double> cen_offset, std::complex<double> shear,
    double aperture, double inner_fake_ap, double outer_fake_ap,
    const ConfigFile& params, long& flag);

#endif
