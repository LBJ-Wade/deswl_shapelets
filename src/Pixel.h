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


#if 0

typedef std::vector<Pixel> PixelList;

#elif 1

#include "boost/shared_ptr.hpp"

//#include "pool_allocator.h"
//#define PIXELLIST_BLOCK 1024*1024*100  // 100 MB per block

#include "boost/pool/pool_alloc.hpp"
#define PIXELLIST_BLOCK 
// Still define it, since used in #ifdef, but value isn't relevant anymore.
// Eventually change this to something like "USE_POOL" or somesuch.

class PixelList
{
  public :

    inline PixelList() : 
      use_block_mem(false), v1(new std::vector<Pixel>()) {}
    inline PixelList(const int n) : 
      use_block_mem(false), v1(new std::vector<Pixel>(n)) {}
    inline ~PixelList() {}

    // My pool allocator is not thread-safe, so this should only be used
    // for PixelLists that won't be modified in a parallel block.
    inline void UseBlockMem() 
    { 
      // This should be done before any elements are added.
      if (v1.get()) Assert(v1->size() == 0);
      if (v2.get()) Assert(v2->size() == 0);
      v1.reset();
      v2.reset(new std::vector<Pixel,pool_alloc>());
      use_block_mem = true; 
    }

    inline size_t size() const
    {
      if (use_block_mem) return v2->size();
      else return v1->size();
    }
    inline void reserve(const int n)
    {
      if (use_block_mem) v2->reserve(n);
      else v1->reserve(n);
    }
    inline void resize(const int n)
    {
      if (use_block_mem) v2->resize(n);
      else v1->resize(n);
    }
    inline void clear()
    {
      if (use_block_mem) v2->clear();
      else v1->clear();
    }
    inline void push_back(const Pixel& p)
    {
      if (use_block_mem) v2->push_back(p);
      else v1->push_back(p);
    }
    inline Pixel& operator[](const int i)
    {
      if (use_block_mem) return (*v2)[i];
      else return (*v1)[i];
    }
    inline const Pixel& operator[](const int i) const
    {
      if (use_block_mem) return (*v2)[i];
      else return (*v1)[i];
    }

  private :

    bool use_block_mem;
    boost::shared_ptr<std::vector<Pixel> > v1;
    //typedef pool_allocator<Pixel,PIXELLIST_BLOCK> pool_alloc;
    typedef boost::pool_allocator<Pixel> pool_alloc;
    boost::shared_ptr<std::vector<Pixel,pool_alloc> > v2;
};

#else // This one never uses the pool_allocator, but does use shared_ptr

#include "boost/shared_ptr.hpp"

class PixelList
{
  public :

    inline PixelList() : v1(new std::vector<Pixel>()) {}
    inline PixelList(const int n) : v1(new std::vector<Pixel>(n)) {}
    inline ~PixelList() {}

    inline size_t size() const
    { return v1->size(); }
    inline void reserve(const int n)
    { v1->reserve(n); }
    inline void resize(const int n)
    { v1->resize(n); }
    inline void clear()
    { v1->clear(); }
    inline void push_back(const Pixel& p)
    { v1->push_back(p); }
    inline Pixel& operator[](const int i)
    { return (*v1)[i]; }
    inline const Pixel& operator[](const int i) const
    { return (*v1)[i]; }

  private :

    boost::shared_ptr<std::vector<Pixel> > v1;
};

#endif

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
