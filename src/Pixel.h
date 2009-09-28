#ifndef PIXEL_H
#define PIXEL_H

#include <complex>
#include <string>
#include "Image.h"
#include "Transformation.h"
#include "ConfigFile.h"
#include "pool_allocator.h"

struct Pixel 
{ 
  Pixel() : z(0.), I(0.), wt(0.) {}
  Pixel(double _u, double _v, double _I, double _wt) :
    z(_u,_v), I(_I), wt(_wt) {}
  ~Pixel() {}

  std::complex<double> z;
  double I,wt;
};

#define PIXELLIST_BLOCK 1024*1024*100  // 100 MB per block

#if 1
typedef std::vector<Pixel> PixelList;
#else
class PixelList
{
  public :

    inline PixelList() : 
      use_block_mem(false), v1(new std::vector<Pixel>()), v2(0) {}
    inline PixelList(const int n) : 
      use_block_mem(false), v1(new std::vector<Pixel>(n)), v2(0) {}
    inline PixelList(const PixelList& p2);
    inline PixelList& operator=(const PixelList& p2);
    inline ~PixelList() {}

    // My pool allocator is not thread-safe, so this should only be used
    // for PixelLists that won't be modified in a parallel block.
    inline void UseBlockMem() 
    { 
      // This should be done before any elements are added.
      if (v1.get()) Assert(v1->size() == 0);
      if (v2.get()) Assert(v2->size() == 0);
      v1.reset(0);
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
    std::auto_ptr<std::vector<Pixel> > v1;
    typedef pool_allocator<Pixel,PIXELLIST_BLOCK> pool_alloc;
    std::auto_ptr<std::vector<Pixel,pool_alloc> > v2;
};

PixelList::PixelList(const PixelList& p2) : 
  use_block_mem(false), v1(new std::vector<Pixel>(p2.size())), v2(0) 
{
  if (p2.use_block_mem) 
  { std::copy(p2.v2->begin(),p2.v2->end(),v1->begin()); }
  else *v1 = *p2.v1;
}

PixelList& PixelList::operator=(const PixelList& p2)
{
  if (use_block_mem)
  {
    if (p2.use_block_mem) *v2 = *p2.v2;
    else
    { std::copy(p2.v1->begin(),p2.v1->end(),v2->begin()); }
  }
  else 
  {
    if (p2.use_block_mem) 
    { std::copy(p2.v2->begin(),p2.v2->end(),v1->begin()); }
    else *v1 = *p2.v1;
  }
  return *this;
}
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
