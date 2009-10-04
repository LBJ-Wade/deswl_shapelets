
#include "Pixel.h"

PixelList::PixelList() :
  use_block_mem(false), v1(new std::vector<Pixel>()) {}

PixelList::PixelList(const int n) :
  use_block_mem(false), v1(new std::vector<Pixel>(n)) {}

PixelList::~PixelList()
{
#ifdef _OPENMP
#pragma omp critical
#endif
  {
    v2.reset();
  }
}

void PixelList::UseBlockMem() 
{
  // This should be done before any elements are added.
  if (v1.get()) Assert(v1->size() == 0);
  if (v2.get()) Assert(v2->size() == 0);
  v1.reset();
#ifdef _OPENMP
#pragma omp critical
#endif
  {
    v2.reset(new std::vector<Pixel,pool_alloc>());
  }
  use_block_mem = true; 
}

size_t PixelList::size() const
{
  if (use_block_mem) return v2->size();
  else return v1->size();
}

void PixelList::reserve(const int n)
{
  if (use_block_mem) 
#ifdef _OPENMP
#pragma omp critical
#endif
  {
    v2->reserve(n);
  }
  else v1->reserve(n);
}

void PixelList::resize(const int n)
{
  if (use_block_mem) 
#ifdef _OPENMP
#pragma omp critical
#endif
  {
    v2->resize(n);
  }
  else v1->resize(n);
}

void PixelList::clear()
{
  if (use_block_mem) 
#ifdef _OPENMP
#pragma omp critical
#endif
  {
    v2->clear();
  }
  else v1->clear();
}

void PixelList::push_back(const Pixel& p)
{
  if (use_block_mem) 
#ifdef _OPENMP
#pragma omp critical
#endif
  {
    v2->push_back(p);
  }
  else v1->push_back(p);
}

Pixel& PixelList::operator[](const int i)
{
  if (use_block_mem) return (*v2)[i];
  else return (*v1)[i];
}

const Pixel& PixelList::operator[](const int i) const
{
  if (use_block_mem) return (*v2)[i];
  else return (*v1)[i];
}

