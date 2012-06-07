
#include "Pixel.h"

PixelList::PixelList() :
    _use_pool(false), _v1(new std::vector<Pixel>()) 
{}

PixelList::PixelList(const int n) :
    _use_pool(false), _v1(new std::vector<Pixel>(n)) 
{}

PixelList::PixelList(const PixelList& rhs) :
    _use_pool(rhs._use_pool), _v1(rhs._v1), _v2(rhs._v2) 
{}

PixelList& PixelList::operator=(const PixelList& rhs)
{
    _use_pool = rhs._use_pool;
    _v1 = rhs._v1;
#ifdef _OPENMP
#pragma omp critical (PixelList)
#endif
    {
        _v2 = rhs._v2;
    }
    return *this;
}

PixelList::~PixelList()
{
#ifdef _OPENMP
#pragma omp critical (PixelList)
#endif
    {
        _v2.reset();
    }
}

void PixelList::usePool() 
{
#ifdef PIXELLIST_USE_POOL
    // This should be done before any elements are added.
    if (_v1.get()) Assert(_v1->size() == 0);
    if (_v2.get()) Assert(_v2->size() == 0);
    _v1.reset();
#ifdef _OPENMP
#pragma omp critical (PixelList)
#endif
    {
        _v2.reset(new std::vector<Pixel,PoolAllocPixel>());
    }
    _use_pool = true; 
#endif
}

void PixelList::dumpPool(std::ostream& os) 
{
#ifdef PIXELLIST_USE_POOL
    os<<"Dump Pool: \n";
    if (XDEBUG) PoolAllocPixel::dump(os);
    else PoolAllocPixel::summary(os);
#else
    os<<"Not using PoolAlloc\n";
#endif
}

void PixelList::reclaimMemory()
{
#ifdef PIXELLIST_USE_POOL
    PoolAllocPixel::reclaimMemory();
#endif
}

int PixelList::size() const
{
    if (_use_pool) return _v2->size();
    else return _v1->size();
}

void PixelList::reserve(const int n)
{
    if (_use_pool) {
#ifdef _OPENMP
#pragma omp critical (PixelList)
#endif
        {
            _v2->reserve(n);
        }
    } else {
        _v1->reserve(n);
    }
}

int PixelList::capacity() const
{ return _use_pool ? _v2->capacity() : _v1->capacity(); }

void PixelList::resize(const int n)
{
    if (_use_pool) {
#ifdef _OPENMP
#pragma omp critical (PixelList)
#endif
        {
            _v2->resize(n);
        }
    } else {
        _v1->resize(n);
    }
}

void PixelList::clear()
{
    if (_use_pool) {
#ifdef _OPENMP
#pragma omp critical (PixelList)
#endif
        {
            _v2->clear();
        }
    } else {
        _v1->clear();
    }
}

void PixelList::push_back(const Pixel& p)
{
    if (_use_pool) {
#ifdef _OPENMP
#pragma omp critical (PixelList)
#endif
        {
            _v2->push_back(p);
        }
    } else {
        _v1->push_back(p);
    }
}

Pixel& PixelList::operator[](const int i)
{
    if (_use_pool) return (*_v2)[i];
    else return (*_v1)[i];
}

const Pixel& PixelList::operator[](const int i) const
{
    if (_use_pool) return (*_v2)[i];
    else return (*_v1)[i];
}

struct PixelListSorter
{
    Position _cen;
    PixelListSorter(const Position& cen) : _cen(cen) {}
    bool operator()(const Pixel& p1, const Pixel& p2) const
    { return std::norm(p1.getPos()-_cen) < std::norm(p2.getPos()-_cen); }
};

void PixelList::sort(const Position& cen) 
{
    PixelListSorter sorter(cen);
    if (_use_pool) std::sort(_v2->begin(),_v2->end(),sorter);
    else std::sort(_v1->begin(),_v1->end(),sorter);
}
