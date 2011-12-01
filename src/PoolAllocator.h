#ifndef POOL_ALLOCATOR_H_INCLUDED_GF
#define POOL_ALLOCATOR_H_INCLUDED_GF

// This is based on a pool_allocator found at:
// http://www.codeguru.com/cpp/cpp/cpp_mfc/stl/article.php/c4079

#include <ostream>
#include "dbg.h"
#include "Pool.h"

// Note: Most of the names of things here are mandated by the standard
// library.  So they don't always conform to the LSST style.
//
template <typename T, int blockSize> class PoolAllocator;

template <int blockSize> class PoolAllocator<void,blockSize>
{
public:

    typedef void* pointer;
    typedef const void* const_pointer;
    // reference to void members are impossible.
    typedef void value_type;

    template <class U>
    struct rebind { typedef PoolAllocator<U,blockSize> other; };
};    

namespace PoolAlloc
{
    inline void destruct(char *) {}
    inline void destruct(wchar_t*) {}

    template <typename T> 
    inline void destruct(T *t) { t->~T(); }
} // namespace

template <typename T, int blockSize>
class PoolAllocator
{

public:

    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef T value_type;

    template <class U> 
    struct rebind { typedef PoolAllocator<U,blockSize> other; };

    PoolAllocator() {}

    pointer address(reference x) const { return &x; }

    const_pointer address(const_reference x) const { return &x; }

    pointer allocate(
        size_type size, 
        typename PoolAllocator<void,blockSize>::const_pointer = 0)
    {
        return static_cast<pointer>(mem.allocate(size*sizeof(T)));
    }

    //for Dinkumware:
    char *_Charalloc(size_type n) 
    { return static_cast<char*>(mem.allocate(n)); }
    // end Dinkumware

    template <class U> PoolAllocator(const PoolAllocator<U,blockSize>&) {}

    void deallocate(pointer p, size_type n)
    {
        mem.deallocate(p, n);
    }

    void deallocate(void *p, size_type n)
    {
        mem.deallocate(p, n);
    }

    size_type max_size() const throw() 
    { return size_t(-1) / sizeof(value_type); }

    void construct(pointer p, const T& val)
    {
        new(static_cast<void*>(p)) T(val);
    }

    void construct(pointer p)
    {
        new(static_cast<void*>(p)) T();
    }

    void destroy(pointer p) { PoolAlloc::destruct(p); }

    static void dump(std::ostream& os) { mem.dump(os); };

    static size_t totalMemoryUsed() { return mem.totalMemoryUsed(); }

private:

    static Pool<blockSize> mem;
};

template <typename T, int blockSize> 
Pool<blockSize> PoolAllocator<T,blockSize>::mem;

template <typename T, typename U, int bs>
inline bool operator==(
    const PoolAllocator<T,bs>&, const PoolAllocator<U,bs>&)
{ return true; }

template <typename T, typename U, int bs>
inline bool operator!=(
    const PoolAllocator<T,bs>&, const PoolAllocator<U,bs>&)
{ return false; }


// For VC6/STLPort 4-5-3 see /stl/_alloc.h, line 464
//"If custom allocators are being used without member template classes support :
// user (on purpose) is forced to define rebind/get operations !!!"
#ifdef _WIN32
#define POOL_ALLOC_CDECL __cdecl
#else
#define POOL_ALLOC_CDECL
#endif

namespace std
{
    template <class _Tp1, class _Tp2, int bs>
    inline PoolAllocator<_Tp2,bs>& POOL_ALLOC_CDECL
    __stl_alloc_rebind(PoolAllocator<_Tp1,bs>& __a, const _Tp2*) 
    {  
        return (PoolAllocator<_Tp2,bs>&)(__a); 
    }


    template <class _Tp1, class _Tp2, int bs>
    inline PoolAllocator<_Tp2,bs> POOL_ALLOC_CDECL
    __stl_alloc_create(const PoolAllocator<_Tp1,bs>&, const _Tp2*) 
    { 
        return PoolAllocator<_Tp2,bs>(); 
    }

} // namespace std
// end STLPort

#endif

