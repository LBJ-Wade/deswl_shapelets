#ifndef POOL_H_INCLUDED_GF
#define POOL_H_INCLUDED_GF


#include <list>
#include <set>

// Define the block structure that we will use in Pool:
struct PoolBlock
{
    PoolBlock* prev;
    PoolBlock* next;
    size_t size;
    bool free;
};

struct PoolBlockPtr
{
    PoolBlock* p;
    PoolBlockPtr(PoolBlock* _p) : p(_p) {}
    operator PoolBlock*() const { return p; }
    operator char*() const { return (char*)p; }
    operator void*() const { return (void*)p; }
    PoolBlock& operator*() const { return *p; }
    PoolBlock* operator->() const { return p; }
};

// Sort by size first, then ptr if same size.
inline bool operator<(const PoolBlockPtr& b1, const PoolBlockPtr& b2)
{
    if ( b1->size == b2->size ) return (void*)b1 < (void*)b2;
    else return (b1->size < b2->size);
}
inline bool operator==(const PoolBlockPtr b1, const PoolBlockPtr b2)
{
    return ( b1->size == b2->size );
}


template <int blockSize>
class Pool
{
public :

    Pool();
    ~Pool();

    size_t totalMemoryUsed() const;
    void* allocate(size_t size);
    void deallocate(void *p, size_t = 0);

    void dump();
    void checkp(char* p, int size);
    void check();

private :

    // Private data:
    std::list<char*> _allBlocks;
    typedef std::list<char*>::iterator ListIt;

    std::set<PoolBlockPtr> _freeBlocks;
    typedef std::set<PoolBlockPtr>::iterator SetIt;

    // This is just used as a reference comparison for lower_bound.
    mutable PoolBlock _comparisonBlock;

    // Copy and assign are not defined:
    Pool(const Pool &);
    Pool& operator=(const Pool&);

    // killer struct used by destructor to delete memory
    struct Killer
    {
        void operator()(char *p){delete [] p;}
    };
    static void kill(char *p){delete [] p;}

    // Private member function to allocate more memory
    PoolBlock* growPool();
};

// Need to instantiate all blockSizes that you want to use:

#include "Pixel.h"
extern template class Pool<PIXELLIST_BLOCK>;

#endif

