#ifndef POOL_H_INCLUDED_GF
#define POOL_H_INCLUDED_GF


#include <list>
#include <set>

// Define the block structure that we will use in pool:
struct block
{
  block* prev;
  block* next;
  size_t size;
  bool free;
};

struct block_ptr
{
  block* p;
  block_ptr(block* _p) : p(_p) {}
  operator block*() const { return p; }
  operator char*() const { return (char*)p; }
  operator void*() const { return (void*)p; }
  block& operator*() const { return *p; }
  block* operator->() const { return p; }
};

// Sort by size first, then ptr if same size.
inline bool operator<(const block_ptr b1, const block_ptr b2)
{
  if ( b1->size == b2->size ) return (void*)b1 < (void*)b2;
  else return (b1->size < b2->size);
}
inline bool operator==(const block_ptr b1, const block_ptr b2)
{
  return ( b1->size == b2->size );
}


template <int block_size>
class pool
{
  public :

    pool();
    ~pool();

    size_t total_memory_used() const;
    void* allocate(size_t size);
    void deallocate(void *p, size_t = 0);

    void dump();
    void checkp(char* p, int size);
    void check();

  private :

    // Private data:
    std::list<char*> pool_mem;
    typedef std::list<char*>::iterator listit;

    std::set<block_ptr> free_blocks;
    typedef std::set<block_ptr>::iterator setit;

    mutable block comp_block;

    // Copy and assign are not defined:
    pool(const pool &);
    pool& operator=(const pool&);

    // killer struct used by destructor to delete memory
    struct killer
    {
      void operator()(char *p){delete [] p;}
    };
    static void kill(char *p){delete [] p;}

    // Private member function to allocate more memory
    block* grow_pool();
};

#endif

