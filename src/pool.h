#ifndef POOL_H_INCLUDED_GF
#define POOL_H_INCLUDED_GF


#ifdef _WIN32
#pragma warning(push)
#pragma warning(disable:4018) //signed/unsigned mismatch
#pragma warning(disable:4290) // exception spec ignored
#endif

#include <exception>
#include <list>
#include <algorithm>
#include <iostream>  //for dump()
#include <valgrind/memcheck.h>

template <int block_size>
class pool
{
  public :

    pool()
    {
      //std::cout<<"Initialize pool\n";
      pool_mem.push_back(new char[block_size]);
      VALGRIND_CREATE_MEMPOOL(&pool_mem,0,false);
      VALGRIND_MAKE_MEM_NOACCESS(pool_mem.back(),block_size);
      std::cout<<"Adding block "<<(void*)pool_mem.back()<<" ... "<<(void*)(pool_mem.back()+block_size)<<std::endl;
      std::cout<<"total memory = "<<total_memory_used()/(1024.*1024.)<<" MB\n";
      blocks = reinterpret_cast<block*>(*(pool_mem.begin()));
      blocks->prev = 0;
      blocks->next = 0;
      blocks->free = true;
      blocks->last = true;
      blocks->size = block_size - sizeof(block);
      check();
    };

    ~pool()
    {
      std::for_each(pool_mem.begin(), pool_mem.end(), killer());
      VALGRIND_DESTROY_MEMPOOL(&pool_mem);
    }

    size_t total_memory_used() const 
    { return pool_mem.size() * block_size; }

    void* allocate(size_t size)
    {
      size_t size2 = size + sizeof(block);
      //std::cout<<"start allocate "<<std::ios::hex<<size<<std::endl;
      //dump();
      check();
      if(size2>block_size) throw std::bad_alloc();
      block *b = blocks;
      while(!b->free || b->size < size)
      {
	if(!b->next) grow_pool(b);
	b = b->next;
      }
      //std::cout<<"b = "<<b<<std::endl;
      //std::cout<<"b->size = "<<b->size<<std::endl;
      if(b->size - size < 2*sizeof(block))
      {
	//std::cout<<"block works, and the rest is small"<<std::endl;
	b->free = 0;
	//std::cout<<"return "<<b<<std::endl;
	checkp((char*)b + sizeof(block),size);
	check();
	char* ret = reinterpret_cast<char *>(b) + sizeof(block);
	VALGRIND_MEMPOOL_ALLOC(&pool_mem,ret,size);
	return ret;
      }
      else
      {
	//std::cout<<"block works, and the rest is not small"<<std::endl;
	block* new_block = reinterpret_cast<block*>(
	    reinterpret_cast<char*>(b) + size2);
	//std::cout<<"new block = "<<b<<" + "<<size2<<std::endl;
	//std::cout<<"new block = "<<new_block<<std::endl;

	if (b->next) b->next->prev = new_block;
	new_block->next = b->next;
	b->next = new_block;
	new_block->prev = b;

	new_block->size = b->size - size2;
	b->size = size;
	//std::cout<<"new_block size = "<<b->size<<" - "<<size2<<std::endl;
	//std::cout<<"new_block size = "<<new_block->size<<std::endl;
	//std::cout<<"b size = "<<b->size<<std::endl;

	b->free = false;
	new_block->free = true;

	if (b->last) new_block->last = true;
	else new_block->last = false;
	b->last = false;

	//std::cout<<"return "<<b<<std::endl;
	checkp((char*)b+sizeof(block),size);
	check();
	char* ret = reinterpret_cast<char *>(b) + sizeof(block);
	VALGRIND_MEMPOOL_ALLOC(&pool_mem,ret,size);
	return ret;
      }
    }

    void deallocate(void *p, size_t = 0)
    {
      //std::cout<<"start deallocate "<<p<<std::endl;
      check();
      if(!p) return;
      block* b = reinterpret_cast<block*>(
	  static_cast<char*>(p) - sizeof(block));
      //std::cout<<"b = "<<b<<std::endl;
      if (b->prev && b->next && !b->prev->last && !b->last && 
	  b->prev->free && b->next->free)
      {
	//std::cout<<"combine with prev and next\n";
	b->prev->size += b->size + b->next->size + 2*sizeof(block);
	if (b->next->last) b->prev->last = true;
	b->prev->next = b->next->next;
	if (b->next->next) b->next->next->prev = b->prev;
      }
      else if (b->prev && !b->prev->last && b->prev->free)
      {
	//std::cout<<"combine with prev\n";
	b->prev->size += b->size + sizeof(block);
	if (b->last) b->prev->last = true;
	b->prev->next = b->next;
	if (b->next) b->next->prev = b->prev;
      }
      else if (!b->last && b->next && b->next->free)
      {
	//std::cout<<"combine with next\n";
	b->size += b->next->size + sizeof(block);
	if (b->next->last) b->last = true;
	b->next = b->next->next;
	if (b->next) b->next->prev = b;
	b->free = 1;
      }
      else 
      {
	//std::cout<<"don't combine\n";
	b->free = 1;
      }
      //std::cout<<"done deallocate"<<std::endl;
      check();
      VALGRIND_MEMPOOL_FREE(&pool_mem,p);
    }
    void dump()
    {
      std::cout<<"dump:\n";
      std::cout<<"pool_mem:\n";
      typename std::list<char*>::iterator it;
      for (it=pool_mem.begin(); it!=pool_mem.end(); ++it)
      {
	std::cout<<(void*)*it<<" ... "<<(void*)(*it+block_size)<<std::endl;
      }
      for (block* b = blocks; b; b=b->next)
      {
	if (b->free) std::cout<<"F ";
	else std::cout<<"  ";
	if (b->last) std::cout<<"L ";
	else std::cout<<"  ";
	std::cout<<b<<"  Size="<<b->size;
	std::cout<<"   prev="<<b->prev<<", next="<<b->next;
	std::cout<<std::endl;
      }
      std::cout<<"end dump\n";
    }
    void checkp(char* p, int size)
    {
      typename std::list<char*>::iterator it;
      for (it=pool_mem.begin(); it!=pool_mem.end(); ++it)
      {
	if (p >= *it && (p+size) <= (*it + block_size)) return;
      }
      std::cout<<"p is not in pool_mem\n";
      std::cout<<"p = "<<(void*)p<<std::endl;
      std::cout<<"p+size = "<<(void*)(p+size)<<std::endl;
      dump();
      exit(1);
    }
    void check()
    {
      for (block* b = blocks; b; b=b->next)
      {
	int size2 = b->size + sizeof(block);
	checkp((char*)b,size2);
	if (b->next) 
	{
	  if ((char*)b+size2 > (char*)b->next && !b->last)
	  {
	    std::cout<<"Bad b size (1):\n";
	    std::cout<<"b = "<<b<<"  "<<b->size<<std::endl;
	    std::cout<<"b->next = "<<b->next<<std::endl;
	    std::cout<<"b+size2 = "<<(void*)((char*)b+size2)<<std::endl;
	    dump();
	    exit(1);
	  }
	  if ((char*)b+size2 != (char*)b->next)
	  {
	    // then b->next should be in pool_mem, and b should be last
	    bool nextisinpool=false;
	    typename std::list<char*>::iterator it;
	    for (it=pool_mem.begin(); it!=pool_mem.end(); ++it)
	      if ((char*)b->next == *it) nextisinpool = true;
	    if (!nextisinpool)
	    { 
	      std::cout<<"Bad b size (2):\n";
	      std::cout<<"b = "<<b<<"  "<<b->size<<std::endl;
	      std::cout<<"b->next = "<<b->next<<std::endl;
	      std::cout<<"b->next is no in pool\n";
	      dump();
	      exit(1);
	    }
	    if (!b->last)
	    {
	      std::cout<<"Bad b size (3):\n";
	      std::cout<<"b = "<<b<<"  "<<b->size<<std::endl;
	      std::cout<<"b->next = "<<b->next<<std::endl;
	      std::cout<<"b is not marked as last\n";
	      dump();
	      exit(1);
	    }
	  }
	}
      }
    }

  private :

    std::list<char*> pool_mem;
    struct block
    {
      block *prev;
      block *next;
      size_t size;
      bool free;
      bool last;

      block(block *prev_, block *next_, size_t size_, bool free_, bool last_) : 
	prev(prev_), next(next_), size(size_), free(free_), last(last_) {}
      ~block() {}
    };
    block* blocks;
    pool(const pool &);
    pool& operator=(const pool&);
    struct killer
    {
      void operator()(char *p){delete [] p;}
    };
    static void kill(char *p){delete [] p;}
    void grow_pool(block *b)
    {
      //std::cout<<"Start grow:"<<std::endl;
      check();
      block* new_block;
      char *p = new char[block_size];
      VALGRIND_MAKE_MEM_NOACCESS(p,block_size);
      pool_mem.push_back(p);
      std::cout<<"Adding block "<<(void*)pool_mem.back()<<" ... "<<(void*)(pool_mem.back()+block_size)<<std::endl;
      std::cout<<"total memory = "<<total_memory_used()/(1024.*1024.)<<" MB\n";
      new_block = reinterpret_cast<block*>(p);
      new_block->prev = b;
      new_block->next = false;
      new_block->last = true;
      new_block->free = 1;
      new_block->size = block_size-sizeof(block);
      b->next = new_block;
      //std::cout<<"done grow"<<std::endl;
      //dump();
      check();
    }
};

#ifdef _WIN32
#pragma warning(pop)
#endif

#endif

