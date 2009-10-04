

#ifdef _WIN32
#pragma warning(push)
#pragma warning(disable:4018) //signed/unsigned mismatch
#pragma warning(disable:4290) // exception spec ignored
#endif

#include <exception>
#include <algorithm>
#include <iostream>  //for dump()

#include "dbg.h"
#include "pool.h"

#ifdef VG
#include <valgrind/memcheck.h>
#endif

template <int block_size>
pool<block_size>::pool()
{
  //std::cerr<<"Initialize pool\n";
#ifdef VG
  VALGRIND_CREATE_MEMPOOL(&pool_mem,0,false);
#endif
  //std::cerr<<"Adding block "<<(void*)pool_mem.back()<<" ... "<<(void*)(pool_mem.back()+block_size)<<std::endl;
  //std::cerr<<"total memory = "<<total_memory_used()/(1024.*1024.)<<" MB\n";
  block* b = grow_pool();
  free_blocks.insert(b);
  //check();
}

template <int block_size>
pool<block_size>::~pool()
{
  std::for_each(pool_mem.begin(), pool_mem.end(), killer());
#ifdef VG
  VALGRIND_DESTROY_MEMPOOL(&pool_mem);
#endif
}

template <int block_size>
size_t pool<block_size>::total_memory_used() const 
{ return pool_mem.size() * block_size; }

template <int block_size>
void* pool<block_size>::allocate(size_t size)
{
  //std::cerr<<"start allocate "<<size<<std::endl;
  //dump();
  //check();

  // Make sure the requested size isn't bigger than we can handle
  // given our block_size.
  size_t size2 = size + sizeof(block);
  if(size2>block_size) throw std::bad_alloc();

  // Get the smallest free block that can handle the requested size:
  comp_block.size = size; // need this for lower_bound to work.
  block* b;
  setit low = free_blocks.lower_bound(&comp_block);
  if (low == free_blocks.end()) 
  {
    b = grow_pool();
    //std::cerr<<"from grow_pool: b.size = "<<b->size<<std::endl;
  }
  else
  {
    b = *low;
    free_blocks.erase(low);
    //std::cerr<<"from lower_bound: b.size = "<<b->size<<std::endl;
  }
  //std::cerr<<"b = "<<(void*)b<<std::endl;
  //std::cerr<<"b->size = "<<b->size<<std::endl;

  // If the block is the right size or only slightly larger,
  // then we can use the whole block for this allocation,
  // and we don't need to do anything with the remainder.
  if(b->size - size < 2*sizeof(block))
  {
    //std::cerr<<"block works, and the rest is small"<<std::endl;
    b->free = 0;
    char* ret = reinterpret_cast<char *>(b) + sizeof(block);
    //std::cerr<<"return "<<(void*)ret<<std::endl;
#ifdef VG
    VALGRIND_MEMPOOL_ALLOC(&pool_mem,ret,size);
#endif
    //checkp(ret,size);
    //check();
    return ret;
  }
  // Otherwise, we need to split the block into two pieces.
  // We return the first piece, and push the other piece back onto
  // the free_blocks list.
  else
  {
    //std::cerr<<"block works, and the rest is not small"<<std::endl;
    block* new_block = reinterpret_cast<block*>(
	reinterpret_cast<char*>(b) + size2);
    //std::cerr<<"new block = "<<(void*)b<<" + "<<size2<<std::endl;
    //std::cerr<<"new block = "<<(void*)new_block<<std::endl;

    // Deal with the pointers:
    if (b->next) b->next->prev = new_block;
    new_block->next = b->next;
    b->next = new_block;
    new_block->prev = b;

    // Set the sizes:
    new_block->size = b->size - size2;
    //std::cerr<<"new_block size = "<<b->size<<" - "<<size2<<std::endl;
    //std::cerr<<"new_block size = "<<new_block->size<<std::endl;
    b->size = size;
    //std::cerr<<"b size = "<<b->size<<std::endl;

    // Set the free flag:
    b->free = false;
    new_block->free = true;

    // Push the new block onto the free_blocks list:
    free_blocks.insert(new_block);

    char* ret = reinterpret_cast<char *>(b) + sizeof(block);
    //std::cerr<<"return "<<(void*)ret<<std::endl;
#ifdef VG
    VALGRIND_MEMPOOL_ALLOC(&pool_mem,ret,size);
#endif
    //checkp(ret,size);
    //check();
    return ret;
  }
}

template <int block_size>
void pool<block_size>::deallocate(void *p, size_t)
{
  //std::cerr<<"start deallocate "<<p<<std::endl;
  //dump();
  //check();

  // Trivial case: p = 0
  if(!p) return;

  // Get the block the corresponds to p.
  block* b = reinterpret_cast<block*>(
      static_cast<char*>(p) - sizeof(block));
  //std::cerr<<"b = "<<(void*)b<<std::endl;

  // See if we need to combine this newly freed block with either 
  // the previous block or the next block:
  // 1) Combine with both:
  if (b->prev && b->next && b->prev->free && b->next->free)
  {
    //std::cerr<<"combine with prev and next\n";
    free_blocks.erase(b->prev);
    free_blocks.erase(b->next);
    b->prev->size += b->size + b->next->size + 2*sizeof(block);
    b->prev->next = b->next->next;
    if (b->next->next) b->next->next->prev = b->prev;
    // Need to erase and insert it to get it in the right place,
    // since free_blocks is sorted by size.
    free_blocks.insert(b->prev);
  }
  // 2) Combine with prev
  else if (b->prev && b->prev->free)
  {
    //std::cerr<<"combine with prev\n";
    free_blocks.erase(b->prev);
    b->prev->size += b->size + sizeof(block);
    b->prev->next = b->next;
    if (b->next) b->next->prev = b->prev;
    free_blocks.insert(b->prev);
  }
  // 3) Combine with next
  else if (b->next && b->next->free)
  {
    //std::cerr<<"combine with next\n";
    free_blocks.erase(b->next);
    b->size += b->next->size + sizeof(block);
    b->next = b->next->next;
    if (b->next) b->next->prev = b;
    b->free = 1;
    free_blocks.insert(b);
  }
  // 4) No combining
  else 
  {
    //std::cerr<<"don't combine\n";
    b->free = 1;
    free_blocks.insert(b);
  }
#ifdef VG
  VALGRIND_MEMPOOL_FREE(&pool_mem,p);
#endif
  //std::cerr<<"done deallocate"<<std::endl;
  //check();
}

template <int block_size>
void pool<block_size>::dump()
{
  std::cerr<<"dump:\n";
  std::cerr<<"pool_mem:\n";
  for (listit it=pool_mem.begin(); it!=pool_mem.end(); ++it)
  {
    std::cerr<<(void*)*it<<" ... "<<(void*)(*it+block_size)<<std::endl;
    for (block* b = reinterpret_cast<block*>(*it); b; b=b->next)
    {
      std::cerr<<(b->free ? "  F " : "    ");
      std::cerr<<(void*)b<<"  Size="<<b->size;
      std::cerr<<"   prev="<<(void*)b->prev<<", next="<<(void*)b->next;
      std::cerr<<std::endl;
    }
  }
  std::cerr<<"free_blocks: \n";
  for (setit it=free_blocks.begin(); it!=free_blocks.end(); ++it)
  {
    std::cerr<<" Size="<<(*it)->size<<"   "<<(void*)(*it)<<std::endl;
  }
  std::cerr<<"end dump\n";
}

template <int block_size>
void pool<block_size>::checkp(char* p, int size)
{
  typename std::list<char*>::iterator it;
  for (it=pool_mem.begin(); it!=pool_mem.end(); ++it)
  {
    if (p >= *it && (p+size) <= (*it + block_size)) return;
  }
  std::cerr<<"p is not in pool_mem\n";
  std::cerr<<"p = "<<(void*)p<<std::endl;
  std::cerr<<"p+size = "<<(void*)(p+size)<<std::endl;
  dump();
  exit(1);
}

template <int block_size>
void pool<block_size>::check()
{
  std::set<block_ptr> free_blocks2;

  for (listit blocks=pool_mem.begin(); blocks!=pool_mem.end(); ++blocks)
  {
    for (block* b = reinterpret_cast<block*>(*blocks); b; b=b->next)
    {
      int size2 = b->size + sizeof(block);
      if (b->free) free_blocks2.insert(b);

      // Check that this block fits into full pool_mem block.
      if ((char*)b < *blocks)
      {
	std::cerr<<"b is off the front of the block\n";
	std::cerr<<"b = "<<(void*)b<<"  "<<b->size<<std::endl;
	std::cerr<<"*blocks = "<<(void*)(*blocks)<<std::endl;
	dump();
	exit(1);
      }
      if ((char*)b + size2 > *blocks + block_size)
      {
	std::cerr<<"b extends off the back of the block\n";
	std::cerr<<"b = "<<(void*)b<<"  "<<b->size<<std::endl;
	std::cerr<<"b+size2 = "<<(void*)((char*)b+size2)<<std::endl;
	std::cerr<<"*blocks = "<<(void*)(*blocks)<<std::endl;
	std::cerr<<"*blocks + block_size = "<<(void*)(*blocks+block_size)<<std::endl;
	dump();
	exit(1);
      }

      // Check that b->next is correct:
      if (b->next) 
      {
	if ((char*)b+size2 != (char*)b->next)
	{
	  std::cerr<<"Bad b size:\n";
	  std::cerr<<"b = "<<(void*)b<<"  "<<b->size<<std::endl;
	  std::cerr<<"b+size2 = "<<(void*)((char*)b+size2)<<std::endl;
	  std::cerr<<"b->next = "<<(void*)(b->next)<<std::endl;
	  dump();
	  exit(1);
	}
	if (b->next->prev != b)
	{
	  std::cerr<<"Pointer error:\n";
	  std::cerr<<"b = "<<(void*)b<<"  "<<b->size<<std::endl;
	  std::cerr<<"b->next->prev = "<<(void*)(b->next->prev)<<std::endl;
	  dump();
	  exit(1);
	}
      }
      // If no b->next, then make sure block fills the rest of pool block:
      else
      {
	if ((char*)b+size2 != *blocks+block_size)
	{
	  std::cerr<<"Last b is wrong size:\n";
	  std::cerr<<"b = "<<(void*)b<<"  "<<b->size<<std::endl;
	  std::cerr<<"b+size2 = "<<(void*)((char*)b+size2)<<std::endl;
	  std::cerr<<"*blocks+block_size = "<<(void*)(*blocks+block_size)<<std::endl;
	  dump();
	  exit(1);
	}
      }
    }
  }

  if (free_blocks != free_blocks2)
  {
    std::cerr<<"free_blocks doesn't match actual set of free blocks\n";
    std::cerr<<"Found:\n";
    for(setit it=free_blocks2.begin();it!=free_blocks2.end();++it)
    {
      std::cerr<<"  "<<(*it)->size<<"  "<<(void*)(*it)<<std::endl;
    }
    std::cerr<<"free_blocks structure has:\n";
    for(setit it=free_blocks.begin();it!=free_blocks.end();++it)
    {
      std::cerr<<"  "<<(*it)->size<<"  "<<(void*)(*it)<<std::endl;
    }
    dump();
    exit(1);
  }
}

template <int block_size>
block* pool<block_size>::grow_pool()
{
  //std::cerr<<"Start grow:"<<std::endl;
  char *p = new char[block_size];
#ifdef VG
  VALGRIND_MAKE_MEM_NOACCESS(p,block_size);
#endif
  pool_mem.push_back(p);
  //std::cerr<<"Adding block "<<(void*)pool_mem.back()<<" ... "<<(void*)(pool_mem.back()+block_size)<<std::endl;
  //std::cerr<<"total memory = "<<total_memory_used()/(1024.*1024.)<<" MB\n";
  block* new_block = reinterpret_cast<block*>(p);
  new_block->prev = 0;
  new_block->next = 0;
  new_block->free = true;
  new_block->size = block_size-sizeof(block);
  //std::cerr<<"done grow"<<std::endl;
  return new_block;
}


// Need to instantiate all block_sizes that you want to use:

#include "Pixel.h"
template class pool<PIXELLIST_BLOCK>;

#ifdef _WIN32
#pragma warning(pop)
#endif


