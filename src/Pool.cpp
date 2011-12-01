

#ifdef _WIN32
#pragma warning(push)
#pragma warning(disable:4018) //signed/unsigned mismatch
#pragma warning(disable:4290) // exception spec ignored
#endif

#include <exception>
#include <algorithm>

#ifdef VG
#include <valgrind/memcheck.h>
#endif

#include "dbg.h"
#include "Pool.h"

template <int blockSize>
Pool<blockSize>::Pool()
{
#ifdef VG
    VALGRIND_CREATE_MEMPOOL(&_allBlocks,0,false);
#endif
    PoolBlock* b = growPool();
    _freeBlocks.insert(b);
}

template <int blockSize>
Pool<blockSize>::~Pool()
{
    std::for_each(_allBlocks.begin(), _allBlocks.end(), Killer());
#ifdef VG
    VALGRIND_DESTROY_MEMPOOL(&_allBlocks);
#endif
}

template <int blockSize>
int Pool<blockSize>::totalMemoryUsed() const 
{ return _allBlocks.size() * blockSize; }

template <int blockSize>
void* Pool<blockSize>::allocate(size_t size)
{
    // Make sure the requested size isn't bigger than we can handle
    // given our blockSize.
    size_t size2 = size + sizeof(PoolBlock);
    if(size2>blockSize) throw std::bad_alloc();

    // Get the smallest free block that can handle the requested size:
    _comparisonBlock.size = size; // need this for lower_bound to work.
    PoolBlock* b;
    SetIt low = _freeBlocks.lower_bound(&_comparisonBlock);
    if (low == _freeBlocks.end()) {
        b = growPool();
    } else {
        b = *low;
        _freeBlocks.erase(low);
    }

    // If the block is the right size or only slightly larger,
    // then we can use the whole block for this allocation,
    // and we don't need to do anything with the remainder.
    if(b->size - size < 2*sizeof(PoolBlock)) {
        b->free = 0;
        char* ret = reinterpret_cast<char *>(b) + sizeof(PoolBlock);
#ifdef VG
        VALGRIND_MEMPOOL_ALLOC(&_allBlocks,ret,size);
#endif
        //checkp(ret,size);
        //check();
        return ret;
    } else {
        // Otherwise, we need to split the block into two pieces.
        // We return the first piece, and push the other piece back onto
        // the _freeBlocks list.
        PoolBlock* newBlock = reinterpret_cast<PoolBlock*>(
            reinterpret_cast<char*>(b) + size2);

        // Deal with the pointers:
        if (b->next) b->next->prev = newBlock;
        newBlock->next = b->next;
        b->next = newBlock;
        newBlock->prev = b;

        // Set the sizes:
        newBlock->size = b->size - size2;
        b->size = size;

        // Set the free flag:
        b->free = false;
        newBlock->free = true;

        // Push the new block onto the _freeBlocks list:
        _freeBlocks.insert(newBlock);

        char* ret = reinterpret_cast<char *>(b) + sizeof(PoolBlock);
#ifdef VG
        VALGRIND_MEMPOOL_ALLOC(&_allBlocks,ret,size);
#endif
        return ret;
    }
}

template <int blockSize>
void Pool<blockSize>::deallocate(void *p, size_t)
{
    // Trivial case: p = 0
    if(!p) return;

    // Get the block the corresponds to p.
    PoolBlock* b = reinterpret_cast<PoolBlock*>(
        static_cast<char*>(p) - sizeof(PoolBlock));

    // See if we need to combine this newly freed block with either 
    // the previous block or the next block:
    // 1) Combine with both:
    if (b->prev && b->next && b->prev->free && b->next->free) {
        _freeBlocks.erase(b->prev);
        _freeBlocks.erase(b->next);
        b->prev->size += b->size + b->next->size + 2*sizeof(PoolBlock);
        b->prev->next = b->next->next;
        if (b->next->next) b->next->next->prev = b->prev;
        // Need to erase and insert it to get it in the right place,
        // since _freeBlocks is sorted by size.
        _freeBlocks.insert(b->prev);
    } else if (b->prev && b->prev->free) {
        // 2) Combine with prev
        _freeBlocks.erase(b->prev);
        b->prev->size += b->size + sizeof(PoolBlock);
        b->prev->next = b->next;
        if (b->next) b->next->prev = b->prev;
        _freeBlocks.insert(b->prev);
    } else if (b->next && b->next->free) {
        // 3) Combine with next
        _freeBlocks.erase(b->next);
        b->size += b->next->size + sizeof(PoolBlock);
        b->next = b->next->next;
        if (b->next) b->next->prev = b;
        b->free = 1;
        _freeBlocks.insert(b);
    } else {
        // 4) No combining
        b->free = 1;
        _freeBlocks.insert(b);
    }
#ifdef VG
    VALGRIND_MEMPOOL_FREE(&_allBlocks,p);
#endif
}

template <int blockSize>
void Pool<blockSize>::dump(std::ostream& os)
{
    size_t tot_alloc=_allBlocks.size() * blockSize;
    size_t tot_used(0), tot_free(0), tot_free2(0);
    os<<"Pool<"<<blockSize<<"> dump:\n";
    os<<"_allBlocks:\n";
    for (ListIt block=_allBlocks.begin(); block!=_allBlocks.end(); ++block) {
        os<<(void*)*block<<" ... "<<(void*)(*block+blockSize)<<std::endl;
        for (PoolBlock* b = reinterpret_cast<PoolBlock*>(*block);
             b; b=b->next) {
            os<<(b->free ? "  F " : "    ");
            os<<(void*)b<<"  Size="<<b->size;
            os<<"   prev="<<(void*)b->prev<<", next="<<(void*)b->next;
            os<<std::endl;
            if (!b->free) tot_used += b->size;
            else tot_free += b->size;
        }
    }
    os<<"freeBlocks: \n";
    for (SetIt block=_freeBlocks.begin(); block!=_freeBlocks.end(); ++block) {
        os<<" Size="<<(*block)->size<<"   "<<(void*)(*block)<<std::endl;
        tot_free2 += (*block)->size;
    }
    os<<"total memory allocated = "<<tot_alloc/(1024*1024)<<" MB\n";
    os<<"total memory used = "<<tot_used/(1024*1024)<<" MB\n";
    os<<"total memory free = "<<tot_free/(1024*1024)<<" MB\n";
    os<<"(should also be)  = "<<tot_free2/(1024*1024)<<" MB\n";
}

template <int blockSize>
void Pool<blockSize>::summary(std::ostream& os)
{
    size_t tot_alloc=_allBlocks.size() * blockSize;
    size_t tot_used(0), tot_free(0);
    os<<"Pool<"<<blockSize<<"> summary:\n";
    for (ListIt block=_allBlocks.begin(); block!=_allBlocks.end(); ++block) {
        os<<(void*)*block<<" ... "<<(void*)(*block+blockSize)<<std::endl;
        for (PoolBlock* b = reinterpret_cast<PoolBlock*>(*block);
             b; b=b->next) {
            if (!b->free) tot_used += b->size;
            else tot_free += b->size;
        }
    }
    os<<"total memory allocated = "<<tot_alloc/(1024*1024)<<" MB\n";
    os<<"total memory used = "<<tot_used/(1024*1024)<<" MB\n";
    os<<"total memory free = "<<tot_free/(1024*1024)<<" MB\n";
}

template <int blockSize>
void Pool<blockSize>::checkp(char* p, int size, std::ostream& os)
{
    for (ListIt block=_allBlocks.begin(); block!=_allBlocks.end(); ++block) {
        if (p >= *block && (p+size) <= (*block + blockSize)) return;
    }
    os<<"p is not in _allBlocks\n";
    os<<"p = "<<(void*)p<<std::endl;
    os<<"p+size = "<<(void*)(p+size)<<std::endl;
    dump(os);
    exit(1);
}

template <int blockSize>
void Pool<blockSize>::check(std::ostream& os)
{
    std::set<PoolBlockPtr> freeBlocks2;

    for (ListIt block=_allBlocks.begin(); block!=_allBlocks.end(); ++block) {
        for (PoolBlock* b = reinterpret_cast<PoolBlock*>(*block);
             b; b=b->next) {
            int size2 = b->size + sizeof(PoolBlock);
            if (b->free) freeBlocks2.insert(b);

            // Check that this block fits into full block.
            if ((char*)b < *block) {
                os<<"b is off the front of the block\n";
                os<<"b = "<<(void*)b<<"  "<<b->size<<std::endl;
                os<<"*block = "<<(void*)(*block)<<std::endl;
                dump(os);
                exit(1);
            }
            if ((char*)b + size2 > *block + blockSize) {
                os<<"b extends off the back of the block\n";
                os<<"b = "<<(void*)b<<"  "<<b->size<<std::endl;
                os<<"b+size2 = "<<(void*)((char*)b+size2)<<std::endl;
                os<<"*block = "<<(void*)(*block)<<std::endl;
                os<<"*block + blockSize = "<<
                    (void*)(*block+blockSize)<<std::endl;
                dump(os);
                exit(1);
            }

            // Check that b->next is correct:
            if (b->next) {
                if ((char*)b+size2 != (char*)b->next) {
                    os<<"Bad b size:\n";
                    os<<"b = "<<(void*)b<<"  "<<b->size<<std::endl;
                    os<<"b+size2 = "<<(void*)((char*)b+size2)<<std::endl;
                    os<<"b->next = "<<(void*)(b->next)<<std::endl;
                    dump(os);
                    exit(1);
                }
                if (b->next->prev != b) {
                    os<<"Pointer error:\n";
                    os<<"b = "<<(void*)b<<"  "<<b->size<<std::endl;
                    os<<"b->next->prev = "<<
                        (void*)(b->next->prev)<<std::endl;
                    dump(os);
                    exit(1);
                }
            } else {
                // If no b->next, then make sure block fills the rest of 
                // pool block:
                if ((char*)b+size2 != *block+blockSize) {
                    os<<"Last b is wrong size:\n";
                    os<<"b = "<<(void*)b<<"  "<<b->size<<std::endl;
                    os<<"b+size2 = "<<(void*)((char*)b+size2)<<std::endl;
                    os<<"*block+blockSize = "<<
                        (void*)(*block+blockSize)<<std::endl;
                    dump(os);
                    exit(1);
                }
            }
        }
    }

    if (_freeBlocks != freeBlocks2) {
        os<<"freeBlocks doesn't match actual set of free blocks\n";
        os<<"Found:\n";
        for(SetIt block=freeBlocks2.begin();block!=freeBlocks2.end();++block) {
            os<<"  "<<(*block)->size<<"  "<<(void*)(*block)<<std::endl;
        }
        os<<"_freeBlocks structure has:\n";
        for(SetIt block=_freeBlocks.begin();block!=_freeBlocks.end();++block) {
            os<<"  "<<(*block)->size<<"  "<<(void*)(*block)<<std::endl;
        }
        dump(os);
        exit(1);
    }
}

template <int blockSize>
PoolBlock* Pool<blockSize>::growPool()
{
    char *p = new char[blockSize];
#ifdef VG
    VALGRIND_MAKE_MEM_NOACCESS(p,blockSize);
#endif
    _allBlocks.push_back(p);
    PoolBlock* newBlock = reinterpret_cast<PoolBlock*>(p);
    newBlock->prev = 0;
    newBlock->next = 0;
    newBlock->free = true;
    newBlock->size = blockSize-sizeof(PoolBlock);
    return newBlock;
}

template class Pool<PIXELLIST_BLOCK>;

#ifdef _WIN32
#pragma warning(pop)
#endif


