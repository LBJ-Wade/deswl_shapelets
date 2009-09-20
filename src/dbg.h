//---------------------------------------------------------------------------
#ifndef dbgH
#define dbgH
//---------------------------------------------------------------------------

/* Put the following in the main program file:

std::ostream* dbgout=0;
bool XDEBUG=false;

*/

//
// Set up debugging stuff
//

#include <iostream>
#include <stdexcept>

extern std::ostream* dbgout;
extern bool XDEBUG;

#ifdef _OPENMP
#pragma omp threadprivate( dbgout , XDEBUG )
#endif

struct AssertFailure :
  public std::runtime_error
{
  AssertFailure(const char* e) : std::runtime_error(
      std::string("Error - Assert ") + e + " failed") {}
};

#ifdef NDEBUG
  #define dbg if (false) (*dbgout)
  #define xdbg if (false) (*dbgout)
  #define xxdbg if (false) (*dbgout)
  #define ompdbg if (false) (*dbgout)
  #define ompxdbg if (false) (*dbgout)
  #define Assert(x)
#else
  #define dbg if (dbgout) (*dbgout)
  #define xdbg if (dbgout && XDEBUG) (*dbgout)
  #define xxdbg if (false) (*dbgout)
#ifdef _OPENMP
  #define ompdbg if (false) (*dbgout)
  #define ompxdbg if (false) (*dbgout)
#else
  #define ompdbg if (dbgout) (*dbgout)
  #define ompxdbg if (dbgout && XDEBUG) (*dbgout)
#endif
  #define Assert(x) \
    do { if(!(x)) { \
      dbg << "Error - Assert " #x " failed"<<std::endl; \
      dbg << "on line "<<__LINE__<<" in file "<<__FILE__<<std::endl; \
      throw AssertFailure(#x); } \
    } while(false)
#endif

#endif
