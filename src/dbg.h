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

#if defined(__GNUC__) && defined(OPENMP_LINK)
extern __thread std::ostream* dbgout;
extern __thread bool XDEBUG;
#else
extern std::ostream* dbgout;
extern bool XDEBUG;
#endif


#ifdef _OPENMP
#pragma omp threadprivate( dbgout , XDEBUG )
#endif

struct AssertFailureException :
    public std::runtime_error
{
    AssertFailureException(const char* e) : std::runtime_error(
        std::string("Error - Assert ") + e + " failed") {}
};

#ifdef NDEBUG
#define dbg if (false) (*dbgout)
#define xdbg if (false) (*dbgout)
#define xxdbg if (false) (*dbgout)
#define Assert(x)
#else
#define dbg if (dbgout) (*dbgout)
#define xdbg if (dbgout && XDEBUG) (*dbgout)
#define xxdbg if (false) (*dbgout)
#define Assert(x) \
    do { \
        if(!(x)) { \
            dbg << "Error - Assert " #x " failed"<<std::endl; \
            dbg << "on line "<<__LINE__<<" in file "<<__FILE__<<std::endl; \
            throw AssertFailureException(#x); \
        } \
    } while(false)
#endif

#endif
