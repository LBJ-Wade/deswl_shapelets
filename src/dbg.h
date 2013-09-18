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
#include <sstream>
#include <fstream>

#ifdef MEM_TEST
// valarray includes somethings that get messed up if it is included after
// mmgr.h.  So make sure it gets included first.
#include <valarray>
#include "mmgr.h"
#endif

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

// Some examples of information that can be obtained from /proc/self/status:
// Pid, VmPeak, VmSize, VmHWM, VmRss, Threads
// ( See http://www.kernel.org/doc/man-pages/online/pages/man5/proc.5.html
//   for more information.)

// Get a single item from /proc/self/status.
// If os is provided, output all the lines.
inline double get_proc_stat(std::string stat, std::ostream* os=0) 
{
    std::ifstream proc("/proc/self/status");
    if (!proc) {
        if (os) *os<<"Could not open /proc/self/status";
        return -1.; // Negative value reports an error.
    }
    if (os) *os << "Reading information from /proc/self/status:\n";
    std::string line;
    double value;
    while(getline(proc, line), !proc.fail()) {
        if (os) *os <<line<<std::endl;
        if(line.substr(0, stat.size()) == stat) {
            std::istringstream iss(line);
            std::string name;
            iss >> name >> value;
            if (!os) return value;
        }
    }
    return value;
}

// Return the memory being used in MBytes.
inline double memory_usage(std::ostream* os=0) 
{ return get_proc_stat("VmSize",os) / 1024.; }

inline double peak_memory_usage(std::ostream* os=0)
{ return get_proc_stat("VmPeak",os) / 1024.; }

#endif
