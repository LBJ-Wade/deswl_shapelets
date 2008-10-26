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

extern std::ostream* dbgout;
extern bool XDEBUG;

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
    do { if(!(x)) { \
      dbg << "Error - Assert " #x " failed"<<std::endl; \
      dbg << "on line "<<__LINE__<<" in file "<<__FILE__<<std::endl; \
      std::cerr << "Error - Assert " #x " failed"<<std::endl; \
      std::cerr << "on line "<<__LINE__<<" in file "<<__FILE__<<std::endl; \
      exit(1);} } while(false)
#endif

inline void myerror(const char *s,const char *s2 = "")
{
  dbg << "Error: " << s << ' ' << s2 << std::endl;
  std::cerr << "Error: " << s << ' ' << s2 << std::endl;
  exit(1);
}

#endif
