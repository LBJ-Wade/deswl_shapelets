#ifndef TESTHELPER_H
#define TESTHELPER_H

#include <stdexcept>

extern bool showtests;
extern bool dontthrow;
extern std::string lastsuccess;

inline void PreTest(std::string s)
{
  if (showtests) {
    std::cout<<"Trying: "<<s;
    std::cout.flush();
  }
}

inline void DoTest(bool x, std::string s)
{
  if (x) {
    if (showtests) std::cout<<"  Passed"<<std::endl;
    lastsuccess = s;
  } else {
    if (showtests) std::cout<<"  Failed"<<std::endl;
    if (dontthrow) std::cout<<"Failed test: "<<s<<std::endl;
    else 
#ifdef NOTHROW
    { std::cout<<"Error in test: "<<s<<std::endl; exit(1); }
#else
      throw std::runtime_error(std::string("Error in test: ") + s);
#endif
  }
}

#define Test(x,s) \
  do { \
    PreTest(s); \
    DoTest(x,s); \
  } while (false)

#endif
