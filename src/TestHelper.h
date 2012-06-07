#ifndef TESTHELPER_H
#define TESTHELPER_H

#include <stdexcept>

extern bool show_tests;
extern bool should_throw;
extern std::string last_success;
extern std::ostream* testout;

inline void PreTest(std::string s)
{
    if (show_tests) {
        (*testout)<<"Trying: "<<s;
        (*testout).flush();
    }
}

inline void DoTest(bool x, std::string s)
{
    if (x) {
        if (show_tests) (*testout)<<"  Passed"<<std::endl;
        last_success = s;
    } else {
        if (show_tests) (*testout)<<"  Failed"<<std::endl;
        if (!should_throw) (*testout)<<"Failed test: "<<s<<std::endl;
        else  {
#ifdef NOTHROW
            (*testout)<<"Error in test: "<<s<<std::endl;
            exit(1); 
#else
            throw std::runtime_error(std::string("Error in test: ") + s);
#endif
        }
    }
}

#define Test(x,s) \
    do { \
        PreTest(s); \
        DoTest(x,s); \
    } while (false)

#endif
