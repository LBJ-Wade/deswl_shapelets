#ifndef TESTHELPER_H
#define TESTHELPER_H

#include <stdexcept>

extern bool shouldShowTests;
extern bool shouldThrow;
extern std::string lastSuccess;
extern std::ostream* testout;

inline void preTest(std::string s)
{
    if (shouldShowTests) {
        (*testout)<<"Trying: "<<s;
        (*testout).flush();
    }
}

inline void doTest(bool x, std::string s)
{
    if (x) {
        if (shouldShowTests) (*testout)<<"  Passed"<<std::endl;
        lastSuccess = s;
    } else {
        if (shouldShowTests) (*testout)<<"  Failed"<<std::endl;
        if (!shouldThrow) (*testout)<<"Failed test: "<<s<<std::endl;
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

#define test(x,s) \
    do { \
        preTest(s); \
        doTest(x,s); \
    } while (false)

#endif
