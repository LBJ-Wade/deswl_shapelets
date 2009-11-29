#ifndef TESTHELPER_H
#define TESTHELPER_H

#include <stdexcept>

extern bool shouldShowTests;
extern bool shouldThrow;
extern std::string lastSuccess;

inline void preTest(std::string s)
{
    if (shouldShowTests) {
        std::cout<<"Trying: "<<s;
        std::cout.flush();
    }
}

inline void test(bool x, std::string s)
{
    if (x) {
        if (shouldShowTests) std::cout<<"  Passed"<<std::endl;
        lastSuccess = s;
    } else {
        if (shouldShowTests) std::cout<<"  Failed"<<std::endl;
        if (!shouldThrow) std::cout<<"Failed test: "<<s<<std::endl;
        else  {
#ifdef NOTHROW
            std::cout<<"Error in test: "<<s<<std::endl;
            exit(1); 
#else
            throw std::runtime_error(std::string("Error in test: ") + s);
#endif
        }
    }
}

#define runTest(x,s) \
    do { \
        preTest(s); \
        test(x,s); \
    } while (false)

#endif
