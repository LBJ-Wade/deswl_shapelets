//---------------------------------------------------------------------------
#ifndef StdH
#define StdH
//---------------------------------------------------------------------------
// 	$Id: Std.h,v 1.11 2003/03/06 03:29:16 garyb Exp $

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <complex>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <new>

namespace laguerre {

    using std::vector;
    using std::string;
    using std::ofstream;
    using std::ifstream;
    using std::ostream;
    using std::istream;
    using std::cout;
    using std::cin;
    using std::cerr;
    using std::endl;

    // Useful conventions:
    template <class T>
    static void SWAP(T& a,T& b) {T temp=a; a=b; b=temp;}

    template <class T>
    static T SQR(const T& x) {return x*x;}

    template <class T>
    static const T& MAX(const T& a, const T& b) {return a>b ? a : b;}

    template <class T>
    static const T& MIN(const T& a, const T& b) {return a<b ? a : b;}

    typedef unsigned int uint;
    typedef unsigned short ushort;
    typedef unsigned long ulong;
    typedef std::complex<double> DComplex;


    // ??? replace here with numeric_limits<double>.infinity()
#ifdef HUGE_VAL
    const double DOUBLE_INFINITY=HUGE_VAL;
#else
    const double DOUBLE_INFINITY=1./0.;
#endif
    const double DOUBLE_NEGATIVE_INFINITY=-DOUBLE_INFINITY;
    //const double DOUBLE_NAN = DOUBLE_NEGATIVE_INFINITY + DOUBLE_INFINITY;
#ifndef PI
#ifdef M_PI
#define PI M_PI
#else
#define PI 3.14159265358979323846
#endif
#endif

    // Debugging and exception classes:

    inline void error(const char *s,const char *s2 = "")
    {
        cerr << "Error: " << s << ' ' << s2 << endl;
        exit(1);
    }

    // Define a base exception class with error message:
    class MyException 
    {
    public:
        MyException(const string &m=""): msg(m) {};
        MyException(const char* c): msg(c) {};
        void set(const string &m="") {msg=m;}
        void set(const char* c) {msg=c;}
        void dump(ostream &os) const {os << "Error: " << msg << endl;}
        void quit(const int exit_code=1) {
            cerr << "Error: " << msg << endl;
            exit(exit_code);
        }
        string msg;
    };

}

#endif
