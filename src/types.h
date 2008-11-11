#if !defined (_types_h)
#define _types_h

#include <vector>
using std::vector;

// this is not generic enough
typedef char                    int8;
typedef unsigned char           uint8;
typedef short int               int16;
typedef unsigned short int      uint16;
typedef int                     int32;
typedef unsigned int            uint32;
typedef float                   float32;
typedef double                  float64;
#ifdef _WIN32
typedef __int64                 int64;
typedef unsigned __int64        uint64;
#else
typedef long long               int64;
typedef unsigned long long      uint64;
#endif

#endif /* _types_h */
