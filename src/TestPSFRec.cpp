#include <iostream>

#include "TMV.h"
#include "TMV_Small.h"

#include "dbg.h"

#include "PSFCatalog.h"


std::ostream* dbgout = 0;
bool XDEBUG = true;

#define xdbgout (XDEBUG ? dbgout : 0)

int main(int argc, char **argv) try
{

	std::cout<<"Hello World\n";
}
catch (...) {
  std::cerr<<"Caught Unknown error\n";
}
