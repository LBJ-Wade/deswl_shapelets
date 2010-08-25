%module cwl
%{
#include "../src/ConfigFile.h"
#include "../src/StarCatalog.h"
#include "cwl.h"
%}

// must you declare with throw (const char *)?
%typemap(throws) const char * %{
    PyErr_SetString(PyExc_RuntimeError, $1);
    SWIG_fail;
%}

%include "std_string.i"

%include "cwl.h"


