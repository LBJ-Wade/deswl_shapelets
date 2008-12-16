/*
 * CCFits is a better more complete package, but this is small and does
 * what we need
 */

#ifndef _FITS_FILE_H
#define _FITS_FILE_H

#include <stdexcept>
#include <string>
#include <sstream>
#include <iostream>
#include "fitsio.h"

#include <vector>
#include <cstring>

// in case we don't have the ConfigFile class available
#define HAVE_CONFIG_FILE

#ifdef HAVE_CONFIG_FILE
#include "ConfigFile.h"
#endif

class FitsException : public std::runtime_error {
  public:
    FitsException(const char* m) : std::runtime_error(m) {};
    FitsException(std::string m) : std::runtime_error(m.c_str()) {};
    FitsException(char* m) : std::runtime_error((const char*) m) {};
};


class FitsFile
{
  public:
    FitsFile() {};
    FitsFile(const char* filename, int mode=READONLY, bool create=false);
    FitsFile(std::string filename, int mode=READONLY, bool create=false);
    ~FitsFile();

    void Open(const char* filename, int mode=READONLY, bool create=false);
    void Open(std::string filename, int mode=READONLY, bool create=false);
    void Close();
    // Read a binary fits column

    void ReadScalarCol(
	char* colname, 
	int coltype, 
	char* dptr,
	long long nrows);
    void ReadCell(
	char* colname, 
	int coltype, 
	char* dptr,
	LONGLONG row,
	LONGLONG nel);

#ifdef HAVE_CONFIG_FILE
    void WriteParKey(const ConfigFile& params, const char* name, int type);
#endif

    template <class T> void 
      WriteColumn(
	  int datatype, 
	  int colnum, 
	  LONGLONG firstrow, 
	  LONGLONG firstel,
	  LONGLONG nel,
	  const T* data)
      {
	int fits_status=0;
	fits_write_col(
	    mFptr, datatype, 
	    colnum, firstrow, firstel, nel, 
	    (void*)data, 
	    &fits_status);
	if (!fits_status==0) {
	  fits_report_error(stderr, fits_status); 
	  std::stringstream serr;
	  serr<<"Error writing to fits column "<<colnum;
	  throw FitsException(serr.str());
	}


      }



    void ReadKey(char const* name, int dtype, char* dptr);

    void WriteKey(
	const char* name,
	int datatype,
	const char* cstr, 
	const char* comment) {

      int fits_status=0;
      fits_write_key(
	  mFptr, 
	  datatype, 
	  (char*)name, 
	  (char*)cstr, 
	  (char*)comment,
	  &fits_status);

      if (!fits_status==0) {
	fits_report_error(stderr, fits_status); 
	std::stringstream serr;
	serr<<"Error writing keyword "<<name;
	throw FitsException(serr.str());
      }

    }

    template <class T> void WriteKey(
	const char* name,
	int datatype,
	T& val, 
	const char* comment) {

      int fits_status=0;
      fits_write_key(
	  mFptr, 
	  datatype, 
	  (char*)name, 
	  &val, 
	  (char*)comment,
	  &fits_status);

      if (!fits_status==0) {
	fits_report_error(stderr, fits_status); 
	std::stringstream serr;
	serr<<"Error writing keyword "<<name;
	throw FitsException(serr.str());
      }


    }


    long ReadLongKey(char const* name);

    void GotoHDU(int hdu);

    fitsfile* get_fptr() { 
      return mFptr;
    }

  protected:
    std::string mFileName;
    fitsfile* mFptr;
};

#endif
