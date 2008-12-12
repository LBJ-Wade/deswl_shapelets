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


    /*
    template <class T> void WriteColumn(
	int datatype, 
	int colnum, 
	LONGLONG firstrow, 
	LONGLONG firstel,
	LONGLONG nel,
	T* data);
	*/

    template <class T> void 
      FitsFile::WriteColumn(
	  int datatype, 
	  int colnum, 
	  LONGLONG firstrow, 
	  LONGLONG firstel,
	  LONGLONG nel,
	  T* data)
      {
	int fits_status=0;
	fits_write_col(
	    mFptr, datatype, 
	    colnum, firstrow, firstel, nel, 
	    data, 
	    &fits_status);
	if (!fits_status==0) {
	  fits_report_error(stderr, fits_status); 
	  std::stringstream serr;
	  serr<<"Error writing to fits column "<<colnum;
	  throw FitsException(serr.str());
	}


      }



    void ReadKey(char const* name, int dtype, char* dptr);
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
