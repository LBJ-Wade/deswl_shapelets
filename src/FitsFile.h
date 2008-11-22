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
    FitsFile(const char* filename);
    FitsFile(std::string filename);
    ~FitsFile();

    void Open(const char* filename);
    void Open(std::string filename);
    void Close();
    // Read a binary fits column

    void ReadScalarCol(
	char* colname, 
	int coltype, 
	char* dptr,
	long long nrows);

    void ReadKey(char const* name, int dtype, char* dptr);
    long ReadLongKey(char const* name);


    void GotoHDU(int hdu);

  protected:
    std::string mFileName;
    fitsfile* mFptr;
};

#endif
