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
#include <vector>
#include "fitsio.h"

class FitsException : public std::runtime_error {
  public:
    FitsException(std::string m) : std::runtime_error(m) {};
};

class FitsFile
{
  public:
    FitsFile() {};
    FitsFile(std::string filename, int mode=READONLY, bool create=false);
    ~FitsFile();

    void Open(std::string filename, int mode=READONLY, bool create=false);
    void Close();
    // Read a binary fits column

    void ReadScalarCol(std::string colname, int coltype, void* dptr,
	LONGLONG nrows);
    void ReadCell(std::string colname, int coltype, void* dptr,
	LONGLONG row, LONGLONG nel);

#ifdef ConfigFile_H
    void WriteParKey(const ConfigFile& params, std::string name, int type);
#endif

    void ReadKey(std::string name, int dtype, void* dptr);
    void WriteKey(std::string name, int datatype, const void* value,
	std::string comment);
    void WriteColumn(int datatype, int colnum,
	LONGLONG firstrow, LONGLONG firstel, LONGLONG nel, const void* data);

    void CreateBinaryTable(
	LONGLONG nrows,
	const std::vector<std::string>& names,
	const std::vector<std::string>& types,
	const std::vector<std::string>& units);
    void CreateBinaryTable(
	LONGLONG nrows, int nfields,
	const std::string* names,
	const std::string* types,
	const std::string* units);

    long ReadLongKey(std::string name);

    void GotoHDU(int hdu);

    fitsfile* get_fptr() { return mFptr; }

  protected:
    std::string mFileName;
    fitsfile* mFptr;
};

#endif
