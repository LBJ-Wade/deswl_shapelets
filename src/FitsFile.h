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

// This encapsulates the fitsio type codes
// I am only including the ones that I think will be used.
// If you need to add more, make sure you add the appropriate
// functionality to the places where this is used.
enum Type_FITS {
  XSTRING = TSTRING,
  XSHORT=TSHORT, 
  XINT = TINT,
  XLONG = TLONG,
  XFLOAT = TFLOAT,
  XDOUBLE = TDOUBLE
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

    void ReadScalarCol(std::string colname, Type_FITS coltype, void* dptr,
	long nrows);
    void ReadCell(std::string colname, Type_FITS coltype, void* dptr,
	long row, long nel);

    void ReadKey(std::string name, Type_FITS dtype, void* dptr);
    void WriteKey(std::string name, Type_FITS datatype, const void* value,
	std::string comment);
    void WriteColumn(Type_FITS datatype, int colnum,
	long firstrow, long firstel, long nel, const void* data);

    // names gives the name of each field.
    // nelem gives the number of elements in each field.
    // types gives the type of data in each field.
    // units gives the unit value for each field.
    void CreateBinaryTable(
	long nrows,
	const std::vector<std::string>& names,
	const std::vector<int>& nelem,
	const std::vector<Type_FITS>& types,
	const std::vector<std::string>& units);
    // The same, but as C-arrays with nfields being the number of fields.
    void CreateBinaryTable(
	long nrows, int nfields,
	const std::string* names,
	const int* nelem,
	const Type_FITS* types,
	const std::string* units);

    long ReadLongKey(std::string name);

    void GotoHDU(int hdu);

    fitsfile* get_fptr() { return mFptr; }

  protected:
    std::string mFileName;
    fitsfile* mFptr;
};

#endif
