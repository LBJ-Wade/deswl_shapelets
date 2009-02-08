
#include <cstring>
#include <sstream>
#include "ConfigFile.h"
#include "FitsFile.h"
#include "dbg.h"


FitsFile::FitsFile(std::string filename, int mode, bool create)
{
  Open(filename, mode, create);
}

void FitsFile::Open(std::string filename, int mode, bool create)
{
  std::string serr;
  int fitserr=0;

  mFileName=filename;
  if (create && READWRITE == mode) {
    //std::cout<<"Creating file: "<<filename<<std::endl;
    std::string f = filename;
    f = "!"+f;
    fits_create_file(&mFptr, f.c_str(), &fitserr);
  } else {
    //std::cout<<"Not Creating file: "<<filename<<" mode,create="<<mode<<","<<create<<std::endl;
    fits_open_file(&mFptr, filename.c_str(), mode, &fitserr);
  }
  if (!fitserr==0) {
    fits_report_error(stderr, fitserr); 
    serr="Error opening FITS file: " + mFileName;
    throw FitsException(serr);
  } else {
    //std::cerr<<"Opened fits FITS file: "<<mFileName<<std::endl;
  }
}

FitsFile::~FitsFile()
{
  Close();
}

void FitsFile::Close()
{
  int fitserr=0;
  fits_close_file(mFptr, &fitserr);
}


/* 
 * this reads the entire column.  Must already be pointed to the correct
 * hdu
 */

void FitsFile::ReadScalarCol(std::string colname, 
    Type_FITS coltype, void* dptr, long nrows)
{
  std::ostringstream err;
  int fitserr=0;

  // always read the first row.
  long frow=1;
  // Assume scalar column
  long felem=1;

  int colnum;
  fits_get_colnum(mFptr, CASEINSEN, (char*) colname.c_str(), &colnum, &fitserr);
  if (!fitserr==0) {
    fits_report_error(stderr, fitserr); 
    err<<"Error reading colnum for \""<<colname<<"\""<<std::endl;
    throw FitsException(err.str());
  }

  // using same dblnull everywhere probably OK.  Just needs to be big
  // enough to take the var I think.

  double dblnull;
  int anynull=0;
  fits_read_col(mFptr, coltype, colnum, frow, felem, nrows, &dblnull, 
      (char*) dptr, &anynull, &fitserr);
  if (!fitserr==0) {
    fits_report_error(stderr, fitserr); 
    err<<"Error reading column \""<<colname<<"\""<<std::endl;
    throw FitsException(err.str());
  }

}

template <class IT1, class IT2, class IT3, class IT4>
static void DoCreateBinaryTable(
    long nrows, long nfields, fitsfile* mFptr,
    IT1 names, IT2 nelem, IT3 types, IT4 units) 
{
  // Create a binary table
  int fits_status=0;
  int tbl_type = BINARY_TBL;

  // Convert to C-style arrays of C-style strings:
  char** cnames = new char*[nfields];
  char** cunits = new char*[nfields];
  for(int i=0;i<nfields;++i, ++names, ++units ) {
    cnames[i] = new char[names->size()+1];
    strcpy(cnames[i],names->c_str());
    cunits[i] = new char[units->size()+1];
    strcpy(cunits[i],units->c_str());
  }
  // Make FITS TFORM string:
  char** cforms = new char*[nfields];
  for(int i=0;i<nfields;++i, ++nelem, ++types ) {
    std::ostringstream form;
    form << *nelem;
    switch (*types) {
      case XSTRING :
	form << "A"; break;
      case XSHORT :
	form << "I"; break;
      case XINT : case XLONG :
	form << "J"; break;
      case XFLOAT :
	form << "E"; break;
      case XDOUBLE :
	form << "D"; break;
      default :
	std::stringstream err;
	err<<"WriteParKey Unsupported type: "<<*types;
	throw std::runtime_error(err.str());
    }
    cforms[i] = new char[form.str().size()+1];
    strcpy(cforms[i],form.str().c_str());
  }
  
  fits_create_tbl(mFptr, tbl_type, nrows, nfields, 
      cnames, cforms, cunits, NULL, &fits_status);

  for(int i=0;i<nfields;++i) {
    delete [] cnames[i];
    delete [] cforms[i];
    delete [] cunits[i];
  }

  delete [] cnames;
  delete [] cforms;
  delete [] cunits;

  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error creating FindStars FITS table: ";
    throw FitsException(serr);
  }
}

void FitsFile::CreateBinaryTable(
    long nrows,
    const std::vector<std::string>& names,
    const std::vector<int>& nelem,
    const std::vector<Type_FITS>& types,
    const std::vector<std::string>& units) 
{
  Assert(nelem.size() == names.size());
  Assert(types.size() == names.size());
  Assert(units.size() == names.size());
  DoCreateBinaryTable(nrows,names.size(),mFptr,
      names.begin(),nelem.begin(),types.begin(),units.begin());
}

void FitsFile::CreateBinaryTable(
    long nrows, int nfields,
    const std::string* names,
    const int* nelem,
    const Type_FITS* types,
    const std::string* units) 
{
  DoCreateBinaryTable(nrows,nfields,mFptr,names,nelem,types,units);
}


void FitsFile::ReadCell(std::string colname, Type_FITS coltype, 
    void* dptr, long row, long nel)
{
  std::ostringstream err;
  int fitserr=0;

  // Start at the first element and read nel
  long felem=1;

  int colnum;
  fits_get_colnum(mFptr, CASEINSEN, (char*) colname.c_str(), 
      &colnum, &fitserr);
  if (!fitserr==0) {
    fits_report_error(stderr, fitserr); 
    err<<"Error reading colnum for \""<<colname<<"\""<<std::endl;
    throw FitsException(err.str());
  }

  // using same dblnull everywhere probably OK.  Just needs to be big
  // enough to take the var I think.

  double dblnull;
  int anynull=0;
  fits_read_col(mFptr, coltype, colnum, row, felem, nel, &dblnull, 
      (char*) dptr, &anynull, &fitserr);
  if (!fitserr==0) {
    fits_report_error(stderr, fitserr); 
    err<<"Error reading column \""<<colname<<"\""<<std::endl;
    throw FitsException(err.str());
  }
}

void FitsFile::ReadKey(std::string name, Type_FITS dtype, void* dptr)
{
  int fitserr=0;
  std::stringstream err;
  fits_read_key(mFptr, dtype, (char *)name.c_str(), (char*)dptr,
      NULL, &fitserr);
  if (!fitserr==0)
  {
    fits_report_error(stderr, fitserr); 
    err<<"Error reading header keyword: "<<name
      <<" dtype="<<dtype
      <<std::endl;
    throw FitsException(err.str());
  }

}


long FitsFile::ReadLongKey(std::string name)
{
  long val;

  int fitserr=0;
  std::stringstream err;
  fits_read_key(mFptr, XLONG, (char *)name.c_str(), &val, NULL, &fitserr);
  if (!fitserr==0)
  {
    fits_report_error(stderr, fitserr); 
    err<<"Error reading header keyword: "<<name<<" fitserr="<<fitserr<<"val="<<val<<std::endl;
    throw FitsException(err.str());
  }

  return val;
}



void FitsFile::GotoHDU(int hdu)
{
  int hdutype;
  int fitserr=0;

  std::stringstream err;

  fits_movabs_hdu(mFptr, hdu, &hdutype, &fitserr);
  if (!fitserr==0)
  {
    fits_report_error(stderr, fitserr); 
    err << "Error moving to extension #" <<hdu<<std::endl;
    throw FitsException(err.str());
  }

}

#if 0
void CreateBinaryTable(
    long nrows, 
    vector<string>& names, 
    vector<string>& types,
    vector<string>& units, 
    std::string extension_name)
{
  //BINARY_TBL

}
#endif

void FitsFile::WriteKey(std::string name, Type_FITS datatype,
    const void* value, std::string comment) 
{
  int fits_status=0;
  fits_write_key(
      mFptr, 
      datatype, 
      (char*)(name.c_str()), 
      (char*)value, 
      (char*)(comment.c_str()),
      &fits_status);

  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::stringstream serr;
    serr<<"Error writing keyword "<<name;
    throw FitsException(serr.str());
  }
}

void FitsFile::WriteColumn(Type_FITS datatype, int colnum,
    long firstrow, long firstel, long nel, const void* data)
{
  int fits_status=0;
  fits_write_col(
      mFptr, datatype, 
      colnum, firstrow, firstel, nel, 
      (char*) data, 
      &fits_status);
  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::stringstream serr;
    serr<<"Error writing to fits column "<<colnum;
    throw FitsException(serr.str());
  }
}
