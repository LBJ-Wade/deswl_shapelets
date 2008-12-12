#include "FitsFile.h"


FitsFile::FitsFile(const char* filename, int mode, bool create)
{
  Open(filename, mode, create);
}
FitsFile::FitsFile(std::string filename, int mode, bool create)
{
  Open(filename, mode, create);
}

void FitsFile::Open(std::string filename, int mode, bool create)
{
  Open(filename.c_str(), mode, create);
}
void FitsFile::Open(const char* filename, int mode, bool create)
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
    fits_open_file(&mFptr, filename, mode, &fitserr);
  }
  if (!fitserr==0) {
    fits_report_error(stderr, fitserr); 
    serr="Error opening FITS file: "+mFileName;
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

void 
FitsFile::ReadScalarCol(
    char* colname, 
    int coltype, 
    char* dptr,
    long long nrows)
{
  std::ostringstream err;
  int fitserr=0;

  // always read the first row.
  LONGLONG frow=1;
  // Assume scalar column
  LONGLONG felem=1;

  int colnum;
  fits_get_colnum(mFptr, CASEINSEN, colname, &colnum, &fitserr);
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
      dptr,
      &anynull, &fitserr);
  if (!fitserr==0) {
    fits_report_error(stderr, fitserr); 
    err<<"Error reading column \""<<colname<<"\""<<std::endl;
    throw FitsException(err.str());
  }

}



void 
FitsFile::ReadCell(
    char* colname, 
    int coltype, 
    char* dptr,
    LONGLONG row,
    LONGLONG nel)
{
  std::ostringstream err;
  int fitserr=0;

  // Start at the first element and read nel
  LONGLONG felem=1;

  int colnum;
  fits_get_colnum(mFptr, CASEINSEN, colname, &colnum, &fitserr);
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
      dptr,
      &anynull, &fitserr);
  if (!fitserr==0) {
    fits_report_error(stderr, fitserr); 
    err<<"Error reading column \""<<colname<<"\""<<std::endl;
    throw FitsException(err.str());
  }

}




void FitsFile::ReadKey(char const* name,int dtype,char* dptr)
{
  int fitserr=0;
  std::stringstream err;
  fits_read_key(mFptr, dtype, (char *)name, dptr, NULL, &fitserr);
  if (!fitserr==0)
  {
    fits_report_error(stderr, fitserr); 
    err<<"Error reading header keyword: "<<name
      <<" dtype="<<dtype
      <<std::endl;
    throw FitsException(err.str());
  }

}

long FitsFile::ReadLongKey(char const* name)
{
  long val;

  int fitserr=0;
  std::stringstream err;
  fits_read_key(mFptr, TLONG, (char *)name, &val, NULL, &fitserr);
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



/*
void CreateBinaryTable(
    LONGLONG nrows, 
    vector<string>& names, 
    vector<string>& types,
    vector<string>& units, 
    const char* extension_name)
{
  //BINARY_TBL

}
*/
