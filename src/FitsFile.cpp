
#include "ConfigFile.h"
#include "FitsFile.h"
#include "dbg.h"
#include <cstring>


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


// only if we have config file available
#ifdef ConfigFile_H
void FitsFile::WriteParKey(const ConfigFile& params, std::string cname,
    int type)
{

  std::string key = cname;
  std::string hname_key = key + "_hname";
  std::string comment_key = key + "_comment";

  std::string hname=params.get(hname_key);
  std::string comment=params.get(comment_key);

  switch (type) {
    case TSTRING:
      {
	std::string val=params.get(key);
	WriteKey(hname, type, val.c_str(), comment);
      }
      break;
    case TDOUBLE:
      {
	double val=params.get(key);
	WriteKey(hname, type, &val, comment);
      }
      break;
    case TINT:
      {
	int val=params.get(key);
	WriteKey(hname, type, &val, comment);
      }
    case TLONG:
      {
	long val=params.get(key);
	WriteKey(hname, type, &val, comment);
      }
      break;
    default:
      std::stringstream err;
      err<<"WriteParKey Unsupported type: "<<type;
      throw std::runtime_error(err.str());
  }

}
#endif


/* 
 * this reads the entire column.  Must already be pointed to the correct
 * hdu
 */

void FitsFile::ReadScalarCol(std::string colname, 
    int coltype, void* dptr, LONGLONG nrows)
{
  std::ostringstream err;
  int fitserr=0;

  // always read the first row.
  LONGLONG frow=1;
  // Assume scalar column
  LONGLONG felem=1;

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

void FitsFile::CreateBinaryTable(
    LONGLONG nrows,
    const std::vector<std::string>& names,
    const std::vector<std::string>& types,
    const std::vector<std::string>& units) 
{
  // Create a binary table
  int fits_status=0;
  int tbl_type = BINARY_TBL;
  Assert(types.size() == names.size());
  Assert(units.size() == names.size());
  int nfields = names.size();

  // Convert to C-style arrays of C-style strings:
  char** cnames = new char*[nfields];
  char** ctypes = new char*[nfields];
  char** cunits = new char*[nfields];
  for(int i=0;i<nfields;++i) {
    cnames[i] = new char[names[i].size()+1];
    strcpy(cnames[i],names[i].c_str());
    ctypes[i] = new char[types[i].size()+1];
    strcpy(ctypes[i],types[i].c_str());
    cunits[i] = new char[units[i].size()+1];
    strcpy(cunits[i],units[i].c_str());
  }
  
  fits_create_tbl(mFptr, tbl_type, nrows, nfields, 
      cnames, ctypes, cunits, NULL, &fits_status);

  for(int i=0;i<nfields;++i) {
    delete [] cnames[i];
    delete [] ctypes[i];
    delete [] cunits[i];
  }

  delete [] cnames;
  delete [] ctypes;
  delete [] cunits;

  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error creating FindStars FITS table: ";
    throw FitsException(serr);
  }
}

void FitsFile::CreateBinaryTable(
    LONGLONG nrows, int nfields,
    const std::string* names,
    const std::string* types,
    const std::string* units) 
{
  // Create a binary table
  int fits_status=0;
  int tbl_type = BINARY_TBL;

  // Convert to C-style arrays of C-style strings:
  char** cnames = new char*[nfields];
  char** ctypes = new char*[nfields];
  char** cunits = new char*[nfields];
  for(int i=0;i<nfields;++i) {
    cnames[i] = new char[names[i].size()+1];
    strcpy(cnames[i],names[i].c_str());
    ctypes[i] = new char[types[i].size()+1];
    strcpy(ctypes[i],types[i].c_str());
    cunits[i] = new char[units[i].size()+1];
    strcpy(cunits[i],units[i].c_str());
  }
  
  fits_create_tbl(mFptr, tbl_type, nrows, nfields, 
      cnames, ctypes, cunits, NULL, &fits_status);

  for(int i=0;i<nfields;++i) {
    delete [] cnames[i];
    delete [] ctypes[i];
    delete [] cunits[i];
  }

  delete [] cnames;
  delete [] ctypes;
  delete [] cunits;

  if (!fits_status==0) {
    fits_report_error(stderr, fits_status); 
    std::string serr="Error creating FindStars FITS table: ";
    throw FitsException(serr);
  }
}


void FitsFile::ReadCell(std::string colname, int coltype, 
    void* dptr, LONGLONG row, LONGLONG nel)
{
  std::ostringstream err;
  int fitserr=0;

  // Start at the first element and read nel
  LONGLONG felem=1;

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

void FitsFile::ReadKey(std::string name, int dtype, void* dptr)
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
  fits_read_key(mFptr, TLONG, (char *)name.c_str(), &val, NULL, &fitserr);
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
    LONGLONG nrows, 
    vector<string>& names, 
    vector<string>& types,
    vector<string>& units, 
    std::string extension_name)
{
  //BINARY_TBL

}
#endif

void FitsFile::WriteKey(std::string name, int datatype, const void* value,
    std::string comment) 
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

void FitsFile::WriteColumn(int datatype, int colnum,
    LONGLONG firstrow, LONGLONG firstel, LONGLONG nel, const void* data)
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
