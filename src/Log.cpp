
#include "Log.h"
#include <iostream>
#include <stdexcept>

#include "Params.h"


ShearLog::ShearLog(std::string _logfile) :
  exitcode(SUCCESS), ngals(0),
  nf_range1(0), nf_range2(0),
  nf_edge1(0), nf_npix1(0), nf_small(0),
  nf_edge2(0), nf_npix2(0),
  nf_tmverror(0), nf_othererror(0),
  ns_native(0), nf_native(0),
  ns_mu(0), nf_mu(0),
  ns_gamma(0), nf_gamma(0)
{
  if (_logfile == "") {
    logout = &std::cout;
  } else {
    logout = new std::ofstream(_logfile.c_str(),std::ios_base::app);
    if (!logout) {
      throw std::runtime_error(
	  std::string("Error: Unable to open logfile ") + _logfile
	  + " for append output");
    }
  }
}

ShearLog::~ShearLog() 
{
  WriteLog();
  NoWriteLog(); // deletes logout
}

void ShearLog::NoWriteLog() 
{ 
  if (logout && logout != &std::cout) {
    delete logout;
  }
}

void ShearLog::WriteLogToFitsHeader(FitsFile& fits) const
{

  // Write to Header:
  // exitcode
  // ngals
  // ns_gamma
  // ns_nativ = ns_native
  // nf_rnge1 = nf_range1
  // nf_rnge2 = nf_range2
  // nf_edge1
  // nf_npix1
  // nf_nativ = nf_native
  // nf_small
  // nf_edge2
  // nf_npix2
  // nf_tmver = nf_tmverror
  // nf_other = nf_othererror
  // nf_mu
  // nf_gamma
}

void ShearLog::WriteLog() const
{
  if (logout) {
    // Emit logging information
    if (exitcode) {
      *logout << 
	"STATUS5BEG "<<
	Text(exitcode)<<" "<<
	extraexitinfo<<" "<<
	"STATUS5END"<<std::endl;
    } else {
      *logout << 
	"STATUS2BEG "<<
	"Name=measureshear & "<<
	"ngals="<<ngals<<" & "<<
	"ns_gamma="<<ns_gamma <<" & "<<
	"ns_native="<<ns_native <<" & "<<
	"nf_range1="<<nf_range1 <<" & "<<
	"nf_range2="<<nf_range2 <<" & "<<
	"nf_edge1="<<nf_edge1 <<" & "<<
	"nf_npix1="<<nf_npix1 <<" & "<<
	"nf_native="<<nf_native <<" & "<<
	"nf_small="<<nf_small <<" & "<<
	"nf_edge2="<<nf_edge2 <<" & "<<
	"nf_npix2="<<nf_npix2 <<" & "<<
	"nf_tmverror="<<nf_tmverror <<" & "<<
	"nf_othererror="<<nf_othererror <<" & "<<
	"nf_mu="<<nf_mu <<" & "<<
	"nf_gamma="<<nf_gamma <<" & "<<
	"STATUS2END"<<std::endl;
    }
  }
}

void ShearLog::Write(std::ostream& os) const
{
  os<<"N_Rejected for range error - distortion = "<<nf_range1<<std::endl;
  os<<"N_Rejected for range error - psf interp = "<<nf_range2<<std::endl;
  os<<"N_Rejected for pixel flag (#1) = "<<nf_edge1<<std::endl;
  os<<"N_Rejected for too few pixels (#1) = "<<nf_npix1<<std::endl;

  os<<"Native fits:\n";
  os<<"  N_Success = "<<ns_native<<std::endl;
  os<<"  N_Fail = "<<nf_native<<std::endl;

  os<<"N_Rejected for being too small = "<<nf_small<<std::endl;
  os<<"N_Rejected for pixel flag (#2) = "<<nf_edge2<<std::endl;
  os<<"N_Rejected for too few pixels (#2) = "<<nf_npix2<<std::endl;

  os<<"Mu fits:\n";
  os<<"  N_Success = "<<ns_mu<<std::endl;
  os<<"  N_Fail = "<<nf_mu<<std::endl;

  os<<"Gamma fits:\n";
  os<<"  N_Success = "<<ns_gamma<<std::endl;
  os<<"  N_Fail = "<<nf_gamma<<std::endl;

  os<<"N_Error: TMV Error caught = "<<nf_tmverror<<std::endl;
  os<<"N_Error: Other caught = "<<nf_othererror<<std::endl;
}

PSFLog::PSFLog(std::string _logfile) :
  exitcode(SUCCESS), nstars(0),
  nf_range(0), nf_edge(0), nf_npix(0),
  nf_tmverror(0), nf_othererror(0),
  ns_psf(0), nf_psf(0)
{
  if (_logfile == "") {
    logout = &std::cout;
  } else {
    logout = new std::ofstream(_logfile.c_str(),std::ios_base::app);
    if (!logout) {
      throw std::runtime_error(
	  std::string("Error: Unable to open logfile ") + _logfile
	  + " for append output");
    }
  }
}

PSFLog::~PSFLog() 
{ 
  WriteLog();
  NoWriteLog();
}

void PSFLog::NoWriteLog()
{
  if (logout && logout != &std::cout) {
    delete logout;
  }
}

void PSFLog::WriteLogToFitsHeader(FitsFile& fits) const
{

  // Write to Header:
  // exitcode
  // nstars
  // ns_psf
  // nf_range
  // nf_edge
  // nf_npix
  // nf_tmver = nf_tmverror
  // nf_other = nf_otherror
  // nf_psf
}

void PSFLog::WriteLog() const
{
  if (logout) {
    // Emit logging information
    if (exitcode) {
      *logout << 
	"STATUS5BEG "<<
	Text(exitcode)<<" "<<
	extraexitinfo<<" "<<
	"STATUS5END"<<std::endl;
    } else {
      *logout << 
	"STATUS2BEG "<<
	"Name=measurepsf & "<<
	"nstars="<<nstars <<" & "<<
	"ns_psf="<<ns_psf <<" & "<<
	"nf_range="<<nf_range <<" & "<<
	"nf_edge="<<nf_edge <<" & "<<
	"nf_npix="<<nf_npix <<" & "<<
	"nf_tmverror="<<nf_tmverror <<" & "<<
	"nf_othererror="<<nf_othererror <<" & "<<
	"nf_psf="<<nf_psf <<" & "<<
	"STATUS2END"<<std::endl;
    }
  }
}

void PSFLog::Write(std::ostream& os) const
{
  os<<"N_Rejected for range error - distortion = "<<nf_range<<std::endl;
  os<<"N_Rejected for edge error = "<<nf_edge<<std::endl;
  os<<"N_Rejected for too few pixels = "<<nf_npix<<std::endl;

  os<<"PSF measurement:\n";
  os<<"  N_Success = "<<ns_psf<<std::endl;
  os<<"  N_Fail = "<<nf_psf<<std::endl;

  os<<"N_Error: TMV Error caught = "<<nf_tmverror<<std::endl;
  os<<"N_Error: Other caught = "<<nf_othererror<<std::endl;
}

FindStarsLog::FindStarsLog(std::string _logfile) :
  exitcode(SUCCESS), nobj(0), nstars(0)
{
  if (_logfile == "") {
    logout = &std::cout;
  } else {
    logout = new std::ofstream(_logfile.c_str(),std::ios_base::app);
    if (!logout) {
      throw std::runtime_error(
	  std::string("Error: Unable to open logfile ") + _logfile
	  + " for append output");
    }
  }
}

FindStarsLog::~FindStarsLog() 
{ 
  WriteLog();
  NoWriteLog();
}

void FindStarsLog::NoWriteLog()
{
  if (logout != &std::cout) {
    delete logout;
  }
}

void FindStarsLog::WriteLogToFitsHeader(FitsFile& fits) const
{
  // Write to Header:
  // exitcode
  // nobj
  // nstars
}

void FindStarsLog::WriteLog() const
{
  if (logout) {
    // Emit logging information
    if (exitcode) {
      *logout << 
	"STATUS5BEG "<<
	Text(exitcode)<<" "<<
	extraexitinfo<<" "<<
	"STATUS5END"<<std::endl;
    } else {
      *logout << 
	"STATUS2BEG "<<
	"Name=findstars & "<<
	"nobj="<<nobj <<" & "<<
	"nobj="<<nstars <<" & "<<
	"STATUS2END"<<std::endl;
    }
  }
}

void FindStarsLog::Write(std::ostream& os) const
{
  os<<"N_Input objects = "<<nobj<<std::endl;
  os<<"N_Stars = "<<nstars<<std::endl;
}
