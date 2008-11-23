
#include "Log.h"
#include <iostream>
#include <stdexcept>

ShearLog::ShearLog(std::string logfile, std::string _delim) :
  exitcode(SUCCESS), ngals(0),
  nf_range1(0), nf_range2(0),
  nf_edge1(0), nf_npix1(0), nf_small(0),
  nf_edge2(0), nf_npix2(0),
  nf_tmverror(0), nf_othererror(0),
  ns_native(0), nf_native(0),
  ns_mu(0), nf_mu(0),
  ns_gamma(0), nf_gamma(0),
  delim(_delim)
{
  if (logfile == "") {
    logout = &std::cout;
  } else {
    logout = new std::ofstream(logfile.c_str(),std::ios_base::app);
    if (!logout) {
      throw std::runtime_error(
	  std::string("Error: Unable to open logfile ") + logfile
	  + " for append output");
    }
  }
}

ShearLog::~ShearLog() 
{
  WriteLog();
  if (logout != &std::cout) {
    delete logout;
  }
}

void ShearLog::WriteLog() const
{
  *logout << 
    exitcode <<delim<<
    ngals <<delim<<
    ns_gamma <<delim<<
    ns_native <<delim<<
    nf_range1 <<delim<<
    nf_range2 <<delim<<
    nf_edge1 <<delim<<
    nf_npix1 <<delim<<
    nf_native <<delim<<
    nf_small <<delim<<
    nf_edge2 <<delim<<
    nf_npix2 <<delim<<
    nf_tmverror <<delim<<
    nf_othererror <<delim<<
    nf_mu <<delim<<
    nf_gamma <<std::endl;
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

PSFLog::PSFLog(std::string logfile, std::string _delim) :
  exitcode(SUCCESS), nstars(0),
  nf_range(0), nf_edge(0), nf_npix(0),
  nf_tmverror(0), nf_othererror(0),
  ns_psf(0), nf_psf(0),
  delim(_delim)
{
  if (logfile == "") {
    logout = &std::cout;
  } else {
    logout = new std::ofstream(logfile.c_str(),std::ios_base::app);
    if (!logout) {
      throw std::runtime_error(
	  std::string("Error: Unable to open logfile ") + logfile
	  + " for append output");
    }
  }
}

PSFLog::~PSFLog() 
{ 
  WriteLog();
  if (logout != &std::cout) {
    delete logout;
  }
}

void PSFLog::WriteLog() const
{
  *logout << 
    exitcode <<delim<<
    nstars <<delim<<
    ns_psf <<delim<<
    nf_range <<delim<<
    nf_edge <<delim<<
    nf_npix <<delim<<
    nf_tmverror <<delim<<
    nf_othererror <<delim<<
    nf_psf <<std::endl;
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

FindStarsLog::FindStarsLog(std::string logfile, std::string _delim) :
  exitcode(SUCCESS), nobj(0), nstars(0),
  delim(_delim)
{
  if (logfile == "") {
    logout = &std::cout;
  } else {
    logout = new std::ofstream(logfile.c_str(),std::ios_base::app);
    if (!logout) {
      throw std::runtime_error(
	  std::string("Error: Unable to open logfile ") + logfile
	  + " for append output");
    }
  }
}

FindStarsLog::~FindStarsLog() 
{ 
  WriteLog();
  if (logout != &std::cout) {
    delete logout;
  }
}

void FindStarsLog::WriteLog() const
{
  *logout << 
    exitcode <<delim<<
    nobj <<delim<<
    nstars <<std::endl;
}

void FindStarsLog::Write(std::ostream& os) const
{
  os<<"N_Input objects = "<<nobj<<std::endl;
  os<<"N_Stars = "<<nstars<<std::endl;
}
