
#include "Log.h"
#include <iostream>
#include <stdexcept>
#include "Params.h"

ShearLog::ShearLog(std::string _logfile, std::string _fits_file) :
  exitcode(SUCCESS), ngals(0),
  nf_range1(0), nf_range2(0),
  nf_edge1(0), nf_npix1(0), nf_small(0),
  nf_edge2(0), nf_npix2(0),
  nf_tmverror(0), nf_othererror(0),
  ns_native(0), nf_native(0),
  ns_mu(0), nf_mu(0),
  ns_gamma(0), nf_gamma(0),
  fits_file(_fits_file)
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
  WriteLogToFitsHeader();
  NoWriteLog(); // deletes logout
}

void ShearLog::NoWriteLog() 
{ 
  if (logout && logout != &std::cout) {
    delete logout;
  }
  logout = 0;
}

void ShearLog::WriteLogToFitsHeader() const
{

  if ("" == fits_file) return;

  try {

    //std::cout<<"Writing log to Shear file: "<<fits_file<<std::endl;
    bool create=false;
    FitsFile fits(fits_file.c_str(), READWRITE, create);
    fits.GotoHDU(2);

    // copying out to avoid g++ errors due to const incorrectness
    long val=0;

    // use msexit to differentiate from exit codes for other codes
    val=exitcode;
    fits.WriteKey("msexit", TLONG, val, "Exit code for MeasureShear");

    val=ngals;
    fits.WriteKey("msnobj", TLONG, val, 
	"# of objects processed by MeasureShear");
    val=ns_gamma;
    fits.WriteKey("ns_gamma", TLONG, val, 
	"# of successful shear measurements");
    val=ns_native;
    fits.WriteKey("ns_nativ", TLONG, val, 
	"# of successful native shapeles measurements");

    val=nf_range1;
    fits.WriteKey("nf_rnge1", TLONG, val, 
	"# of MeasureShear failures range1");
    val=nf_range2;
    fits.WriteKey("nf_rnge2", TLONG, val, 
	"# of MeasureShear failures range2");

    val=nf_edge1;
    fits.WriteKey("nf_edge1", TLONG, val, 
	"# of MeasureShear failures edge1");
    val=nf_edge2;
    fits.WriteKey("nf_edge2", TLONG, val, 
	"# of MeasureShear failures edge2");

    val=nf_npix1;
    fits.WriteKey("nf_npix1", TLONG, val, 
	"# of MeasureShear failures pixel extraction1");
    val=nf_npix2;
    fits.WriteKey("nf_npix2", TLONG, val, 
	"# of MeasureShear failures pixel extraction2");

    val=nf_native;
    fits.WriteKey("nf_nativ", TLONG, val, 
	"# of MeasureShear failures native calculations");
    val=nf_small;
    fits.WriteKey("nf_small", TLONG, val, 
	"# of MeasureShear failures too small");

    // prefix ms to make unique
    val=nf_tmverror;
    fits.WriteKey("nf_mstmv", TLONG, val, 
	"# of MeasureShear failures TMV errors");
    val=nf_othererror;
    fits.WriteKey("nf_msoth", TLONG, val, 
	"# of MeasureShear failures unclassified errors");

    val=nf_mu;
    fits.WriteKey("nf_mu", TLONG, val, 
	"# of MeasureShear failures calculating mu");
    val=nf_gamma;
    fits.WriteKey("nf_gamma", TLONG, val, 
	"# of MeasureShear failures calculating shear");
  } 
  catch(...) 
  { if (exitcode == 0) throw; }

}

inline std::string StatusText1(ExitCode exitcode)
{ 
  return std::string("STATUS") +
    char('0'+Status(exitcode)) +
    std::string("BEG"); 
}
inline std::string StatusText2(ExitCode exitcode)
{
  return std::string("STATUS") +
    char('0'+Status(exitcode)) +
    std::string("END");
}

void ShearLog::WriteLog() const
{
  if (logout) {
    // Emit logging information
    if (exitcode) {
      *logout << 
	StatusText1(exitcode)<<" "<<
	Text(exitcode)<<" "<<
	extraexitinfo<<" ";
      if (fits_file != "")
	*logout << " (Name="<<fits_file<<") ";
      *logout<< StatusText2(exitcode)<<std::endl;
    } else {
      std::string name = "measureshear";
      if (fits_file != "") name = fits_file;
      *logout << 
	"QA2BEG "<<
	"Name="<<name <<" & "<<
	"ngals="<<ngals <<" & "<<
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
	"nf_gamma="<<nf_gamma <<
	" QA2END"<<std::endl;
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

PSFLog::PSFLog(std::string _logfile, std::string _fits_file) :
  exitcode(SUCCESS), nstars(0),
  nf_range(0), nf_edge(0), nf_npix(0),
  nf_tmverror(0), nf_othererror(0),
  ns_psf(0), nf_psf(0),
  fits_file(_fits_file)
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
  WriteLogToFitsHeader();
  NoWriteLog();
}

void PSFLog::NoWriteLog()
{
  if (logout && logout != &std::cout) {
    delete logout;
  }
  logout = 0;
}

void PSFLog::WriteLogToFitsHeader() const
{
  if ("" == fits_file) return;

  try {
    bool create=false;
    //std::cout<<"Writing log to PSF file: "<<fits_file<<std::endl;
    FitsFile fits(fits_file.c_str(), READWRITE, create);
    fits.GotoHDU(2);

    // copying out to avoid g++ errors due to const incorrectness
    long val=0;

    val=exitcode;
    // use fsexit to differentiate from exit codes for other codes
    fits.WriteKey("mpexit", TLONG, val, "Exit code for MeasurePSF");
    val=nstars;
    fits.WriteKey("mpnstars", TLONG, val, 
	"# of PSF star candidates used");
    val=ns_psf;
    fits.WriteKey("ns_psf", TLONG, val, 
	"# of successful PSF decompositions");
    val=nf_range;
    fits.WriteKey("nf_range", TLONG, val, 
	"# PSF failures due to range error");
    val=nf_edge;
    fits.WriteKey("nf_edge", TLONG, val, 
	"# PSF failures due to edge error");
    val=nf_npix;
    fits.WriteKey("nf_npix", TLONG, val, 
	"# PSF failures due bad pix extraction");
    val=nf_tmverror;
    fits.WriteKey("nf_mptmv", TLONG, val, 
	"# PSF failures due to TMV errors");
    val=nf_othererror;
    fits.WriteKey("nf_mpoth", TLONG, val, 
	"# PSF failures due to unclassified errors");
    val=nf_psf;
    fits.WriteKey("nf_psf", TLONG, val, 
	"# PSF failures");

    fits.Close();
  }
  catch (...)
  { if (exitcode == 0) throw; }
}

void PSFLog::WriteLog() const
{
  if (logout) {
    // Emit logging information
    if (exitcode) {
      *logout << 
	StatusText1(exitcode)<<" "<<
	Text(exitcode)<<" "<<
	extraexitinfo<<" ";
      if (fits_file != "")
	*logout << " (Name="<<fits_file<<") ";
      *logout<< StatusText2(exitcode)<<std::endl;
    } else {
      std::string name = "measurepsf";
      if (fits_file != "") name = fits_file;
      *logout << 
	"QA2BEG "<<
	"Name="<<name <<" & "<<
	"nstars="<<nstars <<" & "<<
	"ns_psf="<<ns_psf <<" & "<<
	"nf_range="<<nf_range <<" & "<<
	"nf_edge="<<nf_edge <<" & "<<
	"nf_npix="<<nf_npix <<" & "<<
	"nf_tmverror="<<nf_tmverror <<" & "<<
	"nf_othererror="<<nf_othererror <<" & "<<
	"nf_psf="<<nf_psf <<
	" QA2END"<<std::endl;
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


FindStarsLog::FindStarsLog(std::string _logfile, std::string _fits_file) :
  exitcode(SUCCESS), nobj(0), nstars(0),
  fits_file(_fits_file)
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
  WriteLogToFitsHeader();
  NoWriteLog();
}

void FindStarsLog::NoWriteLog()
{
  if (logout != &std::cout) {
    delete logout;
  }
  logout = 0;
}

void FindStarsLog::WriteLogToFitsHeader() const
{
  if ("" == fits_file) return;

  try {
    bool create=false;
    //std::cout<<"Writing log to FindStars file: "<<fits_file<<std::endl;
    FitsFile fits(fits_file.c_str(), READWRITE, create);
    fits.GotoHDU(2);

    // not efficient since entire file must be resized and all data moved
    // forward each time

    // copying out to avoid g++ errors due to const incorrectness
    long val;
    val=exitcode;
    fits.WriteKey("fsexit", TLONG, val, "Exit code for FindStars");
    val=nobj;
    fits.WriteKey("fsnobj", TLONG, val, 
	"# of objects processed in FindStars");
    // Using fsnstar to avoid name collision with same parameter in the
    // PSF log file
    val=nstars;
    fits.WriteKey("fsnstars", TLONG, val, 
	"# of PSF candidate stars found by FindStars");

    fits.Close();
  }
  catch(...)
  { if (exitcode == 0) throw; }
}

void FindStarsLog::WriteLog() const
{
  if (logout) {
    // Emit logging information
    if (exitcode) {
      *logout << 
	StatusText1(exitcode)<<" "<<
	Text(exitcode)<<" "<<
	extraexitinfo<<" ";
      if (fits_file != "")
	*logout << " (Name="<<fits_file<<") ";
      *logout<< StatusText2(exitcode)<<std::endl;
    } else {
      std::string name = "findstars";
      if (fits_file != "") name = fits_file;
      *logout << 
	"QA2BEG "<<
	"Name="<<name <<" & "<<
	"nobj="<<nobj <<" & "<<
	"nstars="<<nstars <<
	" QA2END"<<std::endl;
    }
  }
}

void FindStarsLog::Write(std::ostream& os) const
{
  os<<"N_Input objects = "<<nobj<<std::endl;
  os<<"N_Stars = "<<nstars<<std::endl;
}
