
#include "Log.h"
#include <iostream>
#include <stdexcept>
#include "Params.h"

Log::Log(std::string _logfile, std::string _fits_file) :
  exitcode(SUCCESS), fits_file(_fits_file)
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

Log::~Log() 
{
  // I thought this would work, but it doesn't. The reason is a bit subtle.  
  // Basically, the destructors get called in order from most derived 
  // up to the base.
  // So by this point, the derived portions have already been deleted,
  // so the derived implementations of these functions are already gone.
  // The solution is to put this same thing in each derived destructor.
#if 0
  this->WriteLog();
  this->WriteLogToFitsHeader();
#endif
  NoWriteLog(); // deletes logout
}

void Log::NoWriteLog() 
{ 
  if (logout && logout != &std::cout) {
    delete logout;
  }
  logout = 0;
}

ShearLog::ShearLog(std::string _logfile, std::string _fits_file) :
  Log(_logfile,_fits_file),
  ngals(0), nf_range1(0), nf_range2(0), nf_small(0),
  nf_tmverror(0), nf_othererror(0), ns_native(0), nf_native(0),
  ns_mu(0), nf_mu(0), ns_gamma(0), nf_gamma(0)
{}

ShearLog::~ShearLog() {
  this->WriteLog();
  this->WriteLogToFitsHeader();
}

void ShearLog::WriteLogToFitsHeader() const
{

  if ("" == fits_file) return;

  try {
    //std::cout<<"Writing log to Shear file: "<<fits_file<<std::endl;
    bool create=false;
    FitsFile fits(fits_file.c_str(), READWRITE, create);
    fits.GotoHDU(2);

    // use msexit to differentiate from exit codes for other codes
    fits.WriteKey("msexit", TLONG, &exitcode, "Exit code for MeasureShear");
    fits.WriteKey("msnobj", TLONG, &ngals, 
	"# of objects processed by MeasureShear");
    fits.WriteKey("ns_gamma", TLONG, &ns_gamma, 
	"# of successful shear measurements");
    fits.WriteKey("ns_nativ", TLONG, &ns_native, 
	"# of successful native shapeles measurements");

    fits.WriteKey("nf_rnge1", TLONG, &nf_range1, 
	"# of MeasureShear failures range1");
    fits.WriteKey("nf_rnge2", TLONG, &nf_range2, 
	"# of MeasureShear failures range2");

    fits.WriteKey("nf_nativ", TLONG, &nf_native, 
	"# of MeasureShear failures native calculations");
    fits.WriteKey("nf_small", TLONG, &nf_small, 
	"# of MeasureShear failures too small");

    // prefix ms to make unique
    fits.WriteKey("nf_mstmv", TLONG, &nf_tmverror, 
	"# of MeasureShear failures TMV errors");
    fits.WriteKey("nf_msoth", TLONG, &nf_othererror, 
	"# of MeasureShear failures unclassified errors");

    fits.WriteKey("nf_mu", TLONG, &nf_mu, 
	"# of MeasureShear failures calculating mu");
    fits.WriteKey("nf_gamma", TLONG, &nf_gamma, 
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
	"nf_native="<<nf_native <<" & "<<
	"nf_small="<<nf_small <<" & "<<
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

  os<<"Native fits:\n";
  os<<"  N_Success = "<<ns_native<<std::endl;
  os<<"  N_Fail = "<<nf_native<<std::endl;

  os<<"N_Rejected for being too small = "<<nf_small<<std::endl;

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
  Log(_logfile,_fits_file),
  nstars(0), nf_range(0), nf_tmverror(0), nf_othererror(0),
  ns_psf(0), nf_psf(0) {}

PSFLog::~PSFLog() {
  this->WriteLog();
  this->WriteLogToFitsHeader();
}

void PSFLog::WriteLogToFitsHeader() const
{
  if ("" == fits_file) return;

  try {
    bool create=false;
    //std::cout<<"Writing log to PSF file: "<<fits_file<<std::endl;
    FitsFile fits(fits_file.c_str(), READWRITE, create);
    fits.GotoHDU(2);

    // use fsexit to differentiate from exit codes for other codes
    fits.WriteKey("mpexit", TLONG, &exitcode, "Exit code for MeasurePSF");
    fits.WriteKey("mpnstars", TLONG, &nstars, 
	"# of PSF star candidates used");
    fits.WriteKey("ns_psf", TLONG, &ns_psf, 
	"# of successful PSF decompositions");
    fits.WriteKey("nf_range", TLONG, &nf_range, 
	"# PSF failures due to range error");
    fits.WriteKey("nf_mptmv", TLONG, &nf_tmverror, 
	"# PSF failures due to TMV errors");
    fits.WriteKey("nf_mpoth", TLONG, &nf_othererror, 
	"# PSF failures due to unclassified errors");
    fits.WriteKey("nf_psf", TLONG, &nf_psf, 
	"# PSF failures");
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

  os<<"PSF measurement:\n";
  os<<"  N_Success = "<<ns_psf<<std::endl;
  os<<"  N_Fail = "<<nf_psf<<std::endl;

  os<<"N_Error: TMV Error caught = "<<nf_tmverror<<std::endl;
  os<<"N_Error: Other caught = "<<nf_othererror<<std::endl;
}


FindStarsLog::FindStarsLog(std::string _logfile, std::string _fits_file) :
  Log(_logfile,_fits_file),
  ntot(0), nr_flag(0), nr_mag(0), nr_size(0),
  nobj(0), nallstars(0), nstars(0) {}

FindStarsLog::~FindStarsLog() {
  this->WriteLog();
  this->WriteLogToFitsHeader();
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
    fits.WriteKey("fsexit", TLONG, &exitcode, "Exit code for FindStars");
    fits.WriteKey("fsntot", TLONG, &ntot, 
	"# of total objects processed in FindStars");
    fits.WriteKey("fsnrflag", TLONG, &nr_flag,
	"# of objects rejected because of input flag");
    fits.WriteKey("fsnrmag", TLONG, &nr_mag,
	"# of objects rejected because of mag limits");
    fits.WriteKey("fsnrsize", TLONG, &nr_size,
	"# of objects rejected because of size limits");
    fits.WriteKey("fsnobj", TLONG, &nobj, 
	"# of objects processed in FindStars");
    fits.WriteKey("fsnall", TLONG, &nallstars,
	"# of total stars found by FindStars");
    fits.WriteKey("fsnstars", TLONG, &nstars, 
	"# of good stars found by FindStars");
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
	"ntot="<<ntot <<" & "<<
	"nr_flag="<<nr_flag <<" & "<<
	"nr_mag="<<nr_mag <<" & "<<
	"nr_size="<<nr_size <<" & "<<
	"nobj="<<nobj <<" & "<<
	"nallstars="<<nallstars <<" & "<<
	"nstars="<<nstars <<
	" QA2END"<<std::endl;
    }
  }
}

void FindStarsLog::Write(std::ostream& os) const
{
  os<<"N total input objects = "<<ntot<<std::endl;
  os<<"N rejected because of input flag"<<nr_flag<<std::endl;
  os<<"N rejected because of mag limits"<<nr_mag<<std::endl;
  os<<"N rejected because of size limits"<<nr_size<<std::endl;
  os<<"N objects processed in FindStars"<<nobj<<std::endl;
  os<<"N total stars found by FindStars"<<nallstars<<std::endl;
  os<<"N good stars found by FindStars"<<nstars<<std::endl;
}


