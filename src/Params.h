#ifndef PARAMS_H
#define PARAMS_H

#include <stdexcept>
#include "ConfigFile.h"

// Default value for 
#define DEFVALPOS 9999
#define DEFVALNEG -9999

//
// Flags 
//
#define INPUT_FLAG              0x1
#define TRANSFORM_EXCEPTION     0x2
#define FITTEDPSF_EXCEPTION     0x4
#define TMV_EXCEPTION           0x8
#define STD_EXCEPTION           0x10
#define UNKNOWN_EXCEPTION       0x20
#define EDGE                    0x40
#define LT10PIX                 0x80
#define MEASURE_PSF_FAILED      0x100
#define NATIVE_FAILED           0x200
#define TOO_SMALL               0x400
#define DECONV_FAILED           0x800
#define SHEAR_FAILED            0x1000
#define SHAPELET_FAILED         0x2000
#define UNKNOWN_FAILURE         0x4000
#define SHAPE_REDUCED_ORDER     0x8000
#define SHEAR_LOCAL_MIN         0x10000
#define SHEAR_POOR_FIT          0x20000
#define SHAPE_LOCAL_MIN         0x40000
#define SHAPE_POOR_FIT          0x80000
#define SHEAR_BAD_COVAR         0x100000
#define NO_SINGLE_EPOCH_IMAGES	0x200000
#define BKG_NOPIX		0x400000
#define PSF_INTERP_OUTLIER	0x800000


// Errors specific to the weak lensing code

struct FileNotFound : 
  public std::runtime_error 
{
  FileNotFound(const std::string& filename) throw() :
    std::runtime_error("Error: file "+filename+" not found") {} 
};

struct ParameterError:
  public std::runtime_error 
{
  ParameterError(const std::string& msg) throw() : std::runtime_error(msg) {}
};

struct ReadError:
  public std::runtime_error 
{
  ReadError(const std::string& msg) throw() : std::runtime_error(msg) {}
};

struct WriteError:
  public std::runtime_error 
{
  WriteError(const std::string& msg) throw() : std::runtime_error(msg) {}
};

struct ProcessingError:
  public std::runtime_error 
{
  ProcessingError(const std::string& msg) throw() : std::runtime_error(msg) {}
};

// Errors that may be thrown by the weak lensing code, but 
// defined in other files

// ConfigFile_FileNotFound    -- Treat as FileNotFound
// ConfigFile_KeyNotFound     -- Treat as ParameterError
// ConfigFile_ParameterError  -- Treat as ParameterError
// StarFinderError            -- Treat as ProcessingError
// AssertFailure              -- Treat as ProcessingError
// tmv::Error                 -- Treat as ProcessingError
// std::exception             -- Treat as ProcessingError



//
// Exit codes
//

enum ExitCode { 
  SUCCESS = 0,
  FAILURE,
  FAILURE_FILE_NOT_FOUND,
  FAILURE_PARAMETER_ERROR,
  FAILURE_READ_ERROR,
  FAILURE_WRITE_ERROR,
  FAILURE_PROCESSING_ERROR
};

inline const char* Text(const ExitCode& code)
{
  switch (code) {
    case SUCCESS : return "SUCCESS";
    case FAILURE : return "FAILURE";
    case FAILURE_FILE_NOT_FOUND : return "FAILURE_FILE_NOT_FOUND";
    case FAILURE_PARAMETER_ERROR : return "FAILURE_PARAMETER_ERROR";
    case FAILURE_READ_ERROR : return "FAILURE_READ_ERROR";
    case FAILURE_WRITE_ERROR : return "FAILURE_WRITE_ERROR";
    case FAILURE_PROCESSING_ERROR : return "FAILURE_PROCESSING_ERROR";
    default : return "UNKNOWN";
  }
}

inline int Status(ExitCode code, const ConfigFile& params)
{
  switch (code) {
    case SUCCESS : return params.read("success_status",2);
    case FAILURE : return params.read("failure_status",4);
    case FAILURE_FILE_NOT_FOUND : return params.read("file_not_found_status",5);
    case FAILURE_PARAMETER_ERROR : return params.read("parameter_error_status",5);
    case FAILURE_READ_ERROR : return params.read("read_error_status",5);
    case FAILURE_WRITE_ERROR : return params.read("write_error_status",4);
    case FAILURE_PROCESSING_ERROR : return params.read("processing_error_status",4);
    default : return 0;
  }
}


// tolerance for testing output files
#define TEST_TOL 1.e-6

#endif
