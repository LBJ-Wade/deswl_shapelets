#ifndef PARAMS_H
#define PARAMS_H

// Default value for 
#define DEFVALPOS 9999
#define DEFVALNEG -9999

// Flags for a given module should increase monotonically in this file in 
// order to keep things straight

//
// Flags for the DoMeasurePSF processing
// In DoMeasurePSF or MeasureSinglePSF
//
#define MPSF_TMV_EXCEPTION            0x1
#define MPSF_UNKNOWN_EXCEPTION        0x2
#define MPSF_TRANSFORM_EXCEPTION      0x4
#define MPSF_EDGE1                    0x8
#define MPSF_LT10PIX1                 0x10
#define MPSF_MEASURE_FAILED           0x20


//
// Flags for the DoMeasureShear processing
// In DoMeasureShear or MeasureSingleShear
//
#define MSH_TRANSFORM_EXCEPTION     0x1
#define MSH_FITTEDPSF_EXCEPTION     0x2
#define MSH_TMV_EXCEPTION           0x4
#define MSH_UNKNOWN_EXCEPTION       0x8
#define MSH_EDGE1                   0x10
#define MSH_LT10PIX1                0x20
#define MSH_NATIVE_FAILED           0x40
#define MSH_TOOSMALL                0x80
#define MSH_EDGE2                   0x100
#define MSH_LT10PIX2                0x200
#define MSH_DECONV_FAILED           0x400
#define MSH_SHEAR_FAILED            0x800
#define MSH_SHAPELET_FAILED         0x1000


//
// Exit codes
//
enum ExitCode { 
  SUCCESS			= 0,
  FAILURE			= 1,
  FAILURE_FILE_NOT_FOUND	= 2,
  FAILURE_TMV_ERROR		= 3,
  FAILURE_CONFIGFILE_ERROR	= 4,
  FAILURE_STD_EXCEPTION		= 5,
  FAILURE_READ_ERROR            = 6,
  FAILURE_FORMAT_ERROR          = 7,
  FAILURE_STARFINDER_ERROR      = 8
};

inline const char* Text(const ExitCode& code)
{
  switch (code) {
    case SUCCESS : return "SUCCESS";
    case FAILURE : return "FAILURE";
    case FAILURE_FILE_NOT_FOUND : return "FAILURE_FILE_NOT_FOUND";
    case FAILURE_TMV_ERROR : return "FAILURE_TMV_ERROR";
    case FAILURE_CONFIGFILE_ERROR : return "FAILURE_CONFIGFILE_ERROR";
    case FAILURE_STD_EXCEPTION : return "FAILURE_STD_EXCEPTION";
    case FAILURE_READ_ERROR : return "FAILURE_READ_ERROR";
    case FAILURE_FORMAT_ERROR : return "FAILURE_FORMAT_ERROR";
    case FAILURE_STARFINDER_ERROR : return "FAILURE_FORMAT_ERROR";
    default : return "UNKNOWN";
  }
}

inline int Status(ExitCode code)
{
  switch (code) {
    case SUCCESS : return 2;
    case FAILURE : return 5;
    case FAILURE_FILE_NOT_FOUND : return 5;
    case FAILURE_TMV_ERROR : return 4;
    case FAILURE_CONFIGFILE_ERROR : return 5;
    case FAILURE_STD_EXCEPTION : return 5;
    case FAILURE_READ_ERROR : return 5;
    case FAILURE_FORMAT_ERROR : return 5;
    case FAILURE_STARFINDER_ERROR : return 4;
    default : return 0;
  }
}


// tolerance for testing output files
#define TEST_TOL 1.e-6

#endif
