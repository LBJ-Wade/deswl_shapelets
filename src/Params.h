#ifndef PARAMS_H
#define PARAMS_H

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
    case FAILURE_STD_EXCEPTION : return 4;
    case FAILURE_READ_ERROR : return 5;
    case FAILURE_FORMAT_ERROR : return 5;
    case FAILURE_STARFINDER_ERROR : return 4;
    default : return 0;
  }
}


// tolerance for testing output files
#define TEST_TOL 1.e-6

#endif
