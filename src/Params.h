// Default value for 
#define DEFVALPOS 9999
#define DEFVALNEG -9999

// Flags for a given module should increase monotonically in this file in 
// order to keep things straight

//
// Flags for the DoMeasurePSF processing
//
// These occur in MeasureSinglePSF
#define MSPSF_GETPIX_EXCEPTION 0x1
#define MSPSF_GETPIX_FAILED 0x2
// Note, no try-catch block around Measure(), should we add one?
#define MSPSF_MEASURE_EXCEPTION 0x4
#define MSPSF_MEASURE_FAILED 0x8
#define DMPSF_MSPSF_UNKOWN_EXCEPTION 0x10



//
// Flags for the DoMeasureShear main program processing, and the
// MeasureSingleShear
//
// these occur in main DoMeasureShear loop
#define DMSH_TRANSFORM_EXCEPTION    0x1
#define DMSH_FITTEDPSF_EXCEPTION    0x2
#define DMSH_MSH_TMV_EXCEPTION      0x4
#define DMSH_MSH_UNKNOWN_EXCEPTION  0x8
#define MSH_GETPIX1_FAILED          0x10
#define MSH_GETPIX1_LT10            0x20
#define MSH_NATIVE_FAILED           0x40
#define MSH_TOOSMALL                0x80
#define MSH_GETPIX2_FAILED          0x100
#define MSH_GETPIX2_LT10            0x200
#define MSH_DECONV_FAILED           0x400
#define MSH_SHEAR_FAILED            0x800
