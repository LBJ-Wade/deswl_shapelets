
#include "Params.h"
#include "ConfigFile.h"

void PrintFlags(const std::vector<long>& flags, std::ostream& os)
{
    const int nObj = flags.size();

    std::vector<long> nFlagCount(NFLAGS,0);
    long nNoFlag = 0;

    for(int i=0;i<nObj;++i) {
        long flag = flags[i];
        if (!flag) ++nNoFlag;
        for(long flagNum = 0; flagNum < NFLAGS; ++flagNum) {
            if (flag & (1 << flagNum)) ++nFlagCount[flagNum];
        }
    }
    os<<"Total N = "<<nObj<<std::endl;
    os<<"N with no flags = "<<nNoFlag<<std::endl;
    for(long flagNum = 0; flagNum < NFLAGS; ++flagNum) {
        if (nFlagCount[flagNum]) {
            os<<"N with "<<flagName[flagNum]<<" = "<<
                nFlagCount[flagNum]<<std::endl;
        }
    }
}

// Errors specific to the weak lensing code
const char* Text(const ExitCode& code)
{
    switch (code) {
      case SUCCESS :
           return "SUCCESS";
      case FAILURE :
           return "FAILURE";
      case FAILURE_FILE_NOT_FOUND :
           return "FAILURE_FILE_NOT_FOUND";
      case FAILURE_PARAMETER_ERROR :
           return "FAILURE_PARAMETER_ERROR";
      case FAILURE_READ_ERROR :
           return "FAILURE_READ_ERROR";
      case FAILURE_WRITE_ERROR :
           return "FAILURE_WRITE_ERROR";
      case FAILURE_PROCESSING_ERROR :
           return "FAILURE_PROCESSING_ERROR";
      default :
           return "UNKNOWN";
    }
}

int Status(ExitCode code, const ConfigFile& params)
{
    switch (code) {
      case SUCCESS :
           return params.read("success_status",2);
      case FAILURE :
           return params.read("failure_status",4);
      case FAILURE_FILE_NOT_FOUND :
           return params.read("file_not_found_status",5);
      case FAILURE_PARAMETER_ERROR :
           return params.read("parameter_error_status",5);
      case FAILURE_READ_ERROR :
           return params.read("read_error_status",5);
      case FAILURE_WRITE_ERROR :
           return params.read("write_error_status",4);
      case FAILURE_PROCESSING_ERROR :
           return params.read("processing_error_status",4);
      default :
           return 0;
    }
}
