
#include "Params.h"
#include "ConfigFile.h"

std::string FlagText(long flag)
{
    std::string ret = "(";
    bool first = true;
    for(long flagnum = 0; flagnum < NFLAGS; ++flagnum) {
        if (flag & (1 << flagnum)) {
            if (first) first = false; 
            else ret += ",";
            ret += flag_name[flagnum];
        }
    }
    ret += ")";
    return ret;
}

void PrintFlags(const std::vector<long>& flags, std::ostream& os)
{
    const int nobj = flags.size();

    std::vector<long> nflagcount(NFLAGS,0);
    long nnoflag = 0;

    for(int i=0;i<nobj;++i) {
        long flag = flags[i];
        if (!flag) ++nnoflag;
        for(long flagnum = 0; flagnum < NFLAGS; ++flagnum) {
            if (flag & (1 << flagnum)) ++nflagcount[flagnum];
        }
    }
    os<<"   Total N = "<<nobj<<std::endl;
    os<<"     # with no flags = "<<nnoflag<<std::endl;
    for(long flagnum = 0; flagnum < NFLAGS; ++flagnum) {
        if (nflagcount[flagnum]) {
            os<<"     # with "<<flag_name[flagnum]<<" = "<<
                nflagcount[flagnum]<<std::endl;
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
