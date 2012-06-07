#include <cstdio>
#include "ExecuteCommand.h"

void executeCommand(
    std::string command, std::string& result,
    bool strip_trailing_newline)
{
    result.erase();
    char buffer[256];
    FILE* stream = popen(command.c_str(),"r");
    while ( fgets(buffer, 256, stream) != NULL ) {
        result.append(buffer);
    }
    pclose(stream);

    if (strip_trailing_newline) {
        int len=result.size();
        if (len > 0) {
            if (result[len-1] == '\n') {
                std::string tmp;
                tmp.resize(len-1);
                for (int i=0; i<len-1; ++i) {
                    tmp[i] = result[i];
                }
                result=tmp;
            }
        }
    }
}


