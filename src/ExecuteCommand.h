
#ifndef EXECUTE_COMMAND_H
#define EXECUTE_COMMAND_H

#include <string>
#include "dbg.h"

void ExecuteCommand(
    std::string command, std::string& result,
    bool strip_trailing_newline = false);


#endif
