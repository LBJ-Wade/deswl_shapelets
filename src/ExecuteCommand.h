
#ifndef EXECUTE_COMMAND_H
#define EXECUTE_COMMAND_H

#include <string>

void executeCommand(
    std::string command, std::string& result,
    bool shouldStripTrailingNewline = false);


#endif
