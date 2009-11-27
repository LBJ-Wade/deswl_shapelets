#ifndef NAME_H
#define NAME_H

#include <string>
#include <vector>
#include "ConfigFile.h"
#include "Params.h"

void setRoot(ConfigFile& params);
void setRoot(ConfigFile& params, const std::string& imageFileName);

bool doesFileExist(const std::string& fileName);

std::string makeName(
    const ConfigFile& params, const std::string& what,
    bool isInputPrefix, bool mustExist);

std::vector<std::string> makeMultiName(
    const ConfigFile& params, const std::string& what);

#endif
