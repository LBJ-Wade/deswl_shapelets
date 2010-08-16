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

// By forcing fits we can avoid dealing with all the cases involved with the
// naming scheme
std::string makeFitsName(
    const ConfigFile& params, 
    const std::string& what);

// for names like name.ext, add an extra string to make name{extra}.ext
std::string addExtraToName(
    const std::string& input_name, 
    const std::string extra);

int getHdu(
    const ConfigFile& params, const std::string& what, 
    const std::string& name, int def);

#endif
