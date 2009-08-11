#ifndef NAME_H
#define NAME_H

#include <string>
#include <vector>
#include "ConfigFile.h"
#include "Params.h"

void SetRoot(ConfigFile& params);

bool FileExists(const std::string& strFilename);

std::string Name(const ConfigFile& params, const std::string& what,
    bool input_prefix=false, bool mustexist=false);

std::vector<std::string> MultiName(const ConfigFile& params,
    const std::string& what);

#endif
