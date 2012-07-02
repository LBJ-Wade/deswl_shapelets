#ifndef NAME_H
#define NAME_H

#include <string>
#include <vector>
#include "ConfigFile.h"
#include "dbg.h"
#include "Params.h"

void SetRoot(ConfigFile& params);
void SetRoot(ConfigFile& params, const std::string& image_file_name);

bool DoesFileExist(const std::string& file_name);

std::string MakeName(
    const ConfigFile& params, const std::string& what,
    bool is_input_prefix, bool must_exist);

std::vector<std::string> MakeMultiName(
    const ConfigFile& params, const std::string& what);

// By forcing fits we can avoid dealing with all the cases involved with the
// naming scheme
std::string MakeFitsName(const ConfigFile& params, const std::string& what);

// for names like name.ext, add an extra string to make name{extra}.ext
std::string AddExtraToName(
    const std::string& input_name, const std::string extra);

int GetHdu(
    const ConfigFile& params, const std::string& what, 
    const std::string& name, int def);

#endif
