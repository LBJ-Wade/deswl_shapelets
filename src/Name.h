#ifndef NAME_H
#define NAME_H

#include <string>
#include <vector>
#include "ConfigFile.h"
#include <CCfits/CCfits>
#include "Params.h"

template <typename T>
void CCfitsWriteParKey(
    const ConfigFile& params, 
    CCfits::Table* table,
    std::string key,
    T& tmpvar)
{
  if (params.keyExists(key))  {
    try {
      tmpvar = params.get(key);
      table->addKey(
	  params.get(key+"_hname"), 
	  tmpvar,
	  params.get(key+"_comment"));
    }
    catch (ConvertibleStringError& e)
    {
      throw ConfigFile_ParameterError(key,e.what());
    }
  }
  else {
    throw ParameterError("Error: the key, " + key +
	", is not in the parameter file");
  }
}


void SetRoot(ConfigFile& params);

bool FileExists(const std::string& strFilename);

std::string Name(const ConfigFile& params, const std::string& what,
    bool input_prefix=false, bool mustexist=false);

std::vector<std::string> MultiName(const ConfigFile& params,
    const std::string& what);

#endif
