#ifndef NAME_H
#define NAME_H

#include <string>
#include <vector>
#include <stdexcept>
#include "ConfigFile.h"
#include "FitsFile.h"

void DoWriteParKey(const ConfigFile& params, FitsFile& fits, std::string key);
#define WriteParKey(key) DoWriteParKey(params,fits,key);

struct file_not_found : public std::runtime_error {
  file_not_found( const std::string& filename = std::string() ) throw() :
    std::runtime_error(std::string("Error: file ") + filename
	+ " not found") {} 
};

void SetRoot(ConfigFile& params);

bool FileExists(const std::string& strFilename);

std::string Name(const ConfigFile& params, const std::string& what,
    bool input_prefix=false, bool mustexist=false);

std::vector<std::string> MultiName(const ConfigFile& params,
    const std::string& what);

#endif
