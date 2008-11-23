#ifndef NAME_H
#define NAME_H


// Taken from:
// http://www.techbytes.ca/techbyte103.html
#include <sys/stat.h> 
inline bool FileExists(const std::string& strFilename) 
{ 
  struct stat stFileInfo; 
  int intStat; 

  // Attempt to get the file attributes 
  intStat = stat(strFilename.c_str(),&stFileInfo); 
  if (intStat == 0) { 
    // We were able to get the file attributes 
    // so the file obviously exists. 
    return true;
  } else { 
    // We were not able to get the file attributes. 
    // This may mean that we don't have permission to 
    // access the folder which contains this file. If you 
    // need to do that level of checking, lookup the 
    // return values of stat which will give you 
    // more details on why stat failed. 
    return false;
  } 
}

struct file_not_found : public std::runtime_error {
  file_not_found( const std::string& filename = std::string() ) throw() :
    std::runtime_error(std::string("Error: file ") + filename
	+ " not found") {} 
};

inline std::string Name(const ConfigFile& params, const std::string& what,
    bool input_prefix=false, bool mustexist=false)
{
  if (input_prefix) mustexist = true;

  xdbg<<"Making name for "<<what<<std::endl;
  std::string name;
  if (params.keyExists((what+"_file"))) {
    name = params[what+"_file"];
    if (mustexist) {
      if (!FileExists(name)) throw file_not_found(name);
    }
  } else {
    Assert(params.keyExists("root"));
    Assert(params.keyExists(what+"_ext"));
    std::string root = params["root"];
    std::vector<std::string> ext = params[what+"_ext"];
    std::string pre = "";
    if (params.keyExists((what+"_prefix"))) pre=params[what+"_prefix"];
    else {
      if (input_prefix) {
	if (params.keyExists("input_prefix")) pre=params["input_prefix"];
      } else {
	if (params.keyExists("output_prefix")) pre=params["output_prefix"];
      }
    }

    for(size_t i=0;i<ext.size();i++) {
      name = pre + root + ext[i];
      if (mustexist) {
	if (FileExists(name)) break;
	if (i == ext.size()-1) throw file_not_found(name);
      } else {
	break;
      }
    }
  }
  xdbg<<"name = "<<name<<std::endl;
  return name;
}

#endif
