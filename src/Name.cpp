
#include <sys/stat.h> 
#include "Name.h"
#include "dbg.h"
#include "Params.h"


// This function is taken from:
// http://www.techbytes.ca/techbyte103.html
bool FileExists(const std::string& strFilename) 
{
  struct stat stFileInfo; 
  int intStat; 

  // Attempt to get the file attributes 
  intStat = stat(strFilename.c_str(),&stFileInfo); 
  if (intStat == 0) 
  {
    // We were able to get the file attributes 
    // so the file obviously exists. 
    return true;
  }
  else 
  {
    // We were not able to get the file attributes. 
    // This may mean that we don't have permission to 
    // access the folder which contains this file. If you 
    // need to do that level of checking, lookup the 
    // return values of stat which will give you 
    // more details on why stat failed. 
    return false;
  } 
}

void SetRoot(ConfigFile& params)
{
  Assert(params.keyExists("root") || params.keyExists("image_file"));

  if (!params.keyExists("root")) 
  {
    // Get image name and extensions from params
    dbg<<"Make root from image_file\n";
    std::string image_file = params["image_file"];
    dbg<<"image_file = "<<image_file<<std::endl;
    Assert(params.keyExists("image_ext"));
    std::vector<std::string> all_image_ext = params["image_ext"];

    // Find which extension is correct
    std::string image_ext;
    size_t ext_pos;
    bool found = false;
    for(size_t i=0;i<all_image_ext.size();i++) 
    {
      image_ext = all_image_ext[i];
      dbg<<"image_ext = "<<image_ext<<std::endl;
      if (image_ext.size() > image_file.size()) continue;
      ext_pos = image_file.size() - image_ext.size();
      dbg<<"ext_pos = "<<ext_pos<<std::endl;
      std::string image_ext2 = std::string(image_file,ext_pos);
      dbg<<"image_ext2 = "<<image_ext2<<std::endl;
      if (image_ext2 == image_ext) 
      {
	found = true;
	break;
      }
    }
    if (!found) 
    {
      throw(ParameterError("image file "+image_file+
	    " does not end with "+params["image_ext"]));
    }

    // Calculate the root
    size_t root_pos = image_file.find_last_of('/',ext_pos);
    if (root_pos == std::string::npos) root_pos = 0;
    else ++root_pos;
    dbg<<"root_pos = "<<root_pos<<std::endl;
    std::string root = std::string(image_file,root_pos,ext_pos-root_pos);
    dbg<<"root = "<<root<<std::endl;

    // Calculate the prefix
    std::string prefix = std::string(image_file,0,root_pos);
    dbg<<"prefix = "<<prefix<<std::endl;

    // Double check
    std::string image_file2 = prefix + root + image_ext;
    dbg<<"image_file2 = "<<image_file2<<std::endl;
    Assert(image_file == image_file2);

    // Assign results back to params
    params["root"] = root;
    params["input_prefix"] = prefix;
    if (!params.keyExists("output_prefix")) params["output_prefix"] = prefix;
  }

  Assert(params.keyExists("root"));
}

std::string Name(const ConfigFile& params, const std::string& what,
    bool input_prefix, bool mustexist)
{
  if (input_prefix) mustexist = true;

  xdbg<<"Making name for "<<what<<std::endl;
  std::string name;
  if (params.keyExists((what+"_file"))) 
  {
    name = params[what+"_file"];
    if (mustexist) 
    {
      if (!FileExists(name)) throw FileNotFound(name);
    }
  }
  else 
  {
    Assert(params.keyExists("root"));
    Assert(params.keyExists(what+"_ext"));
    std::string root = params["root"];
    std::vector<std::string> ext = params[what+"_ext"];
    std::string pre = "";
    if (params.keyExists((what+"_prefix"))) pre=params[what+"_prefix"];
    else 
    {
      if (input_prefix) 
      {
	if (params.keyExists("input_prefix")) pre=params["input_prefix"];
      }
      else 
      {
	if (params.keyExists("output_prefix")) pre=params["output_prefix"];
      }
    }

    for(size_t i=0;i<ext.size();i++) 
    {
      name = pre + root + ext[i];
      if (mustexist) 
      {
	if (FileExists(name)) break;
	if (i == ext.size()-1) 
	{
	  std::string allnames;
	  for(size_t j=0;j<ext.size();j++) 
	  {
	    if (j > 0) allnames += ",";
	    name = pre + root + ext[j];
	    allnames += name;
	  }
	  throw FileNotFound(allnames);
	}
      }
      else 
      {
	break;
      }
    }
  }
  xdbg<<"name = "<<name<<std::endl;
  return name;
}

std::vector<std::string> MultiName(const ConfigFile& params,
    const std::string& what)
{
  xdbg<<"Making multi-name for "<<what<<std::endl;
  std::vector<std::string> names;

  if (params.keyExists((what+"_file"))) 
  {
    names = params[what+"_file"];
  }
  else 
  {
    Assert(params.keyExists("root"));
    Assert(params.keyExists(what+"_ext"));
    std::string root = params["root"];
    std::vector<std::string> ext = params[what+"_ext"];
    std::string pre = "";
    if (params.keyExists((what+"_prefix"))) pre=params[what+"_prefix"];
    else if (params.keyExists("output_prefix")) pre=params["output_prefix"];

    for(size_t i=0;i<ext.size();i++) 
    {
      std::string name = pre + root + ext[i];
      xdbg<<"name = "<<name<<std::endl;
      names.push_back(name);
    }
  }
  return names;
}
