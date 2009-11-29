
#include <sys/stat.h> 
#include "Name.h"
#include "dbg.h"
#include "Params.h"


// This function is taken from:
// http://www.techbytes.ca/techbyte103.html
bool doesFileExist(const std::string& fileName) 
{
    struct stat fileInfo; 
    int statValue; 

    // Attempt to get the file attributes 
    statValue = stat(fileName.c_str(),&fileInfo); 
    if (statValue == 0) {
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

void setRoot(ConfigFile& params, const std::string& imageFileName)
{
    // Get image name and extensions from params
    dbg<<"image_file = "<<imageFileName<<std::endl;
    Assert(params.keyExists("image_ext"));
    std::vector<std::string> allImageExt = params["image_ext"];
    const int nExt = allImageExt.size();

    // Find which extension is correct
    std::string imageExt;
    int extPos;
    bool isFound = false;
    for(int i=0;i<nExt;i++) 
    {
        imageExt = allImageExt[i];
        dbg<<"image_ext = "<<imageExt<<std::endl;
        if (imageExt.size() > imageFileName.size()) continue;
        extPos = imageFileName.size() - imageExt.size();
        dbg<<"ext_pos = "<<extPos<<std::endl;
        std::string imageExt2 = std::string(imageFileName,extPos);
        dbg<<"image_ext2 = "<<imageExt2<<std::endl;
        if (imageExt2 == imageExt) {
            isFound = true;
            break;
        }
    }
    if (!isFound) {
        throw ParameterException(
                "image file "+imageFileName+
                " does not end with "+params["image_ext"]);
    }

    // Calculate the root
    int rootPos = imageFileName.find_last_of('/',extPos);
    if (rootPos == std::string::npos) rootPos = 0;
    else ++rootPos;
    dbg<<"root_pos = "<<rootPos<<std::endl;
    std::string root = std::string(imageFileName,rootPos,extPos-rootPos);
    dbg<<"root = "<<root<<std::endl;

    // Calculate the prefix
    std::string prefix = std::string(imageFileName,0,rootPos);
    dbg<<"prefix = "<<prefix<<std::endl;

    // Double check
    std::string imageFileName2 = prefix + root + imageExt;
    dbg<<"image_file2 = "<<imageFileName2<<std::endl;
    Assert(imageFileName == imageFileName2);

    // Assign results back to params
    params["root"] = root;
    params["input_prefix"] = prefix;
    if (!params.keyExists("output_prefix")) params["output_prefix"] = prefix;
}

void setRoot(ConfigFile& params)
{
    Assert(params.keyExists("root") || params.keyExists("image_file"));

    if (!params.keyExists("root")) 
    {
        // Get image name and extensions from params
        dbg<<"Make root from image_file\n";
        std::string imageFileName = params["image_file"];
        setRoot(params,imageFileName);
    }
    Assert(params.keyExists("root"));
}

std::string makeName(
    const ConfigFile& params, const std::string& what,
    bool isInputPrefix, bool mustExist)
{
    if (isInputPrefix) mustExist = true;

    xdbg<<"Making name for "<<what<<std::endl;
    std::string name;
    if (params.keyExists((what+"_file"))) 
    {
        xdbg<<(what+"_file")<<" parameter exists, so use that.\n";
        name = params[what+"_file"];
        if (mustExist) 
        {
            if (!doesFileExist(name)) throw FileNotFoundException(name);
        }
    } else {
        Assert(params.keyExists("root"));
        Assert(params.keyExists(what+"_ext"));
        std::string root = params["root"];
        xdbg<<"No "<<(what+"_file")<<" parameter, so use root.\n";
        xdbg<<"root = "<<root<<std::endl;
        std::vector<std::string> ext = params[what+"_ext"];
        xdbg<<"ext = "<<params[what+"_ext"]<<std::endl;
        std::string pre = "";
        if (params.keyExists((what+"_prefix"))) {
            pre=params[what+"_prefix"];
        } else {
            if (isInputPrefix) {
                if (params.keyExists("input_prefix")) 
                    pre=params["input_prefix"];
            } else {
                if (params.keyExists("output_prefix")) 
                    pre=params["output_prefix"];
            }
        }
        xdbg<<"pre = "<<pre<<std::endl;

        const int nExt = ext.size();
        for(int i=0;i<nExt;i++) {
            name = pre + root + ext[i];
            if (mustExist) {
                if (doesFileExist(name)) break;
                if (i == nExt-1) {
                    std::string allNames;
                    for(int j=0;j<nExt;j++) {
                        if (j > 0) allNames += ",";
                        name = pre + root + ext[j];
                        allNames += name;
                    }
                    throw FileNotFoundException(allNames);
                }
            } else {
                break;
            }
        }
    }
    xdbg<<"name = "<<name<<std::endl;
    return name;
}

std::vector<std::string> makeMultiName(
    const ConfigFile& params, const std::string& what)
{
    xdbg<<"Making multi-name for "<<what<<std::endl;
    std::vector<std::string> nameList;

    if (params.keyExists((what+"_file"))) {
        nameList = params[what+"_file"];
    } else {
        Assert(params.keyExists("root"));
        Assert(params.keyExists(what+"_ext"));
        std::string root = params["root"];
        std::vector<std::string> ext = params[what+"_ext"];
        std::string pre = "";
        if (params.keyExists((what+"_prefix"))) pre=params[what+"_prefix"];
        else if (params.keyExists("output_prefix")) pre=params["output_prefix"];

        const int nExt = ext.size();
        for(int i=0;i<nExt;i++) {
            std::string name = pre + root + ext[i];
            xdbg<<"name = "<<name<<std::endl;
            nameList.push_back(name);
        }
    }
    return nameList;
}
