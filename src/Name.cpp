
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

static int findExtNum(
    const std::vector<std::string>& extList, const std::string& name)
{
    const int nExt = extList.size();
    for(int i=0;i<nExt;++i) {
        std::string ext = extList[i];
        xdbg<<"ext = "<<ext<<std::endl;
        if (ext.size() > name.size()) continue;
        int extPos = name.size() - ext.size();
        xdbg<<"ext_pos = "<<extPos<<std::endl;
        std::string ext2 = std::string(name,extPos);
        xdbg<<"ext2 = "<<ext2<<std::endl;
        if (ext2 == ext) return i;
    }
    return -1;
}

void setRoot(ConfigFile& params, const std::string& imageFileName)
{
    // Get image name and extensions from params
    dbg<<"image_file = "<<imageFileName<<std::endl;
    Assert(params.keyExists("image_ext"));
    std::vector<std::string> allImageExt = params["image_ext"];

    // Find which extension is correct
    int extNum = findExtNum(allImageExt, imageFileName);
    if (extNum == -1) {
        throw ParameterException(
            "image file "+imageFileName+
            " does not end with "+params["image_ext"]);
    }
    std::string imageExt = allImageExt[extNum];

    // Calculate the root
    int extPos = imageFileName.size() - imageExt.size();
    std::string::size_type rootPos = imageFileName.find_last_of('/',extPos);
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

    if (!params.keyExists("root")) {
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
    xdbg<<"Making name for "<<what<<std::endl;

    std::vector<std::string> extList;
    if (params.keyExists(what+"_ext")) extList = params[what+"_ext"];
    std::string startName;

    if (params.keyExists((what+"_file"))) {
        xdbg<<(what+"_file")<<" parameter exists, so use that.\n";
        xdbg<<"value = "<<params[what+"_file"]<<std::endl;
        std::vector<std::string> multiname = params[what+"_file"];
        xdbg<<"multiname.size = "<<multiname.size()<<std::endl;
        if (multiname.size() == 0) {
            throw ParameterException(
                std::string("Name for ") + what + "_file " +
                "parses as a vector with no elements");
        }
        std::string name = multiname[0];
        xdbg<<"name = "<<name<<std::endl;
        int extNum = findExtNum(extList,name);
        if (!mustExist || doesFileExist(name)) {
            // File exists or doesn't need to.
            xdbg<<"name = "<<name<<std::endl;
            return name;
        } else if (extList.size() > 1 && extNum != -1) {
            // File is supposed to exist and doesn't.
            // But there may be other extensions to try.
            std::string ext = extList[extNum];
            int extPos = name.size() - ext.size();
            startName = std::string(name,0,extPos);
        } else {
            // No other options to try, so throw an error.
            throw FileNotFoundException(name);
        }
    } else {
        // Determine the name from the root, prefix, and ext.
        Assert(params.keyExists("root"));
        Assert(params.keyExists(what+"_ext"));
        Assert(extList.size() > 0);
        std::string root = params["root"];
        xdbg<<"No "<<(what+"_file")<<" parameter, so use root.\n";
        xdbg<<"root = "<<root<<std::endl;
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
        startName = pre + root;
    }

    // Add each extension in turn and check if the file exists.
    const int nExt = extList.size();
    for(int i=0;i<nExt;++i) {
        std::string name = startName + extList[i];
        if (!mustExist || doesFileExist(name)) {
            xdbg<<"name = "<<name<<std::endl;
            return name;
        }
    }

    // None of the extensions worked, so throw an exception.
    std::string allNames;
    for(int j=0;j<nExt;++j) {
        if (j > 0) allNames += " , ";
        allNames += startName + extList[j];
    }
    throw FileNotFoundException(allNames);
}

std::string makeFitsName(
    const ConfigFile& params, 
    const std::string& what) {

    // Since we don't know what type the extension is, .fits, .dat, etc.  we will
    // replace .dat with .fits

    // Just get the first name and modify it if necessary
    std::string file = makeName(params,what,false,false);  

    size_t fits_pos = file.rfind(".fits");
    size_t fit_pos = file.rfind(".fit");
    size_t dat_pos = file.rfind(".dat");

    if (fits_pos==std::string::npos && fit_pos==std::string::npos) {
        if (dat_pos != std::string::npos) {
            // replace the .dat with .fits
            std::string e = ".fits";
            file.replace(dat_pos, 4, e);
        } else {
            // didn't find the .dat;  we will just tac it on the end
            file += ".fits";
        }
    }

    return file;
}

std::string addExtraToName(
    const std::string& input_name, 
    const std::string extra) {

    std::string name = input_name;

    if (extra == "") {
        return name;
    }

    size_t pos = name.rfind(".");
    if (pos != std::string::npos) {
        name.insert(pos, extra);
    } else {
        // didn't find the expected extension ending with a .ext;  we will just tac
        // it on the end
        name+= extra;
    }

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
        std::vector<std::string> extList = params[what+"_ext"];
        std::string pre = "";
        if (params.keyExists((what+"_prefix"))) pre=params[what+"_prefix"];
        else if (params.keyExists("output_prefix")) pre=params["output_prefix"];

        const int nExt = extList.size();
        for(int i=0;i<nExt;++i) {
            std::string name = pre + root + extList[i];
            xdbg<<"name = "<<name<<std::endl;
            nameList.push_back(name);
        }
    }
    return nameList;
}

int getHdu(const ConfigFile& params, const std::string& what, 
           const std::string& name, int def)
{
    if (params.keyExists(what+"_hdu_by_ext")) {
        std::vector<int> hduList = params[what+"_hdu_by_ext"];
        std::vector<std::string> extList = params[what+"_ext"];
        int extNum = findExtNum(extList,name);
        return hduList[extNum];
    } else if (params.keyExists(what+"_hdu")) {
        return params[what+"_hdu"];
    } else {
        return def;
    }
}


