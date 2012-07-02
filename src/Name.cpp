
#include <sys/stat.h> 
#include "Name.h"
#include "Params.h"


// This function is taken from:
// http://www.techbytes.ca/techbyte103.html
bool DoesFileExist(const std::string& file_name) 
{
    struct stat file_info; 
    int stat_value; 

    // Attempt to get the file attributes 
    stat_value = stat(file_name.c_str(),&file_info); 
    if (stat_value == 0) {
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

static int FindExtNum(
    const std::vector<std::string>& ext_list, const std::string& name)
{
    const int next = ext_list.size();
    for(int i=0;i<next;++i) {
        std::string ext = ext_list[i];
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

void SetRoot(ConfigFile& params, const std::string& imageFileName)
{
    // Get image name and extensions from params
    dbg<<"image_file = "<<imageFileName<<std::endl;
    Assert(params.keyExists("image_ext"));
    std::vector<std::string> allImageExt = params["image_ext"];

    // Find which extension is correct
    int extNum = FindExtNum(allImageExt, imageFileName);
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

void SetRoot(ConfigFile& params)
{
    Assert(params.keyExists("root") || params.keyExists("image_file"));

    if (!params.keyExists("root")) {
        // Get image name and extensions from params
        dbg<<"Make root from image_file\n";
        std::string imageFileName = params["image_file"];
        SetRoot(params,imageFileName);
    }
    Assert(params.keyExists("root"));
}

std::string MakeName(
    const ConfigFile& params, const std::string& what,
    bool isInputPrefix, bool must_exist)
{
    xdbg<<"Making name for "<<what<<std::endl;

    std::vector<std::string> ext_list;
    if (params.keyExists(what+"_ext")) ext_list = params[what+"_ext"];
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
        int extNum = FindExtNum(ext_list,name);
        if (!must_exist || DoesFileExist(name)) {
            // File exists or doesn't need to.
            xdbg<<"name = "<<name<<std::endl;
            return name;
        } else if (ext_list.size() > 1 && extNum != -1) {
            // File is supposed to exist and doesn't.
            // But there may be other extensions to try.
            std::string ext = ext_list[extNum];
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
        Assert(ext_list.size() > 0);
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
    const int next = ext_list.size();
    for(int i=0;i<next;++i) {
        std::string name = startName + ext_list[i];
        if (!must_exist || DoesFileExist(name)) {
            xdbg<<"name = "<<name<<std::endl;
            return name;
        }
    }

    // None of the extensions worked, so throw an exception.
    std::string allNames;
    for(int j=0;j<next;++j) {
        if (j > 0) allNames += " , ";
        allNames += startName + ext_list[j];
    }
    throw FileNotFoundException(allNames);
}

std::string MakeFitsName(const ConfigFile& params, const std::string& what) 
{
    // Since we don't know what type the extension is, .fits, .dat, etc.  we will
    // replace .dat with .fits

    // Just get the first name and modify it if necessary
    std::string file = MakeName(params,what,false,false);  

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

std::string AddExtraToName(const std::string& input_name, const std::string extra) 
{

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

std::vector<std::string> MakeMultiName(const ConfigFile& params, const std::string& what)
{
    xdbg<<"Making multi-name for "<<what<<std::endl;
    std::vector<std::string> name_list;

    if (params.keyExists((what+"_file"))) {
        name_list = params[what+"_file"];
    } else {
        Assert(params.keyExists("root"));
        Assert(params.keyExists(what+"_ext"));
        std::string root = params["root"];
        std::vector<std::string> ext_list = params[what+"_ext"];
        std::string pre = "";
        if (params.keyExists((what+"_prefix"))) pre=params[what+"_prefix"];
        else if (params.keyExists("output_prefix")) pre=params["output_prefix"];

        const int next = ext_list.size();
        for(int i=0;i<next;++i) {
            std::string name = pre + root + ext_list[i];
            xdbg<<"name = "<<name<<std::endl;
            name_list.push_back(name);
        }
    }
    return name_list;
}

int GetHdu(const ConfigFile& params, const std::string& what, const std::string& name, int def)
{
    if (params.keyExists(what+"_hdu_by_ext")) {
        std::vector<int> hdu_list = params[what+"_hdu_by_ext"];
        std::vector<std::string> ext_list = params[what+"_ext"];
        int extNum = FindExtNum(ext_list,name);
        return hdu_list[extNum];
    } else if (params.keyExists(what+"_hdu")) {
        return params[what+"_hdu"];
    } else {
        return def;
    }
}


