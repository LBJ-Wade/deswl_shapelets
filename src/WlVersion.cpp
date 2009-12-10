#include "WlVersion.h"
#include <iostream>
#include <vector>
#include <sstream>
// for basename
#include <libgen.h>

std::string getWlVersion()
{
    std::string badVers="NOTAG: unparseable";
    std::string thisName="/src/WlVersion.cpp";

    // this gets automatically updated by svn
    // don't forget to do svn propset svn:keywords "HeadURL" thisfile
    std::string svnUrl=
        "$HeadURL$";

    std::vector<std::string> urlTokens;
    std::istringstream urlStream(svnUrl);
    std::string temp;
    while (urlStream >> temp) {
        urlTokens.push_back(temp);
    }

    if (urlTokens.size() < 3) {
        std::cerr<<"URL string is not in 3 elements: \n\t"<<svnUrl<<std::endl;
        return badVers;
    } 
    std::string url=urlTokens[1];

    std::string::size_type ind=url.find(thisName);

    if (ind == std::string::npos) {
        std::cerr<<"Could not find '"<<thisName<<"' in url"<<std::endl;
        return badVers;
    } else {
        std::string front=url.substr(0,ind);
        // now grab the basename
        std::string version = basename((char*)front.c_str());
        return version;
    }
}
