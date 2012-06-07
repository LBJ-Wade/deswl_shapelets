#include <iostream>
#include <vector>
#include <sstream>
#include <libgen.h> // for basename

#include "WlVersion.h"

std::string GetWlVersion()
{
    std::string bad_vers="NOTAG: unparseable";
    std::string this_name="/src/WlVersion.cpp";

    // this gets automatically updated by svn
    // don't forget to do svn propset svn:keywords "HeadURL" thisfile
    std::string svn_url=
        "$HeadURL$";

    std::vector<std::string> urlTokens;
    std::istringstream urlStream(svn_url);
    std::string temp;
    while (urlStream >> temp) {
        urlTokens.push_back(temp);
    }

    if (urlTokens.size() < 3) {
        std::cerr<<"URL string is not in 3 elements: \n\t"<<svn_url<<std::endl;
        return bad_vers;
    } 
    std::string url=urlTokens[1];

    std::string::size_type ind=url.find(this_name);

    if (ind == std::string::npos) {
        std::cerr<<"Could not find '"<<this_name<<"' in url"<<std::endl;
        return bad_vers;
    } else {
        std::string front=url.substr(0,ind);
        // now grab the basename
        std::string version = basename((char*)front.c_str());
        return version;
    }
}
