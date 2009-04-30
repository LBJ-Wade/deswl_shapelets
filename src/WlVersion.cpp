#include "WlVersion.h"

std::string WlVersion()
{
	std::string badvers="NOTAG: unparseable";
	std::string thisname="/ccode/test/test_svn_version.cpp";

	// this gets automatically updated by svn
	// don't forget to do svn propset svn:keywords "HeadURL" thisfile
	std::string svnURL="$HeadURL$";

	std::vector<std::string> urltokens;
	std::istringstream urlstream(svnURL);
	std::string temp;
	while (urlstream >> temp) {
		urltokens.push_back(temp);
	}

	if (urltokens.size() < 3) {
		std::cerr<<"URL string is not in 3 elements: \n\t"<<svnURL<<std::endl;
		return badvers;
	} 
	std::string url=urltokens[1];

	size_t ind=url.find(thisname);

	if (ind == std::string::npos) {
		std::cerr<<"Could not find '"<<thisname<<"' in url"<<std::endl;
		return badvers;
	} else {
		std::string front=url.substr(0,ind);
		// now grab the basename
		std::string version = basename((char*)front.c_str());
		return version;
	}

}


