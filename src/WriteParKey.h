#ifndef WriteParKey_H
#define WriteParKey_H

template <typename T>
static void CCfitsWriteParKey(
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
#endif
