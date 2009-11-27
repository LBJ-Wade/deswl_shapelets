#ifndef WriteParKey_H
#define WriteParKey_H

template <typename T>
static void writeParamToTable(
    const ConfigFile& params, CCfits::Table* table,
    std::string key, T& tempVar)
{
    // read rather than get, so ConfigFile throws a more 
    // descriptive error on failed conversion or missing key.
    tempVar = params.read<T>(key);
    table->addKey(
        params.get(key+"_hname"), 
        tempVar,
        params.get(key+"_comment"));
}

#endif
