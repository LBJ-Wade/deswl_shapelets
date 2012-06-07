#ifndef WriteParam_H
#define WriteParam_H

template <typename T>
static void WriteParamToTable(
    const ConfigFile& params, CCfits::Table* table,
    std::string key, T& temp_var)
{
    // read rather than get, so ConfigFile throws a more 
    // descriptive error on failed conversion or missing key.
    temp_var = params.read<T>(key);
    table->addKey(
        params.get(key+"_hname"), 
        temp_var,
        params.get(key+"_comment"));
}

#endif
