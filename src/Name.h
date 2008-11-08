#ifndef NAME_H
#define NAME_H

inline std::string Name(ConfigFile& params, const std::string& what)
{
  xdbg<<"Making name for "<<what<<std::endl;
  std::string name;
  if (params.keyExists((what+"_file"))) {
    name = params[what+"_file"];
  } else {
    Assert(params.keyExists("root"));
    Assert(params.keyExists(what+"_ext"));
    name = params["root"] + params[what+"_ext"];
  }
  xdbg<<"name = "<<name<<std::endl;
  return name;
}

#endif
