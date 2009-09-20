// ConfigFile.cpp

#include "ConfigFile.h"

ConfigFile::ConfigFile( std::string filename, std::string delimiter,
    std::string comment, std::string sentry ) : 
  myDelimiter(delimiter), myComment(comment), mySentry(sentry)
{
  // Construct a ConfigFile, getting keys and values from given file

  std::ifstream in( filename.c_str() );

  if( !in ) 
#ifdef NOTHROW
  { std::cerr<<"File not found: "<<filename<<std::endl; exit(1); }
#else
    throw ConfigFile_FileNotFound(filename);
#endif

  in >> (*this);
}

ConfigFile::ConfigFile() : myDelimiter("="), myComment("#")
{
  // Construct a ConfigFile without a file; empty
}

void ConfigFile::Load( std::string filename, std::string delimiter,
    std::string comment, std::string sentry )
{
  // Construct a ConfigFile, getting keys and values from given file

  std::string delimiter1 = myDelimiter;
  std::string comment1 = myComment;
  std::string sentry1 = mySentry;

  if (delimiter != "") myDelimiter = delimiter;
  if (comment != "") myComment = comment;
  if (sentry != "") mySentry = sentry;

  std::ifstream in( filename.c_str() );

  if( !in ) 
#ifdef NOTHROW
  { std::cerr<<"File not found: "<<filename<<std::endl; exit(1); }
#else
    throw ConfigFile_FileNotFound(filename);
#endif

  in >> (*this);

  myDelimiter = delimiter1;
  myComment = comment1;
  mySentry = sentry1;
}

ConvertibleString& ConfigFile::getnocheck( const std::string& key )
{
  std::string key2 = key;
  trim(key2);
  return myContents[key2]; 
}

ConvertibleString ConfigFile::get( const std::string& key ) const
{
  std::string key2 = key;
  trim(key2);
  mapci p = myContents.find(key2);
  if (p == myContents.end()) 
#ifdef NOTHROW
  { std::cerr<<"Key not found: "<<key2<<std::endl; exit(1); return key; }
#else
    throw ConfigFile_KeyNotFound(key2);
#endif
  else return p->second;
}

void ConfigFile::remove( const std::string& key )
{
  // Remove key and its value
  myContents.erase( myContents.find( key ) );
  return;
}


bool ConfigFile::keyExists( const std::string& key ) const
{
  // Indicate whether key is found
  mapci p = myContents.find( key );
  return ( p != myContents.end() );
}


void ConfigFile::trim( std::string& s )
{
  // Remove leading and trailing whitespace
  std::string whitespace = " \n\t\v\r\f";
  s.erase( 0, s.find_first_not_of(whitespace) );
  s.erase( s.find_last_not_of(whitespace) + 1);
}


void ConfigFile::Write(std::ostream& os) const
{
  // Save a ConfigFile to os
  for( ConfigFile::mapci p = myContents.begin();
      p != myContents.end();
      ++p )
  {
    os << p->first << " " << myDelimiter << " ";
    os << p->second << std::endl;
  }
}

void ConfigFile::WriteAsComment(std::ostream& os) const
{
  // Save a ConfigFile to os
  for( ConfigFile::mapci p = myContents.begin();
      p != myContents.end();
      ++p )
  {
    std::string f = p->first;
    std::string s = p->second;
    std::replace(f.begin(),f.end(),'\n',' ');
    std::replace(s.begin(),s.end(),'\n',' ');
    os << myComment << " " << f << " " << myDelimiter << " ";
    os << s << std::endl;
  }
}

void ConfigFile::Read(std::istream& is)
{
  // Load a ConfigFile from is
  // Read in keys and values, keeping internal whitespace
  typedef std::string::size_type pos;
  const std::string& delim  = myDelimiter;  // separator
  const std::string& comm   = myComment;    // comment
  const std::string& sentry = mySentry;     // end of file sentry
  const pos skip = delim.length();        // length of separator

  std::string nextline = "";  
  // might need to read ahead to see where value ends

  while( is || nextline.length() > 0 )
  {
    // Read an entire line at a time
    std::string line;
    if( nextline.length() > 0 )
    {
      line = nextline;  // we read ahead; use it now
      nextline = "";
    }
    else
    {
      std::getline( is, line );
    }

    // Ignore comments
    line = line.substr( 0, line.find(comm) );

    // Check for end of file sentry
    if( sentry != "" && line.find(sentry) != std::string::npos ) return;

    // Parse the line if it contains a delimiter
    pos delimPos = line.find( delim );
    if( delimPos < std::string::npos )
    {
      // Extract the key
      std::string key = line.substr( 0, delimPos );
      line.replace( 0, delimPos+skip, "" );

      // See if value continues on the next line
      // Stop at blank line, next line with a key, end of stream,
      // or end of file sentry
      bool terminate = false;
      while( !terminate && is )
      {
	std::getline( is, nextline );
	terminate = true;

	std::string nlcopy = nextline;
	ConfigFile::trim(nlcopy);
	if( nlcopy == "" ) continue;

	nextline = nextline.substr( 0, nextline.find(comm) );
	if( nextline.find(delim) != std::string::npos )
	  continue;
	if( sentry != "" && nextline.find(sentry) != std::string::npos )
	  continue;

	nlcopy = nextline;
	ConfigFile::trim(nlcopy);
	if( nlcopy != "" ) line += "\n";
	line += nextline;
	terminate = false;
      }

      // Store key and value
      ConfigFile::trim(key);
      ConfigFile::trim(line);
      myContents[key] = line;  // overwrites if key is repeated
    }
  }

}
