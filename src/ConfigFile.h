// ConfigFile.h
// Class for reading named values from configuration files
// Richard J. Wagner  v2.1  24 May 2004  wagnerr@umich.edu

// Copyright (c) 2004 Richard J. Wagner
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

// Typical usage
// -------------
// 
// Given a configuration file "settings.inp":
//   atoms  = 25
//   length = 8.0  # nanometers
//   name = Reece Surcher
// 
// Named values are read in various ways, with or without default values:
//   ConfigFile config( "settings.inp" );
//   int atoms = config.read<int>( "atoms" );
//   double length = config.read( "length", 10.0 );
//   string author, title;
//   config.readInto( author, "name" );
//   config.readInto( title, "title", string("Untitled") );
// 
// See file example.cpp for more examples.

#ifndef CONFIGFILE_H
#define CONFIGFILE_H

#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <vector>
#include <algorithm>
#include <stdexcept>

class ConvertibleString : public std::string
{

  public:
    ConvertibleString() : std::string("") {};
    ConvertibleString(const std::string sin) : std::string(sin) {}

    template <typename T> explicit ConvertibleString(const T& x)
    {
      std::stringstream ss(*this);
      ss << x;
    }
    
    template <typename T> explicit ConvertibleString(const std::vector<T>& x)
    {
      std::stringstream ss(*this);
      ss << "{";
      if (x.size() > 0) ss << x[0];
      for(size_t i=1;i<x.size();++i) ss << ',' << x[i];
      ss << "}";
    }

    ~ConvertibleString() {};

    ConvertibleString& operator=(const std::string& rhs)
    {
      std::string::operator=(rhs); 
      return *this;
    }

    template <typename T> ConvertibleString& operator=(const T& x)
    {
      std::stringstream ss;
      ss << x;
      *this = ss.str();
      return *this;
    }

    template <typename T> ConvertibleString& operator=(
	const std::vector<T>& x)
    {
      std::stringstream ss;
      if (x.size() > 0) ss << x[0];
      for(size_t i=1;i<x.size();++i) ss << ' ' << x[i];
      *this = ss.str();
      return *this;
    }

    template <typename T> operator T() const 
    {
#ifdef Use_Zero_Default
      if (*this == "") return T();
#endif

      std::string err="Could not convert ConvertibleString to input type ";
      err += typeid(T).name();
      err += std::string(": this = ") + *this;
      T temp;
      std::stringstream ss(*this);
      ss >> temp;
      if (!ss) 
#ifndef NOTHROW
      { std::cerr<<err<<std::endl; exit(1); }
#else
	throw std::runtime_error(err);
#endif
      return temp;
    }

    template <typename T> operator std::vector<T>() const 
    {
#ifdef Use_Zero_Default
      if (*this == "") return std::vector<T>();
#endif

      std::string err="Could not convert ConvertibleString to input type ";
      err += std::string("std::vector<")+typeid(T).name()+">";
      err += std::string(": this = ") + *this;

      // Two syntaxes: "{1.,2.,3.}" or "1. 2. 3."
      if ((*this)[0] == '{') {
	// Using "{1.,2.,3.}" syntax
	size_t i1 = this->find_first_not_of(" \t\n\r\f",1);
	if (i1 == npos) 
#ifdef NOTHROW
	{ std::cerr<<err<<std::endl; exit(1); return std::vector<T>(); }
#else
	throw std::runtime_error(err);
#endif
	else if ((*this)[i1] == '}') {
	  // Then "{  }"
	  return std::vector<T>();
	} 
	else {
	  char ch;
	  int ncomma = std::count(this->begin(),this->end(),',');
	  std::vector<T> temp(ncomma+1);
	  std::stringstream ss(*this);
	  ss >> ch;
	  if (!ss || ch != '{') 
#ifdef NOTHROW
	  { std::cerr<<err<<std::endl; exit(1); }
#else
	  throw std::runtime_error(err);
#endif
	  ss >> temp[0];
	  if (!ss) 
#ifdef NOTHROW
	  { std::cerr<<err<<std::endl; exit(1); }
#else
	  throw std::runtime_error(err);
#endif
	  for(size_t i=1;i<temp.size();++i) {
	    ss >> ch;
	    if (!ss || ch != ',') 
#ifdef NOTHROW
	    { std::cerr<<err<<std::endl; exit(1); }
#else
	    throw std::runtime_error(err);
#endif
	    ss >> temp[i];
	    if (!ss) 
#ifdef NOTHROW
	    { std::cerr<<err<<std::endl; exit(1); }
#else
	    throw std::runtime_error(err);
#endif
	  }
	  ss >> ch;
	  if (!ss || ch != '}') 
#ifdef NOTHROW
	    { std::cerr<<err<<std::endl; exit(1); }
#else
	    throw std::runtime_error(err);
#endif
	  return temp;
	}
      } else {
	// Using "1. 2. 3." syntax
	std::stringstream ss(*this);
	std::vector<T> temp;
	T x;
	while (ss >> x) temp.push_back(x);
	return temp;
      }
    }

};

template <> inline ConvertibleString::operator bool() const
{
#ifdef Use_Zero_Default
  if (*this == "") return false;
#endif

  // make string all caps
  std::string sup = *this;
  for ( std::string::iterator p = sup.begin(); p != sup.end(); ++p )
    *p = toupper(*p); 

  if ( sup=="FALSE" || sup=="F" || sup=="NO" || sup=="N" ||
      sup=="0" || sup=="NONE" )
    return false;
  else if ( sup=="TRUE" || sup=="T" || sup=="YES" || sup=="Y" ||
      sup=="1" )
    return true;
  else {
    std::string err="Could not convert ConvertibleString to input type bool";
    err += std::string(": this = ") + *this;
#ifdef NOTHROW
    std::cerr<<err<<std::endl; exit(1);
    return false;
#else
    throw std::runtime_error(err);
#endif
  }
}

class ConfigFile 
{
  // Data
  protected:
    std::string myDelimiter;  // separator between key and value
    std::string myComment;    // separator between value and comments
    std::string mySentry;     // optional string to signal end of file
    std::map<std::string,ConvertibleString> myContents;   
    // extracted keys and values

    typedef std::map<std::string,ConvertibleString>::iterator mapi;
    typedef std::map<std::string,ConvertibleString>::const_iterator mapci;

    // Methods
  public:
    ConfigFile( std::string filename,
	std::string delimiter = "=",
	std::string comment = "#",
	std::string sentry = "EndConfigFile" );
    ConfigFile();

    void Load( std::string filename,
	std::string delimiter = "",
	std::string comment = "",
	std::string sentry = "" );

    // Search for key and read value or optional default value
    ConvertibleString& get( const std::string& key );
    ConvertibleString get( const std::string& key ) const;
    inline ConvertibleString& operator[]( const std::string& key )
    { return get(key); }
    inline ConvertibleString operator[]( const std::string& key ) const
    { return get(key); }

    template<class T> inline T read( const std::string& key ) const;
    template<class T> inline T read(
	const std::string& key, const T& value ) const;
    template<class T> inline bool readInto( 
	T& var, const std::string& key ) const;
    template<class T> inline bool readInto( 
	T& var, const std::string& key, const T& value ) const;

    // Modify keys and values
    template<class T> inline void add( std::string key, const T& value );
    void remove( const std::string& key );

    // Check whether key exists in configuration
    bool keyExists( const std::string& key ) const;

    // Check or change configuration syntax
    std::string getDelimiter() const { return myDelimiter; }
    std::string getComment() const { return myComment; }
    std::string getSentry() const { return mySentry; }
    std::string setDelimiter( const std::string& s )
    { std::string old = myDelimiter;  myDelimiter = s;  return old; }  
    std::string setComment( const std::string& s )
    { std::string old = myComment;  myComment = s;  return old; }

    // Write or read configuration
    void Write(std::ostream& os) const;
    void WriteAsComment(std::ostream& os) const;
    void Read(std::istream& is);
    void Append(const std::string& s)
    { std::istringstream ss(s); Read(ss); }

  protected:
    static void trim( std::string& s );

    // Exception types
  public:
    struct file_not_found : public std::runtime_error {
      file_not_found( const std::string& filename = std::string() ) throw() :
	std::runtime_error(std::string("ConfigFile error: file ") + filename
	    + " not found") {} 
    };

    struct key_not_found : public std::runtime_error { 
      // thrown only by T read(key) variant of read()
      key_not_found( const std::string& keyname = std::string() ) throw() :
	std::runtime_error(std::string("ConfigFile error: key ") + keyname
	    + " not found") {} 
    };
};

inline std::ostream& operator<<( std::ostream& os, const ConfigFile& cf )
{ cf.Write(os); return os; }
inline std::istream& operator>>( std::istream& is, ConfigFile& cf )
{ cf.Read(is); return is; }

template<class T>
T ConfigFile::read( const std::string& key ) const
{
  // Read the value corresponding to key
  mapci p = myContents.find(key);
  if( p == myContents.end() ) 
#ifdef NOTHROW
  { std::cerr<<"Key not found: "<<kay<<std::endl; exit(1); }
#else
    throw key_not_found(key);
#endif
  return p->second;
}


template<class T>
T ConfigFile::read( const std::string& key, const T& value ) const
{
  // Return the value corresponding to key or given default value
  // if key is not found
  mapci p = myContents.find(key);
  if( p == myContents.end() ) return value;
  return p->second;
}


template<class T>
bool ConfigFile::readInto( T& var, const std::string& key ) const
{
  // Get the value corresponding to key and store in var
  // Return true if key is found
  // Otherwise leave var untouched
  mapci p = myContents.find(key);
  bool found = ( p != myContents.end() );
  if( found ) var = static_cast<T>( p->second );
  return found;
}


template<class T>
bool ConfigFile::readInto( T& var, const std::string& key, const T& value ) const
{
  // Get the value corresponding to key and store in var
  // Return true if key is found
  // Otherwise set var to given default
  mapci p = myContents.find(key);
  bool found = ( p != myContents.end() );
  if( found ) var = p->second;
  else var = value;
  return found;
}

template<class T>
void ConfigFile::add( std::string key, const T& value )
{
  // Add a key with given value
  trim(key);
  myContents[key] = value;
}

#endif  // CONFIGFILE_H

// Release notes:
// v1.0  21 May 1999
//   + First release
//   + Template read() access only through non-member readConfigFile()
//   + ConfigurationFileBool is only built-in helper class
// 
// v2.0  3 May 2002
//   + Shortened name from ConfigurationFile to ConfigFile
//   + Implemented template member functions
//   + Changed default comment separator from % to #
//   + Enabled reading of multiple-line values
// 
// v2.1  24 May 2004
//   + Made template specializations inline to avoid compiler-dependent linkage
//   + Allowed comments within multiple-line values
//   + Enabled blank line termination for multiple-line values
//   + Added optional sentry to detect end of configuration file
//   + Rewrote messy trimWhitespace() function as elegant trim()
//
// v3.0  27 October 2006 by Mike Jarvis and Erin Sheldon
//   + Added operator[] access notation using ConvertibleString's
//   + Old read and add versions kept for backwards compatibility
//   + Added possibility of vector configuration values using either
//     whitespace delimited values or {,,,} notation.
