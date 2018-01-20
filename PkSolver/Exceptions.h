#ifndef __Exceptions_h
#define __Exceptions_h

#include <stdexcept>
#include <string>


class FileNotFoundException : public std::runtime_error
{
public:
  FileNotFoundException(const std::string& fileName) 
    : std::runtime_error("File \"" + fileName + "\" could not be opened.")
  {}
};

class WrongFileFormatException : public std::runtime_error
{
public:
  WrongFileFormatException(const std::string& fileName)
    : std::runtime_error("Content of file \"" + fileName + "\" has the wrong format.")
  {}
};

class ImageNullException : public std::runtime_error
{
public:
  ImageNullException(const std::string& imageName)
    : std::runtime_error(imageName + " does not exist.")
  {}
};

class FailedDictionaryLookup : public std::runtime_error
{
public:
  FailedDictionaryLookup(const std::string& key)
    : std::runtime_error("Failed lookup up key \"" + key + "\" in dictionary.")
  {}
};

#endif
