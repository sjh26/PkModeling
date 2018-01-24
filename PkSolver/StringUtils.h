#ifndef __StringUtils_h
#define __StringUtils_h

#include <string>
#include <sstream>
#include <vector>
#include <iterator>

namespace StringUtils
{

  std::vector<std::string> split(const std::string &s, char delim);

}
#endif
