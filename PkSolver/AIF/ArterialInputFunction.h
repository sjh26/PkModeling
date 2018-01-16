#ifndef __ArterialInputFunction_h
#define __ArterialInputFunction_h

#include <vector>

class ArterialInputFunction
{
public:
  ArterialInputFunction() {}
  
  virtual ~ArterialInputFunction() {}

  virtual std::vector<float> getAIF() const = 0;
};

#endif
