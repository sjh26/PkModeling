#ifndef __ArterialInputFunction_h
#define __ArterialInputFunction_h

#include <vector>

class ArterialInputFunction
{
public:
  ArterialInputFunction() {}
  
  virtual ~ArterialInputFunction() {}

  virtual std::vector<float> getSignalValues() const = 0;
  virtual unsigned int getSignalSize() const = 0;
};

#endif
