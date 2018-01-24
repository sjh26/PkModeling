#ifndef __BolusArrivalTimeEstimator_h
#define __BolusArrivalTimeEstimator_h

#include <stddef.h>

namespace BolusArrivalTime 
{

  class BolusArrivalTimeEstimator
  {
  public:
    BolusArrivalTimeEstimator() {}
    
    virtual ~BolusArrivalTimeEstimator() {}
    
    virtual int getBATIndex(int signalSize, const float* signal, float* optRet_maxSlope = NULL) const = 0;

  };

}
#endif
