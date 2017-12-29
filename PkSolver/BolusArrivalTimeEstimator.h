#ifndef __BolusArrivalTimeEstimator_h
#define __BolusArrivalTimeEstimator_h

namespace BolusArrivalTime 
{

  class BolusArrivalTimeEstimator
  {
  public:
    BolusArrivalTimeEstimator() {}
    
    virtual ~BolusArrivalTimeEstimator() {}
    
    virtual int getBATIndex(int signalSize, const float* signal, float* optRet_maxSlope = nullptr) const = 0;

  };

}
#endif
