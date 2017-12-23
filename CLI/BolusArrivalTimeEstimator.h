#ifndef __BolusArrivalTimeEstimator_h
#define __BolusArrivalTimeEstimator_h

namespace BolusArrivalTime 
{

  class BolusArrivalTimeEstimator
  {
  public:
    BolusArrivalTimeEstimator(int defaultBolusArrivalTimeIndex);
    
    virtual ~BolusArrivalTimeEstimator() {}
    
    virtual int getBATIndex(int signalSize, const float* signal) = 0;

  protected:
    const int m_defaultBATIndex;

  };

}
#endif
