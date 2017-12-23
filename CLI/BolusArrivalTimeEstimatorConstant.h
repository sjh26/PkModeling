#ifndef __BolusArrivalTimeEstimatorConstant_h
#define __BolusArrivalTimeEstimatorConstant_h

#include "BolusArrivalTimeEstimator.h"

namespace BolusArrivalTime
{

  class BolusArrivalTimeEstimatorConstant : public BolusArrivalTimeEstimator
  {
  public:
    BolusArrivalTimeEstimatorConstant(int defaultBolusArrivalTimeIndex) : BolusArrivalTimeEstimator(defaultBolusArrivalTimeIndex) {}

    virtual ~BolusArrivalTimeEstimatorConstant() {}

    virtual int getBATIndex(int signalSize, const float* signal)
    {
      return m_defaultBATIndex;
    }

  };

}
#endif
