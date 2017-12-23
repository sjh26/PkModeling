#ifndef __BolusArrivalTimeEstimatorPeakGradient_h
#define __BolusArrivalTimeEstimatorPeakGradient_h

#include "BolusArrivalTimeEstimator.h"

namespace BolusArrivalTime
{

  class BolusArrivalTimeEstimatorPeakGradient : public BolusArrivalTimeEstimator
  {
  public:
    BolusArrivalTimeEstimatorPeakGradient(int defaultBolusArrivalTimeIndex);

    virtual ~BolusArrivalTimeEstimatorPeakGradient() {}

    virtual int getBATIndex(int signalSize, const float* signal);

  };

}
#endif
