#ifndef __BolusArrivalTimeEstimatorPeakGradient_h
#define __BolusArrivalTimeEstimatorPeakGradient_h

#include "BolusArrivalTimeEstimator.h"

namespace BolusArrivalTime
{

  class BolusArrivalTimeEstimatorPeakGradient : public BolusArrivalTimeEstimator
  {
  public:
    BolusArrivalTimeEstimatorPeakGradient() {}

    virtual ~BolusArrivalTimeEstimatorPeakGradient() {}

    virtual int getBATIndex(int signalSize, const float* signal, float* optRet_maxSlope = nullptr) const;

  private:
    virtual int getArrivalIndex(int start, int maxSlopeIdx, const float* signal) const;

  };

}
#endif
