#ifndef __BolusArrivalTimeEstimatorConstant_h
#define __BolusArrivalTimeEstimatorConstant_h

#include "BolusArrivalTimeEstimator.h"

namespace BolusArrivalTime
{

  class BolusArrivalTimeEstimatorConstant : public BolusArrivalTimeEstimator
  {
  public:
    BolusArrivalTimeEstimatorConstant(int defaultBolusArrivalTimeIndex) : m_defaultBATIndex(defaultBolusArrivalTimeIndex) {}

    virtual ~BolusArrivalTimeEstimatorConstant() {}

    virtual int getBATIndex(int signalSize, const float* signal, float* optRet_maxSlope = NULL) const
    {
      if (optRet_maxSlope) {
        *optRet_maxSlope = 0.0;
      }
      return m_defaultBATIndex;
    }

  private:
    const int m_defaultBATIndex;

  };

}
#endif
