#include "BolusArrivalTimeEstimatorPeakGradient.h"
#include "PkSolver.h"

namespace BolusArrivalTime
{
  BolusArrivalTimeEstimatorPeakGradient::BolusArrivalTimeEstimatorPeakGradient(int defaultBolusArrivalTimeIndex) : BolusArrivalTimeEstimator(defaultBolusArrivalTimeIndex) {}

  int BolusArrivalTimeEstimatorPeakGradient::getBATIndex(int signalSize, const float* signal)
  {
    return 0;
  }


}
