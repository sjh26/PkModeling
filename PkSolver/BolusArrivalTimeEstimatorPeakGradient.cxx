#include "BolusArrivalTimeEstimatorPeakGradient.h"
#include "PkSolver.h"
#include "SignalComputationUtils.h"

namespace BolusArrivalTime
{
  using namespace PkSolver;

  int BolusArrivalTimeEstimatorPeakGradient::getBATIndex(int signalSize, const float* signal, float* optRet_maxSlope /*= nullptr*/) const
  {
    if (signalSize <= 0) {
      throw NoSignalException();
    }

    int skipFront = 0;                  // Leading points to ignore
    int skipBack = 2;                   // Trailing points to ignore
    
    float* signalDerivative = new float[signalSize];
    memcpy(signalDerivative, signal, signalSize * sizeof(float));
    itk::compute_derivative(signalSize, signal, signalDerivative);

    int maxSlopeIdx = getMaxPositionInRange(skipFront, signalSize - skipBack, signalDerivative);
    int arrivalIdx = getArrivalIndex(skipFront, maxSlopeIdx, signalDerivative);
    
    if (optRet_maxSlope) {
      *optRet_maxSlope = signalDerivative[maxSlopeIdx];
    }
    delete[] signalDerivative;
    return arrivalIdx;
  }

  int BolusArrivalTimeEstimatorPeakGradient::getArrivalIndex(int start, int maxSlopeIdx, const float* signal) const
  {
    float thresh = signal[maxSlopeIdx] / 10.0;
    int arrivalTime = maxSlopeIdx;
    for (arrivalTime = maxSlopeIdx; arrivalTime >= start; arrivalTime--)
    {
      if (signal[arrivalTime] < thresh) {
        break;
      }
    }
    return arrivalTime + 1;
  }


}
