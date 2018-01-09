#include "SignalComputationUtils.h"

namespace PkSolver
{

  int getMaxPosition(int signalSize, const float* signal)
  {
    return getMaxPositionInRange(0, signalSize, signal);
  }

  int getMaxPositionInRange(int start, int stop, const float* signal)
  {
    float max = signal[start];
    int maxIndex = start;
    for (int i = start; i < stop; i++)
      if (signal[i] > max)
      {
        max = signal[i];
        maxIndex = i;
      }
    return maxIndex;
  }

}
