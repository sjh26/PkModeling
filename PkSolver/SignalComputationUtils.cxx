#include "SignalComputationUtils.h"

namespace SignalUtils
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

  std::vector<float> resampleSignal(std::vector<float> signalTime, std::vector<float> signal, std::vector<float> referenceTime)
  {
    std::vector<float>::size_type timeSize = referenceTime.size();
    std::vector<float> resampledSignal(timeSize);

    std::vector<float>::iterator resampledSignalIter = resampledSignal.begin();
    std::vector<float>::iterator referenceTimeIter = referenceTime.begin();

    std::vector<float>::iterator signalIter = signal.begin();
    std::vector<float>::iterator signalTimeIter = signalTime.begin();

    std::vector<float>::iterator signalIterNext = signalIter;
    signalIterNext++;
    std::vector<float>::iterator signalTimeIterNext = signalTimeIter;
    signalTimeIterNext++;

    for (; referenceTimeIter != referenceTime.end(); ++referenceTimeIter, ++resampledSignalIter)
    {
      // Three cases
      // (1) extrapolate the signal on the low end of the range of prescribed timings
      // (2) interpolate the signal
      // (3) extrapolate the signal on the high end of the range of prescribed timings
      //
      // Case (1) is handled implictly by the initialization and conditionals.
      if (*signalTimeIter <= *referenceTimeIter)
      {
        // Case (2) from above)
        // find the prescribed times that straddle the current time to interpolate
        while (signalTimeIterNext != signalTime.end() && *signalTimeIterNext < *referenceTimeIter)
        {
          ++signalTimeIter;
          ++signalTimeIterNext;
          ++signalIter;
          ++signalIterNext;
        }
      }
      if (signalTimeIterNext == signalTime.end())
      {
        // we'll need to extrapolate (Case (3) from above)
        signalTimeIterNext = signalTimeIter;
        --signalTimeIter;
        signalIterNext = signalIter;
        --signalIter;
      }

      // interpolate signal;
      float a = *signalIter + ((*referenceTimeIter - *signalTimeIter) / (*signalTimeIterNext - *signalTimeIter)) * (*signalIterNext - *signalIter);
      *resampledSignalIter = a;
    }

    return resampledSignal;
  }

}
