#ifndef __SignalComputationUtils_h
#define __SignalComputationUtils_h

#include <vector>

namespace SignalUtils {
  int getMaxPosition(int signalSize, const float* signal);
  int getMaxPositionInRange(int start, int stop, const float* signal);
  std::vector<float> resampleSignal(std::vector<float> signalTime, std::vector<float> signal, std::vector<float> referenceTime);
}

#endif
