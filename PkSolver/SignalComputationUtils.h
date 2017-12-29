#ifndef __SignalComputationUtils_h
#define __SignalComputationUtils_h

namespace PkSolver {
  int getMaxPosition(int signalSize, const float* signal);
  int getMaxPositionInRange(int start, int stop, const float* signal);
}

#endif
