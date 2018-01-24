#ifndef __ArterialInputFunctionPopulation_h
#define __ArterialInputFunctionPopulation_h

#include "ArterialInputFunction.h"

//! Provides a generated population AIF.
//
//! See "Experimentally-Derived Functional Form for a Population-Averaged High-
//! Temporal-Resolution Arterial Input Function for Dynamic Contrast-Enhanced
//! MRI" - Parker, Robers, Macdonald, Buonaccorsi, Cheung, Buckley, Jackson,
//! Watson, Davies, Jayson.  Magnetic Resonance in Medicine 56:993-1000 (2006)
class ArterialInputFunctionPopulation : public ArterialInputFunction
{
public:
  //! Initialize the AIF generator.
  //! referenceSignalTime : sequence time, presumed in units of seconds.
  //! bolusArrivalTimeFraction : fractional point between 0 and 1 when the bolus is
  //!     desired to arrive.  Choose 0.0 to have it at the very beginning,
  //!     1.0 to have it at the end.
  ArterialInputFunctionPopulation(const std::vector<float>& referenceSignalTime, const float bolusArrivalTimeFraction = 0.1);

  virtual ~ArterialInputFunctionPopulation() {}

  virtual std::vector<float> getSignalValues() const;
  virtual unsigned int getSignalSize() const;

private:
  std::vector<float> computeAIF() const;

  const std::vector<float> m_referenceSignalTime;
  const float m_bolusArrivalTimeFraction;
  std::vector<float> m_aif;
};

#endif
