#ifndef __ArterialInputFunctionPrescribed_h
#define __ArterialInputFunctionPrescribed_h

#include "ArterialInputFunction.h"

#include <string>


class ArterialInputFunctionPrescribed : public ArterialInputFunction
{
public:
  //! Creates an AIF from a CSV file, columns should be
  //! SignalTime,SignalValue
  //! Loaded AIF is resampled to provided referenceTime.
  //! Throws Exceptions if errors on loading.
  ArterialInputFunctionPrescribed(const std::string& aifFileName, const std::vector<float>& referenceTime);

  virtual ~ArterialInputFunctionPrescribed() {}

  virtual std::vector<float> getSignalValues() const;
  virtual unsigned int getSignalSize() const;

private:
  std::vector<float> loadResampledAIF() const;
  void loadAifFromCsv(const std::string& fileName, std::vector<float>& timing, std::vector<float>& aif) const;

  const std::string m_aifFileName;
  const std::vector<float> m_referenceTime;
  std::vector<float> m_aif;
};


#endif
