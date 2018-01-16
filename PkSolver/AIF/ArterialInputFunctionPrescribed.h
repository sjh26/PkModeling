#ifndef __ArterialInputFunctionPrescribed_h
#define __ArterialInputFunctionPrescribed_h

#include "ArterialInputFunction.h"

class ArterialInputFunctionPrescribed : public ArterialInputFunction
{
public:
  ArterialInputFunctionPrescribed(const std::string& aifFileName, const std::vector<float>& referenceTime);

  virtual ~ArterialInputFunctionPrescribed() {}

  virtual std::vector<float> getAIF() const;

private:
  std::vector<float> loadResampledAIF() const;
  void loadAifFromCsv(const std::string& fileName, std::vector<float>& timing, std::vector<float>& aif) const;

  const std::string m_aifFileName;
  const std::vector<float> m_referenceTime;
  std::vector<float> m_aif;
};


#endif
