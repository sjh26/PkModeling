#include "ArterialInputFunctionPrescribed.h"

#include <sstream>
#include "StringUtils.h"
#include "Exceptions.h"
#include "SignalComputationUtils.h"
#include "IO/CSVReader.h"


ArterialInputFunctionPrescribed::ArterialInputFunctionPrescribed(const std::string& aifFileName,
                                                                const std::vector<float>& referenceTime)
                                                                : m_aifFileName(aifFileName),
                                                                  m_referenceTime(referenceTime)
{
  m_aif = loadResampledAIF();
}

std::vector<float> ArterialInputFunctionPrescribed::getSignalValues() const
{
  return m_aif;
}

unsigned int ArterialInputFunctionPrescribed::getSignalSize() const
{
  return m_aif.size();
}

std::vector<float> ArterialInputFunctionPrescribed::loadResampledAIF() const
{
  std::vector<float> aifTime;
  std::vector<float> aif;
  loadAifFromCsv(m_aifFileName, aifTime, aif);
  return SignalUtils::resampleSignal(aifTime, aif, m_referenceTime);
}

void ArterialInputFunctionPrescribed::loadAifFromCsv(const std::string& fileName, std::vector<float>& timing, std::vector<float>& aif) const
{
  timing.clear();
  aif.clear();

  CSVReader csvReader(fileName);
  while (csvReader.hasMoreRows())
  {
    std::vector<std::string> sValues = csvReader.nextRow();
    if (sValues.size() < 2) {
      continue;
    }
    
    float time = -1.0, value = -1.0;
    try {
      time = std::stof(sValues[0]);
      value = std::stof(sValues[1]);
    }
    catch (const std::invalid_argument& ia) {
      // not a float, probably the column labels, skip the row
      continue;
    }
    timing.push_back(time);
    aif.push_back(value);
  }
  if (timing.size() < 2) {
    throw WrongFileFormatException(fileName);
  }
}

