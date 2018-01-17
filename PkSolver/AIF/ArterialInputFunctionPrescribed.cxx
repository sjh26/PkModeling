#include "ArterialInputFunctionPrescribed.h"

#include <sstream>
#include <fstream>
#include "StringUtils.h"
#include "Exceptions.h"
#include "SignalComputationUtils.h"


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

  std::string line;
  std::ifstream csv;
  csv.open(fileName.c_str());
  if (csv.fail())
  {
    throw FileNotFoundException(fileName);
  }

  while (!csv.eof())
  {
    getline(csv, line);

    if (line[0] == '#')
    {
      continue;
    }

    std::vector<std::string> svalues = StringUtils::split(line, ',');

    if (svalues.size() < 2)
    {
      // not enough values on the line
      continue;
    }

    float time = -1.0, value = -1.0;
    try {
      time = std::stof(svalues[0]);
      value = std::stof(svalues[1]);
    }
    catch (const std::invalid_argument& ia) {
      // not a float, probably the column labels, skip the row
      continue;
    }

    timing.push_back(time);
    aif.push_back(value);
  }

  csv.close();

  if (timing.size() == 0)
  {
    throw WrongFileFormatException(fileName);
  }
}

