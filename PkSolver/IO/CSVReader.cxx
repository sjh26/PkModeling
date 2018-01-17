#include "CSVReader.h"

#include <sstream>
#include "StringUtils.h"
#include "Exceptions.h"

CSVReader::CSVReader(const std::string& fileName)
{
  m_inputStream = openFileStream(fileName);
  moveToNextValidRow();
}

bool CSVReader::hasMoreRows() const
{
  return m_hasMoreRows;
}

std::vector<std::string> CSVReader::nextRow()
{
  std::vector<std::string> rowEntries = StringUtils::split(m_nextRow, ',');
  moveToNextValidRow();
  return rowEntries;
}

std::ifstream CSVReader::openFileStream(const std::string& fileName) const
{
  std::ifstream fileStream;
  fileStream.open(fileName.c_str());
  if (fileStream.fail())
  {
    throw FileNotFoundException(fileName);
  }
  return fileStream;
}

void CSVReader::moveToNextValidRow()
{
  m_hasMoreRows = false;
  m_nextRow = "";
  while (!m_inputStream.eof() && !m_hasMoreRows)
  {
    std::string line;
    std::getline(m_inputStream, line);
    if (line[0] == '#') {
      continue;
    }
    m_hasMoreRows = true;
    m_nextRow = line;
  }
  if (!m_hasMoreRows) {
    m_inputStream.close();
  }
}

