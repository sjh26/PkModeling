#ifndef __CSVReader_h
#define __CSVReader_h

#include <string>
#include <vector>
#include <fstream>

class CSVReader
{
public:
  //! Init the CSV Reader.
  //! Throws Exception if file cannot be opened.
  CSVReader(const std::string& fileName);

  virtual ~CSVReader() {};

  bool hasMoreRows() const;
  std::vector<std::string> nextRow();

private:
  std::ifstream openFileStream(const std::string& fileName) const;
  void moveToNextValidRow();

  std::ifstream m_inputStream;
  std::string m_nextRow;
  bool m_hasMoreRows;
};

#endif
