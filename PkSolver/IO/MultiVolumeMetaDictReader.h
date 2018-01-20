#ifndef __MultiVolumeMetaDictReader_h
#define __MultiVolumeMetaDictReader_h

#include "itkMetaDataDictionary.h"

class MultiVolumeMetaDictReader
{
public:
  MultiVolumeMetaDictReader(const itk::MetaDataDictionary& dictionary);

  virtual ~MultiVolumeMetaDictReader() {}

  float get(const std::string& key);
  std::vector<float> getTiming();

private:
  const itk::MetaDataDictionary& m_dictionary;
};

#endif
