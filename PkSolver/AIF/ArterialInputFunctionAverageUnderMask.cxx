#include "ArterialInputFunctionAverageUnderMask.h"

ArterialInputFunctionAverageUnderMask::ArterialInputFunctionAverageUnderMask(const itk::VectorImage<float, 3>* inputVectorVolume,
                                                                             const itk::Image<unsigned short, 3>* maskVolume)
                                                                             : m_inputVectorVolume(inputVectorVolume), 
                                                                               m_maskVolume(maskVolume)
{
  m_aif = computeAIF();
}

std::vector<float> ArterialInputFunctionAverageUnderMask::getAIF() const
{
  return m_aif;
}

std::vector<float> ArterialInputFunctionAverageUnderMask::computeAIF() const
{
  VectorVolumeConstIterator inputVectorVolumeIter(m_inputVectorVolume, m_inputVectorVolume->GetRequestedRegion());
  MaskVolumeConstIterator maskVolumeIter(m_maskVolume, m_maskVolume->GetRequestedRegion());
  inputVectorVolumeIter.GoToBegin();
  maskVolumeIter.GoToBegin();

  long numberVoxels = 0;
  std::vector<float> aif(m_inputVectorVolume->GetNumberOfComponentsPerPixel(), 0.0);

  while (!inputVectorVolumeIter.IsAtEnd())
  {
    if (maskVolumeIter.Get())
    {
      numberVoxels++;
      VectorVoxel vectorVoxel = inputVectorVolumeIter.Get();
      aif = add(aif, vectorVoxel);
    }
    ++maskVolumeIter;
    ++inputVectorVolumeIter;
  }

  return divide(aif, numberVoxels);
}

inline
std::vector<float> ArterialInputFunctionAverageUnderMask::add(std::vector<float> vec1, const VectorVoxel& vec2) const
{
  for (long i = 0; i < vec1.size(); i++)
  {
    vec1[i] += vec2[i];
  }
  return vec1;
}

inline
std::vector<float> ArterialInputFunctionAverageUnderMask::divide(std::vector<float> vec, long div) const
{
  for (long i = 0; i < vec.size(); i++)
  {
    vec[i] /= div;
  }
  return vec;
}



