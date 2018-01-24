#ifndef __ArterialInputFunctionAverageUnderMask_h
#define __ArterialInputFunctionAverageUnderMask_h

#include "ArterialInputFunction.h"

#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageRegionConstIterator.h"


class ArterialInputFunctionAverageUnderMask : public ArterialInputFunction
{
public:
  typedef itk::VectorImage<float, 3> VectorVolume;
  typedef itk::Image<unsigned short, 3> MaskVolume;

  ArterialInputFunctionAverageUnderMask(const VectorVolume* inputVectorVolume, const MaskVolume* maskVolume);

  virtual ~ArterialInputFunctionAverageUnderMask() {}

  virtual std::vector<float> getSignalValues() const;
  virtual unsigned int getSignalSize() const;

private:
  typedef itk::ImageRegionConstIterator<VectorVolume> VectorVolumeConstIterator;
  typedef itk::ImageRegionConstIterator<MaskVolume> MaskVolumeConstIterator;
  typedef itk::VariableLengthVector<float> VectorVoxel;

  std::vector<float> computeAIF() const;
  inline std::vector<float> add(std::vector<float> vec1, const VectorVoxel& vec2) const;
  inline std::vector<float> divide(std::vector<float> vec, long div) const;

  const itk::VectorImage<float, 3>* const m_inputVectorVolume;
  const itk::Image<unsigned short, 3>* const m_maskVolume;
  std::vector<float> m_aif;
};

#endif
