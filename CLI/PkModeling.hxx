#ifndef __PkModeling_h
#define __PkModeling_h

#include "Configuration.h"

#include "itkMetaDataObject.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkMultiThreader.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "itkPluginUtilities.h"

#include "itkSignalIntensityToConcentrationImageFilter.h"
#include "itkConcentrationToQuantitativeImageFilter.h"

#include "AIF/ArterialInputFunctionPrescribed.h"
#include "AIF/ArterialInputFunctionPopulation.h"
#include "AIF/ArterialInputFunctionAverageUnderMask.h"

#include "BAT/BolusArrivalTimeEstimator.h"
#include "BAT/BolusArrivalTimeEstimatorConstant.h"
#include "BAT/BolusArrivalTimeEstimatorPeakGradient.h"

#include "IO/MultiVolumeMetaDictReader.h"

#include <sstream>
#include <fstream>
#include <memory>



class PkModeling {
// Helper typedefs to make the code easier to read with all the templates
private:
  static const unsigned int VectorVolumeDimension = 3;
  typedef itk::VectorImage<float, VectorVolumeDimension> VectorVolumeType;
  typedef itk::ImageFileReader<VectorVolumeType>         VectorVolumeReaderType;

  static const unsigned int MaskVolumeDimension = 3;
  typedef itk::Image<unsigned short, MaskVolumeDimension> MaskVolumeType;
  typedef itk::ImageFileReader<MaskVolumeType>            MaskVolumeReaderType;

  typedef itk::Image<float, VectorVolumeDimension> OutputVolumeType;

  typedef itk::ResampleImageFilter<MaskVolumeType, MaskVolumeType>     ResamplerType;
  typedef itk::NearestNeighborInterpolateImageFunction<MaskVolumeType> InterpolatorType;

  typedef itk::SignalIntensityToConcentrationImageFilter<VectorVolumeType, MaskVolumeType, VectorVolumeType> ConvertFilterType;
  typedef itk::ConcentrationToQuantitativeImageFilter<VectorVolumeType, MaskVolumeType, OutputVolumeType>    QuantifierType;

// Member Variables
private:
  const Configuration m_config;

  // Input Data
  VectorVolumeType::Pointer m_inputVectorVolume;
  MaskVolumeType::Pointer m_aifMaskVolume;
  MaskVolumeType::Pointer m_T1MapVolume;
  MaskVolumeType::Pointer m_roiMaskVolume;
  std::unique_ptr<MultiVolumeMetaDictReader> m_imageMetaDict;

  // Filters
  ConvertFilterType::Pointer m_signalToConcentrationsConverter;
  QuantifierType::Pointer m_concentrationsToQuantitativeImageFilter;

  // Computation Strategies for Filters
  std::unique_ptr<BolusArrivalTime::BolusArrivalTimeEstimator> m_batEstimator;
  std::unique_ptr<ArterialInputFunction> m_aif;

  // Progress Watchers
  std::vector<itk::PluginFilterWatcher> m_progressWatchers;


// Public interface
public:
  PkModeling(Configuration config) : m_config(config) {}
  virtual ~PkModeling() {}

  int execute()
  {
    initialize();
    setupProcessingPipeline();
    runProcessingPipeline();
    writeResults();
    return EXIT_SUCCESS;
  }

// Internal functions
private:
  void initialize()
  {
    m_inputVectorVolume = getVectorVolume(m_config.InputFourDImageFileName);
    m_batEstimator = getBatEstimator();

    m_imageMetaDict.reset(new MultiVolumeMetaDictReader(m_inputVectorVolume->GetMetaDataDictionary()));

    m_aifMaskVolume = getMaskVolumeOrNull(m_config.AIFMaskFileName);
    m_T1MapVolume = getMaskVolumeOrNull(m_config.T1MapFileName);
    m_roiMaskVolume = getResampledMaskVolumeOrNull(m_config.ROIMaskFileName, m_inputVectorVolume);
  }

  void setupProcessingPipeline()
  {
    setupSignalToConcentrationsConverter();
    setupAIF();
    setupConcentrationsToQuantitativeImageFilter();
  }

  void runProcessingPipeline()
  {
    m_concentrationsToQuantitativeImageFilter->Update();
  }

  void writeResults()
  {
    writeMultiVolumeIfFileNameValid(m_config.OutputConcentrationsImageFileName, m_signalToConcentrationsConverter->GetOutput(), m_inputVectorVolume);
    writeMultiVolumeIfFileNameValid(m_config.OutputFittedDataImageFileName, m_concentrationsToQuantitativeImageFilter->GetFittedDataOutput(), m_inputVectorVolume);

    writeVolumeIfFileNameValid(m_config.OutputKtransFileName, m_concentrationsToQuantitativeImageFilter->GetKTransOutput());
    writeVolumeIfFileNameValid(m_config.OutputVeFileName, m_concentrationsToQuantitativeImageFilter->GetVEOutput());
    writeVolumeIfFileNameValid(m_config.OutputMaxSlopeFileName, m_concentrationsToQuantitativeImageFilter->GetMaxSlopeOutput());
    writeVolumeIfFileNameValid(m_config.OutputAUCFileName, m_concentrationsToQuantitativeImageFilter->GetAUCOutput());
    writeVolumeIfFileNameValid(m_config.OutputRSquaredFileName, m_concentrationsToQuantitativeImageFilter->GetRSquaredOutput());
    writeVolumeIfFileNameValid(m_config.OutputBolusArrivalTimeImageFileName, m_concentrationsToQuantitativeImageFilter->GetBATOutput());
    writeVolumeIfFileNameValid(m_config.OutputOptimizerDiagnosticsImageFileName, m_concentrationsToQuantitativeImageFilter->GetOptimizerDiagnosticsOutput());

    if (m_config.ComputeFpv) {
      writeVolumeIfFileNameValid(m_config.OutputFpvFileName, m_concentrationsToQuantitativeImageFilter->GetFPVOutput());
    }
  }

  void setupSignalToConcentrationsConverter()
  {
    m_signalToConcentrationsConverter = ConvertFilterType::New();
    m_signalToConcentrationsConverter->SetInput(m_inputVectorVolume);
    m_signalToConcentrationsConverter->SetROIMask(m_roiMaskVolume);
    m_signalToConcentrationsConverter->SetT1PreBlood(m_config.T1PreBloodValue);
    m_signalToConcentrationsConverter->SetT1PreTissue(m_config.T1PreTissueValue);
    m_signalToConcentrationsConverter->SetTR(m_imageMetaDict->get("MultiVolume.DICOM.RepetitionTime"));
    m_signalToConcentrationsConverter->SetFA(m_imageMetaDict->get("MultiVolume.DICOM.FlipAngle"));
    m_signalToConcentrationsConverter->SetBatEstimator(m_batEstimator.get());
    m_signalToConcentrationsConverter->SetRGD_relaxivity(m_config.RelaxivityValue);
    m_signalToConcentrationsConverter->SetS0GradThresh(m_config.S0GradValue);
    m_signalToConcentrationsConverter->SetT1Map(m_T1MapVolume);
    if (m_config.AIFMode == "AverageUnderAIFMask") {
      m_signalToConcentrationsConverter->SetAIFMask(m_aifMaskVolume);
    }

    m_progressWatchers.push_back(itk::PluginFilterWatcher(m_signalToConcentrationsConverter, "Concentrations", m_config.CLPProcessInformation, 1.0 / 20.0, 0.0));
  }

  void setupAIF()
  {
    if (m_config.AIFMode == "Prescribed")
    {
      m_aif.reset(new ArterialInputFunctionPrescribed(m_config.PrescribedAIFFileName, m_imageMetaDict->getTiming()));
    }
    else if (m_config.AIFMode == "Population")
    {
      m_aif.reset(new ArterialInputFunctionPopulation(m_imageMetaDict->getTiming()));
    }
    else
    {
      m_aif.reset(new ArterialInputFunctionAverageUnderMask(m_signalToConcentrationsConverter->GetOutput(), m_aifMaskVolume));
    }
  }

  void setupConcentrationsToQuantitativeImageFilter()
  {
    m_concentrationsToQuantitativeImageFilter = QuantifierType::New();
    m_concentrationsToQuantitativeImageFilter->SetInput(m_signalToConcentrationsConverter->GetOutput());
    m_concentrationsToQuantitativeImageFilter->SetAIF(m_aif.get());
    m_concentrationsToQuantitativeImageFilter->SetAUCTimeInterval(m_config.AUCTimeInterval);
    m_concentrationsToQuantitativeImageFilter->SetTiming(m_imageMetaDict->getTiming());
    m_concentrationsToQuantitativeImageFilter->SetfTol(m_config.FTolerance);
    m_concentrationsToQuantitativeImageFilter->SetgTol(m_config.GTolerance);
    m_concentrationsToQuantitativeImageFilter->SetxTol(m_config.XTolerance);
    m_concentrationsToQuantitativeImageFilter->Setepsilon(m_config.Epsilon);
    m_concentrationsToQuantitativeImageFilter->SetmaxIter(m_config.MaxIter);
    m_concentrationsToQuantitativeImageFilter->Sethematocrit(m_config.Hematocrit);
    m_concentrationsToQuantitativeImageFilter->SetBatEstimator(m_batEstimator.get());
    m_concentrationsToQuantitativeImageFilter->SetROIMask(m_roiMaskVolume);
    if (m_config.ComputeFpv) {
      m_concentrationsToQuantitativeImageFilter->SetModelType(itk::LMCostFunction::TOFTS_3_PARAMETER);
    }
    else {
      m_concentrationsToQuantitativeImageFilter->SetModelType(itk::LMCostFunction::TOFTS_2_PARAMETER);
    }

    m_progressWatchers.push_back(itk::PluginFilterWatcher(m_concentrationsToQuantitativeImageFilter, "Quantifying", m_config.CLPProcessInformation, 19.0 / 20.0, 1.0 / 20.0));
  }

  MaskVolumeType::Pointer getMaskVolumeOrNull(const std::string& maskFileName)
  {
    MaskVolumeReaderType::Pointer maskVolumeReader = MaskVolumeReaderType::New();
    MaskVolumeType::Pointer maskVolume = NULL;
    try
    {
      maskVolumeReader->SetFileName(maskFileName.c_str());
      maskVolumeReader->Update();
      maskVolume = maskVolumeReader->GetOutput();
    }
    catch (...) 
    { }
    return maskVolume;
  }

  MaskVolumeType::Pointer getResampledMaskVolumeOrNull(const std::string& maskFileName, const VectorVolumeType::Pointer referenceVolume)
  {
    MaskVolumeType::Pointer maskVolume = getMaskVolumeOrNull(maskFileName);
    if (maskVolume.IsNotNull()) {
      ResamplerType::Pointer resampler = ResamplerType::New();
      InterpolatorType::Pointer interpolator = InterpolatorType::New();

      resampler->SetOutputDirection(referenceVolume->GetDirection());
      resampler->SetOutputSpacing(referenceVolume->GetSpacing());
      resampler->SetOutputStartIndex(referenceVolume->GetBufferedRegion().GetIndex());
      resampler->SetSize(referenceVolume->GetBufferedRegion().GetSize());
      resampler->SetOutputOrigin(referenceVolume->GetOrigin());
      resampler->SetInput(maskVolume);
      resampler->SetInterpolator(interpolator);
      resampler->Update();

      maskVolume = resampler->GetOutput();
    }
    return maskVolume;
  }

  VectorVolumeType::Pointer getVectorVolume(const std::string& volumeFileName)
  {
    VectorVolumeReaderType::Pointer multiVolumeReader = VectorVolumeReaderType::New();
    multiVolumeReader->SetFileName(volumeFileName.c_str());
    multiVolumeReader->Update();
    return multiVolumeReader->GetOutput();
  }

  std::unique_ptr<BolusArrivalTime::BolusArrivalTimeEstimator> getBatEstimator()
  {
    std::unique_ptr<BolusArrivalTime::BolusArrivalTimeEstimator> batEstimator(new BolusArrivalTime::BolusArrivalTimeEstimatorConstant(m_config.ConstantBAT));
    if (m_config.BATCalculationMode == "PeakGradient") {
      batEstimator.reset(new BolusArrivalTime::BolusArrivalTimeEstimatorPeakGradient());
    }
    return batEstimator;
  }

  template <typename TOutVolume>
  void writeVolumeIfFileNameValid(std::string fileName, const TOutVolume* outVolume)
  {
    try
    {
      itk::ImageFileWriter<TOutVolume>::Pointer volumeWriter = itk::ImageFileWriter<TOutVolume>::New();
      volumeWriter->SetFileName(fileName.c_str());
      volumeWriter->SetInput(outVolume);
      volumeWriter->SetUseCompression(1);
      volumeWriter->Update();
    }
    catch (...)
    { }
  }

  void writeMultiVolumeIfFileNameValid(std::string fileName, const VectorVolumeType::Pointer outVolume, const VectorVolumeType::Pointer referenceVolume)
  {
    // this line is needed to make Slicer recognize this as a VectorVolume and not a MultiVolume
    outVolume->SetMetaDataDictionary(referenceVolume->GetMetaDataDictionary());
    writeVolumeIfFileNameValid(fileName, outVolume.GetPointer());
  }
};


#endif
