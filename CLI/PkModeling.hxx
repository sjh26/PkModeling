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

#include <sstream>
#include <fstream>
#include <memory>



class PkModeling {
private:
  const Configuration m_config;

public:
  static const unsigned int VectorVolumeDimension = 3;
  typedef itk::VectorImage<float, VectorVolumeDimension>     VectorVolumeType;
  typedef itk::VectorImage<float, VectorVolumeDimension>     FloatVectorVolumeType;
  typedef VectorVolumeType::RegionType              VectorVolumeRegionType;
  typedef itk::ImageFileReader<VectorVolumeType>             VectorVolumeReaderType;
  typedef itk::ImageFileWriter<FloatVectorVolumeType>        VectorVolumeWriterType;

  static const unsigned int MaskVolumeDimension = 3;
  typedef itk::Image<unsigned short, MaskVolumeDimension>      MaskVolumeType;
  typedef itk::ImageFileReader<MaskVolumeType>                 MaskVolumeReaderType;

  typedef itk::Image<float, VectorVolumeDimension> OutputVolumeType;
  typedef itk::ImageFileWriter< OutputVolumeType> OutputVolumeWriterType;

  typedef itk::ResampleImageFilter<MaskVolumeType, MaskVolumeType> ResamplerType;
  typedef itk::NearestNeighborInterpolateImageFunction<MaskVolumeType> InterpolatorType;

  typedef itk::SignalIntensityToConcentrationImageFilter<VectorVolumeType, MaskVolumeType, FloatVectorVolumeType> ConvertFilterType;
  typedef itk::ConcentrationToQuantitativeImageFilter<FloatVectorVolumeType, MaskVolumeType, OutputVolumeType> QuantifierType;


  PkModeling(Configuration config) : m_config(config) {}
  virtual ~PkModeling() {}

#define SimpleAttributeGetMethodMacro(name, key, type)     \
type Get##name(itk::MetaDataDictionary& dictionary)           \
  {\
  type value = type(); \
  if (dictionary.HasKey(key))\
        {\
    /* attributes stored as strings */ \
    std::string valueString; \
    itk::ExposeMetaData(dictionary, key, valueString);  \
    std::stringstream valueStream(valueString); \
    valueStream >> value; \
        }\
        else\
      {\
    itkGenericExceptionMacro("Missing attribute '" key "'.");\
      }\
  return value;\
  }

  //SimpleAttributeGetMethodMacro(EchoTime, "MultiVolume.DICOM.EchoTime", float);
  SimpleAttributeGetMethodMacro(RepetitionTime, "MultiVolume.DICOM.RepetitionTime", float);
  SimpleAttributeGetMethodMacro(FlipAngle, "MultiVolume.DICOM.FlipAngle", float);

  std::vector<float> GetTiming(itk::MetaDataDictionary& dictionary)
  {
    std::vector<float> triggerTimes;

    if (dictionary.HasKey("MultiVolume.FrameIdentifyingDICOMTagName"))
    {
      std::string tag;
      itk::ExposeMetaData(dictionary, "MultiVolume.FrameIdentifyingDICOMTagName", tag);
      if (dictionary.HasKey("MultiVolume.FrameLabels"))
      {
        // Acquisition parameters stored as text, FrameLabels are comma separated
        std::string frameLabelsString;
        itk::ExposeMetaData(dictionary, "MultiVolume.FrameLabels", frameLabelsString);
        std::stringstream frameLabelsStream(frameLabelsString);
        if (tag == "TriggerTime" || tag == "AcquisitionTime" || tag == "SeriesTime" || tag == "ContentTime")
        {
          float t;
          float t0 = 0.0;
          bool first = true;
          while (frameLabelsStream >> t)
          {
            t /= 1000.0;  // convert to seconds
            if (first)
            {
              t0 = t;
              first = false;
            }
            t = t - t0;

            triggerTimes.push_back(t);
            frameLabelsStream.ignore(1); // skip the comma
          }
        }
        else
        {
          itkGenericExceptionMacro("Unrecognized frame identifying DICOM tag name " << tag);
        }
        // what other frame identification methods are there?
      }
      else
      {
        itkGenericExceptionMacro("Missing attribute 'MultiVolume.FrameLabels'.")
      }
    }
    else
    {
      itkGenericExceptionMacro("Missing attribute 'MultiVolume.FrameIdentifyingDICOMTagName'.");
    }

    return triggerTimes;
  }


  MaskVolumeType::Pointer getMaskVolumeOrNull(const std::string& maskFileName)
  {
    MaskVolumeReaderType::Pointer maskVolumeReader = MaskVolumeReaderType::New();
    MaskVolumeType::Pointer maskVolume = NULL;
    if (maskFileName != "")
    {
      maskVolumeReader->SetFileName(maskFileName.c_str());
      maskVolumeReader->Update();
      maskVolume = maskVolumeReader->GetOutput();
    }
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

  std::unique_ptr<ArterialInputFunction> getAIF(const std::vector<float>& Timing, 
                                                const VectorVolumeType::Pointer inputVectorVolume, 
                                                const MaskVolumeType::Pointer aifMaskVolume)
  {
    std::unique_ptr<ArterialInputFunction> aif;
    if (m_config.AIFMode == "Prescribed")
    {
      aif.reset(new ArterialInputFunctionPrescribed(m_config.PrescribedAIFFileName, Timing));
    }
    else if (m_config.AIFMode == "Population")
    {
      aif.reset(new ArterialInputFunctionPopulation(Timing));
    }
    else
    {
      aif.reset(new ArterialInputFunctionAverageUnderMask(inputVectorVolume, aifMaskVolume));
    }
    return aif;
  }

  void writeVectorVolumeIfFileNameValid(std::string fileName, const FloatVectorVolumeType::Pointer outVolume, const VectorVolumeType::Pointer referenceVolume)
  {
    if (!fileName.empty())
    {
      // this line is needed to make Slicer recognize this as a VectorVolume and not a MultiVolume
      outVolume->SetMetaDataDictionary(referenceVolume->GetMetaDataDictionary());

      VectorVolumeWriterType::Pointer multiVolumeWriter = VectorVolumeWriterType::New();
      multiVolumeWriter->SetFileName(fileName.c_str());
      multiVolumeWriter->SetInput(outVolume);
      multiVolumeWriter->SetUseCompression(1);
      multiVolumeWriter->Update();
    }
  }

  void writeVolumeIfFileNameValid(std::string fileName, const OutputVolumeType::Pointer outVolume)
  {
    if (!fileName.empty())
    {
      OutputVolumeWriterType::Pointer volumeWriter = OutputVolumeWriterType::New();
      volumeWriter->SetInput(outVolume);
      volumeWriter->SetFileName(fileName.c_str());
      volumeWriter->SetUseCompression(1);
      volumeWriter->Update();
    }
  }

  void saveRequestedOutputs(const ConvertFilterType::Pointer signalToConcentrationsConverter, 
                            const QuantifierType::Pointer concentrationsToQuantitativeImageFilter,
                            const VectorVolumeType::Pointer inputVectorVolume)
  {
    writeVectorVolumeIfFileNameValid(m_config.OutputConcentrationsImageFileName, signalToConcentrationsConverter->GetOutput(), inputVectorVolume);
    writeVectorVolumeIfFileNameValid(m_config.OutputFittedDataImageFileName, concentrationsToQuantitativeImageFilter->GetFittedDataOutput(), inputVectorVolume);

    writeVolumeIfFileNameValid(m_config.OutputKtransFileName, concentrationsToQuantitativeImageFilter->GetKTransOutput());
    writeVolumeIfFileNameValid(m_config.OutputVeFileName, concentrationsToQuantitativeImageFilter->GetVEOutput());
    writeVolumeIfFileNameValid(m_config.OutputMaxSlopeFileName, concentrationsToQuantitativeImageFilter->GetMaxSlopeOutput());
    writeVolumeIfFileNameValid(m_config.OutputAUCFileName, concentrationsToQuantitativeImageFilter->GetAUCOutput());
    writeVolumeIfFileNameValid(m_config.OutputRSquaredFileName, concentrationsToQuantitativeImageFilter->GetRSquaredOutput());
    writeVolumeIfFileNameValid(m_config.OutputBolusArrivalTimeImageFileName, concentrationsToQuantitativeImageFilter->GetBATOutput());
    writeVolumeIfFileNameValid(m_config.OutputOptimizerDiagnosticsImageFileName, concentrationsToQuantitativeImageFilter->GetOptimizerDiagnosticsOutput());

    if (m_config.ComputeFpv) {
      writeVolumeIfFileNameValid(m_config.OutputFpvFileName, concentrationsToQuantitativeImageFilter->GetFPVOutput());
    }
  }


  int execute()
  {
    VectorVolumeType::Pointer inputVectorVolume = getVectorVolume(m_config.InputFourDImageFileName);
    std::unique_ptr<BolusArrivalTime::BolusArrivalTimeEstimator> batEstimator = getBatEstimator();

    std::vector<float> Timing;
    float FAValue = 0.0;
    float TRValue = 0.0;
    try
    {
      Timing = GetTiming(inputVectorVolume->GetMetaDataDictionary());
      FAValue = GetFlipAngle(inputVectorVolume->GetMetaDataDictionary());
      TRValue = GetRepetitionTime(inputVectorVolume->GetMetaDataDictionary());
    }
    catch (itk::ExceptionObject &exc)
    {
      itkGenericExceptionMacro(<< exc.GetDescription()
        << " Image " << m_config.InputFourDImageFileName.c_str()
        << " does not contain sufficient attributes to support algorithms.");
    }

    //Read masks
    MaskVolumeType::Pointer aifMaskVolume = getMaskVolumeOrNull(m_config.AIFMaskFileName);
    MaskVolumeType::Pointer T1MapVolume = getMaskVolumeOrNull(m_config.T1MapFileName);
    MaskVolumeType::Pointer roiMaskVolume = getResampledMaskVolumeOrNull(m_config.ROIMaskFileName, inputVectorVolume);


    /////////////////////////////// PROCESSING /////////////////////

    //Convert to concentration values
    ConvertFilterType::Pointer signalToConcentrationsConverter = ConvertFilterType::New();
    signalToConcentrationsConverter->SetInput(inputVectorVolume);

    if (m_config.AIFMode == "AverageUnderAIFMask")
    {
      signalToConcentrationsConverter->SetAIFMask(aifMaskVolume);
    }

    if (m_config.ROIMaskFileName != "")
    {
      signalToConcentrationsConverter->SetROIMask(roiMaskVolume);
    }

    signalToConcentrationsConverter->SetT1PreBlood(m_config.T1PreBloodValue);
    signalToConcentrationsConverter->SetT1PreTissue(m_config.T1PreTissueValue);
    signalToConcentrationsConverter->SetTR(TRValue);
    signalToConcentrationsConverter->SetFA(FAValue);
    signalToConcentrationsConverter->SetBatEstimator(batEstimator.get());
    signalToConcentrationsConverter->SetRGD_relaxivity(m_config.RelaxivityValue);
    signalToConcentrationsConverter->SetS0GradThresh(m_config.S0GradValue);

    if (m_config.T1MapFileName != "")
    {
      signalToConcentrationsConverter->SetT1Map(T1MapVolume);
    }

    itk::PluginFilterWatcher watchConverter(signalToConcentrationsConverter, "Concentrations", m_config.CLPProcessInformation, 1.0 / 20.0, 0.0);
    signalToConcentrationsConverter->Update();

    std::unique_ptr<ArterialInputFunction> aif = getAIF(Timing, signalToConcentrationsConverter->GetOutput(), aifMaskVolume);

    //Calculate parameters
    QuantifierType::Pointer concentrationsToQuantitativeImageFilter = QuantifierType::New();
    concentrationsToQuantitativeImageFilter->SetInput(signalToConcentrationsConverter->GetOutput());
    concentrationsToQuantitativeImageFilter->SetAIF(aif.get());
    concentrationsToQuantitativeImageFilter->SetAUCTimeInterval(m_config.AUCTimeInterval);
    concentrationsToQuantitativeImageFilter->SetTiming(Timing);
    concentrationsToQuantitativeImageFilter->SetfTol(m_config.FTolerance);
    concentrationsToQuantitativeImageFilter->SetgTol(m_config.GTolerance);
    concentrationsToQuantitativeImageFilter->SetxTol(m_config.XTolerance);
    concentrationsToQuantitativeImageFilter->Setepsilon(m_config.Epsilon);
    concentrationsToQuantitativeImageFilter->SetmaxIter(m_config.MaxIter);
    concentrationsToQuantitativeImageFilter->Sethematocrit(m_config.Hematocrit);
    concentrationsToQuantitativeImageFilter->SetBatEstimator(batEstimator.get());
    if (m_config.ROIMaskFileName != "")
    {
      concentrationsToQuantitativeImageFilter->SetROIMask(roiMaskVolume);
    }

    if (m_config.ComputeFpv)
    {
      concentrationsToQuantitativeImageFilter->SetModelType(itk::LMCostFunction::TOFTS_3_PARAMETER);
    }
    else
    {
      concentrationsToQuantitativeImageFilter->SetModelType(itk::LMCostFunction::TOFTS_2_PARAMETER);
    }

    itk::PluginFilterWatcher watchQuantifier(concentrationsToQuantitativeImageFilter, "Quantifying", m_config.CLPProcessInformation, 19.0 / 20.0, 1.0 / 20.0);
    concentrationsToQuantitativeImageFilter->Update();

    ///////////////////////////////////// OUTPUT ////////////////////

    saveRequestedOutputs(signalToConcentrationsConverter, concentrationsToQuantitativeImageFilter, inputVectorVolume);

    return EXIT_SUCCESS;
  }

};


#endif
