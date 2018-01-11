/*=========================================================================

  Program:   PkModeling module
  Module:    $HeadURL: http://svn.slicer.org/Slicer4/trunk/Modules/CLI/PkModeling/PkModeling.cxx $
  Language:  C++
  Date:      $Date: 2012-06-06 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.
  =========================================================================*/

#include "PkModeling.hxx"

#define TESTMODE_ERROR_TOLERANCE 0.1


std::unique_ptr<PkModelingBase> getNewPkModelingInstanceOrNull(std::string inputImageFileName)
{
  std::unique_ptr<PkModelingBase> pkModeling;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;
  itk::GetImageType(inputImageFileName, pixelType, componentType);

  switch (componentType)
  {
  case itk::ImageIOBase::CHAR:
  case itk::ImageIOBase::UCHAR:
  case itk::ImageIOBase::SHORT:
    pkModeling.reset(new PkModeling<short, short>());
    break;
  case itk::ImageIOBase::USHORT:
  case itk::ImageIOBase::INT:
    pkModeling.reset(new PkModeling<int, short>());
    break;
  case itk::ImageIOBase::UINT:
  case itk::ImageIOBase::ULONG:
    pkModeling.reset(new PkModeling<unsigned long, short>());
    break;
  case itk::ImageIOBase::LONG:
    pkModeling.reset(new PkModeling<long, short>());
    break;
  case itk::ImageIOBase::FLOAT:
    pkModeling.reset(new PkModeling<float, short>());
    break;
  case itk::ImageIOBase::DOUBLE:
    pkModeling.reset(new PkModeling<float, short>());
    break;
  case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
  default:
    std::cout << "unknown component type" << std::endl;
    break;
  }
  return pkModeling;
}


int main(int argc, char * argv[])
{
  PARSE_ARGS;


  // this line is here to be able to see the full output on the dashboard even
  // when the test succeeds (to see the reproducibility error measure)
  std::cout << "ctest needs: CTEST_FULL_OUTPUT" << std::endl;


  std::unique_ptr<PkModelingBase> pkModeling = getNewPkModelingInstanceOrNull(InputFourDImageFileName);
  if (!pkModeling) {
    return EXIT_FAILURE;
  }

  try
  {
    pkModeling->DoIt(argc, argv);
  }
  catch (itk::ExceptionObject & excep)
  {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
