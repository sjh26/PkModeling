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


int main(int argc, char * argv[])
{
  PARSE_ARGS_INTO_CONFIG;

  // this line is here to be able to see the full output on the dashboard even
  // when the test succeeds (to see the reproducibility error measure)
  std::cout << "ctest needs: CTEST_FULL_OUTPUT" << std::endl;

  PkModeling pkModeling(configuration);

  try
  {
    pkModeling.execute();
  }
  catch (std::exception& excep)
  {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
