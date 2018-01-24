#-----------------------------------------------------------------------------
# Testing DROs expecting to fail
#-----------------------------------------------------------------------------
set(inputDataBaseName ${CMAKE_CURRENT_SOURCE_DIR}/../../../Data/RegressionTests/DROs/Input/DRO)
set(referenceDataBaseDir ${CMAKE_CURRENT_SOURCE_DIR}/../../../Data/RegressionTests/DROs/Reference/)

set(testName DRO3min5secinf_FailOnMissingPrescribedAIFArg)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ModuleEntryPointExpectFail
    --aifMode Prescribed
    ${inputDataBaseName}3min5secinf.nrrd                   
)
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

set(testName DRO3min5secinf_FailOnInvalidPrescribedAIFFile)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ModuleEntryPointExpectFail
    --aifMode Prescribed
    --prescribedAIF INVALID-FILENAME
    ${inputDataBaseName}3min5secinf.nrrd                   
)
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

set(testName DRO3min5secinf_FailOnMissingPrescribedAIFArgButHavingAIFMaskFile)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ModuleEntryPointExpectFail
    --aifMode Prescribed
    --aifMask ${inputDataBaseName}-AIF.nrrd
    ${inputDataBaseName}3min5secinf.nrrd                   
)
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

set(testName DRO3min5secinf_FailOnMissingAIFMaskArg)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ModuleEntryPointExpectFail
    --aifMode AverageUnderAIFMask
    ${inputDataBaseName}3min5secinf.nrrd                   
)
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

set(testName DRO3min5secinf_FailOnInvalidAIFMaskFile)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ModuleEntryPointExpectFail
    --aifMode AverageUnderAIFMask
    --aifMask INVALID-FILENAME
    ${inputDataBaseName}3min5secinf.nrrd                   
)
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

set(testName DRO3min5secinf_FailOnMissingAIFMaskArgButHavingPrescribedAIFFile)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ModuleEntryPointExpectFail
    --aifMode AverageUnderAIFMask
    --prescribedAIF ${inputDataBaseName}3min5secinf_PrescribedAIF.csv
    ${inputDataBaseName}3min5secinf.nrrd                   
)
set_property(TEST ${testName} PROPERTY LABELS ${CLP})



