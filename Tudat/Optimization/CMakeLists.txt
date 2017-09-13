 #    Copyright (c) 2010-2017, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Add source files.
set(OPTIMIZATION_SOURCES
  "${SRCROOT}${OPTIMIZATIONDIR}/decisionVariableSettings.cpp"
  "${SRCROOT}${OPTIMIZATIONDIR}/missionLinking.cpp"
  "${SRCROOT}${OPTIMIZATIONDIR}/objectiveFunction.cpp"
)

# Add header files.
set(OPTIMIZATION_HEADERS 
  "${SRCROOT}${OPTIMIZATIONDIR}/decisionVariableSettings.h"
  "${SRCROOT}${OPTIMIZATIONDIR}/missionLinking.h"
  "${SRCROOT}${OPTIMIZATIONDIR}/objectiveFunction.h"
)

# Add unit test files.
set(OPTIMIZATION_UNITTESTS
  "${SRCROOT}${OPTIMIZATIONDIR}/UnitTests/unitTestDecisionVariableSettings.cpp"
)

# Add static libraries.
add_library(tudat_optimization STATIC ${OPTIMIZATION_SOURCES} ${OPTIMIZATION_HEADERS})
setup_tudat_library_target(tudat_optimization "${SRCROOT}${OPTIMIZATIONDIR}")

# Add unit tests.
add_executable(test_decision_variable_settings "${SRCROOT}${OPTIMIZATIONDIR}/UnitTests/unitTestDecisionVariableSettings.cpp")
setup_custom_test_program(test_decision_variable_settings "${SRCROOT}${OPTIMIZATIONDIR}")
target_link_libraries(test_decision_variable_settings ${TUDAT_PROPAGATION_LIBRARIES}
    tudat_basic_astrodynamics pthread ${Boost_LIBRARIES} tudat_optimization )

