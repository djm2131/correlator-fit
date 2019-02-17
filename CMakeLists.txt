cmake_minimum_required( VERSION 3.0 )

project( CorrelatorFit )

find_package( LibXml2 REQUIRED )
find_package( GSL REQUIRED )

add_compile_definitions( USE_OMP )

include_directories( ${PROJECT_SOURCE_DIR}/include ${LIBXML2_INCLUDE_DIR} )

set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -fPIC -fopenmp" )

set( SOURCE_FILES include/xml_parser.h include/jack_stats.h include/fitter.h include/fitter_controls.h include/fit_results.h include/fit_functions.h include/eff_mass.h include/correlator.h src/correlator.cpp src/fitter.cpp include/gsl_pinv.h src/gsl_pinv.cpp )
add_library( cf STATIC ${SOURCE_FILES} )
target_link_libraries( cf GSL::gsl GSL::gslcblas ${LIBXML2_LIBRARIES} )
install(TARGETS cf DESTINATION ${PROJECT_SOURCE_DIR}/lib)

add_subdirectory(main)