add_executable( correlator-fit.x main.cpp )
target_link_libraries( correlator-fit.x cf )
install(TARGETS correlator-fit.x DESTINATION ${PROJECT_SOURCE_DIR}/bin)