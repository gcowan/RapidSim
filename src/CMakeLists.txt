
file(GLOB RapidSim_sources ${PROJECT_SOURCE_DIR}/src/*.cc) 

ADD_EXECUTABLE ( RapidSim.exe ${RapidSim_sources} ${PROJECT_SOURCE_DIR}/src/RapidSim.C )

TARGET_LINK_LIBRARIES( RapidSim.exe Core RIO RooFit RooFitCore RooStats Hist Tree Matrix Physics MathCore )

# install target
 
install(TARGETS RapidSim.exe DESTINATION ${CMAKE_INSTALL_BINDIR} RUNTIME DESTINATION bin) 

