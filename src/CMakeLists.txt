
file(GLOB RapidSim_sources ${PROJECT_SOURCE_DIR}/src/*.cc)

ADD_EXECUTABLE ( RapidSim.exe ${RapidSim_sources} ${PROJECT_SOURCE_DIR}/src/RapidSim.C )

if(EvtGen_FOUND)
  TARGET_LINK_LIBRARIES( RapidSim.exe ${ROOT_LIBRARIES} ${EVTGEN} ${EVTGENEXT} )
ELSE()
  TARGET_LINK_LIBRARIES( RapidSim.exe ${ROOT_LIBRARIES} )
ENDIF()

# install target

install(TARGETS RapidSim.exe DESTINATION ${CMAKE_INSTALL_BINDIR} RUNTIME DESTINATION bin)

