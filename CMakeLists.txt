#  CMAKE CONFIG
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 FATAL_ERROR )
project(RapidSim)

set(CMAKE_VERBOSE_MAKEFILE ON)

#This module respects HFS, e.g. defines lib or lib64 when it is needed. 
include("GNUInstallDirs") 


SET( ${PROJECT_NAME}_MAJOR_VERSION 1 )
SET( ${PROJECT_NAME}_MINOR_VERSION 1 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 1 )

SET( CMAKE_BUILD_TYPE Release )

EXECUTE_PROCESS(
            COMMAND git rev-parse HEAD
        OUTPUT_VARIABLE GIT_REV )
MESSAGE( "GIT VERSION:  ${GIT_REV}" )

EXECUTE_PROCESS(
        COMMAND date +%Y_%m_%d
        OUTPUT_VARIABLE BUILD_DATE )

SET ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} )

set(CMAKE_CXX_FLAGS    "${CMAKE_CXX_FLAGS} -msse -msse2 -msse3 -m3dnow -fPIC -fmerge-all-constants -D__ROOFIT_NOBANNER -Wall -Wextra -Werror -Wsign-compare -Wmissing-noreturn -Wno-non-virtual-dtor")

set(CMAKE_CXX_FLAGS_DEBUG          "-g3 -Wall") 
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG") 
set(CMAKE_CXX_FLAGS_RELEASE        "-O2 -DNDEBUG") 
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g3 -O2")


# setting linker flags

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 
  SET(CMAKE_EXE_LINKER_FLAGS "-Wl,-dead_strip_dylibs")   
  SET(CMAKE_SHARED_LINKER_FLAGS "-Wl,-dead_strip_dylibs")   
  SET(STATIC_LIBRARY_FLAGS "-Wl,-dead_strip_dylibs")
ELSEIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  SET(CMAKE_EXE_LINKER_FLAGS "-Wl,--as-needed")
  SET(CMAKE_SHARED_LINKER_FLAGS "-Wl,--as-needed")
  SET(STATIC_LIBRARY_FLAGS "-Wl,--as-needed")	
ENDIF()


set (ROOT_FIND_COMPONENTS Core RIO Tree)
  find_file(ROOT_CONFIG ROOTConfig.cmake HINTS $ENV{ROOTSYS} ${ROOTSYS} ${ROOT_DIR} NO_DEFAULT_PATH PATH_SUFFIXES cmake)
if (ROOT_CONFIG)
  message(STATUS "ROOTConfig.cmake found and will be used: ${ROOT_CONFIG}")
  include(${ROOT_CONFIG})
  message(STATUS "ROOT use file: ${ROOT_USE_FILE}")
  include(${ROOT_USE_FILE})

  if(ROOT_LIBRARIES AND ROOT_INCLUDE_DIRS)
    set(ROOT_FOUND TRUE)
  endif()

else()      
  message(STATUS "ROOTConfig.cmake not found, trying to use FindROOT.cmake...")      
  find_package(ROOT)   
  INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIR} )
  LINK_DIRECTORIES( ${ROOT_LIBRARY_DIR} )
  set(CMAKE_CXX_FLAGS    "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}" ) 
endif()

message(STATUS "ROOT includes: ${ROOT_INCLUDE_DIR}")   
message(STATUS "ROOT libraries: ${ROOT_LIBRARIES}")   
message(STATUS "ROOT library directory: ${ROOT_LIBRARY_DIR}")
message(STATUS "ROOT_FOUND: ${ROOT_FOUND}")

if(ROOT_FOUND)     
  include_directories(${ROOT_INCLUDE_DIRS})     
else()     
  message(ERROR "ROOT needed for RapidSim not found")     
endif()

add_subdirectory(src)

set(bindir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

if (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  configure_file(${CMAKE_SOURCE_DIR}/validation/runValidation.sh.in ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/runValidation.sh   @ONLY)
endif()

install(DIRECTORY validation DESTINATION ${CMAKE_INSTALL_DATADIR} )
install(DIRECTORY rootfiles DESTINATION ${CMAKE_INSTALL_DATADIR} )
install(DIRECTORY utils DESTINATION ${CMAKE_INSTALL_DATADIR} )
install(DIRECTORY config DESTINATION ${CMAKE_INSTALL_DATADIR} )



