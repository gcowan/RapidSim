# include ROOT library and include files

SET( CMAKE_CXX_FLAGS "-Weffc++ -Os -msse -msse2 -msse3 -m3dnow -ansi -fPIC -fmerge-all-constants -funroll-all-loops -g -D__ROOFIT_NOBANNER -Wconversion -Weffc++ -Wextra -Wsign-compare -Wfloat-equal -Wmissing-noreturn -Wall -Wno-non-virtual-dtor -std=c++11" )

SET (LIBRARIES functions RooGounarisSakurai RooRelBreitWigner FastDecay)
FOREACH ( lib ${LIBRARIES} ) 
    ADD_LIBRARY ( ${lib} ${PROJECT_SOURCE_DIR}/src/${lib}.cc )
TARGET_LINK_LIBRARIES( ${lib} dl Core RooFit RooFitCore RooStats RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread Gui pthread m)
ENDFOREACH ( lib )

ADD_EXECUTABLE ( FastSim.exe ${PROJECT_SOURCE_DIR}/src/FastSim.C )
TARGET_LINK_LIBRARIES( FastSim.exe
    FastDecay functions RooGounarisSakurai RooRelBreitWigner dl Core RooFit RooFitCore RooStats RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread Gui pthread m)

