# include ROOT library and include files

SET( CMAKE_CXX_FLAGS "-Os -msse -msse2 -msse3 -m3dnow -ansi -fPIC -fmerge-all-constants -funroll-all-loops -g -D__ROOFIT_NOBANNER -Wall -Wextra -Werror -Wsign-compare -Wmissing-noreturn -Wno-non-virtual-dtor -std=c++11" )

SET (FILES functions RooGounarisSakurai RooRelBreitWigner RapidParticle RapidParticleData RapidDecay RapidMomentumSmearHisto RapidMomentumSmearGauss)
FOREACH ( lib ${FILES} ) 
	SET (LIBPATH ${LIBPATH} ${PROJECT_SOURCE_DIR}/src/${lib}.cc)
ENDFOREACH ( lib )

ADD_LIBRARY ( RapidSim ${LIBPATH} )
TARGET_LINK_LIBRARIES( RapidSim dl Core RooFit RooFitCore RooStats RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread Gui pthread m)

ADD_EXECUTABLE ( RapidSim.exe ${PROJECT_SOURCE_DIR}/src/RapidSim.C )
TARGET_LINK_LIBRARIES( RapidSim.exe
    RapidSim dl Core RooFit RooFitCore RooStats RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread Gui pthread m) 

