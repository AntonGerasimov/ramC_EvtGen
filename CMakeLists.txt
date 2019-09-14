##Test
MESSAGE("Hello, world")
IF(DEFINED APPLE)
    MESSAGE("APPLE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ENDIF(DEFINED APPLE)

set(CMAKE_OSX_DEPLOYMENT_TARGET "")
INCLUDE(SetCompileFlag.cmake)

cmake_minimum_required(VERSION 2.6)
project(directPsi CXX Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

## include CERN ROOT
include(FindROOT.cmake)
MESSAGE(${ROOT_LIBRARIES})
include(FindEvtGen.cmake)
include(FindHepMC.cmake)
include(FindRAMC.cmake)

## set fortran flags
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-fno-automatic")



include_directories(~/local/include)


include_directories(src)
include_directories(${ROOT_INCLUDE_DIR})
include_directories(${EVTGEN_INCLUDE_DIR}) 
include_directories(${RAMC_INCLUDE_DIR}) 

include_directories(./)
include_directories(src/)




#adding sources for prog.exe
add_executable (prog.exe "src/prog.c++" "src/EvtWnPi2.cpp")
target_link_libraries(prog.exe ${EVTGEN_LIBRARY_DIR}/archive/libEvtGen.a)
target_link_libraries(prog.exe ${HEPMC_LIBRARY_DIR}/libHepMC.a)
target_link_libraries(prog.exe ${ROOT_LIBRARIES}) 
target_link_libraries(prog.exe ${RAMC_LIBRARY_DIR}/libramC.a) 


#adding sources for pion.exe
add_executable (pion.exe "src/pion.c++" "src/EvtWnPi2.cpp")
target_link_libraries(pion.exe ${EVTGEN_LIBRARY_DIR}/archive/libEvtGen.a)
target_link_libraries(pion.exe ${HEPMC_LIBRARY_DIR}/libHepMC.a)
target_link_libraries(pion.exe ${ROOT_LIBRARIES}) 
target_link_libraries(pion.exe ${RAMC_LIBRARY_DIR}/libramC.a) 
