################################################################################
 #    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    #
 #                                                                              #
 #              This software is distributed under the terms of the             # 
 #         GNU Lesser General Public Licence version 3 (LGPL) version 3,        #  
 #                  copied verbatim in the file "LICENSE"                       #
 ################################################################################
# - Try to find EvtGen instalation
# Once done this will define
#

MESSAGE(STATUS "Looking for EvtGen ...")

FIND_PATH(RAMC_INCLUDE_DIR NAMES ramC.h PATHS
   /Users/luchinsky/Work/ramC/src/
	/home/anton/programs/ramc/src/
  NO_DEFAULT_PATH
)


FIND_PATH(RAMC_LIBRARY_DIR NAMES libramC.a PATHS
  ${SIMPATH}/lib
 /Users/luchinsky/Work/ramC/build/
/home/anton/programs/ramc/build/
  NO_DEFAULT_PATH
)


if (RAMC_INCLUDE_DIR AND RAMC_LIBRARY_DIR)
   set(RAMC_FOUND TRUE)
endif (RAMC_INCLUDE_DIR AND RAMC_LIBRARY_DIR)

if (RAMC_FOUND)
  if (NOT RAMC_FOUND_QUIETLY)
    MESSAGE(STATUS "Looking for RAMC... - found ${RAMC_LIBRARY_DIR}")
    SET(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${RAMC_LIBRARY_DIR})
  endif (NOT RAMC_FOUND_QUIETLY)
else (RAMC_FOUND)
  if (RAMC_FOUND_REQUIRED)
    message(FATAL_ERROR "Looking for RAMC... - Not found")
  endif (RAMC_FOUND_REQUIRED)
endif (RAMC_FOUND)

