# CMAKE version of Ubuntu 16.04 is 3.5.
cmake_minimum_required(VERSION 3.5)
project(test)

set(LINKLIBS)
set(COMPILE_OPTION)
set(SRC)
set(INCL)
  
if(NOT MSVC)
	list(APPEND LINKLIBS
	# <math.h>
	m
	# openMP
    -fopenmp
	)
else()
	list(APPEND LINKLIBS
    
	)	
endif()


list(APPEND SRC
	RtBase.h
	RtInput.h
	RtOutput.h
	WAV.h
	MLDR.h
)
	
list(APPEND INCL
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/STFT/cpp

list(APPEND COMPILE_OPTION
  -g
)


if(NOT MSVC)
	list(APPEND COMPILE_OPTION
	# openMP
	  -fopenmp
	# debugging
	  -g
	# -O3 : maximize performace
	# -O2 : maximize performance as long as not increasing size of executable
	# -O1 : minimize compile time and others are secondary
	# -Os : minimize size of executable
	  -O3
	)
	
else()
list(APPEND COMPILE_OPTION
# Visual Studio uses '/' for compile option
	  /openmp
	  # /MP - compiled in mult-process : faster build time
	  /MP
	)
endif()


include(RtAudio.cmake)
include(Libsamplerate.cmake)

add_executable(${PROJECT_NAME}_output test/test_output.cpp ${SRC})
target_include_directories(${PROJECT_NAME}_output	PRIVATE	${INCL}  )
target_link_libraries(${PROJECT_NAME}_output ${LINKLIBS})
target_compile_options(${PROJECT_NAME}_output PRIVATE ${COMPILE_OPTION} )
set_target_properties(${PROJECT_NAME}_output PROPERTIES  CXX_STANDARD 14 )