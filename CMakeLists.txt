# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0 FATAL_ERROR)

set(PROJECT_NAME_STR RooUnfold)
set(PROJECT_LIB_NAME ${PROJECT_NAME_STR})
PROJECT(${PROJECT_NAME_STR})
include_directories(include)

set(CINTFILE ${PROJECT_NAME_STR}Dict.cxx)
file(GLOB LINK_DEF_FILE include/*_LinkDef.h)
file(GLOB HEADER_FILES include/*.h)
list(REMOVE_ITEM HEADER_FILES ${LINK_DEF_FILE})


if((CMAKE_SYSTEM_PROCESSOR MATCHES "i386") AND (CMAKE_SIZEOF_VOID_P EQUAL 8) AND (APPLE))
	set(CMAKE_OSX_ARCHITECTURES "x86_64")
	MESSAGE(STATUS "Building ${PROJECT_NAME_STR} for ${CMAKE_OSX_ARCHITECTURES} architecture on ${CMAKE_SYSTEM_NAME}")
else()
	MESSAGE(STATUS "Building ${PROJECT_NAME_STR} on ${CMAKE_SYSTEM_NAME}")
endif()
MESSAGE(STATUS "Using compiler ${CMAKE_CXX_COMPILER_ID}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
MESSAGE(STATUS "Using CMAKE module path ${CMAKE_SOURCE_DIR}")

# get ROOT
find_package(ROOT "5.34" REQUIRED)
MESSAGE(STATUS "** ROOT Include path: ${ROOT_INCLUDE_DIRS}")
MESSAGE(STATUS "** ROOT Library path: ${ROOT_LIBRARY_DIRS}")
MESSAGE(STATUS "** ROOT Libraries: ${ROOT_LIBRARIES}")

include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIRS})

# load helper macros for using rootcint 
INCLUDE( ${ROOT_DICT_MACROS_FILE} )
set(ROOT_DICT_INPUT_HEADERS ${HEADER_FILES})
set(ROOT_DICT_OUTPUT_DIR ${PROJECT_SOURCE_DIR})
# and add LinkDef file to the end
LIST( APPEND ROOT_DICT_INPUT_HEADERS ${LINK_DEF_FILE} )

# generate dictionary sources using rootcint
GEN_ROOT_DICT_SOURCE( ${CINTFILE} )

# get boost
set(Boost_USE_MULTITHREADED ON)
find_package(Boost COMPONENTS
				unit_test_framework
             REQUIRED)
MESSAGE(STATUS "** Boost Include path: ${Boost_INCLUDE_DIR}")
MESSAGE(STATUS "** Boost Library path: ${Boost_LIBRARY_DIRS}")
MESSAGE(STATUS "** Boost Libraries: ${Boost_LIBRARIES}")
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

#Check the compiler and set the compile and link flags
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
# properly set the C++ standard according to the standard used for the root build
# Enable C++11 by default if found in ROOT
message(STATUS "Determining C++ standard ...")
if(ROOT_HAS_CXX11 AND NOT DISABLE_CXX11)
  message(STATUS "Enabling C++11")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

# Enable C++14 by default if found in ROOT
if(ROOT_HAS_CXX14 AND NOT DISABLE_CXX14)
  message(STATUS "Enabling C++14")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif()
# first attempt to make cmake work again on OS X
if((CMAKE_CXX_COMPILER_ID STREQUAL "Clang") AND (APPLE))
	MESSAGE(STATUS "** std library for clang: libc++")
	set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -stdlib=libc++")
endif()

# adding the library
aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC)
LIST(APPEND SRC ${ROOT_DICT_OUTPUT_SOURCES} )
add_library(${PROJECT_LIB_NAME} SHARED ${SRC})
target_link_libraries(${PROJECT_LIB_NAME}
	 ${Boost_LIBRARIES}
	 ${ROOT_LIBRARIES}
	 )

# tests
enable_testing()
set(PROJECT_TEST_NAME ${PROJECT_NAME_STR}_test)
include_directories(${COMMON_INCLUDES})
file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp)
add_executable(${PROJECT_TEST_NAME} ${TEST_SRC_FILES})
target_link_libraries(${PROJECT_TEST_NAME}
	 ${Boost_LIBRARIES}
	 ${ROOT_LIBRARIES}
	 ${PROJECT_LIB_NAME}
	 )

add_test(test1 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${PROJECT_TEST_NAME})

# install
install(TARGETS ${PROJECT_LIB_NAME} ${PROJECT_TEST_NAME}
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)
