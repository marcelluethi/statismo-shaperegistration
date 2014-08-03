CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(statismo-shaperegistration)

FIND_PACKAGE(STATISMO REQUIRED)

FIND_PACKAGE(ITK REQUIRED)
include (${ITK_USE_FILE})
	

# ITK comes with its own version of HDF5. As we want to use the version that comes with ITK,
# we have to find the directories where the include files are located and set the inlcude dir #accordingly
find_path(ITK_HDF5_INCLUDE_DIR hdf5.h PATHS ${ITK_INCLUDE_DIRS} PATH_SUFFIXES itkhdf5 itkhdf5/cpp )
find_path(ITK_HDF5_INCLUDE_DIR_CPP H5Cpp.h PATHS ${ITK_INCLUDE_DIRS} PATH_SUFFIXES itkhdf5 itkhdf5/cpp )

set(ITK_HDF5_INCLUDE_DIRS ${ITK_HDF5_INCLUDE_DIR} ${ITK_HDF5_INCLUDE_DIR_CPP})

include_directories(${STATISMO_INCLUDE_DIRS} ${ITK_HDF5_INCLUDE_DIRS}) 

link_directories( ${ITK_LIBRARY_DIRS})


add_executable (shaperegistration src/shaperegistration.cpp)
target_link_libraries (shaperegistration ${ITK_LIBRARIES})

add_executable (buildgpmodel src/buildgpmodel.cpp)
target_link_libraries (buildgpmodel ${ITK_LIBRARIES})
