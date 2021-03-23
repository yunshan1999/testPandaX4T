#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.6 FATAL_ERROR)
project(main)

#----------------------------------------------------------------------------
# Find Root package, activating all available UI and Vis drivers by default
#

find_package(ROOT 6.0 REQUIRED COMPONENTS)


#----------------------------------------------------------------------------
# Setup Root include directories and compile definitions
# Setup include directory for this project
#
include(${ROOT_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.c)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Root libraries
#
add_executable(main main.C ${sources} ${headers})
target_link_libraries(main ${ROOT_LIBRARIES})




#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS main DESTINATION bin)