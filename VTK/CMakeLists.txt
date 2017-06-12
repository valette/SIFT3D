
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
add_executable(SIFT3D SIFT3D.cxx)
target_link_libraries(SIFT3D PUBLIC sift3D imutil ${VTK_LIBRARIES})
