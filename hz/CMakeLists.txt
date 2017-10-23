cmake_minimum_required(VERSION 2.6)
project(haptic_test)


# By default, build programs with compiler optimizations and debugging symbols
if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()


## libbarrett
find_package(Barrett REQUIRED)
include_directories(${BARRETT_INCLUDE_DIRS})
link_directories(${BARRETT_LIBRARY_DIRS})
add_definitions(${BARRETT_DEFINITIONS})

function(add_programs)
	foreach (prog ${ARGV})
		add_executable(${prog} ${prog}.cpp)
		target_link_libraries(${prog} ${BARRETT_LIBRARIES})
	endforeach()
endfunction()

add_programs(
	haptic_test
)
