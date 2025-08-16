cmake_minimum_required(VERSION 3.15)
project(WindTunnelSim)

set(CMAKE_CXX_STANDARD 20)

# SFML paths
set(SFML_INCLUDE_DIR "C:/dev/SFML-2.5.0/include")
set(SFML_LIBRARY_DIR "C:/dev/SFML-2.5.0/lib")

# Add executable
add_executable(WindTunnelSim simulation.h simulation.cpp main.cpp)

# Include headers
target_include_directories(WindTunnelSim PRIVATE ${SFML_INCLUDE_DIR})

# Link libraries manually
target_link_libraries(WindTunnelSim PRIVATE
        ${SFML_LIBRARY_DIR}/sfml-graphics.lib
        ${SFML_LIBRARY_DIR}/sfml-window.lib
        ${SFML_LIBRARY_DIR}/sfml-system.lib
)
