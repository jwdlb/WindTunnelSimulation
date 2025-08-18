cmake_minimum_required(VERSION 3.18)
project(WindTunnelSim LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 20)
enable_language(CUDA)

# SFML paths
set(SFML_INCLUDE_DIR "C:/dev/SFML-2.5.0/include")
set(SFML_LIBRARY_DIR "C:/dev/SFML-2.5.0/lib")

# GLEW paths
set(GLEW_INCLUDE_DIR "C:/dev/glew-2.1.0/include")
set(GLEW_LIBRARY_DIR "C:/dev/glew-2.1.0/lib/Release/x64")

# ImGui-SFML paths
set(IMGUI_SFML_INCLUDE_DIR "C:/dev/vcpkg/installed/x64-windows/include")
set(IMGUI_SFML_LIBRARY_DIR "C:/dev/vcpkg/installed/x64-windows/lib")

# Add executable
add_executable(WindTunnelSim main.cpp simulation.cpp simulation.h simulationGPU.h simulationGPU.cu)

# Include headers
target_include_directories(WindTunnelSim PRIVATE ${SFML_INCLUDE_DIR} ${GLEW_INCLUDE_DIR} ${IMGUI_SFML_INCLUDE_DIR})

find_package(CUDAToolkit REQUIRED)

# Link libraries
target_link_libraries(WindTunnelSim PRIVATE
        ${SFML_LIBRARY_DIR}/sfml-graphics.lib
        ${SFML_LIBRARY_DIR}/sfml-window.lib
        ${SFML_LIBRARY_DIR}/sfml-system.lib
        ${GLEW_LIBRARY_DIR}/glew32.lib
        ${IMGUI_SFML_LIBRARY_DIR}/imgui-sfml.lib
        ${IMGUI_SFML_LIBRARY_DIR}/imgui.lib
        opengl32.lib
        CUDA::cudart
        CUDA::cublas
        CUDA::cufft
)


# CUDA architecture for RTX 4060
# set_target_properties(WindTunnelSim PROPERTIES CUDA_ARCHITECTURES "89")