# Set the project name
project (solver)

# Add a library with the above sources
add_library(${PROJECT_NAME} INTERFACE)
add_library(sub::solver ALIAS ${PROJECT_NAME})

target_include_directories( ${PROJECT_NAME} INTERFACE ${PROJECT_SOURCE_DIR}/include)


