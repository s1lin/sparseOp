# Set the project name
project (structure)

# Add a library with the above sources

add_library(${PROJECT_NAME} INTERFACE)
add_library(sub::structure ALIAS ${PROJECT_NAME})

target_include_directories( ${PROJECT_NAME} INTERFACE ${PROJECT_SOURCE_DIR}/include)
