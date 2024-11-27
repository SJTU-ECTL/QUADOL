set(util_include ${PROJECT_SOURCE_DIR}/util)
set(util_include ${PROJECT_SOURCE_DIR}/util PARENT_SCOPE)

include_directories(${ectl_include})
include_directories(${ext_include})

aux_source_directory(. DIR_SRC)
add_library(util_libraries ${DIR_SRC})
