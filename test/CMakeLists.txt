project(fccphat_tests)

add_subdirectory(googletest)

include_directories("${PROJECT_SOURCE_DIR}/googletest/googletest/include")
include_directories("${PROJECT_SOURCE_DIR}/../include")

add_executable(fccphat_tests
  src/main.cpp
        src/system_tests.cpp)
target_link_libraries(fccphat_tests
  fccphat
  gtest
)

add_custom_command(TARGET fccphat_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${PROJECT_SOURCE_DIR}/resources $<TARGET_FILE_DIR:fccphat_tests>/resources)
