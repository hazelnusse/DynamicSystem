add_library(foo foo.cc)

add_executable(fooexec main.cc)
target_link_libraries(fooexec foo)

if (PROJECT_BUILD_TESTS)
  add_subdirectory(tests)
endif()

install(TARGETS foo
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

