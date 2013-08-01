add_executable(test_disc test_disc.cc)
target_link_libraries(test_disc gtest_main)
add_test(NAME test_disc
         COMMAND test_disc)

