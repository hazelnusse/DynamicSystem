add_executable(test_foo test_foo.cc)
target_link_libraries(test_foo gtest_main foo)
add_test(NAME test_foo
         COMMAND test_foo)

