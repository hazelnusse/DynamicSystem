CMake + Googletest C++ project template
=======================================

This project provides the following structure

    .
    |-cmake
    |---Modules
    |-googletest
    |-include
    |-source
    |---tests

I found myself rewriting a lot of boilerplate code so I thought I would put
this all together as a starting point for C++ projects that use CMake and
Googletest.

After cloning this repo, you need to do:

    $ git submodule init
    $ git submodule update

Which pulls the latest googletest submodule. After that, you can build the
project, which consists of a minimal library and a minimal test. They dont'
really do anything other than exercise the build and test system and serve as
an example of how to add your own real code.

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make
    $ make test

