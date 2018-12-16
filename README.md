Dynamic System simulator for C++
================================

As of Boost 1.53.0, Boost.odeint is available to integrate ODE's in C++. This
package provides a simple template class `dynamics::DynamicSystem` that is
designed to be subclassed using the curiously recurring template pattern
(CTRP). The class allows you to easily simulate the ODE's for a given set of
initial conditions, and also makes it very easy to compute output quantities
(i.e., functions of the states and possibly time) during simulation (as opposed
to after the simulation).

To see an example of how to subclass from `DynamicSystem`, see the test in
`./source/tests/test_disc.cc` which is for the rolling disc equations of
motion.

After cloning this repo, you need to do:

    $ git submodule init
    $ git submodule update

Which pulls the latest googletest submodule. After that, you can build the
project with

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make
    $ make test

