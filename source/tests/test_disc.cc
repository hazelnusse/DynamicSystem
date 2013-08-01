#include "gtest/gtest.h"
#include "foo.h"

#include <array>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <utility>
#include "dynamic_system.hpp"

namespace dynamics {
class RollingDisc;

template <> struct Derived_traits<RollingDisc> {
  using state_type = std::array<double, 6>;
  using output_type = std::array<double, 2>;  // To hold KE and PE
};

class RollingDisc : public DynamicSystem<RollingDisc> {
 public:
  using state_type = Derived_traits<RollingDisc>::state_type;
  using output_type = Derived_traits<RollingDisc>::output_type;

  RollingDisc(double m = 1.0, double g = 1.0, double r = 1.0)
    : m_{m}, g_{g}, r_{r} {}

  /* State ordering
    q1 = x[0] // yaw 
    q2 = x[1] // lean
    q3 = x[2] // spin
    u1 = x[3] // b.x (in disc plane, forward) angular velocity measure number
    u2 = x[4] // b.y (perpendicular to disc plane) angular velocity measure number
    u3 = x[5] // b.z (in disc plane, down) angular velocity measure number
  */
  void operator()(const state_type & x, state_type & dxdt, const double)
  {
    // Kinematic differential equations
    dxdt[0] = x[5]/cos(x[1]);
    dxdt[1] = x[3];
    dxdt[2] = x[4] - x[5]*tan(x[1]);
    // Dynamic differential equations
    dxdt[3] = (4*g_*sin(x[1])/r_ + 6*x[4]*x[5] - x[5]*x[5]*tan(x[1]))/5;
    dxdt[4] = -2*x[3]*x[5]/3;
    dxdt[5] = -(2*x[4] - x[5]*tan(x[1]))*x[3];
  }

  output_type outputs(const state_type & x, const double t)
  {
    return {{1.0*t, 2.0*t*t}};
  }

 private:
  double m_, g_, r_;
};

} // namespace dynamics

template <class T>
void print_trajectory(const T & x, const std::vector<double> & t)
{
  std::cout << "Number of steps: " << x.size() << std::endl;
  const int rows = x[0].size();
  Eigen::Map<Eigen::MatrixXd> x_map(const_cast<double *>(x[0].data()), rows, 1);
  for (size_t i = 0; i < t.size(); ++i) {
    new (&x_map) Eigen::Map<Eigen::MatrixXd>(const_cast<double *>(x[i].data()), rows, 1);
    std::cout << t[i] << " " << x_map.transpose() << std::endl;
  }
}

TEST(DynamicSystem, TrueEqualsTrue)
{
  using namespace dynamics;
  RollingDisc kane;
  RollingDisc::state_type xi_kane = {{0.0, -0.01, 0.0, 0.0, 0.0, 0.4116}};
  std::vector<RollingDisc::state_type> x_kane;
  std::vector<double> t_kane;
  size_t steps;
  std::tie(x_kane, t_kane, steps) = kane.simulate(xi_kane, 0.0, 10.0, 0.01);
  print_trajectory(x_kane, t_kane);
  std::vector<RollingDisc::output_type> y_kane;
  std::tie(x_kane, t_kane, steps, y_kane) = kane.simulate_outputs(xi_kane, 0.0, 10.0, 0.01);
  print_trajectory(y_kane, t_kane);
}

