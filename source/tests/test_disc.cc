#include <array>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <utility>
#include "gtest/gtest.h"
#include "dynamic_system.hpp"

namespace dynamics {

// Foward declaration, REQUIRED for template code to work
class RollingDisc;

// Template specialization of traits class to define the state type and the
// output type
template <> struct DynamicSystem_traits<RollingDisc> {
  using state_type = std::array<double, 6>;   // Container for states
  using output_type = std::array<double, 3>;  // Container for KE, PE, KE + PE
};

// Rolling disc class that uses curiously recurring template patter to inherit
// from DynamicSystem<RollingDisc>.  This is a way to get static polymorphism.
class RollingDisc : public DynamicSystem<RollingDisc> {
 public:
  // Convenience typedefs for the state type and the output type
  using state_type = DynamicSystem_traits<RollingDisc>::state_type;
  using output_type = DynamicSystem_traits<RollingDisc>::output_type;

  RollingDisc(double m = 1.0, double g = 1.0, double r = 1.0)
    : m_{m}, g_{g}, r_{r} {}

  /* Function which computes state derivatives. State ordering is:
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

  // Function which computes outputs
  output_type outputs(const state_type & x, const double t)
  {
    double ke, pe;
    ke = (1.0/8.0)*m_*pow(r_, 2)*(5*pow(x[3], 2) + 6*pow(x[4], 2) + pow(x[5], 2));
    pe = g_*m_*r_*(cos(x[1]) - 1);
    return {{ke, pe, ke+pe}};
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
  RollingDisc::state_type xi_kane = {{0.0, 0.0, 0.0, 0.01, 0.0, 0.4116}};
  std::vector<RollingDisc::state_type> x_kane;
  std::vector<double> t_kane;
  size_t steps;
  std::tie(x_kane, t_kane, steps) = kane.simulate(xi_kane, 0.0, 1.0, 0.1);
  print_trajectory(x_kane, t_kane);
  std::vector<RollingDisc::output_type> y_kane;
  std::tie(x_kane, t_kane, steps, y_kane) = kane.simulate_outputs(xi_kane, 0.0, 1.0, 0.1);
  print_trajectory(y_kane, t_kane);
}

