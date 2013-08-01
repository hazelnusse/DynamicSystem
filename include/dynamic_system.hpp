/*
 * ============================================================================
 *
 *         Author:  Dale Lukas Peterson (dlp), hazelnusse@gmail.com
 *
 *    Description:  DynamicSystem class for simulating ODE's using Boost.odeint
 *
 * ============================================================================
 */
#include <utility>
#include <vector>
namespace dynamics {

// Declare templated traits class for Derived
template<typename Derived> struct DynamicSystem_traits {};

template <class Derived>
class DynamicSystem {
public:
  using st = typename DynamicSystem_traits<Derived>::state_type;
  using ot = typename DynamicSystem_traits<Derived>::output_type;

  std::tuple<std::vector<st>, std::vector<double>, size_t>
  simulate(const st & xi, double ti, double tf, double ts)
  {
    using namespace boost::numeric::odeint;
    std::vector<st> x;
    std::vector<double> t;
    st xi_ = xi;
    size_t steps = integrate_const(
            make_controlled(1.0e-8, 1.0e-8, runge_kutta_fehlberg78<st>()),
            *static_cast<Derived *>(this), xi_, ti, tf, ts, state_observer(x, t));

    return std::make_tuple(x, t, steps);
  }
  
  std::tuple<std::vector<st>, std::vector<double>, size_t, std::vector<ot>>
  simulate_outputs(const st & xi, double ti, double tf, double ts)
  {
    using namespace boost::numeric::odeint;
    std::vector<st> x;
    std::vector<double> t;
    std::vector<ot> y;
    st xi_ = xi;

    size_t steps = integrate_const(
            make_controlled(1.0e-8, 1.0e-8, runge_kutta_fehlberg78<st>()),
            *static_cast<Derived *>(this), xi_, ti, tf, ts,
            output_observer(x, t, y, *static_cast<Derived *>(this)));

    return std::make_tuple(x, t, steps, y);
  }

private:
  class state_observer {
  public:
    state_observer(std::vector<st> & x,
                   std::vector<double> & t)
      : x_(x), t_(t) {}

    void operator()(const st& x, double t)
    {
      x_.push_back(x);
      t_.push_back(t);
    }

  private:
    std::vector<st> & x_;
    std::vector<double> & t_;
  };

  class output_observer {
  public:
    output_observer(std::vector<st> & x,
                    std::vector<double> & t,
                    std::vector<ot> & y, Derived & ds)
      : x_(x), t_(t), y_(y), ds_(ds) {}

    void operator()(const st& x, double t)
    {
      x_.push_back(x);
      t_.push_back(t);
      y_.push_back(ds_.outputs(x, t));
    }

  private:
    std::vector<st> & x_;
    std::vector<double> & t_;
    std::vector<ot> & y_;
    Derived & ds_;
  };
};

} // namespace dynamics

