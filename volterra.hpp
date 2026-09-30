// declarations

#ifndef FN_VOLTERRA_HPP
#define FN_VOLTERRA_HPP

#include <vector>

namespace fn {

struct UserPoint
{
  double x;
  double y;
};

struct InternalPoint
{
  double xrel;
  double yrel;
  double H;
};

struct OutputPoint
{
  double x;
  double y;
  double H;
};

struct Parameters
{
  double A;
  double B;
  double C;
  double D;
  size_t N;
  double delta_t;
};

class Simulation
{
  std::vector<InternalPoint> points_;
  Parameters pars_;
  UserPoint p0_;

 public:
  size_t size() const;

  void addParameters(Parameters const&);

  void addUserPoint(UserPoint const&);

  void evolve();

  const std::vector<InternalPoint>&
  points() const; // note: not strictly necessary, but one might want to know
                  // the internal representation vector for later analysis

  std::vector<OutputPoint> result() const;
};

std::vector<OutputPoint> simulation(Simulation const&);

} // namespace fn

#endif