// main executable. Functions defined elsewhere.

#include "volterra.hpp"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>

int main()
{
  try {
    fn::Simulation sim;

    std::ifstream infile{"./input.txt"};

    if (!infile) {
      throw std::runtime_error{"Impossible to open input file"};
    }

    double A, B, C, D;

    double delta_t;
    double x0, y0;

    double input_N;

    if ((!(infile >> A >> B >> C >> D >> input_N >> delta_t >> x0 >> y0))

    ) {
      throw std::runtime_error{"Invalid or incomplete data in input.txt"};
    }

    if (input_N != std::floor(input_N)) {
      throw std::runtime_error{"N must be an integer"};
    }

    size_t N = static_cast<size_t>(input_N);

    sim.addParameters({A, B, C, D, N, delta_t});
    sim.addUserPoint({x0, y0});
    sim.evolve();

    fn::Parameters pars = {A, B, C, D, N, delta_t};
    auto res            = sim.result();

    std::ofstream outfile{"results.txt"};

    if (!outfile) {
      throw std::runtime_error{"Impossible to open output file"};
    }

    outfile << std::left << std::setw(30) << "Time, t" << std::setw(30)
            << "Density of preys, x(t)" << std::setw(30)
            << "Density of predators, y(t)" << std::setw(30)
            << "Integral H(x,y)" << "\n";

    auto s = res.size();
    for (size_t i = 0; i < s; ++i) {
      const auto t_i = pars.delta_t * static_cast<double>(i);
      outfile << std::setw(30) << t_i << std::setw(30) << res[i].x
              << std::setw(30) << res[i].y << std::setw(30) << res[i].H << "\n";
    }

  } catch (std::exception const& e) {
    std::cerr << "Caught exception: '" << e.what() << "'\n";
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Caught unknown exception\n";
    return EXIT_FAILURE;
  }
}