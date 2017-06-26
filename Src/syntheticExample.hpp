// Generate article synthetic examples
#pragma once

#include "matrix.hpp"

namespace JointDeconv
{
  struct SyntheticExample
  {
    RegularMatrix peaks;  // True peaks, each row contains {h,μ,σ} of the Gaussian peak

    Vector y, y_baseline, y_peak, y_noise;

    double y_first, y_last;

    friend std::ostream& operator<<(std::ostream& out, const SyntheticExample& toPrint);
  };

  SyntheticExample create_syntheticExample(const int c, const double sigma_noise);
}
