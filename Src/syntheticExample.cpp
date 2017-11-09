#include "syntheticExample.hpp"

#include <cmath>
#include <iomanip>
#include <random>

namespace JointDeconv
{
  SyntheticExample create_syntheticExample(const int c,
                                           const double sigma_noise,
                                           const int randomGenerator_seed)
  {
    assert(sigma_noise >= 0);
    assert((c == 0) || (c == -1) || (c == 1));

    const Size_t n = 500;
    SyntheticExample to_return = {.peaks = RegularMatrix({{1, 50, 10},
                                                          {0.5, 90, 10},
                                                          {0.5, 170, 10},
                                                          {3, 200, 10},
                                                          {2, 230, 10},
                                                          {1, 260, 10},
                                                          {0.5, 350, 10},
                                                          {3, 370, 10},
                                                          {2, 390, 10},
                                                          {1, 410, 10}}),

                                  .y = Vector(n),
                                  .y_baseline = Vector(n),
                                  .y_peak = Vector(n),
                                  .y_noise = Vector(n),
                                  .y_first = 0,
                                  .y_last = 0};

    // Creates baseline

    const double s = (c == 1) ? 5 : 2;
    for (Index_t i = 0; i < n; ++i)
    {
      to_return.y_baseline[i] =
          s + c * exp(-3. * i / ((double)n)) - 2 * i / ((double)n);
    }

    // Creates peaks
    const Size_t n_peak = to_return.peaks.I_size();

    for (Index_t peak_idx = 0; peak_idx < n_peak; peak_idx++)
    {
      const auto peak = [&](const Index_t i) {
        return to_return.peaks(peak_idx, 0) *
               exp(-0.5 * pow((i - to_return.peaks(peak_idx, 1)) /
                                  to_return.peaks(peak_idx, 2),
                              2));
      };

      for (Index_t i = 0; i < n; ++i)
      {
        to_return.y_peak[i] += peak(i);
      }
    }

    // Add noise
    std::random_device rd;
    std::mt19937 gen(rd());
    gen.seed(randomGenerator_seed);
    std::normal_distribution<> d(0, sigma_noise);
    for (Index_t i = 0; i < n; ++i)
    {
      to_return.y_noise[i] = d(gen);
    }
    // No noise at boundary
    // because the method use them as hard constraints.
    // (with real spectra we use smoothed ones).
    //
    to_return.y_noise[0] = 0;
    to_return.y_noise[n - 1] = 0;

    // Sum all contributions
    for (Index_t i = 0; i < n; ++i)
    {
      to_return.y[i] =
          to_return.y_peak[i] + to_return.y_baseline[i] + to_return.y_noise[i];
    }

    // Finalize structure initialization
    to_return.y_first = to_return.y_baseline[0];
    to_return.y_last = to_return.y_baseline[n - 1];

    return to_return;
  }

  std::ostream& operator<<(std::ostream& out, const SyntheticExample& toPrint)
  {
    const auto old_settings = out.flags();
    const auto old_precision = out.precision();
    out << std::setprecision(8);

    const Size_t n = toPrint.y.size();

    for (Index_t i = 0; i < n; ++i)
    {
      out << "\n"
          << std::setw(5) << i << std::setw(20) << toPrint.y[i] << std::setw(20)
          << toPrint.y_baseline[i] << std::setw(20) << toPrint.y_peak[i]
          << std::setw(20) << toPrint.y_noise[i];
    }

    out.flags(old_settings);
    out.precision(old_precision);

    return out;
  }
}
