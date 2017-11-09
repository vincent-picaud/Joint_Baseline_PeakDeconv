#include "snip.hpp"

namespace JointDeconv
{
  void snip(const Vector& spectrum, Vector& baseline, const Size_t halfWindowSize)
  {
    const Size_t n = spectrum.size();
    Vector buffer(n);
    baseline = spectrum;
    buffer = spectrum;

    for (Index_t i = halfWindowSize; i > 0; --i)
    {
      for (Index_t j = i; j < n - i; ++j)
      {
        double a = buffer[j];
        double b = (buffer[j - i] + buffer[j + i]) / 2;
        if (b < a)
        {
          a = b;
        }
        baseline[j] = a;
      }

      for (Index_t j = i; j < n - i; ++j)
      {
        buffer[j] = baseline[j];
      }
    }
  }
}  // JointDeconv
