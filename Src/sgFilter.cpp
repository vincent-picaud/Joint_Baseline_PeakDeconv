#include "sgFilter.hpp"
#include "vector.hpp"

namespace JointDeconv
{
  // Only used for sequential/joint approache comparison no option
  // to modify sg windows size (that would requires some sg
  // coefficients recomputing)
  //
  void sgFilter(const Vector& spectrum, Vector& filtered)
  {
    const Size_t sg_halfWidth = 3;
    const double sg_coefficient[2 * sg_halfWidth + 1] = {1/7., 1/7., 1/7., 1/7., 1/7., 1/7., 1/7.};

    const Size_t n = spectrum.size();
    assert((n > sg_halfWidth) && "input spectrum is too short for filtering");
    assert((n == filtered.size()));

    //////////////////////////////////////////////
    // Boundary extension (requires input copy) //
    //////////////////////////////////////////////

    Vector buffer(n + 2 * sg_halfWidth);

    for (Index_t i = 0; i < sg_halfWidth; i++)
    {
      buffer[i] = spectrum[sg_halfWidth - i];
    }

    for (Index_t i = 0; i < n; i++)
    {
      buffer[i + sg_halfWidth] = spectrum[i];
    }

    for (Index_t i = 0; i < sg_halfWidth; i++)
    {
      buffer[i + sg_halfWidth + n] = spectrum[n - 2 - i];
    }

    //////////////////
    // Apply filter //
    //////////////////

    for (Index_t i = 0; i < n; i++)
    {
      double sum = 0;
      for (Index_t k = -sg_halfWidth; k <= sg_halfWidth; k++)
      {
        sum += sg_coefficient[k + sg_halfWidth] * buffer[i + k + sg_halfWidth];
      }
      filtered[i] = sum;
    }
  }

}  // JointDeconv
