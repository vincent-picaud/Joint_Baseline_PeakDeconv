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
    const Size_t sg_halfWidth = 19; // as in the article sg_width = 2*19+1 = 39
    const double sg_coefficient[2 * sg_halfWidth + 1] = {
        -0.03377110694183864,  -0.02439024390243901,  -0.015516454540844778,
        -0.007149738857055922, 0.0007099031489275417, 0.008062471477105627,
        0.014907966127478321,  0.021246387100045637,  0.027077734394807564,
        0.03240200801176411,   0.03721920795091527,   0.041529334212261036,
        0.04533238679580143,   0.048628365701536426,  0.05141727092946605,
        0.053699102479590276,  0.05547386035190913,   0.056741544546422586,
        0.05750215506313067,   0.057755691902033356,  0.05750215506313067,
        0.056741544546422586,  0.05547386035190913,   0.053699102479590276,
        0.05141727092946605,   0.048628365701536426,  0.04533238679580143,
        0.041529334212261036,  0.03721920795091527,   0.03240200801176411,
        0.027077734394807564,  0.021246387100045637,  0.014907966127478321,
        0.008062471477105627,  0.0007099031489275417, -0.007149738857055922,
        -0.015516454540844778, -0.02439024390243901,  -0.03377110694183864};

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
