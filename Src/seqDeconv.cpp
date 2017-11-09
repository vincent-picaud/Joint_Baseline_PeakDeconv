#include "seqDeconv.hpp"
#include "jointDeconv_details.hpp"
#include "matrix.hpp"
#include "qpDenseExplicit.hpp"
#include "qpSolver_projectedGradient.hpp"

#include <cmath>
#include <limits>

#include <lapacke.h>

namespace JointDeconv
{
  void seqDeconv_GaussianPeaks(const Vector& x,
                               const Vector& y_unscaled,
                               const double peak_minimal_height_unscaled,
                               const double sigma_first,
                               const double sigma_last,
                               Vector& deconvolvedPeak,
                               Vector& convolvedPeak,
                               const SeqDeconv_InputParameters& inputParameters)
  {
    // Sanity check
    //
    assert(x.size() >= 2);
    assert(x.size() == y_unscaled.size());
    assert(x.size() == deconvolvedPeak.size());
    assert(x.size() == convolvedPeak.size());
    assert(sigma_first > 0);
    assert(sigma_last > 0);

    ////////////////////
    // Initialization //
    ////////////////////

    const Size_t n = x.size();

    // scale y
    const double y_max = max(y_unscaled);
    const double y_min = min(y_unscaled);

    if (y_max == y_min)
    {
      // Constant signal
      deconvolvedPeak = 0;
      convolvedPeak = 0;
      return;
    }

    Vector y(n);
    y = y_unscaled;
    //    y -= y_min; <- CAVEAT we assumed that y has no baseline, but
    //    its noise is centered: do not do =-y_min
    y /= (y_max - y_min);
    const double peak_minimal_height =
        peak_minimal_height_unscaled / (y_max - y_min);

    // Compute L (once for all)
    //
    const RegularMatrix L = Details::create_L(x, sigma_first, sigma_last);

    // Vector/Matrices that depend on λ1 or λ2 need rebuild between
    // the two steps
    //
    QpDenseExplicit qp;

    const auto build = [&](const double lambda_1, const double lambda_2) {

      SymmetricMatrix Q;
      Vector q;

      Q = Details::create_Q(L, lambda_2);
      q = Details::create_q(L, lambda_1, y);

      qp = QpDenseExplicit(std::move(Q), std::move(q));
    };

    // Trigger build
    //
    build(inputParameters.lambda_1, inputParameters.lambda_2);

    // Bounds constraints
    //
    Vector deconvolvedPeak_lb(n);
    deconvolvedPeak_lb = 0;
    Vector deconvolvedPeak_ub(n);
    deconvolvedPeak_ub = std::numeric_limits<double>::max();

    // Solve it!
    //
    qpSolver_projectedGradient_outputParameters outputParam;

    try
    {
      outputParam =
          qpSolver_projectedGradient(qp,
                                     deconvolvedPeak_lb,
                                     deconvolvedPeak_ub,
                                     deconvolvedPeak,
                                     &inputParameters.solver_inputParameters);
    }
    catch (qpSolver_projectedGradient_outputParameters& e_outputParam)
    {
      std::cerr << "\n#Error during step 1: " << e_outputParam << std::endl;
      outputParam = e_outputParam;
    }

    ////////////////////////////
    // Support regularization //
    ////////////////////////////

    deconvolvedPeak_ub = 0;

    for (Index_t i = 1; i + 1 < n; ++i)
    {
      if ((deconvolvedPeak[i] > peak_minimal_height) &&
          (((deconvolvedPeak[i] > deconvolvedPeak[i - 1]) &&
            (deconvolvedPeak[i] >= deconvolvedPeak[i + 1])) ||
           ((deconvolvedPeak[i] >= deconvolvedPeak[i - 1]) &&
            (deconvolvedPeak[i] > deconvolvedPeak[i + 1]))))
      {
        deconvolvedPeak_ub[i] = std::numeric_limits<double>::max();
      }
    }

    // Recreate matrix with λ1=0 and λ2=0
    build(0, 0);

    // Resolve the pb (step2)
    //
    try
    {
      outputParam =
          qpSolver_projectedGradient(qp,
                                     deconvolvedPeak_lb,
                                     deconvolvedPeak_ub,
                                     deconvolvedPeak,
                                     &inputParameters.solver_inputParameters);
    }
    catch (qpSolver_projectedGradient_outputParameters& e_outputParam)
    {
      std::cerr << "\n#Error during step 2:" << e_outputParam << std::endl;
      outputParam = e_outputParam;
    }

    ////////////////////////////
    // Compute convolvedPeak  //
    ////////////////////////////

    // rescale solution:
    deconvolvedPeak *= (y_max - y_min);
    // convolved peak L.xp
    Mv(1, Identity_c, L, deconvolvedPeak, 0, convolvedPeak);
  }

  void seqDeconv_GaussianPeaks(const Vector& y_unscaled,
                               const double peak_minimal_height,
                               const double sigma_first,
                               const double sigma_last,
                               Vector& deconvolvedPeak,
                               Vector& convolvedPeak,
                               const SeqDeconv_InputParameters& inputParameters)
  {
    const Size_t n = y_unscaled.size();
    Vector x(n);
    for (Index_t i = 0; i < n; ++i)
    {
      x[i] = i;
    }

    return seqDeconv_GaussianPeaks(x,
                                   y_unscaled,
                                   peak_minimal_height,
                                   sigma_first,
                                   sigma_last,
                                   deconvolvedPeak,
                                   convolvedPeak,
                                   inputParameters);
  }

}  // JointDeconv
