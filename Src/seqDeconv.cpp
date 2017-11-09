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
                               Vector& baseline,
                               const seqDeconv_InputParameters& inputParameters)
  {
    // Sanity check
    //
    assert(x.size() >= 2);
    assert(x.size() == y_unscaled.size());
    assert(x.size() == deconvolvedPeak.size());
    assert(x.size() == convolvedPeak.size());
    assert(x.size() == baseline.size());
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
      baseline = y_max;
      return;
    }

    Vector y(n);
    y = y_unscaled;
    y -= y_min;
    y /= (y_max - y_min);
    const double peak_minimal_height =
        peak_minimal_height_unscaled / (y_max - y_min);

    // Compute L (once for all)
    //
    const RegularMatrix L = Details::create_L(x, sigma_first, sigma_last);

    // // Vector/Matrices that depend on μ, λ1 or λ2 need rebuild between the
    // two
    // // steps
    // //
    // SymmetricMatrix tildeBmu;
    // SymmetricMatrix inv_tildeBmu;
    // SymmetricMatrix tildeAmu;

    // // SymmetricMatrix Q;
    // // Vector          q;
    // QpDenseExplicit qp;

    // const auto build =
    //     [&](const double mu, const double lambda_1, const double lambda_2) {
    //       tildeBmu = create_tildeBmu_semiParametric(n, mu);
    //       inv_tildeBmu = inverse_tildeBmu(tildeBmu);
    //       tildeAmu = create_tildeAmu(inv_tildeBmu);

    //       SymmetricMatrix Q;
    //       Vector q;

    //       Q = Details::create_Q(L, tildeAmu, lambda_2);
    //       q = Details::create_q(
    //           L, tildeAmu, mu, lambda_1, y, yBaseline_left, yBaseline_right);

    //       qp = QpDenseExplicit(std::move(Q), std::move(q));
    //     };

    // // Trigger build
    // //
    // build(
    //     inputParameters.mu, inputParameters.lambda_1,
    //     inputParameters.lambda_2);

    // // Bounds constraints
    // //
    // Vector deconvolvedPeak_lb(n);
    // deconvolvedPeak_lb = 0;
    // Vector deconvolvedPeak_ub(n);
    // deconvolvedPeak_ub = std::numeric_limits<double>::max();

    // // Solve it!
    // //
    // qpSolver_projectedGradient_outputParameters outputParam;

    // try
    // {
    //   outputParam =
    //       qpSolver_projectedGradient(qp,
    //                                  deconvolvedPeak_lb,
    //                                  deconvolvedPeak_ub,
    //                                  deconvolvedPeak,
    //                                  &inputParameters.solver_inputParameters);
    // }
    // catch (qpSolver_projectedGradient_outputParameters& e_outputParam)
    // {
    //   std::cerr << "\n#Error during step 1: " << e_outputParam << std::endl;
    //   outputParam = e_outputParam;
    // }

    // ////////////////////////////
    // // Support regularization //
    // ////////////////////////////

    // deconvolvedPeak_ub = 0;

    // for (Index_t i = 1; i + 1 < n; ++i)
    // {
    //   if ((deconvolvedPeak[i] > peak_minimal_height) &&
    //       (((deconvolvedPeak[i] > deconvolvedPeak[i - 1]) &&
    //         (deconvolvedPeak[i] >= deconvolvedPeak[i + 1])) ||
    //        ((deconvolvedPeak[i] >= deconvolvedPeak[i - 1]) &&
    //         (deconvolvedPeak[i] > deconvolvedPeak[i + 1]))))
    //   {
    //     deconvolvedPeak_ub[i] = std::numeric_limits<double>::max();
    //   }
    // }

    // // Recreate matrix with lambda_1 & 2 = 0
    // build(inputParameters.mu, 0, 0);

    // // Resolve the pb (step2)
    // //
    // try
    // {
    //   outputParam =
    //       qpSolver_projectedGradient(qp,
    //                                  deconvolvedPeak_lb,
    //                                  deconvolvedPeak_ub,
    //                                  deconvolvedPeak,
    //                                  &inputParameters.solver_inputParameters);
    // }
    // catch (qpSolver_projectedGradient_outputParameters& e_outputParam)
    // {
    //   std::cerr << "\n#Error during step 2:" << e_outputParam << std::endl;
    //   outputParam = e_outputParam;
    // }

    // //////////////////////////////////////
    // // Compute convolvedPeak & baseline //
    // //////////////////////////////////////

    // // rescale solution: CAVEAT do NOT translate by += y_min
    // deconvolvedPeak *= (y_max - y_min);
    // // convolved peak L.xp
    // Mv(1, Identity_c, L, deconvolvedPeak, 0, convolvedPeak);
    // // baseline B^{-1}(tilde_y − L.xp)
    // Vector buffer(n);
    // Details::compute_ytilde(inputParameters.mu,
    //                         y_unscaled,
    //                         yBaseline_left_unscaled,
    //                         yBaseline_right_unscaled,
    //                         buffer);
    // buffer -= convolvedPeak;
    // Mv(1, inv_tildeBmu, buffer, 0, baseline);
  }

  void seqDeconv_GaussianPeaks(const Vector& y_unscaled,
                               const double peak_minimal_height,
                               const double sigma_first,
                               const double sigma_last,
                               Vector& deconvolvedPeak,
                               Vector& convolvedPeak,
                               Vector& baseline,
                               const seqDeconv_InputParameters& inputParameters)
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
                                   baseline,
                                   inputParameters);
  }

}  // JointDeconv
