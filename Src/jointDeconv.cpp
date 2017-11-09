#include "jointDeconv.hpp"
#include "jointDeconv_details.hpp"
#include "matrix.hpp"
#include "qpDenseExplicit.hpp"
#include "qpSolver_projectedGradient.hpp"

#include <cmath>
#include <limits>

#include <lapacke.h>

namespace JointDeconv
{
  namespace
  {
    // Creates L
    //
    // -> see jointDeconv_details.hpp
    //

    // Creates tilde_B_mu, smooth semiparametric baseline (TODO: add B-spline)
    //
    SymmetricMatrix create_tildeBmu_semiParametric(const int n, const double mu)
    {
      assert(n >= 2);
      assert(mu > 0);

      SymmetricMatrix tildeBmu(n, n);
      tildeBmu = 0;

      // Generic terme
      for (Index_t i = 0; i < n; ++i)
      {
        tildeBmu(i, i) = 1 + 2 * mu;
      }
      for (Index_t i = 1; i < n; ++i)
      {
        tildeBmu(i, i - 1) = -mu;
      }

      // Boundaries
      tildeBmu(0, 0) = 1 + mu;
      tildeBmu(1, 0) = 0;

      tildeBmu(n - 1, n - 2) = 0;
      tildeBmu(n - 1, n - 1) = 1 + mu;

      return tildeBmu;
    }

    struct Lapack_exception : std::exception
    {
      const char* what() const throw() final
      {
        static const char error_message_buffer[] = "Lapack exception";

        return error_message_buffer;
      }
    };

    SymmetricMatrix inverse_tildeBmu(const SymmetricMatrix& tildeBmu)
    {
      const Size_t n = tildeBmu.I_size();

      // Cholesky L^t.L <- B
      //
      SymmetricMatrix B_LtL = tildeBmu;

      lapack_int lapack_info =
          LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', n, B_LtL.data(), B_LtL.ld());
      if (lapack_info != 0)
      {
        std::cerr << "\n#Error Q matrix is too badly conditioned (lambda_1 or "
                     "mu penalty too high?)"
                  << std::endl;
        throw Lapack_exception();
      };

      // B <- invB (inplace)
      //
      RegularMatrix invB(n, n);
      invB = 0;
      for (Index_t i = 0; i < n; i++)
      {
        invB(i, i) = 1;
      };
      lapack_info = LAPACKE_dpotrs(LAPACK_COL_MAJOR,
                                   'L',
                                   n,
                                   n,
                                   B_LtL.data(),
                                   B_LtL.ld(),
                                   invB.data(),
                                   invB.ld());
      if (lapack_info != 0)
      {
        std::cerr << "\n#Error Q matrix is too badly conditioned (lambda_1 or "
                     "mu penalty too high?)"
                  << std::endl;
        throw Lapack_exception();
      };

      return invB;  // note: SymmetricMatrix(RegularMatrix&&) do the job
    }

    // A=I-tildeBmu^-1
    SymmetricMatrix create_tildeAmu(const SymmetricMatrix& inverse_tildeBmu)
    {
      SymmetricMatrix tildeAmu(inverse_tildeBmu);
      scal(-1, tildeAmu);
      const Size_t n = tildeAmu.I_size();
      for (Index_t i = 0; i < n; i++)
      {
        tildeAmu(i, i) += 1;
      };

      return tildeAmu;  // note: SymmetricMatrix(RegularMatrix&&) do the job
    }

    // Creates Q = λ2.I + L^t.Aμ.L
    //
    // -> see jointDeconv_details.hpp
    //

    // Creates q = λ1.I - L^t.Aμ.~y - L^t(y-~y)
    //
    // -> see jointDeconv_details.hpp
    //
  }  // anonymous

  void jointDeconv_GaussianPeaks(
      const Vector& x,
      const Vector& y_unscaled,
      const double yBaseline_left_unscaled,
      const double yBaseline_right_unscaled,
      const double peak_minimal_height_unscaled,
      const double sigma_first,
      const double sigma_last,
      Vector& deconvolvedPeak,
      Vector& convolvedPeak,
      Vector& baseline,
      const JointDeconv_InputParameters& inputParameters)
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
    const double yBaseline_left =
        (yBaseline_left_unscaled - y_min) / (y_max - y_min);
    const double yBaseline_right =
        (yBaseline_right_unscaled - y_min) / (y_max - y_min);
    const double peak_minimal_height =
        peak_minimal_height_unscaled / (y_max - y_min);

    // Compute L (once for all)
    //
    const RegularMatrix L = Details::create_L(x, sigma_first, sigma_last);

    // Vector/Matrices that depend on μ, λ1 or λ2 need rebuild between the two
    // steps
    //
    SymmetricMatrix tildeBmu;
    SymmetricMatrix inv_tildeBmu;
    SymmetricMatrix tildeAmu;

    // SymmetricMatrix Q;
    // Vector          q;
    QpDenseExplicit qp;

    const auto build =
        [&](const double mu, const double lambda_1, const double lambda_2) {
          tildeBmu = create_tildeBmu_semiParametric(n, mu);
          inv_tildeBmu = inverse_tildeBmu(tildeBmu);
          tildeAmu = create_tildeAmu(inv_tildeBmu);

          SymmetricMatrix Q;
          Vector q;

          Q = Details::create_Q(L, tildeAmu, lambda_2);
          q = Details::create_q(
              L, tildeAmu, mu, lambda_1, y, yBaseline_left, yBaseline_right);

          qp = QpDenseExplicit(std::move(Q), std::move(q));
        };

    // Trigger build
    //
    build(
        inputParameters.mu, inputParameters.lambda_1, inputParameters.lambda_2);

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

    // Recreate matrix with lambda_1 & 2 = 0
    build(inputParameters.mu, 0, 0);

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

    //////////////////////////////////////
    // Compute convolvedPeak & baseline //
    //////////////////////////////////////

    // rescale solution: CAVEAT do NOT translate by += y_min
    deconvolvedPeak *= (y_max - y_min);
    // convolved peak L.xp
    Mv(1, Identity_c, L, deconvolvedPeak, 0, convolvedPeak);
    // baseline B^{-1}(tilde_y − L.xp)
    Vector buffer(n);
    Details::compute_ytilde(inputParameters.mu,
                            y_unscaled,
                            yBaseline_left_unscaled,
                            yBaseline_right_unscaled,
                            buffer);
    buffer -= convolvedPeak;
    Mv(1, inv_tildeBmu, buffer, 0, baseline);
  }

  void jointDeconv_GaussianPeaks(
      const Vector& y_unscaled,
      const double y_unscaled_first,
      const double y_unscaled_last,
      const double peak_minimal_height,
      const double sigma_first,
      const double sigma_last,
      Vector& deconvolvedPeak,
      Vector& convolvedPeak,
      Vector& baseline,
      const JointDeconv_InputParameters& inputParameters)
  {
    const Size_t n = y_unscaled.size();
    Vector x(n);
    for (Index_t i = 0; i < n; ++i)
    {
      x[i] = i;
    }

    return jointDeconv_GaussianPeaks(x,
                                     y_unscaled,
                                     y_unscaled_first,
                                     y_unscaled_last,
                                     peak_minimal_height,
                                     sigma_first,
                                     sigma_last,
                                     deconvolvedPeak,
                                     convolvedPeak,
                                     baseline,
                                     inputParameters);
  }
}
