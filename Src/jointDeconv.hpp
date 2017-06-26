// Deconvolution Algorithm entry point
#pragma once

#include "matrix.hpp"
#include "qpSolver_projectedGradient.hpp"

namespace JointDeconv
{
  struct JointDeconv_InputParameters
  {
    double lambda_1 = 0.1, lambda_2 = 0, mu = 100;

    //================================================================

    // CAVEAT: for the moment only one solver,
    //         -> in case of several solvers use std::variant (gcc-7) to store input param
    //
    // enum class SolverEnum
    // {
    //   ProjectedGradient
    //   ....
    // }
    //
    qpSolver_projectedGradient_inputParameters solver_inputParameters;
  };

  void jointDeconv_GaussianPeaks(const Vector&                      x,
                                 const Vector&                      y,
                                 const double                       yBaseline_left,
                                 const double                       yBaseline_right,
                                 const double                       peak_minimal_height,
                                 const double                       sigma_left,
                                 const double                       sigma_right,
                                 Vector&                            deconvolvedPeak,
                                 Vector&                            convolvedPeak,
                                 Vector&                            baseline,
                                 const JointDeconv_InputParameters& inputParameters);

  void jointDeconv_GaussianPeaks(const Vector&                      y,
                                 const double                       yBaseline_left,
                                 const double                       yBaseline_right,
                                 const double                       peak_minimal_height,
                                 const double                       sigma_left,
                                 const double                       sigma_right,
                                 Vector&                            deconvolvedPeak,
                                 Vector&                            convolvedPeak,
                                 Vector&                            baseline,
                                 const JointDeconv_InputParameters& inputParameters);
}
