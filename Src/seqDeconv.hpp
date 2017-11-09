// Deconvolution Algorithm entry point
#pragma once

#include "matrix.hpp"
#include "qpSolver_projectedGradient.hpp"

namespace JointDeconv
{
  struct seqDeconv_InputParameters
  {
    double lambda_1 = 0.1, lambda_2 = 0;
    // here snip width

    qpSolver_projectedGradient_inputParameters solver_inputParameters;
  };

  void seqDeconv_GaussianPeaks(
      const Vector& x,
      const Vector& y,
      const double peak_minimal_height,
      const double sigma_left,
      const double sigma_right,
      Vector& deconvolvedPeak,
      Vector& convolvedPeak,
      Vector& baseline,
      const seqDeconv_InputParameters& inputParameters);

  void seqDeconv_GaussianPeaks(
      const Vector& y,
      const double peak_minimal_height,
      const double sigma_left,
      const double sigma_right,
      Vector& deconvolvedPeak,
      Vector& convolvedPeak,
      Vector& baseline,
      const seqDeconv_InputParameters& inputParameters);
}
