// Solves a bound constrained problem using a projected gradient approach
#pragma once

#include "qp.hpp"

#include <iostream>
#include <limits>

namespace JointDeconv
{
  struct qpSolver_projectedGradient_inputParameters
  {
    Size_t iter_max = 5000;
    double eps_goal = 1e-12;
  };

  struct qpSolver_projectedGradient_outputParameters
  {
    Size_t iter = 0;
    double eps  = std::numeric_limits<double>::max();

    enum ExitEnum
    {
      Unitialized,
      Ok_VanishingProjectedGradient,
      Ok_VanishingStep,
      Failed_MaxIter,
      Failed_CauchyStepIsNaNorNeg
    };

    bool
    converged() const noexcept
    {
      return (status == Ok_VanishingProjectedGradient) || (status == Ok_VanishingStep);
    }

    ExitEnum status = Unitialized;
  };

  static inline std::ostream&
  operator<<(std::ostream& out, qpSolver_projectedGradient_outputParameters& toPrint)
  {
    out << "\niter: " << toPrint.iter << ", eps: " << toPrint.eps << ", status: ";

    switch (toPrint.status)
    {
      case qpSolver_projectedGradient_outputParameters::Unitialized:
        out << "Unitializedn";
        break;
      case qpSolver_projectedGradient_outputParameters::Ok_VanishingProjectedGradient:
        out << "Ok_VanishingProjectedGradient";
        break;
      case qpSolver_projectedGradient_outputParameters::Ok_VanishingStep:
        out << "Ok_VanishingStep";
        break;
      case qpSolver_projectedGradient_outputParameters::Failed_MaxIter:
        out << "Failed_MaxIter";
        break;
      case qpSolver_projectedGradient_outputParameters::Failed_CauchyStepIsNaNorNeg:
        out << "Failed_CauchyStepIsNaNorNeg";
        break;
      default:
        out << "INTERNAL ERROR";
        break;
    }

    return out;
  }

  //================================================================

  qpSolver_projectedGradient_outputParameters qpSolver_projectedGradient(
      const QP&                                         qp,
      const Vector&                                     x_lb,
      const Vector&                                     x_ub,
      Vector&                                           x,
      const qpSolver_projectedGradient_inputParameters* inputParam = nullptr);
}
