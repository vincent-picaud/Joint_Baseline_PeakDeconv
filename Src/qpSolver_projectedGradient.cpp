#include "qpSolver_projectedGradient.hpp"

#include <cassert>
#include <cmath>

namespace JointDeconv
{
  namespace
  {

    bool
    check_dimensions(const QP& qp, const Vector& x_lb, const Vector& x_ub, const Vector& x) noexcept
    {
      bool         ok = true;
      const Size_t n  = qp.size();

      ok &= n == x_lb.size();
      ok &= n == x_ub.size();
      ok &= n == x.size();

      return ok;
    }

    bool
    check_init_x(const Vector& x) noexcept
    {
      bool         ok = true;
      const Size_t n  = x.size();
      Index_t      i  = 0;
      for (; i < n; i++)
      {
        if (std::isnan(x[i]))
        {
          ok = false;
          break;
        }
      }
      ok &= (i == n);
      return ok;
    }

    bool
    check_bounds(const Vector& x_lb, const Vector& x_ub) noexcept
    {
      bool ok = true;
      ok &= x_lb.size() == x_ub.size();

      if (ok)
      {
        const Size_t n = x_lb.size();
        Index_t      i = 0;
        for (; i < n; i++)
        {
          if (std::isnan(x_lb[i]) || std::isnan(x_ub[i]))
          {
            break;
          }
          if (x_lb[i] > x_ub[i])
          {
            break;
          }
        }
        ok &= (i == n);
      }
      return ok;
    }

    bool
    check_arguments(const QP& qp, const Vector& x_lb, const Vector& x_ub, const Vector& x) noexcept
    {
      bool ok = true;
      ok &= check_dimensions(qp, x_lb, x_ub, x);

      if (ok)
      {
        ok &= check_init_x(x);
      }

      if (ok)
      {
        ok &= check_bounds(x_lb, x_ub);
      }

      return ok;
    }

    //================================================================)

    void
    x_projection(const Vector& x_lb, const Vector& x_ub, Vector& x) noexcept
    {
      const Size_t n = x_lb.size();

      for (Index_t i = 0; i < n; i++)
      {
        x[i] = std::min(std::max(x[i], x_lb[i]), x_ub[i]);
      }
    }

    void
    grad_projection(const Vector& x_lb, const Vector& x_ub, const Vector& x, Vector& grad) noexcept
    {
      const Size_t n = x_lb.size();

      for (int i = 0; i < n; i++)
      {
        if (((x[i] <= x_lb[i]) && (grad[i] > 0)) ||  // active lower bound
            ((x_ub[i] <= x[i]) && (grad[i] < 0)))    // active upper bound
        {
          grad[i] = 0;
        }
      }
    }
  }  // anonymous

  qpSolver_projectedGradient_outputParameters
  qpSolver_projectedGradient(const QP&                                         qp,
                             const Vector&                                     x_lb,
                             const Vector&                                     x_ub,
                             Vector&                                           x,
                             const qpSolver_projectedGradient_inputParameters* inputParam)
  {
    // Sanity check
    //
    assert(check_arguments(qp, x_lb, x_ub, x));

    // I/O Parameters
    //
    // CAVEAT ok to copy inputParam if it's a lightweight object like here...
    qpSolver_projectedGradient_inputParameters inputParam_actual =
        (inputParam != nullptr) ? *inputParam : qpSolver_projectedGradient_inputParameters();
    inputParam = &inputParam_actual;

    qpSolver_projectedGradient_outputParameters toReturn;

    // Algo begins
    //
    const Size_t n = qp.size();
    Vector       delta_x(n), grad(n), projectedGradient(n);

    x_projection(x_lb, x_ub, x);

    for (toReturn.iter = 0; toReturn.iter < inputParam->iter_max; ++toReturn.iter)
    {
      // Compute projected gradient
      qp.eval_grad(x, grad);
      projectedGradient = grad;
      grad_projection(x_lb, x_ub, x, projectedGradient);
      toReturn.eps    = norm2(projectedGradient);
      const double ng = toReturn.eps * toReturn.eps;
      if (toReturn.eps < inputParam->eps_goal)
      {
        toReturn.status = toReturn.Ok_VanishingProjectedGradient;
        return toReturn;
      }
      // Cauchy step
      qp.eval_Qx(projectedGradient, delta_x);
      const double ng2   = dot(projectedGradient, delta_x);
      const double alpha = ng / ng2;
      if (!std::isfinite(alpha) || (alpha < 0))
      {
        std::cerr << "\n#Error: α is not finite or negative, α=" << alpha;
        toReturn.status = toReturn.Failed_CauchyStepIsNaNorNeg;
        throw toReturn;
      }
      // delta_x = -alpha.projGrad
      delta_x = projectedGradient;
      scal(-alpha, delta_x);
      toReturn.eps = norm2(delta_x);
      // x += proj(x+delta_x)
      axpy(1, delta_x, x);
      x_projection(x_lb, x_ub, x);
      if (toReturn.eps < inputParam->eps_goal)
      {
        toReturn.status = toReturn.Ok_VanishingStep;
        return toReturn;
      }
    }
    toReturn.status = toReturn.Failed_MaxIter;
    throw toReturn;
  }
}
