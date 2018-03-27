#include <cmath>
#include <limits>

#include "matrix.hpp"

#include <lapacke.h>

namespace JointDeconv
{
  namespace Details
  {
    // Creates L matrix
    //
    RegularMatrix create_L(const Vector& x,
                           const double sigma_first_slot,
                           const double sigma_last_slot) noexcept
    {
      assert(x.size() >= 2);
      assert(sigma_first_slot > 0);
      assert(sigma_last_slot > 0);

      const Size_t n = x.size();
      RegularMatrix L(n, n);

      const auto sigma_interpolation = [&](const double xi) {
        const double w0 = (x[n - 1] - xi) / (x[n - 1] - x[0]);
        const double w1 = (xi - x[0]) / (x[n - 1] - x[0]);
        return w0 * sigma_first_slot + w1 * sigma_last_slot;
      };

      for (Index_t i = 0; i < n; ++i)
      {
        const double sigma_i = sigma_interpolation(x[i]);

        for (Index_t j = 0; j < n; ++j)
        {
          double gaussian_arg = (x[j] - x[i]) / sigma_i;
          gaussian_arg *= -0.5 * gaussian_arg;

          L(i, j) = exp(gaussian_arg);
        }
      }
      return L;
    }

    SymmetricMatrix create_Q(const RegularMatrix& L,
                             const SymmetricMatrix& tildeAmu,
                             const double lambda_2) noexcept
    {
      const Size_t n = tildeAmu.I_size();

      RegularMatrix buffer(n, n);
      RegularMatrix toReturn(n, n);

      MM(1, tildeAmu, Identity_c, L, 0, buffer);
      MM(1, Transpose_c, L, Identity_c, buffer, 0, toReturn);

      for (Index_t i = 0; i < n; i++)
      {
        toReturn(i, i) += lambda_2;
      }

      return toReturn;
    }

    // Creates Q = λ2.I + L^t.L
    //
    SymmetricMatrix create_Q(const RegularMatrix& L,
                             const double lambda_2) noexcept
    {
      const Size_t n = L.J_size();

      RegularMatrix toReturn(n, n);

      MM(1, Transpose_c, L, Identity_c, L, 0, toReturn);

      for (Index_t i = 0; i < n; i++)
      {
        toReturn(i, i) += lambda_2;
      }

      return toReturn;
    }

    // Compute \tilde{y} (see publication)
    //
    void compute_ytilde(const double mu,
                        const Vector& y,
                        const double y_first,
                        const double y_last,
                        Vector& tilde_y) noexcept
    {
      assert(y.size() >= 2);

      const size_t n = y.size();
      tilde_y = y;
      tilde_y[0] = (1 + mu) * y_first;
      tilde_y[1] += mu * y_first;
      tilde_y[n - 2] += mu * y_last;
      tilde_y[n - 1] = (1 + mu) * y_last;
    }

    // Creates q = λ1.I - L^t.Aμ.~y - L^t(y-~y)
    //           = λ1.I - L^t.( Aμ.~y + (y-~y) )
    //
    Vector create_q(const RegularMatrix& L,
                    const SymmetricMatrix& tildeAmu,
                    const double mu,
                    const double lambda_1,
                    const Vector& y,
                    const double y_first,
                    const double y_last) noexcept
    {
      // Sanity...
      //
      assert(y.size() >= 2);

      const Size_t n = y.size();
      Vector ty(n), y_ty(n), q(n);

      // ~y
      compute_ytilde(mu, y, y_first, y_last, ty);

      // y-~y
      y_ty = y;
      y_ty -= ty;

      // y_ty = Aμ.~y + (y-~y)
      Mv(1, tildeAmu, ty, 1, y_ty);
      q = lambda_1;
      // λ1.I - L^t.( Aμ.~y + (y-~y) )
      Mv(-1, Transpose_c, L, y_ty, 1, q);

      return q;
    }

    // Creates q = λ1.I - L^t.y
    //
    Vector create_q(const RegularMatrix& L,
                    const double lambda_1,
                    const Vector& y) noexcept
    {
      // Sanity...
      //
      assert(y.size() >= 2);

      Vector q(y.size());
      q = lambda_1;
      Mv(-1, Transpose_c, L, y, 1, q);

      return q;
    }

  }  // Details
}
