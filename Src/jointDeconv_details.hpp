#pragma once

namespace JointDeconv
{
  namespace Details
  {
    // Create L, PSF convolution operator with a Gaussian
    // \sigma_shape = first at position x[0]
    // \sigma_shape = last  at position x[end-1]
    //
    RegularMatrix create_L(const Vector& x,
                           const double sigma_first_slot,
                           const double sigma_last_slot) noexcept;

    // Creates Q = λ2.I + L^t.Aμ.L
    //
    SymmetricMatrix create_Q(const RegularMatrix& L,
                             const SymmetricMatrix& tildeAmu,
                             const double lambda_2) noexcept;

    // Creates Q = λ2.I + L^t.L
    //
    SymmetricMatrix create_Q(const RegularMatrix& L,
                             const double lambda_2) noexcept;

    // Compute \tilde{y} (see publication)
    //
    void compute_ytilde(const double mu,
                        const Vector& y,
                        const double y_first,
                        const double y_last,
                        Vector& tilde_y) noexcept;

    // Creates q = λ1.I - L^t.Aμ.~y - L^t(y-~y)
    //
    Vector create_q(const RegularMatrix& L,
                    const SymmetricMatrix& tildeAmu,
                    const double mu,
                    const double lambda_1,
                    const Vector& y,
                    const double y_first,
                    const double y_last) noexcept;

  }  // Details

}  // JointDeconv
