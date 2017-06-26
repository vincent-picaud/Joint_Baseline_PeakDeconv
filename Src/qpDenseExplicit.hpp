// Defines a qp specialization in case of explicit dense Q,q
#pragma once

#include "matrix.hpp"
#include "qp.hpp"

namespace JointDeconv
{
  class QpDenseExplicit final : public QP
  {
   public:
    QpDenseExplicit() noexcept = default;
    QpDenseExplicit(SymmetricMatrix&& Q, Vector&& q) noexcept;

    Size_t size() const noexcept;

   protected:
    void interface_eval_grad(const Vector& x, Vector& grad) const noexcept;
    void interface_eval_Qx(const Vector& x, Vector& Qx) const noexcept;

   private:
    SymmetricMatrix Q_;
    Vector          q_;
  };
}
