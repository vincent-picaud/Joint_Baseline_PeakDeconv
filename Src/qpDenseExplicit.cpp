#include "qpDenseExplicit.hpp"

namespace JointDeconv
{
  QpDenseExplicit::QpDenseExplicit(SymmetricMatrix&& Q, Vector&& q) noexcept : Q_(std::move(Q)), q_(std::move(q))
  {
    assert(q_.size() == Q_.I_size());
    assert(q_.size() == Q_.J_size());
  }

  Size_t
  QpDenseExplicit::size() const noexcept
  {
    return q_.size();
  }

  void
  QpDenseExplicit::interface_eval_grad(const Vector& x, Vector& grad) const noexcept
  {
    grad = q_;
    Mv(1, Q_, x, 1, grad);
  }
  void
  QpDenseExplicit::interface_eval_Qx(const Vector& x, Vector& Qx) const noexcept
  {
    Mv(1, Q_, x, 0, Qx);
  }
}
