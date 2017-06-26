// Defines quadratic problem interface
#pragma once

#include "vector.hpp"

namespace JointDeconv
{
  class QP
  {
   public:
    virtual ~QP(){};
    virtual Size_t size() const noexcept = 0;

    void
    eval_grad(const Vector& x, Vector& grad) const noexcept
    {
      assert(x.size() == size());
      assert(x.size() == grad.size());
      interface_eval_grad(x, grad);
    }
    void
    eval_Qx(const Vector& x, Vector& Qx) const noexcept
    {
      assert(x.size() == size());
      assert(x.size() == Qx.size());
      interface_eval_Qx(x, Qx);
    }

   protected:
    virtual void interface_eval_grad(const Vector& x, Vector& grad) const noexcept = 0;
    virtual void interface_eval_Qx(const Vector& x, Vector& Qx) const noexcept     = 0;
  };
}
