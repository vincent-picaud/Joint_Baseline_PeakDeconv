#include "vector.hpp"

#include <cblas.h>
#include <limits>

namespace JointDeconv
{
  Vector::Vector() noexcept : data_(), n_(0)
  {
  }
  Vector::Vector(const Size_t n) : data_(new double[n]), n_(n)
  {
  }
  Vector::Vector(Vector&& toMove) noexcept : data_(std::move(toMove.data_)), n_(toMove.n_)
  {
    toMove.n_ = 0;
  }

  //================================================================

  Vector&
  Vector::operator=(const Vector& toCopy) noexcept
  {
    assert(size() == toCopy.size());
    assert((inc() == 1) && (toCopy.inc() == 1));  // otherwise strided copy must be used
    std::copy(toCopy.data(), toCopy.data() + size(), data());
    return *this;
  }
  Vector&
  Vector::operator=(Vector&& toCopy) noexcept
  {
    using std::swap;
    std::swap(data_, toCopy.data_);
    std::swap(n_, toCopy.n_);
    return *this;
  }

  // Vector&
  // Vector::operator+=(const Vector& toCopy) noexcept
  // {
  //   assert(size() == toCopy.size());

  //   for (Index_t i = 0; i < n_; ++i)
  //   {
  //     (*this)[i] += toCopy[i];
  //   }
  //   return *this;
  // }
  // Vector&
  // Vector::operator-=(const Vector& toCopy) noexcept
  // {
  //   assert(size() == toCopy.size());

  //   for (Index_t i = 0; i < n_; ++i)
  //   {
  //     (*this)[i] -= toCopy[i];
  //   }
  //   return *this;
  // }

//________________________________________________________________

#define IMPL(OP)                                             \
  Vector& Vector::operator OP(const Vector& toCopy) noexcept \
  {                                                          \
    assert(size() == toCopy.size());                         \
                                                             \
    for (Index_t i = 0; i < n_; ++i)                         \
    {                                                        \
      (*this)[i] OP toCopy[i];                               \
    }                                                        \
    return *this;                                            \
  }                                                          \
                                                             \
  Vector& Vector::operator OP(const double scalar) noexcept  \
  {                                                          \
    for (Index_t i = 0; i < n_; ++i)                         \
    {                                                        \
      (*this)[i] OP scalar;                                  \
    }                                                        \
    return *this;                                            \
  }

  IMPL(+=);
  IMPL(-=);
  IMPL(*=);
  IMPL(/=);

  // Vector&
  // Vector::operator+=(const double div) noexcept
  // {
  //   assert(div != 0);

  //   for (Index_t i = 0; i < n_; ++i)
  //   {
  //     (*this)[i] += div;
  //   }
  //   return *this;
  // }
  // Vector&
  // Vector::operator/=(const double div) noexcept
  // {
  //   assert(div != 0);

  //   for (Index_t i = 0; i < n_; ++i)
  //   {
  //     (*this)[i] /= div;
  //   }
  //   return *this;
  // }

  //////////////////////////////////////////////////////////////////

  double
  min(const Vector& v) noexcept
  {
    const Size_t n         = v.size();
    double       to_return = std::numeric_limits<double>::max();
    for (Index_t i = 0; i < n; ++i)
    {
      if (to_return > v[i])
        to_return = v[i];
    }
    return to_return;
  }
  double
  max(const Vector& v) noexcept
  {
    const Size_t n         = v.size();
    double       to_return = std::numeric_limits<double>::lowest();

    for (Index_t i = 0; i < n; ++i)
    {
      if (to_return < v[i])
        to_return = v[i];
    }
    return to_return;
  }

  double
  norm2(const Vector& v) noexcept
  {
    const double norm2 = cblas_dnrm2(v.size(), v.data(), v.inc());
    return norm2;
  }

  double
  dot(const Vector& v, const Vector& w) noexcept
  {
    assert(v.size() == w.size());

    return cblas_ddot(v.size(), v.data(), v.inc(), w.data(), w.inc());
  }

  void
  scal(const double alpha, Vector& v) noexcept
  {
    cblas_dscal(v.size(), alpha, v.data(), v.inc());
  }

  void
  axpy(const double alpha, const Vector& x, Vector& y) noexcept
  {
    cblas_daxpy(x.size(), alpha, x.data(), x.inc(), y.data(), y.inc());
  }
}
