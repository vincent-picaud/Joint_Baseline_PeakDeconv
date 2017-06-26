// Defines a minimalist vector class
#pragma once

#include "index.hpp"

#include <cassert>
#include <iomanip>
#include <iostream>
#include <memory>

namespace JointDeconv
{
  template <typename LAMBDA>
  void
  indices_pattern(LAMBDA lambda, const Size_t n)
  {
    for (Index_t i = 0; i < n; i++)
    {
      lambda(i);
    }
  }

  class Vector
  {
   public:
    Vector() noexcept;
    explicit Vector(const Size_t n);
    Vector(Vector&& toMove) noexcept;

    Size_t
    size() const noexcept
    {
      return n_;
    };
    Size_t
    inc() const noexcept
    {
      return 1;
    }
    Size_t
    mem_offset(const Index_t i) const noexcept
    {
      return i * inc();
    }
    const double*
    data() const noexcept
    {
      return data_.get();
    };
    double*
    data() noexcept
    {
      return data_.get();
    };
    bool
    check_index(const Index_t i) const noexcept
    {
      return (i >= 0) && (i < n_);
    }
    Vector&
    operator=(const double toCopy) noexcept
    {
      indices_pattern([&](const Index_t i) { (*this)[i] = toCopy; }, size());
      return *this;
    }
    Vector& operator=(const Vector& toCopy) noexcept;
    Vector& operator=(Vector&& toCopy) noexcept;
    double& operator[](const Index_t i) noexcept
    {
      assert(check_index(i));
      return data_[mem_offset(i)];
    };
    const double& operator[](const Index_t i) const noexcept
    {
      assert(check_index(i));
      return data_[mem_offset(i)];
    };

    Vector& operator+=(const Vector& vector) noexcept;
    Vector& operator+=(const double scalar) noexcept;
    Vector& operator-=(const Vector& vector) noexcept;
    Vector& operator-=(const double scalar) noexcept;
    Vector& operator*=(const Vector& vector) noexcept;
    Vector& operator*=(const double scalar) noexcept;
    Vector& operator/=(const Vector& vector) noexcept;
    Vector& operator/=(const double scalar) noexcept;

   private:
    std::unique_ptr<double[]> data_;
    Size_t                    n_;

   private:
    friend std::ostream&
    operator<<(std::ostream& out, const Vector& toPrint)
    {
      const auto old_settings  = out.flags();
      const auto old_precision = out.precision();
      out << std::setprecision(8);
      out << std::endl;
      for (Index_t i = 0; i < toPrint.size(); i++)
      {
        out << std::setw(16) << toPrint[i] << " ";
      }
      out.flags(old_settings);
      out.precision(old_precision);

      return out;
    }
  };

  // Required math functions
  //
  double min(const Vector& v) noexcept;
  double max(const Vector& v) noexcept;
  double norm2(const Vector& v) noexcept;
  double dot(const Vector& v, const Vector& w) noexcept;
  void scal(const double alpha, Vector& v) noexcept;
  void axpy(const double alpha, const Vector& x, Vector& y) noexcept;
}
