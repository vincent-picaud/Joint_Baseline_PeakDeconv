// Defines a minimalist matrix class, no subtleties
//
// Note: ColMajor
// Note: UpLow = Low for symmetric ones
//
#pragma once

#include "vector.hpp"

#include <cassert>
#include <initializer_list>
#include <iostream>
#include <memory>

namespace JointDeconv
{
  enum class MatrixTypeEnum
  {
    Regular,
    Symmetric
  };

  template <MatrixTypeEnum TYPE>
  class Matrix;
  template <MatrixTypeEnum TYPE>
  std::ostream& operator<<(std::ostream& out, const Matrix<TYPE>& toPrint);

  template <MatrixTypeEnum TYPE>
  class Matrix
  {
    template <MatrixTypeEnum OTHER_TYPE>
    friend class Matrix;

   public:
    using Type = std::integral_constant<MatrixTypeEnum, TYPE>;

   public:
    Matrix() noexcept;
    Matrix(const Size_t n, const Size_t m);
    Matrix(std::initializer_list<std::initializer_list<double>> data);
    Matrix(const Matrix<TYPE>& toCopy);
    Matrix(Matrix<TYPE>&&) noexcept;
    template <MatrixTypeEnum OTHER_TYPE>
    Matrix(Matrix<OTHER_TYPE>&& toMove) noexcept : data_(std::move(toMove.data_)), n_(toMove.n_), m_(toMove.m_)
    {
      toMove.n_ = 0;
      toMove.m_ = 0;
    };

    Size_t
    I_size() const noexcept
    {
      return n_;
    };
    Size_t
    J_size() const noexcept
    {
      return m_;
    };
    Size_t
    ld() const noexcept
    {
      return n_;
    }
    Size_t
    mem_offset(const Index_t i, const Index_t j) const noexcept
    {
      return i + j * ld();
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
    check_index(const Index_t i, const Index_t j) const noexcept
    {
      return (i >= 0) && (j >= 0) && (i < n_) && (j < m_) &&
             ((TYPE == MatrixTypeEnum::Regular) || ((TYPE == MatrixTypeEnum::Symmetric) && (i >= j)));
    }
    Matrix& operator=(const double toCopy) noexcept;
    Matrix& operator=(const Matrix& toCopy) noexcept;
    Matrix& operator=(Matrix&& toCopy) noexcept;
    double&
    operator()(const Index_t i, const Index_t j) noexcept
    {
      assert(check_index(i, j));

      return data_[mem_offset(i, j)];
    };
    const double&
    operator()(const Index_t i, const Index_t j) const noexcept
    {
      assert(check_index(i, j));

      return data_[mem_offset(i, j)];
    };

   private:
    std::unique_ptr<double[]> data_;
    Size_t                    n_, m_;

   private:
    friend std::ostream& operator<<<TYPE>(std::ostream& out, const Matrix<TYPE>& toPrint);
  };

  using RegularMatrix   = Matrix<MatrixTypeEnum::Regular>;
  using SymmetricMatrix = Matrix<MatrixTypeEnum::Symmetric>;

  // Some MatrixOp Modifiers
  //
  enum class MatrixOpEnum
  {
    Transpose,
    Identity
  };

  template <MatrixOpEnum OP>
  using MatrixOp_t  = std::integral_constant<MatrixOpEnum, OP>;
  using Transpose_t = MatrixOp_t<MatrixOpEnum::Transpose>;
  using Identity_t  = MatrixOp_t<MatrixOpEnum::Identity>;

  template <MatrixOpEnum OP>
  constexpr auto         MatrixOp_c  = MatrixOp_t<OP>();
  constexpr auto         Transpose_c = Transpose_t();
  constexpr auto         Identity_c  = Identity_t();

  template <MatrixOpEnum OP>
  using MatrixOp_TransposedDimensions_t = std::integral_constant<bool, (OP == MatrixOpEnum::Transpose)>;

  // Required math functions
  //
  void scal(const double alpha, SymmetricMatrix& S) noexcept;
  void scal(const double alpha, RegularMatrix& M) noexcept;

  void Mv(const double alpha, const SymmetricMatrix& S, const Vector& x, const double beta, Vector& y) noexcept;

  template <MatrixOpEnum OP>
  void Mv(const double alpha,
          const MatrixOp_t<OP>,
          const RegularMatrix& M,
          const Vector&        x,
          const double         beta,
          Vector&              y) noexcept;

  template <MatrixOpEnum OP>
  void MM(const double           alpha,
          const SymmetricMatrix& M1,
          const MatrixOp_t<OP>,
          const RegularMatrix& M2,
          const double         beta,
          RegularMatrix&       M3) noexcept;

  template <MatrixOpEnum OP>
  void MM(const double alpha,
          const MatrixOp_t<OP>,
          const RegularMatrix&   M1,
          const SymmetricMatrix& M2,
          const double           beta,
          RegularMatrix&         M3) noexcept;

  template <MatrixOpEnum OP_1, MatrixOpEnum OP_2>
  void MM(const double alpha,
          const MatrixOp_t<OP_1>,
          const RegularMatrix& M1,
          const MatrixOp_t<OP_2>,
          const RegularMatrix& M2,
          const double         beta,
          RegularMatrix&       M3) noexcept;
}
