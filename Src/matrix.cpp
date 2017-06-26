#include "matrix.hpp"

#include <cblas.h>
#include <algorithm>
#include <iomanip>

namespace JointDeconv
{
  namespace
  {
    template <MatrixTypeEnum TYPE, typename LAMBDA>
    void
    indices_pattern(LAMBDA lambda, const Size_t n, const Size_t m)
    {
      if (TYPE == MatrixTypeEnum::Regular)
      {
        for (Index_t j = 0; j < m; j++)
        {
          for (Index_t i = 0; i < n; i++)
          {
            lambda(i, j);
          }
        }
      }
      else
      {
        assert(TYPE == MatrixTypeEnum::Symmetric);
        // CAVEAT UpLow = Low
        for (Index_t j = 0; j < m; j++)
        {
          for (Index_t i = j; i < n; i++)
          {
            lambda(i, j);
          }
        }
      }
    }
  }  // anonymous

  template <MatrixTypeEnum TYPE>
  Matrix<TYPE>::Matrix() noexcept : data_(), n_(0), m_(0)
  {
  }

  template <MatrixTypeEnum TYPE>
  Matrix<TYPE>::Matrix(const Matrix<TYPE>& toCopy) : Matrix(toCopy.I_size(), toCopy.J_size())
  {
    *this = toCopy;
  }

  template <MatrixTypeEnum TYPE>
  Matrix<TYPE>::Matrix(const Size_t n, const Size_t m) : data_(new double[n * m]), n_(n), m_(m)
  {
    assert((Type::value == MatrixTypeEnum::Regular) || ((Type::value == MatrixTypeEnum::Symmetric) && (n == m)));
  }

  template <MatrixTypeEnum TYPE>
  Matrix<TYPE>::Matrix(std::initializer_list<std::initializer_list<double>> data)
  {
    // Any data?
    //
    n_ = data.size();

    if (n_ == 0)
    {
      m_ = 0;
      return;
    }

    // Defines column size
    //
    const Size_t m = data.begin()->size();  // lambda capture does not like [&this]
    // Check that all col have the same size
    assert(std::all_of(data.begin(), data.end(), [m](const auto& row) { return m == (Size_t)row.size(); }));
    m_ = m;  // ok same J_size for all rows

    data_     = std::move(decltype(data_)(new double[n_ * m_]));
    Index_t i = 0, j;
    for (const auto& row : data)
    {
      j = 0;
      for (const auto col : row)
      {
        (*this)(i, j) = col;
        ++j;
      }
      ++i;
    };
  }

  template <MatrixTypeEnum TYPE>
  Matrix<TYPE>::Matrix(Matrix<TYPE>&& toMove) noexcept : data_(std::move(toMove.data_)), n_(toMove.n_), m_(toMove.m_)
  {
    toMove.n_ = 0;
    toMove.m_ = 0;
  }

  template <MatrixTypeEnum TYPE>
  Matrix<TYPE>&
  Matrix<TYPE>::operator=(const double toCopy) noexcept
  {
    indices_pattern<Type::value>([&](const Index_t i, const Index_t j) { (*this)(i, j) = toCopy; }, I_size(), J_size());
    return *this;
  }
  template <MatrixTypeEnum TYPE>
  Matrix<TYPE>&
  Matrix<TYPE>::operator=(const Matrix<TYPE>& toCopy) noexcept
  {
    assert(I_size() == toCopy.I_size());
    assert(J_size() == toCopy.J_size());

    indices_pattern<Type::value>(
        [&](const Index_t i, const Index_t j) { (*this)(i, j) = toCopy(i, j); }, I_size(), J_size());
    return *this;
  }
  template <MatrixTypeEnum TYPE>
  Matrix<TYPE>&
  Matrix<TYPE>::operator=(Matrix<TYPE>&& toCopy) noexcept
  {
    using std::swap;
    std::swap(data_, toCopy.data_);
    std::swap(n_, toCopy.n_);
    std::swap(m_, toCopy.m_);
    return *this;
  }

  template <MatrixTypeEnum TYPE>
  std::ostream&
  operator<<(std::ostream& out, const Matrix<TYPE>& toPrint)
  {
    const auto old_settings  = out.flags();
    const auto old_precision = out.precision();
    out << std::setprecision(8);
    if (Matrix<TYPE>::Type::value == MatrixTypeEnum::Regular)
    {
      for (Index_t i = 0; i < toPrint.I_size(); i++)
      {
        out << std::endl;
        for (Index_t j = 0; j < toPrint.J_size(); j++)
        {
          out << std::setw(16) << toPrint(i, j) << " ";
        }
      }
    }
    else
    {
      assert(Matrix<TYPE>::Type::value == MatrixTypeEnum::Symmetric);
      for (Index_t i = 0; i < toPrint.I_size(); i++)
      {
        out << std::endl;
        for (Index_t j = 0; j <= i; j++)
        {
          out << std::setw(12) << toPrint(i, j) << " ";
        }
        for (Index_t j = i + 1; j < toPrint.J_size(); j++)
        {
          out << std::setw(12) << "*"
              << " ";
        }
      }
    }
    out.flags(old_settings);
    out.precision(old_precision);

    return out;
  }

#define IMPL(TYPE)             \
  template class Matrix<TYPE>; \
  template std::ostream& operator<<<TYPE>(std::ostream& out, const Matrix<TYPE>& toPrint);

  IMPL(MatrixTypeEnum::Regular);
  IMPL(MatrixTypeEnum::Symmetric);

  //////////////////////////////////////////////////////////////////

  namespace
  {
    template <MatrixTypeEnum TYPE>
    void
    internal_scal(const double alpha, Matrix<TYPE>& S) noexcept
    {
      indices_pattern<TYPE>([&](const Index_t i, const Index_t j) { S(i, j) *= alpha; }, S.I_size(), S.J_size());
    }
  }  // anon

  void
  scal(const double alpha, SymmetricMatrix& S) noexcept
  {
    internal_scal(alpha, S);
  }
  void
  scal(const double alpha, RegularMatrix& S) noexcept
  {
    internal_scal(alpha, S);
  }

  namespace
  {
    void
    check_dimensions_MM_fix_dimensions(const std::true_type, Size_t& I_size, Size_t& J_size) noexcept
    {
      std::swap(I_size, J_size);
    }
    void
    check_dimensions_MM_fix_dimensions(const std::false_type, Size_t& I_size, Size_t& J_size) noexcept
    {
      (void)I_size;
      (void)J_size;
    }

    // Op(M).x+y pattern
    //
    template <MatrixOpEnum OP>
    bool
    check_dimensions_Mv(MatrixOp_t<OP>, Size_t M_I_size, Size_t M_J_size, Size_t x_size, Size_t y_size) noexcept
    {
      check_dimensions_MM_fix_dimensions(typename MatrixOp_TransposedDimensions_t<OP>::type(), M_I_size, M_J_size);

      bool ok = true;

      ok &= M_I_size == y_size;
      ok &= M_J_size == x_size;

      return ok;
    }

    //________________

    constexpr auto
    toBlasMatrixOp(const Transpose_t) noexcept
    {
      return CblasTrans;
    };

    constexpr auto
    toBlasMatrixOp(const Identity_t) noexcept
    {
      return CblasNoTrans;
    };
  }  // anon

  void
  Mv(const double alpha, const SymmetricMatrix& S, const Vector& x, const double beta, Vector& y) noexcept
  {
    assert(check_dimensions_Mv(Identity_c, S.I_size(), S.J_size(), x.size(), y.size()));

    // Call blas-dsymv
    //
    // CAVEAT: assume UpLow = Low
    cblas_dsymv(
        CblasColMajor, CblasLower, S.I_size(), alpha, S.data(), S.ld(), x.data(), x.inc(), beta, y.data(), y.inc());
  }

  template <MatrixOpEnum OP>
  void
  Mv(const double         alpha,
     const MatrixOp_t<OP> op,
     const RegularMatrix& M,
     const Vector&        x,
     const double         beta,
     Vector&              y) noexcept
  {
    assert(check_dimensions_Mv(op, M.I_size(), M.J_size(), x.size(), y.size()));

    cblas_dgemv(CblasColMajor,
                toBlasMatrixOp(op),
                M.I_size(),
                M.J_size(),
                alpha,
                M.data(),
                M.ld(),
                x.data(),
                x.inc(),
                beta,
                y.data(),
                y.inc());
  }

#define IMP(OP)                                    \
  template void Mv<OP>(const double         alpha, \
                       const MatrixOp_t<OP> op,    \
                       const RegularMatrix& M,     \
                       const Vector&        x,     \
                       const double         beta,  \
                       Vector&              y) noexcept;

  IMP(MatrixOpEnum::Identity);
  IMP(MatrixOpEnum::Transpose);

#undef IMP

  namespace
  {

    // Op(M1).Op(M2)+M3 pattern
    //
    template <MatrixOpEnum OP_1, MatrixOpEnum OP_2>
    bool
    check_dimensions_MM(MatrixOp_t<OP_1>,
                        Size_t M_1_I_size,
                        Size_t M_1_J_size,
                        MatrixOp_t<OP_2>,
                        Size_t M_2_I_size,
                        Size_t M_2_J_size,
                        Size_t M_3_I_size,
                        Size_t M_3_J_size) noexcept
    {
      check_dimensions_MM_fix_dimensions(
          typename MatrixOp_TransposedDimensions_t<OP_1>::type(), M_1_I_size, M_1_J_size);
      check_dimensions_MM_fix_dimensions(
          typename MatrixOp_TransposedDimensions_t<OP_2>::type(), M_2_I_size, M_2_J_size);

      bool ok = true;

      ok &= M_1_I_size == M_3_I_size;
      ok &= M_2_J_size == M_3_J_size;
      ok &= M_1_J_size == M_2_I_size;

      return ok;
    }

  }  // Anon

  template <MatrixOpEnum OP>
  void
  MM(const double           alpha,
     const SymmetricMatrix& M1,
     const MatrixOp_t<OP>   op,
     const RegularMatrix&   M2,
     const double           beta,
     RegularMatrix&         M3) noexcept
  {
    // Sanity check
    //
    assert(check_dimensions_MM(
        Identity_c, M1.I_size(), M1.J_size(), op, M2.I_size(), M2.J_size(), M3.I_size(), M3.J_size()));

    cblas_dsymm(CblasColMajor,
                CblasLeft,
                CblasLower,
                M3.I_size(),
                M3.J_size(),
                alpha,
                M1.data(),
                M1.ld(),
                M2.data(),
                M2.ld(),
                beta,
                M3.data(),
                M3.ld());
  }

  template <MatrixOpEnum OP>
  void
  MM(const double           alpha,
     const MatrixOp_t<OP>   op,
     const RegularMatrix&   M1,
     const SymmetricMatrix& M2,
     const double           beta,
     RegularMatrix&         M3) noexcept
  {
    // Sanity check
    //
    assert(check_dimensions_MM(
        Identity_c, M1.I_size(), M1.J_size(), op, M2.I_size(), M2.J_size(), M3.I_size(), M3.J_size()));

    cblas_dsymm(CblasColMajor,
                CblasRight,
                CblasLower,
                M3.I_size(),
                M3.J_size(),
                alpha,
                M2.data(),
                M2.ld(),
                M1.data(),
                M1.ld(),
                beta,
                M3.data(),
                M3.ld());
  }

#define IMP(OP)                                             \
  template void MM<OP>(const double           alpha,        \
                       const SymmetricMatrix& M1,           \
                       const MatrixOp_t<OP>   op,           \
                       const RegularMatrix&   M2,           \
                       const double           beta,         \
                       RegularMatrix&         M3) noexcept; \
                                                            \
  template void MM<OP>(const double           alpha,        \
                       const MatrixOp_t<OP>   op,           \
                       const RegularMatrix&   M1,           \
                       const SymmetricMatrix& M2,           \
                       const double           beta,         \
                       RegularMatrix&         M3) noexcept;

  IMP(MatrixOpEnum::Identity);
  IMP(MatrixOpEnum::Transpose);

#undef IMP

  template <MatrixOpEnum OP_1, MatrixOpEnum OP_2>
  void
  MM(const double           alpha,
     const MatrixOp_t<OP_1> op1,
     const RegularMatrix&   M1,
     const MatrixOp_t<OP_2> op2,
     const RegularMatrix&   M2,
     const double           beta,
     RegularMatrix&         M3) noexcept
  {
    // Sanity check
    //
    assert(check_dimensions_MM(op1, M1.I_size(), M1.J_size(), op2, M2.I_size(), M2.J_size(), M3.I_size(), M3.J_size()));

    cblas_dgemm(CblasColMajor,
                toBlasMatrixOp(op1),
                toBlasMatrixOp(op2),
                M3.I_size(),
                M3.J_size(),
                (MatrixOp_TransposedDimensions_t<OP_1>::value ? M1.I_size() : M1.J_size()),
                alpha,
                M2.data(),
                M2.ld(),
                M1.data(),
                M1.ld(),
                beta,
                M3.data(),
                M3.ld());
  }

#define IMP(OP_1, OP_2)                                     \
  template void MM<OP_1, OP_2>(const double alpha,          \
                               const MatrixOp_t<OP_1> op1,  \
                               const RegularMatrix&   M1,   \
                               const MatrixOp_t<OP_2> op2,  \
                               const RegularMatrix&   M2,   \
                               const double           beta, \
                               RegularMatrix&         M3) noexcept;

  IMP(MatrixOpEnum::Identity, MatrixOpEnum::Identity);
  IMP(MatrixOpEnum::Identity, MatrixOpEnum::Transpose);
  IMP(MatrixOpEnum::Transpose, MatrixOpEnum::Identity);
  IMP(MatrixOpEnum::Transpose, MatrixOpEnum::Transpose);
}
