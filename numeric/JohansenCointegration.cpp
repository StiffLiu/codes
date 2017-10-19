/*
Center(M)  : Center the columns of matrix M at their means
Diff(M): M[1:r,] - M[0:(r-1),]
Lagged(M,nlags): a matrix C, where C[,k*(nlags + 1) + i] = M[(nlags-i):(r-i),k] where 0 <= i < nlags + 1
M/val: element-wise division of M by val
M': transpose of M
A*B: mutlplication of A and B
M`: Inverse matrix of M
Divide(X,Y): (X'*X)`*X'*Y


Let M = Center(M)
Let DM = Diff(M)
C[r-nlags-1, c*(nlags + 1)]:   Center(Lagged(DM,nlags))
A[r-nlags-1, c]:               Center(DM[nlags+1:,])
B[r-nlags-1, c]:               Center(M[1:r-nlags,])

Performance of above is low due to matrix multiplication order
  hat(C)[r-nlags, r-nlags]:                      C*(C'*C)`*C'             #This matrix is large, this will reduce the performance greatly, and take lots of memory
  ra[r-nlags-1, c]:                                A - hat(C)*A
  rb[r-nlags-1, c]:                                B - hat(C)*B

Below is more efficient way:
  inv(C)[c*(nlags+1), r-nlags]:                  (C'*C)` * C'

  ra[r-nlags-1, c]:                                A - C*inv(C)*A
  rb[r-nlags-1, c]:                                B - C*inv(C)*B

Skk[c, c]:                                       rb'*rb / (r - nlags - 1)
Sk0[c, c]:                                       rb'*ra / (r - nlags - 1)
S00[c, c]:                                       ra'*ra / (r - nlags - 1)

eigenInputMat: Skk`*(Sk0*S00`)*Sk0'
*/

// #define USE_GSL
// The above macro definition can improve performance a lot.
#define BOOST_UBLAS_NDEBUG
#include "JohansenCointegration.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <typeinfo>
#include <iostream>
#include <array>
#ifdef USE_GSL
#include <gsl/gsl_eigen.h>
#else
#include <Eigen/Eigenvalues>
#endif

namespace era{
  namespace numeric{
    using MatrixD = boost::numeric::ublas::matrix<double>;
    using VectorD = boost::numeric::ublas::vector<double>;

    MatrixD dblArray2Matrix(double *val_, size_t row_, size_t col_){
      MatrixD data(row_, col_);
      for(decltype(row_) i = 0;i < row_;++ i)
        for(decltype(col_) j = 0;j < col_;++ j)
          data(i, j) = val_[i * col_ + j];
      return data;
    }

    class JohansenCointegration::Impl{
    public:
      static double MaxEigenCriticalValues[][3];
      static const int NumEigenCriticalValues;
      friend class JohansenCointegration;

      template<class Matrix>
      Impl(Matrix m_, unsigned int nlags_){
        auto col = m_.size2();
        if(col > NumEigenCriticalValues) throw std::invalid_argument("number of variables too large");
        using namespace boost::numeric::ublas;
        center(m_);

        MatrixD eigenInput;
        size_t residualRows;
        {
          MatrixD s_kk, s_k0, s_00;
          {
            MatrixD ra, rb;
            {
              MatrixD dm = diff(m_);
              size_t row = dm.size1();
              auto lags = nlags_ + 1;
              auto c = lagged(dm, nlags_);
              center(c);
              residualRows = row - lags;

              auto invC = inv(c);
              ra = resid(center(project(dm, {lags, row}, {0, col})), invC, c);
              rb = resid(center(project(m_, {1, row - nlags_}, {0, col})), invC, c);
            }
            s_00 = prod(trans(ra), ra);
            s_kk = prod(trans(rb), rb);
            s_k0 = prod(trans(rb), ra);

            auto factor = 1.0 / residualRows;
            s_00 *= factor;
            s_kk *= factor;
            s_k0 *= factor;
          }
          MatrixD m1 = prod(inverse(s_kk), s_k0);
          MatrixD m2 = prod(inverse(s_00), trans(s_k0));
          eigenInput = prod(m1, m2);
        }
        VectorD eigValues;
        auto eigVecSorted = eigenVectorsSorted(eigenInput, eigValues);
        _tests = statistics(eigVecSorted, eigValues, residualRows, col);
        for(size_t i = 0;i < eigValues.size();++ i) _eig.push_back(eigValues(i));
        for(size_t i = 0;i < eigVecSorted.size2();++ i)
          for(size_t j = 0;j < eigVecSorted.size1();++ j)
            _eig.push_back(eigVecSorted(j, i));
        _numVars = col;
      }

      const auto& tests(){ return _tests;}

      size_t count() const {
        for (size_t i = 0; i < _tests.size(); i++){
          const auto& t = _tests[i];
          if (t._s < t._90 && t._s < t._95 && t._s < t._99) {
            return i;
          }
        }
        return _tests.size();
      }
    private:
      template<class Matrix>
      static auto diff(const Matrix& m_){
        using namespace boost::numeric::ublas;
        auto row = m_.size1();
        auto col = m_.size2();
        auto m1 = project(m_, {0, row - 1}, {0, col});
        auto m2 = project(m_, {1, row}, {0, col});
        return MatrixD(m2 - m1);
      }
      template<class Matrix>
      static auto inverse(const Matrix& m_){
        using namespace boost::numeric::ublas;
        MatrixD m(m_);
        permutation_matrix<size_t> pm(m.size1());
        int res = lu_factorize(m, pm);
        if (res != 0) throw std::invalid_argument("matrix is singular");
        MatrixD inv = identity_matrix<double>(m.size1());
        lu_substitute(m, pm, inv);
        return inv;
      }
      template<class Matrix>
      static auto lagged(const Matrix& m_, unsigned int nlags_){
        using namespace boost::numeric::ublas;
        auto row = m_.size1();
        auto col = m_.size2();
        auto lags = nlags_ + 1;
        auto laggedCols = col * lags;
        auto newRow = row - lags;
        MatrixD c(newRow, laggedCols);
        for(decltype(col) i = 0;i < col;++ i){
          auto base = i * lags;
          for(decltype(lags) j = 0;j < lags;++ j){
            // Lagged(M,nlags): a matrix C, 
            //   where C[,i*(nlags_ + 1) + j] = M[(nlags_-j):(r-j),k] where 0 <= j < nlags_ + 1
            column(c, base + j) = subslice(column(m_, i), nlags_ - j, 1, newRow);
          }

        }
        return c;
      }
      template<class Matrix>
      static auto inv(const Matrix& c_){
        using namespace boost::numeric::ublas;
        return MatrixD(prod(inverse(prod(trans(c_), c_)), trans(c_)));
      }
      template<class Matrix>
      static auto hat(const Matrix& m_, unsigned int nlags_){
        using namespace boost::numeric::ublas;
        auto c = lagged(m_, nlags_);
        return MatrixD(prod(c, inv(c)));
      }

      template<class M1, class M2, class M3>
      static auto resid(const M1& v_, const M2& inv_, const M3& c_){
        using namespace boost::numeric::ublas;
        return MatrixD(v_ - prod(c_, MatrixD(prod(inv_, v_))));
      }

      static auto center(MatrixD& m_){
        auto col = m_.size2();
        auto row = m_.size1();
        for(decltype(col) i = 0;i < col;++ i){
          double avg = 0;
          for(const auto& v : column(m_, i)) avg += v;
          avg /= row;
          for(auto& v : column(m_, i)) v -= avg;
        }
      }
      template<class Matrix>
      static auto center(const Matrix& m_){
        using namespace boost::numeric::ublas;
        MatrixD m(m_);
        center(m);
        return m;
      }
#ifdef USE_GSL
      template<class M, class V>
      static auto eigenVectorsSorted(const M& m_, V& eigens_){
        auto col = m_.size1();

        gsl_matrix* eigenInputMat = gsl_matrix_alloc(m_.size1(), m_.size2());
        gsl_vector_complex* evalPtr = gsl_vector_complex_alloc(col);
        gsl_matrix_complex* ematPtr = gsl_matrix_complex_alloc(col, col);
        gsl_eigen_nonsymmv_workspace* worspacePtr = gsl_eigen_nonsymmv_alloc(col);

        for(decltype(col) i = 0;i < m_.size1();++ i){
          for(size_t j = 0;j < m_.size2();++ j){
            gsl_matrix_set(eigenInputMat, i, j, m_(i, j));
          }
        }

        gsl_eigen_nonsymmv(eigenInputMat, evalPtr, ematPtr, worspacePtr);

        gsl_matrix_free(eigenInputMat);
        gsl_eigen_nonsymmv_free(worspacePtr);

        gsl_eigen_nonsymmv_sort(evalPtr, ematPtr, GSL_EIGEN_SORT_ABS_DESC);

        MatrixD vectors(ematPtr->size1, ematPtr->size2);
        for(int i = 0; i < ematPtr->size1; i++){
          for(int j = 0; j < ematPtr->size2; j++){
            auto comp = gsl_matrix_complex_get(ematPtr, i, j);
            auto real = comp.dat[0];
            auto imag = comp.dat[1];
            vectors(i, j) = real;
          }
        }

        VectorD values(evalPtr->size);
        for (int i = 0; i < evalPtr->size; i++){
          gsl_complex comp = gsl_vector_complex_get(evalPtr, i);
          auto real = comp.dat[0];
          auto imag = comp.dat[1];
          values[i] = real;
        }
        eigens_ = std::move(values);

        gsl_vector_complex_free(evalPtr);
        gsl_matrix_complex_free(ematPtr);

        return vectors;
      }
#else
      template<class M, class V>
      static auto eigenVectorsSorted(const M& m_, V& eigens_){
        using namespace Eigen;
        MatrixXd eigenInput(m_.size1(), m_.size2());
        for(size_t i = 0;i < m_.size1();++ i){
          for(size_t j = 0;j < m_.size2();++ j){
            eigenInput(i, j) = m_(i, j);
          }
        }

        EigenSolver<MatrixXd>  eigenSolver(eigenInput);
        const auto& eigValues = eigenSolver.eigenvalues();
        const auto& eigVectors = eigenSolver.eigenvectors();
        std::vector<size_t> sortedEigIndices;
        for(int i = 0;i < eigValues.size();++ i) sortedEigIndices.push_back(i);
        std::sort(sortedEigIndices.begin(), sortedEigIndices.end(),
          [&eigValues](auto i_, auto j_){return fabs(eigValues[i_].real()) > fabs(eigValues[j_].real());});

        MatrixD vectors(eigVectors.rows(), eigVectors.cols());
        VectorD values(eigValues.size());
        for(size_t j = 0;j < sortedEigIndices.size();++ j){
          auto index = sortedEigIndices[j];
          const auto& real = eigValues[index].real();
          int factor = (real < 0 ? -1 : 1);

          values[j] = fabs(real);
          for(int i = 0;i < eigVectors.rows();++ i)
            vectors(i, j) = factor * eigVectors(i, index).real();
        }
        eigens_ = std::move(values);
        return vectors;
      }
#endif

      template<class M, class V>
      static auto statistics(const M& sortedEigVec_, V& eigValues_, double numSamples_, size_t numVars_){
        Tests tests(eigValues_.size());
        for (size_t i = 0; i < eigValues_.size(); i++) {
          auto& t = tests[i];
          auto index = numVars_ - i - 1;
          // LR Max Eigen Value
          t._s = -numSamples_ * log(1 - eigValues_(i));
          t._90 = MaxEigenCriticalValues[index][0];
          t._95 = MaxEigenCriticalValues[index][1];
          t._99 = MaxEigenCriticalValues[index][2];
        }

        return tests;
      }

      using Tests = std::vector<Test>;
      std::vector<double> _eig;
      size_t _numVars;
      Tests _tests;
    };
    double JohansenCointegration::Impl::MaxEigenCriticalValues[][3] = {
      {2.7055 ,    3.8415,     6.6349},
      {12.2971,   14.2639,      18.52},
      {18.8928,   21.1314,     25.865},
      {25.1236,   27.5858,    32.7172},
      {31.2379,   33.8777,    39.3693},
      {37.2786,   40.0763,    45.8662},
      {43.2947,   46.2299,    52.3069},
      {49.2855,   52.3622,    58.6634},
      {55.2412,   58.4332,     64.996},
      {61.2041,    64.504,    71.2525},
      {67.1307,   70.5392,    77.4877},
      {73.0563,   76.5734,    83.7105}};
    const int JohansenCointegration::Impl::NumEigenCriticalValues =
      sizeof(MaxEigenCriticalValues) / sizeof(*MaxEigenCriticalValues);

    JohansenCointegration::JohansenCointegration(double* val_, unsigned int row_,
      unsigned int col_, unsigned int lags_) : _impl(new Impl(dblArray2Matrix(val_, row_, col_), lags_)){
    }

    JohansenCointegration::~JohansenCointegration(){
      delete _impl;
    }

    unsigned int JohansenCointegration::count() const{
      return _impl->count();
    }

    const double* JohansenCointegration::vector(unsigned int index_) const{
      auto index = (index_ + 1) * _impl->_numVars;
      return _impl->_eig.empty() || index >= _impl->_eig.size() ? nullptr : &_impl->_eig[index];
    }

    const double* JohansenCointegration::values() const{
      return _impl->_eig.empty() ? nullptr : &_impl->_eig[0];
    }

    unsigned int JohansenCointegration::numVars() const{
      return _impl->_numVars;
    }

    const JohansenCointegration::Test* JohansenCointegration::tests() const{
      const auto& t = _impl->tests();
      return t.empty() ? nullptr : &t[0];
    }

    unsigned int JohansenCointegration::numTests() const{
      return _impl->tests().size();
    }
  }
}

