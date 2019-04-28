#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <iostream>

static inline auto norm_cdf(double v_){
  static const double SQRT2 = std::sqrt(2);
  return 0.5 * (1 + std::erf(v_ / SQRT2));
}
/*
  Given N random variables, each with the same distribution and are independent:
    P(x=i) = 1/M, where i = 1, 2, ..., M
  Let Sn denote the sum of these N random variables, then Sn is also a random variable
  This class calculates the probability density function for Sn.   
 */
template<unsigned int N, unsigned int M> class MultiNomialDistribution{
public:
  static_assert(N > 0 && M > 0, "N and M should all be positive.");
  static constexpr auto total(unsigned int n_) {
    return M * n_ * (n_ + 1) / 2 - n_ * (n_ - 1) / 2;
  }
  static constexpr auto lbound(unsigned int n_) { return n_; }
  static constexpr auto ubound(unsigned int n_) { return n_ * M; }
  static constexpr auto idx_no_check(unsigned int n_, unsigned int v_) {
    return total(n_ - 1) + v_ - n_;
  }
  static constexpr auto idx(unsigned int n_, unsigned int v_) {
    return v_ < lbound(n_) || v_ > ubound(n_) ? LEN : idx_no_check(n_, v_);
  }
  static constexpr unsigned long LEN = total(N);
  double _probs[LEN + 1];
  constexpr auto p(unsigned int n_, unsigned int v_) {
    return n_ <= N ? _probs[idx(n_, v_)] : throw std::invalid_argument("n_ greater than N");
  }
  constexpr auto p(unsigned int v_) { return p(N, v_); }
  // cumulative density function.
  // probability of the range: (-inf, v_)
  constexpr double cdf(unsigned int n_, unsigned int v_) {
    if (v_ == 0) return 0.0;
    -- v_;
    auto l = lbound(n_);
    auto u = (v_ > ubound(n_) ? ubound(n_) : v_);
    if (l > u) return 0.0;

    auto begin = idx(n_, l);
    auto r = _probs[begin];
    for(++ begin, ++ l; l <= u; ++ begin, ++l) r += _probs[begin];
    return r;
  }
  constexpr auto cdf(unsigned int v_) { return cdf(N, v_); }
  constexpr MultiNomialDistribution(){
    for(int i = 0;i < M;++ i) _probs[i] = 1.0 / M;
    for(int i = 2;i <= N;++ i){
      for(int v = lbound(i);v <= ubound(i);++ v){
        auto& p = _probs[idx_no_check(i, v)];
        p = 0;
        if (M <= v) for(int k = 1;k <= M;++ k) p += _probs[idx(i - 1, v - k)];
        else for(int k = 1;k <= v;++ k) p += _probs[idx(i - 1, v - k)];
        p /= M;
      }
    }
    _probs[LEN] = 0;
  }

  // By center limit theory, when N approaches infinity:
  //   (Sn - E(Sn)) / sqrt(Var(Sn)) 
  // should be normal distributed.
  // This function estimates the cdf using normal distribution.
  static auto ecdf(unsigned int v_){
    static const auto mean = (M + 1) / 2.0;
    static const auto var = (M + 1) * (2 * M + 1) / 6.0 - mean * mean;
    static const auto nmean = N * mean;
    static const auto sqrtnvar = std::sqrt(N * var);
    return norm_cdf((v_ - nmean) / sqrtnvar);
  }
  template<class Stream> Stream& print(Stream& os_){
    for(const auto& v : _probs) os_ << v << ' ';
    return os_;
  }
};

int test(int argc_, char* argv_[]){
  static const int N = 100, M = 6;
  MultiNomialDistribution<N, M> mnd;
  for (unsigned int n = N; n < N * M; ++ n){
    auto p1 = mnd.cdf(n);
    auto p2 = mnd.ecdf(n);
    std::cout << "p(sum < " << n << ")[real, estimate using normal dist, delta]: [" << p1 << ',' << p2 << ',' << (p2 - p1) << "]\n";
  }
  return 0;
}
int main(int argc_, char* argv_[]){
  return test(argc_, argv_);
}