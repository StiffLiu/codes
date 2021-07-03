#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <deque>


template<class V> class OSTree{
  using Vector = std::deque<V>;

  template<class Vec> class IndexCompareT{
    const Vec& _vec;
  public:
    IndexCompareT(const Vec& vec_) : _vec(vec_){}
    bool operator()(size_t i_, size_t j_) const { return _vec(i_, j_); }
  };
  using IndexCompare = IndexCompareT<OSTree>;

  // A red-black tree table storing values and their
  // order statistics. Note that since the tree uses
  // tree_order_statistics_node_update as its update
  // policy, then it includes its methods by_order
  // and order_of_key.
  using ostree = __gnu_pbds::tree<size_t, __gnu_pbds::null_type, IndexCompare,
    // This policy updates nodes' metadata for order statistics.
    __gnu_pbds::rb_tree_tag, __gnu_pbds::tree_order_statistics_node_update>;

  Vector _vec;
  ostree _tree{IndexCompare{*this}};
  V _default{};
  mutable const V* _key{};
  static const size_t npos = static_cast<size_t>(-1);
  size_t _mincnt{1}, _maxcnt{}, _erased{};
  void real_pop(){
    _tree.erase(_erased);
    _vec.pop_front();
    // _tree.erase will use _erase indiectly in operator[], 
    //   under this case the unincremented _erase is expected.
    ++ _erased;
  }
  bool operator()(size_t i_, size_t j_) const {

    if(i_ == npos) return j_ != npos && *_key < (*this)[j_];
    else if(j_ == npos) return !(*_key < (*this)[i_]);

    const auto& v1 = (*this)[i_];
    const auto& v2 = (*this)[j_];
    if(v1 < v2) return true;
    if(v2 < v1) return false;
    return i_ < j_;
  }
public:
  OSTree(V default_={}, size_t mincnt_=1, size_t maxcnt_={})
    : _default(std::move(default_)), _mincnt{mincnt_}, _maxcnt{maxcnt_}{}
  void push(const V& v_){
    auto cur = _vec.size();
    auto total = cur + _erased;
    _vec.push_back(v_);
    _tree.insert(total);
    if(_maxcnt != 0 && cur >= _maxcnt) real_pop();
  }
  bool pop() {
    return !_tree.empty() && (real_pop(), true);
  }
  auto order_of_key(const V& value_) const {
    _key = &value_;
    return _tree.order_of_key(npos);
  }
  auto percentile_of_key(const V& value_) const {
    return order_of_key(value_) / static_cast<double>(_tree.size());
  }
  auto find_by_order(size_t rank_) const {
    const auto& cnt = _tree.size();
    if(cnt <= _mincnt) return _default;
    if(rank_ >= cnt) return (*this)[*_tree.rbegin()];/*We're not empty here*/
    auto pos = _tree.find_by_order(rank_);
    if(pos == _tree.end()) return _default;
    return (*this)[*pos];
  }
  
  auto find_by_percentile(double p_) const {
    return find_by_order(std::round(p_ * _tree.size()));
  }
  auto& operator[](size_t i_) const {
    return _vec[i_ - _erased];
  }
  auto& mincnt() const { return _mincnt; }
  auto& maxcnt() const { return _maxcnt; }
  auto& vec() const { return _vec; }
  auto size() const { return _tree.size(); }
  void clear() {
    _tree.clear();
    _vec.clear();
    _erased = 0;
  }
};
