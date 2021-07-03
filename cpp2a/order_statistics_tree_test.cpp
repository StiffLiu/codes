#include "order_statistics_tree.h"
#include <random>
#include <iostream>
#include <algorithm>
#include <cassert>

template<class Tree, class Values>
void check_ranks(const Tree& tree_, Values values_, size_t mincnt_){
  std::sort(values_.begin(), values_.end());
  // for(const auto& v : values_) std::cout << v << ' ';
  // std::cout << ':';
  const auto& cnt = values_.size();
  for(const auto& p : {0.0, 0.0333, 0.1, 0.5, 0.8, 0.9, 1.0}){
    auto rank = std::round(p * cnt);
    auto value = (cnt < mincnt_ ? 0 : (rank >= cnt ? *values_.rbegin() : values_[rank]));
    // std::cout << " (" << p << ", " << rank << ", " << tree_.find_by_order(rank) << ", " << value << ")";
    // std::cout.flush();
    auto value1 = tree_.find_by_order(rank);
    auto value2 = tree_.find_by_percentile(p);
    if (value1 != value){
      for(const auto& v : values_) std::cout << v << ' ';
      std::cout << " (" << p << ", " << rank << ", " << tree_.find_by_order(rank) << ", " << value << ")";
      std::cout << "\nsize: " << tree_.size() << ' ' << value1 << "!=" << value << std::endl;
    }
    assert(value1 == value);
    assert(value2 == value);
  }
  for(size_t i = 0, j = 1;j < cnt; ++ i, ++ j){
    if(values_[i] == values_[j]) continue;
    const auto oi = tree_.order_of_key(values_[i]);
    if(oi != j){
      for(const auto& v : values_) std::cout << v << ' ';
      std::cout << cnt << "==" << tree_.size() << " order mismatch: " << values_[i] << ' ' << oi << ' ' << j << std::endl;
    }
    assert(oi == j);
    assert(std::abs(tree_.percentile_of_key(values_[i]) - j / static_cast<double>(cnt)) < 1e-6);
  }
  const auto oi = tree_.order_of_key(values_[cnt - 1]);
  if(oi != cnt){
    for(const auto& v : values_) std::cout << v << ' ';
    std::cout << cnt << "==" << tree_.size() << " order mismatch: " << values_[cnt - 1] << ' ' << oi << ' ' << cnt << std::endl;
  }
  assert(oi == cnt);
  if(cnt > 0){
    assert(std::abs(tree_.percentile_of_key(values_[cnt - 1]) - 1.0) < 1e-6);
  }
  //std::cout << '\n';
}

template<bool Case> void test(){
  std::cout << (Case ? "Unlimited" : "Limited") << std::endl;
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_int_distribution<int> uni(1, 50);//(Case ? 50 : 250));
  auto tree = (Case ? OSTree<int>() : OSTree<int>(0, 20, 30000));
  std::vector<int> values;
  auto count = 30, iteration = (Case ? 20 : 4000000 * 5), pops=0;

  for(int i = 0;i < iteration;++ i){
    if constexpr (Case){
      tree.clear();
      for(int j = 0;j < count;++ j){
        auto v = uni(rng);
        tree.push(v);
        values.push_back(v);
      }
      while(!values.empty()){
        check_ranks(tree, values, tree.mincnt() + 1);
        values.erase(values.begin());
        tree.pop();
      }
      check_ranks(tree, std::move(values), tree.mincnt() + 1);
    } else {
      if (i % int(iteration * 0.001) == 1){
        std::cout << "\b\b." << int(i / double(iteration)*100);
        std::cout.flush();
      }
      const auto& v = uni(rng);
      values.push_back(v);
      tree.push(v);
      if(values.size() > tree.maxcnt()) values.erase(values.begin());
      check_ranks(tree, values, tree.mincnt() + 1);
      if(uni(rng) % 10 == 0){
        for(int j = 0;j < 3 && !values.empty();++ j){
          values.erase(values.begin());
          tree.pop();
          check_ranks(tree, values, tree.mincnt() + 1);
          ++ pops;
        }
      }
    }
  }

  if(!Case) std::cout << " pops: " << pops << std::endl;
}

void benchmark(){
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_real_distribution<double> uni(1.0, 50.0), percentile(0.0, 1.00001);
  auto tree = OSTree<double>(0, 100, 30000);
  auto iteration = 400000;
  auto start = clock();
  for(int i = 0;i < iteration;++ i){
    const auto& v = uni(rng);
    tree.push(v);
    tree.find_by_percentile(percentile(rng));
  }
  auto duration = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC) * 1e6;
  std::cout << "Time: " << duration << "us, average: " << duration / iteration << "us" << std::endl;

  start = clock();
  for(int i = 0;i < iteration;++ i){
    tree.find_by_percentile(percentile(rng));
  }
  duration = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC) * 1e6;
  std::cout << "Time: " << duration << "us, average: " << duration / iteration << "us" << std::endl;

  start = clock();
  for(int i = 0;i < iteration;++ i){
    const auto& v = uni(rng);
    tree.order_of_key(v);
    tree.percentile_of_key(v);
  }
  duration = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC) * 1e6;
  std::cout << "Time: " << duration << "us, average: " << duration / iteration << "us" << std::endl;
}

int main(int argc_, char* argv_[]){
  benchmark();
  test<true>();
  test<false>();
  return 0;
}
