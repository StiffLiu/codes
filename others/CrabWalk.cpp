#include <map>
#include <set>
#include <list>
#include <cmath>
#include <ctime>
#include <deque>
#include <queue>
#include <stack>
#include <string>
#include <bitset>
#include <cstdio>
#include <limits>
#include <vector>
#include <climits>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <numeric>
#include <sstream>
#include <iostream>
#include <algorithm>
using namespace std;


/*
 * Complete the function below.
 * Let:
 *  1. N denote the total number of steps.
 *  2. n1 denote the number of steps the crab moves left
 *  3. n2 denote the number of steps the crab moves right
 *  4. n3 denote the number of steps the crab moves up
 *  5. n4 denote the number of steps the crab moves down
 *  6. C(N - 1, n1, n2, n3, n4) denote the binomial combination with: n1+n2+n3+n4=N-1
 * Since the crab stops at the N-th move, the total number of all possible moves is then
 *   C(N - 1, n1, n2, n3, n4). 
 * With each move, the crab ends at distance sqrt((n1-n2) * (n1-n2) + (n3-n4)*(n3-n4))
 * When N is large C(N - 1, n1, n2, n3, n4) is too large to compute except some approximation
 *  is used(for example Stirling's formula).
 * There're two approaches to this problem
 *  1. Analytical way, which needs mathematics(Combinatory theories and probability theories)
 *  2. Simulated way(Monto Carlo), we simulate crab moves using computer programs.
 */

float ExpectedCrabWalkSimulated(int numsteps){
  static const int TotalIterations{100000};
  // The brute force way, since pow(0.8, 100)=2.037e-10, which is the probability that 
  //  the crab moves more than 100 steps
  std::map<std::pair<int, int>, int> counts;
  for(int i = 0;i < TotalIterations;++ i){
      int h{}, v{};
      for(int j = 0;j < numsteps;++ j){
          int seed = rand()%5;
          if(seed == 0) break;
          switch(seed){
              case 1: ++h;break;
              case 2: --h;break;
              case 3: ++v;break;
              case 4: --v;break;
              default:break;
          }
      }
      ++ counts[std::pair<int, int>(h, v)];
  }
  
  float edist{};
  for(const auto& kvp : counts){
      const auto& p = kvp.first;
      const auto& frequency = kvp.second;
      const auto& dist = sqrt(p.first * p.first + p.second * p.second) * frequency / TotalIterations;
      edist += dist;
  }
  // There's another simulation way, 
  //  which is by choosing splitting positions from 0,1,2,...,numsteps.
  // This might work, however we omit it's implementaion here.
  return edist;
}
/*
 * Since the crab will stop at any step before numsteps
 * So we will have numsteps mutually exclusive situations,
 *  which is: the crab with stop exactly at the i-th move(for 1 <= i < numsteps)
 *        and the crab does not stop at (numsteps - 1)-th step
 * The number of all possible positions for numsteps move is as following:
 *  0. numsteps=0: 1
 *  1. numsteps=1: 1 + 4
 *  2. numsteps=2: 1 + 4 + 8
 *  3. numsteps=3: 1 + 4 + 8 + 12
 *  ......
 *  n. numsteps=n: 1 + 2 * n * (n + 1)
 * So we can setup an one-to-one mapping from the positions to the integers in the range [0, 2*n*(n+1)]
 */
class Mappings{
    int _lastnum{0};
    std::vector<std::pair<int, int> > _i2p;
    std::vector<double> _dist;
public:
    Mappings(){
        _i2p.push_back({0,0});
    }
    void ensure(int num_){
        while(_lastnum < num_){
            ++ _lastnum;
            int i = _lastnum;
            for(;i >= 0;-- i) _i2p.push_back({i, _lastnum - i});
            for(;i >= -_lastnum;--i) _i2p.push_back({i, _lastnum + i});
            i = -_lastnum + 1;
            for(;i <= 0;++ i) _i2p.push_back({i, -i - _lastnum});
            for(;i < _lastnum;++ i) _i2p.push_back({i, -_lastnum + i});
        }
        size_t start = _dist.size();
        size_t count = _i2p.size();
        for(;start < count;++ start){
            const auto& p = _i2p[start];
            _dist.push_back(sqrt(p.first * p.first + p.second * p.second));
        }
    }
    static int total(int n_){
        return 1 + 2 * n_ * (n_ + 1);
    }
    /*
     * Calculate expected distance given probability p_
     * Assuming n_ <= _dist.size()
     */
    float edist(std::vector<double>& p_, int n_) const{
        float sum{};
        for(size_t i = 0;i < n_;++ i) sum += p_[i] * _dist[i];
        return sum;
    }
    /*
     * Assume that n_ < _i2p.size()
     */
    // requires c++17: const auto& n2pos(int n_) const { return _i2p[n_]; }
    const std::pair<int,int>& n2pos(int n_) const { return _i2p[n_]; }
    int pos2n(int x_, int y_){
        int s_ = (x_ < 0 ? (-x_): x_) + (y_ < 0 ? (-y_) : y_);
        if(s_ == 0) return 0;
        int base = total(s_ - 1) - 1, offset{};
        if(x_ > 0 && y_ >= 0) offset = y_ + 1;
        if(x_ <= 0 && y_ > 0) offset = s_ - x_;
        if(x_ < 0 && y_ <= 0) offset = 2 * s_ - y_ + 1;
        if(x_ >= 0 && y_ < 0) offset = 3 * s_ + x_;
        return base + offset;
    }
};
static Mappings instance;
float ExpectedCrabWalkIterative(int numsteps){
    instance.ensure(numsteps + 1);

    // probabilities for each position.
    std::vector<double> current, next;
    int total = instance.total(numsteps + 1);
    double lastp = 0.2, accump{};
    float edist{};
    current.resize(total);
    next.resize(total);
    current[0] = 1;

    for(int i = 1;i <= numsteps;++ i){
        edist += instance.edist(current, instance.total(i - 1)) * lastp;
        // lastp *= 0.8;

        int ntotal = instance.total(i);
        for(int j = 0;j < ntotal;++ j){
            const auto& pos = instance.n2pos(j);
            // update the possbilities.
            next[j] = 0.2 * (current[instance.pos2n(pos.first - 1, pos.second)] +
              current[instance.pos2n(pos.first + 1, pos.second)] +
              current[instance.pos2n(pos.first, pos.second + 1)] +
              current[instance.pos2n(pos.first, pos.second - 1)]);
        }
        current.swap(next);
    }
    edist += instance.edist(current, total);
    return edist;
}

float ExpectedCrabWalk(int numsteps) {
    //return ExpectedCrabWalkSimulated(numsteps);
    return ExpectedCrabWalkIterative(numsteps);
}

int main() {
   std::cout << ExpectedCrabWalkIterative(2) << std::endl;
   return 0;
}
