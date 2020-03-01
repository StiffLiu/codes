// Example taken from: https://kirit.com/How%20C%2B%2B%20coroutines%20work
#include <coroutine>
#include <iostream>
#include <typeinfo>
#include <memory>
#include <atomic>
#include <thread>

struct suspend_never_all {
  using s_init = std::suspend_never;
  using s_return = std::suspend_never;
  using s_final = std::suspend_never;
  using s_yield = std::suspend_never;
};
struct init_suspend_always : public suspend_never_all {
  using s_init = std::suspend_always;
};
struct final_suspend_always : public suspend_never_all {
  using s_final = std::suspend_always;
};
struct init_final_suspend_always :  public suspend_never_all {
  using s_init = std::suspend_always;
  using s_final = std::suspend_always;
};
struct init_yield_final_suspend_always :  public suspend_never_all {
  using s_init = std::suspend_always;
  using s_final = std::suspend_always;
  using s_yield = std::suspend_always;
};
struct yield_final_suspend_always :  public suspend_never_all {
  using s_final = std::suspend_always;
  using s_yield = std::suspend_always;
};

template<typename T, typename sync_traits> struct sync{
  struct promise_type;
  static constexpr bool IsLazy =
    !std::is_same_v<typename sync_traits::s_init, 
      std::suspend_never>;
  static constexpr bool PromiseOutLives =
    std::is_same_v<typename sync_traits::s_final, 
      std::suspend_always>;
  using handle_type = std::conditional_t<!PromiseOutLives,
    std::unique_ptr<T>, std::coroutine_handle<promise_type> >;
  explicit sync(handle_type coro_) : _coro(std::move(coro_)) {
    std::cout << "Created a sync object [" << this << ']' << std::endl;
  }
  sync(const sync& s_) = delete;
  sync(sync&& s_): _coro(s_._coro) {
    std::cout << "Sync moved leaving behind a husk " << std::endl;
    s_._coro = nullptr;
  }
  sync& operator=(const sync& s_) = delete;
  sync& operator=(sync&& s_){
    if(this != &s_){
      std::cout << "Sync copy moved leaving behind a husk " << std::endl;
      _coro = s_._coro;
      s_._coro = nullptr;
    }
    return *this;
  }
  /*sync(const sync& s_){
    std::cout << "Copied a sync objet" << std::endl;
  }*/
  ~sync(){
    std::cout << "Sync gone " << this << std::endl;
    if constexpr(PromiseOutLives){
      if(_coro) _coro.destroy();
    }
  }
  T get(){
    std::cout << this << ": we got asked for the return value..." << std::endl;
    if constexpr(PromiseOutLives){
      if constexpr(IsLazy){
        if(!_coro.done()) _coro.resume();
      }
      return _coro.promise()._v;
    }else{
      if constexpr(IsLazy){
        std::cout << "We got unitialized value!!!" << std::endl;
      }
      return *_coro;
    }
  }
  struct promise_type;
  struct awaitable_type{
    std::coroutine_handle<promise_type> _coro;
    bool is_ready() const {
      if constexpr(PromiseOutLives){
        if constexpr(IsLazy){
          return _coro.done();
        }
        return true;
      }else{
        return !IsLazy;
      }
    }
    bool await_ready(){
      const auto ready = is_ready();
      std::cout << "Await " << (ready ? "is ready" : "not ready") << std::endl;
      return ready;
    }
    auto await_suspend(auto awaiting_){
      std::cout << "About to resume the lazy" << std::endl;
      _coro.resume();
      std::cout << "About to resume the awaiter" << std::endl;
      // awaiting_.resume();
      return awaiting_;
    }
    auto await_resume(){
      const auto& r = _coro.promise()._v;
      std::cout << "Await value is returned: " << r << std::endl;
      return r;
    }
  };
  auto operator co_await(){
    return awaitable_type{_coro};
  }
  struct promise_type{
    promise_type(){
      std::cout << "Promise created " << this << std::endl;
    }
    ~promise_type(){
      std::cout << "Promise died " << this << std::endl;
    }
    auto initial_suspend(){
      std::cout << "Started the coroutine, don't stop now!" << std::endl;
      return typename sync_traits::s_init{};
    }
    auto return_value(T v_){
      std::cout << "Got an answer of " << v_ << std::endl;
      if constexpr(PromiseOutLives){
        _v = std::move(v_);
      }else{
        *_v = std::move(v_);
      }
      return typename sync_traits::s_return{};
    }
    auto final_suspend(){
      std::cout << "Finished the coro" << std::endl;
      return typename sync_traits::s_final{};
    }
    auto get_return_object(){
      std::cout << "Send back a sync" << std::endl;
      if constexpr(PromiseOutLives){
        return sync{handle_type::from_promise(*this)};
      }else{
        auto v = std::make_unique<T>();
        _v = v.get();
        return sync{std::move(v)};
      }
    }
    void unhandled_exception(){
      std::exit(1);
    }
    std::conditional_t<!PromiseOutLives, T*, T> _v{};
  };
  handle_type _coro{};
};

template<typename traits>
sync<int, traits> answer(){
  std::cout << "Thinking deep thoughts..." << std::endl;
  co_return 42;
}

template<typename traits>
sync<int, traits> await_answer(){
  std::cout << "Started await_answer" << std::endl;
  auto a = answer<traits>();
  std::cout << "Got a coroutine, let's get a value" << std::endl;
  auto v = co_await (a);
  std::cout << "And the coroutine values is: " << v << std::endl;
  v = co_await (a);
  std::cout << "And the coroutine values is still: " << v << std::endl;
  co_return 0;
}

template<typename traits>
int test_basic(int argc_, char* argv_[]){
  std::cout << "========traits:" << typeid(traits).name()
            << "========\n";
  auto a = answer<traits>();
  std::cout << "Type: " << typeid(a).name() << std::endl;
  std::cout << "Got a coroutine, let's get a value" << std::endl;
  auto v = a.get();
  std::cout << "And the coroutine value is: " << v << std::endl;
  v = a.get();
  std::cout << "And the coroutine value is still: " << v << std::endl;
  return 0;
}

template<typename T, class sus_traits=init_yield_final_suspend_always>
struct generator{
  struct promise_type;
  using handle_type = std::coroutine_handle<promise_type>;
  generator(handle_type coro_) : _coro(std::move(coro_)) {
    std::cout << "Created a generator object" << std::endl;
  }
  auto move_next() {
    std::cout << "Moving to next" << std::endl;
    _coro.resume();
    auto still_going = !_coro.done();
    std::cout << (still_going ? "There's another!": "We're done!") << std::endl;
    return still_going;
  }
  auto& current_value() const { return _coro.promise()._v; }
  auto begin() { return iterator{*this, move_next()}; }
  auto end() { return iterator{*this, true}; }

  ~generator(){
    std::cout << "~generator "
      << (!_coro ? "(empty)" : "(contains a coroutine)") << std::endl;
    if(_coro) _coro.destroy();
  }
  generator(const generator& s_) = delete;
  generator& operator=(const generator& s_) = delete;
  generator(generator&& s_): _coro(s_._coro) {
    std::cout << "Generator moved leaving behind a husk " << std::endl;
    s_._coro = nullptr;
  }
  generator& operator=(generator&& s_){
    if(this != &s_){
      std::cout << "Generator copy moved leaving behind a husk " << std::endl;
      _coro = s_._coro;
      s_._coro = nullptr;
    }
    return *this;
  }
  struct iterator{
    generator& _owner;
    bool _done{};
    bool operator!=(const iterator& itor_) const{
      return _done == itor_._done;
    }
    iterator& operator++(){
      _done =_owner.move_next();
      return *this;
    }
    auto& operator*() const { return _owner.current_value(); }
  };
  struct promise_type{
    promise_type(){
      std::cout << "Promise created" << std::endl;
    }
    ~promise_type(){
      std::cout << "Promise died" << std::endl;
    }
    auto initial_suspend(){
      std::cout << "Started the coroutine, let's pause!" << std::endl;
      return typename sus_traits::s_init{};
    }
    auto final_suspend(){
      std::cout << "Finished the coro" << std::endl;
      return typename sus_traits::s_final{};
    }
    auto yield_value(T v_){
      std::cout << "Got an answer of " << v_ << std::endl;
      _v = std::move(v_);
      return typename sus_traits::s_yield{};
    }
    auto return_void(){
      std::cout << "return_void" << std::endl;
      return typename sus_traits::s_return{};
    }
    auto get_return_object(){
      std::cout << "Send back a generator" << std::endl;
      return generator{handle_type::from_promise(*this)};
    }
    void unhandled_exception(){
      std::cout << "Exception..." << std::endl;
      std::exit(1);
    }
    T _v{};
  };
  handle_type _coro;
};

template<class sus_traits> generator<int, sus_traits> count(){
  int v{1};
  while(v <=3 ){
  // for(auto v : {1, 2, 3}){ // it's now fails to parse in gcc master.
    std::cout << "Going to yeild " << v << std::endl;
    co_yield v;
    ++ v;
  }
}

static constexpr auto is_prime(auto n_){
  if(n_ <= 1) return false;
  if(n_ <= 3) return true;
  if(n_ % 2 == 0) return false;
  for(decltype(n_) i = 3;i * i <= n_;i += 2)
    if(n_ % i == 0) return false;
  return true;
}

generator<int> primes(int n_){
  if(n_ >= 2) co_yield 2;
  for(int n = 3;n <= n_;++ n){
    if(is_prime(n)) co_yield n;
  }
}

generator<long long> primes(){
  long long n{3};
  co_yield 2;
  while(n > 2){
    if(is_prime(n)) co_yield n;
    ++ n;
  }
}

template<class traits>
int test_generator(int argc_, char* argv_[]){
  std::cout << "========traits:" << typeid(traits).name()
            << "========\n";
  {
    auto g = count<traits>();
    while(g.move_next())
      std::cout << ">> " << g.current_value() << std::endl;
  }
  if constexpr (std::is_same_v<traits, yield_final_suspend_always>) {
    std::cout << "========do first========" << std::endl;
    auto g = count<traits>();
    do{
      std::cout << ">> " << g.current_value() << std::endl;
    }while(g.move_next());
  }
  return 0;
}

int test_iterator(int argc_, char* argv_[]){
  std::cout << "========" << __func__ << "========\n";
  for(auto v : count<init_yield_final_suspend_always>())
    std::cout << ">> " << v << std::endl;
  return 0;
}

int test_primes(int argc_, char* argv_[]){
  std::cout << "========" << __func__ << "========\n";
  for(auto v : primes(10))
    std::cout << ">> " << v << std::endl;
  return 0;
}

int test_primes_infinity(int argc_, char* argv_[]){
  std::cout << "========" << __func__ << "========\n";
  std::atomic_bool stopped{false};
  std::jthread t([&stopped](std::stop_token t_){
    for(auto v : primes()){
      std::cout << ">> " << v << std::endl;
      if(t_.stop_requested()) break;
    }
  });
  std::cin.get();
  return 0;
}

int test_await(int argc_, char* argv_[]){
  std::cout << "========" << __func__ << "========\n";
  return await_answer<init_yield_final_suspend_always>().get();
}

int main(int argc_, char* argv_[]){
  test_basic<suspend_never_all>(argc_, argv_);
  test_basic<final_suspend_always>(argc_, argv_);
  test_basic<init_suspend_always>(argc_, argv_);
  test_basic<init_final_suspend_always>(argc_, argv_);
  test_generator<init_yield_final_suspend_always>(argc_, argv_);
  test_generator<yield_final_suspend_always>(argc_, argv_);
  test_iterator(argc_, argv_);
  test_await(argc_, argv_);
  test_primes(argc_, argv_);
  test_primes_infinity(argc_, argv_);
  return 0;
}
