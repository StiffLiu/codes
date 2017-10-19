#ifndef _era_numeric_JohansenCointegration_H_
#define _era_numeric_JohansenCointegration_H_

#include "../Settings.h"
namespace era{
  namespace numeric{
    class ERA_LIB_EXPORT JohansenCointegration{
    public:
      JohansenCointegration(double* val_, unsigned int row_, unsigned int col_, unsigned int lags_);
      JohansenCointegration(const JohansenCointegration&) = delete;
      JohansenCointegration& operator=(const JohansenCointegration&) = delete;
      JohansenCointegration(JohansenCointegration&& jc_){
        _impl = jc_._impl;
        jc_._impl = nullptr;
      };
      ~JohansenCointegration();
      struct Test {
        double _s, _90, _95, _99;
      };
      class Impl;

      unsigned int count() const;
      const double* vector(unsigned int index_) const;
      const double* values() const;
      unsigned int numVars() const;
      const Test* tests() const;
      unsigned int numTests() const;

    private:
      Impl* _impl{};
    };
  }
}

#endif
