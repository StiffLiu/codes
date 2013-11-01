#include "Constant.h"
#include <cmath>
using namespace tmlg;
template<>
const ConstantD::type ConstantD::pi = asin(1.0)*2;
template<>
const ConstantD::type ConstantD::tolerance = 1.e-10;
template<>
const ConstantD::type ConstantD::cameraNearLimit = 0.001;
template<>
const ConstantF::type ConstantF::pi =  ConstantD::pi;
template<>
const ConstantF::type  ConstantF::tolerance = 1.e-6;
template<>
const ConstantF::type  ConstantF::cameraNearLimit = 0.001;

