#ifndef _TMLG_CONSTANT_H_
#define _TMLG_CONSTANT_H_
namespace tmlg{
	template<class T>
	struct Constant{
		typedef T type;
		const static T pi;
		const static T tolerance;
		const static T cameraNearLimit;
	};
	typedef Constant<double> ConstantD;
	typedef Constant<float> ConstantF;
}
#endif //_TMLG_CONSTANT_H_
