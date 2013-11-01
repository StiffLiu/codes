#ifndef _TMLG_COLOR_H_
#define _TMLG_COLOR_H_
namespace tmlg{
	void hsl2rgb(double h, double sl, double l, double& r, double& g, double& b);
	void rgb2hsl (double r, double g, double b, double& h, double& s, double& l);
	template<class T>
	class Color{
		T xC;
		T yC;
		T zC;
	public:
		operator const T*() const{return (const T*)this;}
		Color() : xC(0), yC(0), zC(0){}
		Color(T x, T y, T z) : xC(x), yC(y), zC(z){}
		const T& x()const{ return xC;}
		const T& y()const{ return yC;}
		const T& z()const{ return zC;}
		T& x(){ return xC;}
		T& y(){ return yC;}
		T& z(){ return zC;}
		Color toRGB() const{
			double x, y, z;
			hsl2rgb(xC, yC, zC, x, y, z);
			return Color((T)x, (T)y, (T)z);
		}
		Color toHSL() const{
			double x, y, z;
			rgb2hsl(xC, yC, zC, x, y, z);
			return Color((T)x, (T)y, (T)z);
		}
	};
	template<class T>
	class ColorA : public Color<T>{
		T aC;
	public:
		ColorA() : aC(1){}
		ColorA(T r, T g, T b, T a) : Color<T>(r, g, b), aC(a){}
		const T& a()const{ return aC;}
		T& a(){ return aC;}
	};
	typedef Color<float> ColorF;
	typedef Color<double> ColorD;	
	typedef ColorA<float> ColorAF;
	typedef ColorA<double> ColorAD;	
}
#endif
