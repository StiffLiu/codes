#ifndef _TMLG_POINT3_H_
#define _TMLG_POINT3_H_
#include "Constant.h"
#include <ostream>
namespace tmlg{
	template<class T>
	struct Point3{
		T coord[3];
		Point3(){
			coord[0] = coord[1] = coord[2] = 0;
		}
		Point3(T x, T y, T z){
			this->coord[0] = x;
			this->coord[1] = y;
			this->coord[2] = z;
		}
		Point3(const T coord[3]){
			this->coord[0] = coord[0];
			this->coord[1] = coord[1];
			this->coord[2] = coord[2];
		}
		template<class U>
		Point3& operator=(const U coord[3]){
			this->coord[0] = coord[0];
			this->coord[1] = coord[1];
			this->coord[2] = coord[2];
			return *this;
		}
		const T& x() const{
			return coord[0];
		}
		const T& y() const{
			return coord[1];
		}
		const T& z() const{
			return coord[2];
		}
		T& x() {
			return coord[0];
		}
		T& y() {
			return coord[1];
		}
		T& z() {
			return coord[2];
		}
		operator T*(){
			return coord;
		}
		operator const T*() const{
			return coord;
		}
		T dist(const Point3& pt) const{
			Point3 tmp(*this);
			return tmp.subtract(pt).norm();
		}
		T norm() const{
			return sqrt(coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2]);
		}
		T normSquare() const{
			return coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2];
		}
		void normalize(){
			double norm = this->norm();
			coord[0] /= norm;
			coord[1] /= norm;
			coord[2] /= norm;
		}
		Point3& operator*=(T val){
			return multiply(val);
		}
		Point3 operator-()const{
			Point3 tmp = *this;
			tmp *= (T)-1;
			return tmp;
		}
		Point3& multiply(T val){
			coord[0] *= val;
			coord[1] *= val;
			coord[2] *= val;
			return *this;
		}
		Point3& add(const Point3& pt){
			coord[0] += pt.coord[0];
			coord[1] += pt.coord[1];
			coord[2] += pt.coord[2];
			return *this;
		}
		Point3 proj(const Point3& vec){
			T normS = normSquare();
			if(normS <= Constant<T>::tolerance)
				return vec;
			Point3 vecProj = vec;
			vecProj *=  vec.dotProduct(vec) / normS;
			return vecProj;
		}
		Point3 orth(const Point3& vec){
			Point3 vecOrth = vec;
			vecOrth.subtract(proj(vec));
			return vecOrth;
		}
		Point3& subtract(const Point3& pt){
			coord[0] -= pt.coord[0];
			coord[1] -= pt.coord[1];
			coord[2] -= pt.coord[2];
			return *this;
		}
		Point3& operator-=(const Point3& pt){
			this->subtract(pt);
			return *this;
		}
		static Point3 rotateBy(const Point3& vec, T angle, const Point3& axis){
			T vecNormSquare = vec.normSquare(), axisNormSquare = axis.normSquare();
			if(vecNormSquare <= Constant<T>::tolerance || 
				axisNormSquare <= Constant<T>::tolerance || 
				(angle <= Constant<T>::tolerance && angle >= -Constant<T>::tolerance))
				return vec;
			Point3 vecProj = axis;
			Point3 vecOrtho = vec;
			Point3 tmp;
			vecProj *=  vec.dotProduct(axis) / axisNormSquare;
			vecOrtho.subtract(vecProj);
			tmp = axis.crossProduct(vecOrtho);
			if(tmp.normSquare() <= Constant<T>::tolerance){
				tmp = axis;
				if(vec.dotProduct(axis) < 0)
					tmp *= -1;
				tmp.normalize();
				return tmp;
			}
			else{
				Point3 rotatedVecOrth = vecOrtho;
				T vecOrthoNorm = vecOrtho.norm();
				vecOrtho.normalize();
				tmp.normalize();
				tmp *= sin(angle) * vecOrthoNorm;
				rotatedVecOrth *= cos(angle);
				rotatedVecOrth.add(tmp);
				rotatedVecOrth.add(vecProj);
				return rotatedVecOrth;
			}
		}
		double dotProduct(const Point3& pt) const{
			return coord[0] * pt.coord[0] + coord[1] * pt.coord[1] + coord[2] * pt.coord[2];
		}
		friend Point3 operator-(const Point3& a, const Point3& b){
			Point3 tmp = a;
			tmp.subtract(b);
			return tmp;
		}
		friend Point3 operator+(const Point3& a, const Point3& b){
			Point3 tmp = a;
			tmp.add(b);
			return tmp;
		}
		Point3 crossProduct(const Point3& pt) const{
			return Point3(coord[1] * pt.coord[2] - coord[2] * pt.coord[1],
				coord[2] * pt.coord[0] - coord[0] * pt.coord[2],
				coord[0] * pt.coord[1] - coord[1] * pt.coord[0]);
		}
		bool operator==(const Point3& pt) const{
			return coord[0] == pt.coord[0] && coord[1] = pt.coord[1] && coord[2] = pt.coord[2];
		}
		friend std::ostream& operator<<(std::ostream& os, const Point3& t){
			using namespace std;
			os << '(' << t.coord[0] << ',' << t.coord[1] << ',' << t.coord[2] << ')';
			return os;
		}
	};
	typedef Point3<float> Point3F;
	typedef Point3<double> Point3D;
}
#endif //_TMLG_POINT3_H_
