#ifndef _TMLG_TRANSFORM_H_
#define _TMLG_TRANSFORM_H_
#include "Point3.h"
namespace tmlg{
	template<class T>
	class Transform{
		T m[16];
	public:
		Transform(){
			setIdentity();
		}
		Transform(T angle, const Point3<T>& axis){
			setIdentity();
			Point3<T> rotatedX = Point3<T>::rotateBy(Point3<T>(1, 0, 0), angle, axis);
			Point3<T> rotatedY = Point3<T>::rotateBy(Point3<T>(0, 1, 0), angle, axis);
			Point3<T> rotatedZ = Point3<T>::rotateBy(Point3<T>(0, 0, 1), angle, axis);
			m[0] = rotatedX.x();m[4] = rotatedX.y();m[8]=rotatedX.z();
			m[1] = rotatedY.x();m[5] = rotatedY.y();m[9]=rotatedY.z();
			m[2] = rotatedZ.x();m[6] = rotatedZ.y();m[10]=rotatedZ.z();
		}
		Transform(T x, T y, T z){
			setIdentity();
			m[3] = x;
			m[7] = y;
			m[11] = z;
		}
		//(*this) t
		Transform& multiply(const Transform& t){
			for(unsigned int i = 0;i < 16;i += 4){
				T tmp0 = m[i] * t.m[0] + m[i+1] * t.m[4] + m[i+2] * t.m[8] + m[i+3] * t.m[12];
				T tmp1 = m[i] * t.m[1] + m[i+1] * t.m[5] + m[i+2] * t.m[9] + m[i+3] * t.m[13];
				T tmp2 = m[i] * t.m[2] + m[i+1] * t.m[6] + m[i+2] * t.m[10] + m[i+3] * t.m[14];
				T tmp3 = m[i] * t.m[3] + m[i+1] * t.m[7] + m[i+2] * t.m[11] + m[i+3] * t.m[15];
				m[i] = tmp0;m[i+1] = tmp1;m[i+2] = tmp2;m[i+3] = tmp3;
			}
			return *this;
		}
		bool operator==(const Transform& t){
			T sum = 0;
			T tolerance = Constant<T>::tolerance * Constant<T>::tolerance;
			for(int i = 0;i < 16;++ i){
				double delta = m[i] - t.m[i];
				sum += delta * delta;
				if(sum >= tolerance)
					return false;
			}
			return true;
		}
		void setIdentity(){
			m[1] = m[2] = m[3] = 0;
			m[4] = m[6] = m[7] = 0;
			m[8] = m[9] = m[11] = 0;
			m[12] = m[13] = m[14] = 0;
			m[0] = m[5] = m[10] = m[15] = 1;
		}
		Transform(const T* m){
			for(unsigned int i = 0;i < 16;++ i)
				this->m[i] = m[i];
		}
		void transpose(){
			using std::swap;
			swap(m[1], m[4]);
			swap(m[2], m[8]);
			swap(m[3], m[12]);
			swap(m[6], m[9]);
			swap(m[7], m[13]);
			swap(m[11], m[14]);
		}
		template<class U>
		Point3<U> compute(const Point3<U>& pt) const{
			return Point3<U>(m[0] * pt.x() + m[1] * pt.y() + m[2] * pt.z() + m[3], 
				m[4] * pt.x() + m[5] * pt.y() + m[6] * pt.z() + m[7], 
				m[8] * pt.x() + m[9] * pt.y() + m[10] * pt.z() + m[11]);
		}
		template<class U>
		Point3<U> compute(U x, U y, U z) const{
			return Point3<U>(m[0] * x + m[1] * y + m[2] * z + m[3], 
				m[4] * x + m[5] * y + m[6] * z + m[7], 
				m[8] * x + m[9] * y + m[10] * z + m[11]);
		}
		friend std::ostream& operator<<(std::ostream& os, const Transform& t){
			using namespace std;
			const char *sep = "\t  ";
			os << t.m[0] << sep << t.m[1] << sep << t.m[2] << sep << t.m[3] << endl;
			os << t.m[4] << sep << t.m[5] << sep << t.m[6] << sep << t.m[7] << endl;
			os << t.m[8] << sep << t.m[9] << sep << t.m[10] << sep << t.m[11] << endl;
			os << t.m[12] << sep << t.m[13] << sep << t.m[14] << sep << t.m[15] << endl;
			return os;
		}
		static Transform computeTransform(const Point3<T>& ptBottom, const Point3<T>& ptTop, const Point3<T>& axis){
			Point3<T> tmp(ptTop);
			T rotateAngle = 0;
			Point3<T> rotateAxis = axis.crossProduct(tmp.subtract(ptBottom));
			Transform tempMatrix;
			if(ptBottom.norm() > Constant<T>::tolerance){
				tempMatrix = Transform(ptBottom.x(), ptBottom.y(), ptBottom.z());
			}
			if(rotateAxis.norm() > Constant<T>::tolerance){
				rotateAngle = acos(tmp.dotProduct(axis) / tmp.norm());
				rotateAxis.normalize();
			}else if(tmp.dotProduct(axis) < 0){
				rotateAngle = Constant<T>::pi / 2;
				rotateAxis = Point3<T>(0, 1, 0);
			}
			//rotateAngle = rotateAngle / Constant<T>::pi * 180;
			if(rotateAngle  > Constant<T>::tolerance && rotateAxis.norm() > Constant<T>::tolerance)
				tempMatrix.multiply(Transform(rotateAngle, rotateAxis));
			return tempMatrix;
		}
	};
	typedef Transform<float> TransformF;
	typedef Transform<double> TransformD;
}
#endif //_TMLG_TRANSFORM_H_
