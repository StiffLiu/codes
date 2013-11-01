#ifndef _TMLG_CONE_H_
#define _TMLG_CONE_H_
#include "Transform.h"
#include "Cylinder.h"
namespace tmlg{
	template<class T>
	class Cone{	
		Point3<T> ptBottom;
		Point3<T> ptTop;
		T radius;
		Transform<T> modelMatrix;
	public:	
		Cone(const T ptB[3], const T ptT[3], T radius){
			ptBottom = ptB;
			ptTop = ptT;
			this->radius = radius;			
			modelMatrix = Transform<T>::computeTransform(ptBottom, ptTop, Point3<T>(0, 1, 0));
		}
		Point3<T> getAxis() const{
			return ptTop - ptBottom;
		}
		T getHeight() const {
			return ptBottom.dist(ptTop);
		}
		bool isValid(unsigned int slice)const{
			T height = ptBottom.dist(ptTop);
			if(height <= Constant<T>::tolerance) return false;			
			if(slice < minSlice) return false;
			if(radius <= Constant<T>::tolerance) return false;
			return true;
		}
		bool calcData(T *body, T *bodyNorm, T *cap, unsigned int slice)const{
			if(!isValid(slice)) return false;
			T inc = 2 * Constant<T>::pi / slice;
			T height = ptBottom.dist(ptTop);
			if(body != NULL){
				Point3<T>* pBody = reinterpret_cast<Point3<T>*>(body), 
					*pNorm = reinterpret_cast<Point3<T>*>(bodyNorm);				
				Point3<T> projAxis = getAxis();
				*pNorm = projAxis;pNorm->normalize();++ pNorm;
				projAxis *= radius / sqrt(radius * radius + height * height);
				*pBody = ptTop;++ pBody;
				for(unsigned int i = 0;i <= slice;++ i, ++ pBody, ++ pNorm){
					T sinInc = sin(i * inc), cosInc = cos(i * inc);
					*pBody = modelMatrix.compute((T)cosInc * radius, (T)0, (T)sinInc * radius);
					*pNorm = *pBody - ptBottom;
					pNorm->add(projAxis);
					pNorm->normalize();
				}
			}
			if(cap != NULL){
				Point3<T>* pCap = reinterpret_cast<Point3<T>*>(cap);
				*pCap = ptBottom; ++ pCap;
				for(unsigned int i = 0;i <= slice;++ i, ++ pCap){
					*pCap = modelMatrix.compute((T)cos(i * inc)* radius, (T)0, (T)sin(i * inc)* radius);
				}
			}
			return true;
		}
		static const unsigned int minSlice = 3;
	};
	typedef Cone<float> ConeF;
	typedef Cone<double> ConeD;
}
#endif // _TMLG_CONE_H_
