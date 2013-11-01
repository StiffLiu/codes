#ifndef _TMLG_CYLINDER_H_
#define _TMLG_CYLINDER_H_
#include "Transform.h"
#include <cmath>
namespace tmlg{
	template<class T>
	class Cylinder{	
		Point3<T> ptBottom;
		Point3<T> ptTop;
		T topRadius;
		T botRadius;
		Transform<T> modelMatrix;
	public:		
		Cylinder(const T ptB[3], const T ptT[3], T topRadius, T botRadius){
			ptBottom = ptB;
			ptTop = ptT;
			this->topRadius = topRadius;
			this->botRadius = botRadius;
			modelMatrix = Transform<T>::computeTransform(ptBottom, ptTop, Point3<T>(0, 1, 0));
		}
		const Transform<T>& getModelMatrix() const{
			return modelMatrix;
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
			if(std::min(topRadius, botRadius) <= Constant<T>::tolerance) return false;
			return true;
		}
		bool calcData(T *body, T *bodyNorm, T *topCap, T *botCap, unsigned int slice)const{
			if(!isValid(slice)) return false;
			T inc = 2 * Constant<T>::pi / slice;
			T height = ptBottom.dist(ptTop);
			if(body != NULL){
				Point3<T>* pBody = reinterpret_cast<Point3<T>*>(body), 
					*pNorm = reinterpret_cast<Point3<T>*>(bodyNorm);
				T diff = topRadius - botRadius;
				Point3<T> projAxis = ptTop - ptBottom;
				projAxis *= diff / sqrt(diff * diff + height * height);
				for(unsigned int i = 0;i <= slice;++ i, pBody += 2, pNorm += 2){
					T sinInc = sin(i * inc), cosInc = cos(i * inc);
					*pBody = modelMatrix.compute((T)cosInc * botRadius, (T)0, (T)sinInc * botRadius);
					*(pBody + 1) = modelMatrix.compute((T)cosInc * topRadius, height, (T)sinInc * topRadius);	
					*pNorm = *pBody - ptBottom;
					pNorm->add(projAxis);
					pNorm->normalize();
					*(pNorm + 1) = *pNorm;
				}
			}
			if(topCap != NULL){
				Point3<T>* pTopCap = reinterpret_cast<Point3<T>*>(topCap);
				*pTopCap = ptTop; ++ pTopCap;
				for(unsigned int i = 0;i <= slice;++ i, ++ pTopCap){
					*pTopCap = modelMatrix.compute((T)cos(i * inc)* topRadius, height, (T)sin(i * inc)* topRadius);
				}
			}
			if(botCap != NULL){
				Point3<T>* pBotCap = reinterpret_cast<Point3<T>*>(botCap);
				*pBotCap = ptBottom; ++ pBotCap;
				for(unsigned int i = 0;i <= slice;++ i, ++ pBotCap){
					*pBotCap = modelMatrix.compute((T)cos(i * inc)* botRadius, (T)0, (T)sin(i * inc)* botRadius);
				}
			}
			return true;
		}
		static const unsigned int minSlice = 3;
	};
	typedef Cylinder<float> CylinderF;
	typedef Cylinder<double> CylinderD;
}
#endif //_TMLG_CYLINDER_H_
