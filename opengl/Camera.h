#ifndef _TMLG_CAMERA_H_
#define _TMLG_CAMERA_H_
#include "Transform.h"
#include "Point3.h"
namespace tmlg{
	template<class T>
	class Camera{
		Point3<T> eyePos;
		Point3<T> target;
		Point3<T> upVector;
	public:
		Camera(const T *eyePos, const T *target, const T *upVector)
			: eyePos(eyePos), target(target), upVector(upVector){
		}
		bool isValid(){
			if(eyePos.dist(target) < Constant<T>::cameraNearLimit)
				return false;
			if(upVector.crossProduct(target - eyePos).norm() < Constant<T>::tolerance)
				return false;
			return true;
		}
		const Point3<T>& getTarget() const{
			return target;
		}
		const Point3<T>& getEyePos() const{
			return eyePos;
		}
		const Point3<T>& getUpVector() const{
			return upVector;
		}
		Point3<T> directionVector(){
			return target - eyePos;
		}
		bool setTarget(const Point3<T>& newTarget){
			if(eyePos.dist(newTarget) < Constant<T>::cameraNearLimit)
				return false;
			if(upVector.crossProduct(newTarget - eyePos).norm() < Constant<T>::tolerance)
				return false;
			target = newTarget;
			return true;
		}
		bool setEyePos(const Point3<T>& newEyePos){
			if(newEyePos.dist(target) < Constant<T>::cameraNearLimit)
				return false;
			if(upVector.crossProduct(target - newEyePos).norm() < Constant<T>::tolerance)
				return false;
			eyePos = newEyePos;
			return true;
		}
		bool setUpVector(const Point3<T>& newUpVector){
			if(newUpVector.crossProduct(target - eyePos).norm() < Constant<T>::tolerance)
				return false;
			upVector = newUpVector;
			return true;
		}
		void orbit(T angle){
			Point3<T> direction = directionVector();
			Transform<T> tmp(angle, upVector);
			Point3<T> distDirection = tmp.compute(direction);
			eyePos = target - distDirection;
		}
		bool zoom(T dist){
			if(eyePos.dist(target) - dist < Constant<T>::cameraNearLimit)
				return false;
			Point3<T> direction = directionVector();
			direction.normalize();
			direction *= dist;
			eyePos.add(direction);
			return true;
		}
		bool movePositionX(T dist){
			Point3<T> newEyePos = eyePos;
			newEyePos.x() += dist;
			return setEyePos(newEyePos);
		}
		bool movePositionY(T dist){
			Point3<T> newEyePos = eyePos;
			newEyePos.y() += dist;
			return setEyePos(newEyePos);
		}
		bool movePositionZ(T dist){
			Point3<T> newEyePos = eyePos;
			newEyePos.z() += dist;
			return setEyePos(newEyePos);
		}
		void roll(T angle){
			Point3<T> direction = directionVector();
			Transform<T> tmp(angle, upVector);
			Point3<T> distDirection = tmp.compute(direction);
			target = eyePos + distDirection;
		}
		friend std::ostream& operator<<(std::ostream& os, const Camera& camera){
			using namespace std;
			const char *sep = "  ";
			os << "eye_position=" << camera.getEyePos() << sep;
			os << "target=" << camera.getTarget() << sep;
			os << "up_vector=" << camera.getUpVector();
			return os;
		}
	};
	typedef Camera<float> CameraF;
	typedef Camera<double> CameraD;	
}
#endif //_TMLG_CAMERA_H_
