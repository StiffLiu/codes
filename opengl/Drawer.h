#ifndef _TMLG_DRAWER_H_
#define _TMLG_DRAWER_H_
#include "Evaluator.h"
#include "Cylinder.h"
#include "Cone.h"
namespace tmlg{
	template<class T>
	struct CylinderGeometry{
		Cylinder<T> *c;
		int slice;
		bool drawBotCap;
		bool drawTopCap;
		bool hasTexture;
		CylinderGeometry() : c(NULL), slice(100), drawBotCap(false), drawTopCap(false), hasTexture(false){}
	};
	typedef CylinderGeometry<float> CylinderGeometryF;
	typedef CylinderGeometry<double> CylinderGeometryD;
	template<class T>
	struct ConeGeometry{
		Cone<T> *c;
		bool drawCap;
		bool hasTexture;
		int slice;
		ConeGeometry() : c(NULL), slice(100), drawCap(false), hasTexture(false){}
	};
	typedef ConeGeometry<float> ConeGeometryF;
	typedef ConeGeometry<double> ConeGeometryD;
	template<class T>
	struct BoxGeometry{
		Point3<T> minPt;
		Point3<T> maxPt;
		bool hasTexture;
		BoxGeometry() : hasTexture(false){}
	};
	typedef BoxGeometry<float> BoxGeometryF;
	typedef BoxGeometry<double> BoxGeometryD;
	template<class T>
	struct SurfaceGeometry{
		T startX, endX, startY, endY;
		int rows, cols;
		bool hasTexture;
		bool calculateNorm;
		T (*evalFunc)(T,T);
		Point3<T> (*normFunc)(T, T);
		SurfaceGeometry() : hasTexture(false), evalFunc(NULL), normFunc(NULL), calculateNorm(true){}
	};
	typedef SurfaceGeometry<float> SurfaceGeometryF;
	typedef SurfaceGeometry<double> SurfaceGeometryD;
	template<class T>
	struct ParametricSurfaceGeometry{
		T startU, endU, startV, endV;
		int rows, cols;
		bool hasTexture;
		bool calculateNorm;
		Point3<T> (*evalFunc)(T,T);
		Point3<T> (*normFunc)(T, T);
		ParametricSurfaceGeometry() : hasTexture(false), evalFunc(NULL), normFunc(NULL), calculateNorm(true){}
	};
	typedef ParametricSurfaceGeometry<float> ParametricSurfaceGeometryF;
	typedef ParametricSurfaceGeometry<double> ParametricSurfaceGeometryD;
	class Drawer{
	public:
		void draw(const CylinderGeometryF&);
		void draw(const ConeGeometryF&);
		void draw(const BoxGeometryF&);
		void draw(const SurfaceGeometryF&);
		void draw(const ParametricSurfaceGeometryF&);
	};
}
#endif
