#ifndef _TMLG_EVALUATOR_H_
#define _TMLG_EVALUATOR_H_
#include "Point3.h"
namespace tmlg{
	template<class T>
	class Evaluator{
	public:
		//(rows + 1) * (cols + 1) points
		//2 * rows * cols tristrip;
		//2 * cols tristrip every row;
		
		void eval(T startX, T endX, T startY, T endY,
			int rows, int cols, 
			Point3<T>* points, T (*evalFunc)(T,T),
			Point3<T>* norms = NULL, Point3<T> (*normFunc)(T, T) = NULL,
			Point3<T>* texCoord = NULL){
			T stepX = (endX - startX) / cols;
			T stepY = (endY - startY) / rows;
			rows += 1;
			cols += 1;
			for(int i = 0;i < rows;++ i){
				T y = startY + i * stepY;
				T texCoordY = i / (double)(rows - 1);
				int index = i * cols;
				for(int j = 0;j < cols;++ j){
					T x = startX + j * stepX;
					points[index + j].x() = x;
					points[index + j].y() = y;
					points[index + j].z() = evalFunc(x, y);
					if(norms != NULL){
						*(norms + index + j) = normFunc(x, y);
					}
					if(texCoord != NULL){
						texCoord[index + j].x() = (j / (double)(cols - 1));
						texCoord[index + j].y() = texCoordY;
					}
				}
			}			
		}
		void eval(T startU, T endU, T startV, T endV,
			int rows, int cols, 
			Point3<T>* points, Point3<T> (*evalFunc)(T,T),
			Point3<T>* norms = NULL, Point3<T> (*normFunc)(T, T) = NULL,
			Point3<T>* texCoord = NULL){
			T stepU = (endU - startU) / cols;
			T stepV = (endV - startV) / rows;
			rows += 1;
			cols += 1;
			for(int i = 0;i < rows;++ i){
				T v = startV + i * stepV;
				T texCoordY = i / (double)(rows - 1);
				int index = i * cols;
				for(int j = 0;j < cols;++ j){
					T u = startU + j * stepU;
					*(points + index + j) = evalFunc(u, v);
					if(norms != NULL){
						*(norms + index + j) = normFunc(u, v);
					}
					if(texCoord != NULL){
						texCoord[index + j].x() = (j / (double)(cols - 1));
						texCoord[index + j].y() = texCoordY;
					}
				}
			}			
		}
	};
	typedef Evaluator<float> EvaluatorF;
	typedef Evaluator<double> EvaluatorD;
}
#endif
