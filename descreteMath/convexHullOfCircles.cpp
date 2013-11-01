#include <vector>
const double tolerance = 1e-8;

template<class T>
int whichSide(T *pStart, T *pEnd, T* pt, bool up = true){
	T deltaX1 = pEnd[0] - pStart[0];
	T deltaY1 = pEnd[1] - pStart[1];
	T deltaX2 = pt[0] - pEnd[0];
	T deltaY2 = pt[1] - pEnd[1];
	T ret = deltaX1 * deltaY2 - deltaX2 * deltaY1;
	if(ret > tolerance)
		return 1;
	if(ret < tolerance)
		return - 1;
	return 0;

}
template<class T>
bool isIn(T *c1, T *c2){
	T deltaX1 = c2[0] - c1[0];
	T deltaY1 = c2[1] - c1[1];
	T deltaR = c2[2] - c1[2];
	return deltaR < tolerance && (deltaX1 * deltaX1 + deltaY1 * deltaY1 < deltaR * deltaR);
}
template<class T, class Seq>
void convexHull(T* pts, int ptCount, Seq& ret, int stride = 3){
	if(ptCount < 1)
		return;
	if(ptCount <= 2){
		ret.push_back(0);
		if(ptCount > 1)
			ret.push_back(1);
		return;
	}
	ret.push_back(0);
	for(int i = 1;i < ptCount;++ i){
		ret.push_back(i);
		while(ret.size() > 1 && isIn(pts + ret[ret.size() - 2] * stride, pts + ret[ret.size() - 1] * stride))
			ret.pop_back();
		while(ret.size() > 2 && whichSide(pts + ret[ret.size() - 3] * stride,
			pts + ret[ret.size() - 2] * stride, pts + ret[ret.size() - 1] * stride) != -1){
			ret[ret.size() - 2] = ret[ret.size() - 1];
			ret.pop_back();
		}
	}
	for(int i = ptCount - 1;i >= 0;-- i){
		ret.push_back(i);
		while(ret.size() > 1 && isIn(pts + ret[ret.size() - 2] * stride, pts + ret[ret.size() - 1] * stride))
			ret.pop_back();
		while(ret.size() > 2 && whichSide(pts + ret[ret.size() - 3] * stride,
			pts + ret[ret.size() - 2] * stride, pts + ret[ret.size() - 1] * stride, false) != -1){
			ret[ret.size() - 2] = ret[ret.size() - 1];
			ret.pop_back();
		}
	}
	assert(ret[0] == ret.back());
}
