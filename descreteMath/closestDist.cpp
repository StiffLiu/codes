#include <iostream>
#include <cfloat>
#include <cmath>
#include <set>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cassert>

using namespace std;
/*inline double calcNormal(const pair<double, double>& p1, const pair<double, double>& p2){
	double deltaX = p1.first - p2.first;
	double deltaY = p1.second - p2.second;
	return sqrt(deltaX * deltaX + deltaY * deltaY);
}*/
inline double calcNormal(const pair<double, double>& p1, const pair<double, double>& p2){
	double deltaX = fabs(p1.first - p2.first);
	double deltaY = fabs(p1.second - p2.second);
	//return sqrt(deltaX * deltaX + deltaY * deltaY);
	if(deltaX > deltaY)
		return deltaX;
	return deltaY;
}
double closestDist0(const pair<pair<double, double>, int>* sortedByX, 
	const pair<pair<double, double>, int>* sortedByY, const unsigned int n, pair<int, int>* ret){
	if(n <= 1)
		return DBL_MAX;
	if(n == 2){
		if(ret != NULL){
			ret->first = sortedByY[0].second;
			ret->second = sortedByY[1].second;
		}
		return calcNormal(sortedByX[0].first, sortedByX[1].first);
	}
	unsigned int half = n / 2;
	pair<pair<double, double>, int>* firstHalf = new pair<pair<double, double>, int>[half + 1];
	pair<pair<double, double>, int>* secondHalf = new pair<pair<double, double>, int>[n - 1 - half];
	{
		int s = 0, t = 0;
		for(int i = 0;i < n;++ i)
			if(sortedByY[i].second <= sortedByX[half].second){
				//cout << i << ',' << s << ',' << ',' << sortedByY[i].second << ',' << sortedByX[half].second << ',' << half << endl;
				//assert(s < half + 1);
				firstHalf[s] = sortedByY[i];
				++ s;
			}else{
				//cout << i << ',' << t << ',' << ',' << sortedByY[i].second << ',' << sortedByX[half].second << ',' << half << endl;
				//cin.get();
				//assert(t < n - 1 - half);
				secondHalf[t] = sortedByY[i];
				++ t;
			}
		//cout << sortedByX[half].second << ':' <<  s << ',' << t << "-----" << (half + 1) << ',' << (n - 1 - half) << endl;
		//cin.get();
	}
	double mid = sortedByX[half].first.first;
	pair<int, int> ret1;
	double dist = closestDist0(sortedByX, firstHalf, half + 1, ret);
	double dist2 = closestDist0(sortedByX + half + 1, secondHalf, n - 1 - half, &ret1);
	if(dist > dist2){
		dist = dist2;
		if(ret != NULL)
			*ret = ret1;
	}
	for(unsigned int i = 0;i < n;++ i)
		if(fabs(sortedByY[i].first.first - mid) <= dist){
			int end = min(n - i - 1, (unsigned int)7);
			for(int j = 1;j <= end;++ j){
				//double deltaX = sortedByY[i].first.first - sortedByY[i + j].first.first;
				//double deltaY = sortedByY[i].first.second - sortedByY[i + j].first.second;
				//double newDist = sqrt(deltaX * deltaX + deltaY * deltaY);
				double newDist = calcNormal(sortedByY[i].first, sortedByY[i + j].first);
				if(dist > newDist){
					dist = newDist;
					if(ret != NULL){
						ret->first = sortedByY[i].second;
						ret->second = sortedByY[i + j].second;
					}
				}
			}
		}
	delete []firstHalf;
	delete []secondHalf;
	return dist;			
}
double closestDist1(const pair<double, double>* points, unsigned int n, pair<int, int>* ret = NULL){
	double dist = DBL_MAX;
	for(unsigned int i = 0;i < n;++ i)
		for(unsigned int j = 0;j < i;++ j){
			//double deltaX = points[j].first - points[i].first;
			//double deltaY = points[j].second - points[i].second;
			//double newDist = sqrt(deltaX * deltaX + deltaY * deltaY);
			double newDist = calcNormal(points[j], points[i]);
			if(newDist < dist){
				dist = newDist;
				if(ret != NULL){
					ret->first = i;
					ret->second = j;
				}
			}
		}
	return dist;
}
int main(int argc, char *argv[]){
	vector<pair<double, double> > points;
	vector<pair<pair<double, double>, int> > sortedByX;
	vector<pair<pair<double, double>, int> > sortedByY;
	{
		const int pointCount = 10000000;
		set<pair<double ,double> > pts;
		srand(time(0));
		for(int i = 0;i < pointCount;++ i)
			pts.insert(make_pair((double)rand(), (double)rand()));
		set<pair<double ,double> >::iterator begin = pts.begin(), end = pts.end();
		while(begin != end){
			points.push_back(*begin);
			++begin;
		}
		cout << "points count: " << points.size() << endl;
	}
	{
		set<pair<pair<double ,double>, int> > pts;
		for(int i = 0;i < points.size();++ i){
			sortedByX.push_back(make_pair(make_pair(points[i].first, points[i].second), i));
			pts.insert(make_pair(make_pair(points[i].second, points[i].first), i));
		}
		set<pair<pair<double ,double>, int> >::iterator begin = pts.begin(), end = pts.end();
		while(begin != end){
			sortedByY.push_back(make_pair(make_pair(begin->first.second, begin->first.first), begin->second));
			++begin;
		}
		assert(sortedByY.size() == sortedByX.size());
	}
	pair<int, int> ret;
	clock_t start = clock();
	/*double dist = closestDist1(&points[0], points.size(), &ret);
	cout << "time consumed:" << (clock() - start) / (double)CLOCKS_PER_SEC << 
		", minimum distance: " << dist << ", points: (" << sortedByX[ret.first].first.first << ',' << 
		sortedByX[ret.first].first.second << ") (" << sortedByX[ret.second].first.first << ',' << 
		sortedByX[ret.second].first.second << ")" << endl;*/
	start = clock();
	double dist = closestDist0(&sortedByX[0], &sortedByY[0], sortedByX.size(), &ret);
	cout << "time consumed:" << (clock() - start) / (double)CLOCKS_PER_SEC << 
		", minimum distance: " << dist << ", points: (" << sortedByX[ret.first].first.first << ',' << 
		sortedByX[ret.first].first.second << ") (" << sortedByX[ret.second].first.first << ',' << 
		sortedByX[ret.second].first.second << ")" << endl;
}
