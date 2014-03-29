#ifndef MY_LIB_MY_TEST_BASE_H
#define MY_LIB_MY_TEST_BASE_H
#include <vector>
#include <thread>
#include <mutex>

namespace my_lib{
class TwoDPlot{
public:
	virtual void show() = 0;
	virtual int run(int argc, char *argv[]);
	virtual ~TwoDPlotter() = 0;
	static void drawPoints(double *points, double *colors, unsigned int n);
	static void drawPoints(double *points, unsigned int n, double *color);
	static void drawPath(double *vertices, double *colors, unsigned int n);
	static void drawPath(double *vertices, unsigned int n, double *color);
	static void drawAxis(double minX, double maxX, double minY, double maxY);
protected:
	virtual void init(){
	}
};
template<class Collector>
class StatPlot: public TwoDPlot{
	std::vector<std::vector<double> > graphs;
	std::vector<double> colors;
	Collector collector;
	std::thread collectingThread;
	std::mutex m;
	unsigned int interval = 1000;
	static void collecting(StatPlot *plot);	
	typedef TwoDPlot Super;
public:
	StatPlot(unsigned int numGraphs, const Collector& collector, 
		double *colors = nullptr) : graphs(numGraphs), 
		colors(colors, colors + 3 * numGraphs), collector(collector){
	}
	size_t getNumGraphs(){
		return graphs.size();
	}
	override void show();
	override int run(int argc, char *argv[]);
protected:
	override void init();
}
};
#endif //MY_LIB_MY_TEST_BASE_H     
