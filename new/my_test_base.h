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
	virtual ~TwoDPlot(){
	}
	static void drawPoints(double *points, double *colors, unsigned int n);
	static void drawPoints(double *points, unsigned int n, double *color);
	static void drawPath(double *vertices, double *colors, unsigned int n);
	static void drawPath(double *vertices, unsigned int n, double *color);
	static void drawAxis(double minX, double maxX, double minY, double maxY);
protected:
	virtual void init(){
	}
};
class StatPlotBase : public TwoDPlot{
	typedef TwoDPlot Super;
public:
	StatPlotBase(unsigned int numGraphs, double *colors = nullptr)
	       : graphs(numGraphs), colors(colors, 
		(colors == nullptr ? colors : colors + 3 * numGraphs)){
	}
	size_t getNumGraphs(){
		return graphs.size();
	}
	void show() override;
	~StatPlotBase() ;
protected:
	std::vector<std::vector<double> > graphs;
	std::vector<double> colors;
	std::thread collectingThread;
	std::mutex m;
	unsigned int maxNumPoints = 1000;
	unsigned int interval = 10;
	static void collecting(StatPlotBase *plot);	
	virtual bool collect(double *values){
		return false;
	}
};
template<class Collector>
class StatPlot : public StatPlotBase{
	Collector collector;
protected:
	StatPlot(unsigned int numGraphs, const Collector& collector, 
		double *colors = nullptr) : StatPlotBase(numGraphs, colors), collector(collector){
	}
	void init() override{
		collectingThread = std::thread(collecting, this);
	}
	bool collect(double *values) override{
		return collector(values);
	}
};
}
#endif //MY_LIB_MY_TEST_BASE_H     
