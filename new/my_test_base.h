#ifndef MY_LIB_MY_TEST_BASE_H
#define MY_LIB_MY_TEST_BASE_H
#include <vector>
#include <thread>
#include <mutex>
#include <utility>
#include <string>
#include <iostream>
#include <cfloat>
#include <unordered_map>

namespace my_lib{
	class TwoDPlot{
		public:
			virtual void show() = 0;
			virtual void redisplay();
			virtual int run(int argc, char *argv[]);
			virtual void keyboard(unsigned char key, int x, int y){
			}
			virtual ~TwoDPlot(){
			}
			static void drawPoints(double *points, double *colors, unsigned int n);
			static void drawPoints(double *points, unsigned int n, double *color);
			static void drawPath(double *vertices, double *colors, unsigned int n);
			static void drawPath(double *vertices, unsigned int n, double *color);
			static void drawAxis(double minX, double maxX, double minY, double maxY);
			static void drawString(double x, double y, const char *str);
			static void drawEdges(double *vertices, unsigned int* edges, double *colors, unsigned int n);
			static void drawEdges(double *vertices, unsigned int* edges, unsigned int n, double *color);
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
		struct RenderInfo{
			double xMin = DBL_MAX, xMax = -DBL_MAX, yMin = DBL_MAX, yMax = -DBL_MAX; 
			double titleX = 0.02, titleY = 0.02;
			double charSize = 0.0003;
			bool drawAxis = true;
			std::string title;
		};
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
		virtual bool getBounds(double& xMin, double& xMax, double& yMin, double& yMax) const;
		bool getBounds(RenderInfo& renderInfo) const{
			return getBounds(renderInfo.xMin, renderInfo.xMax, renderInfo.yMin, renderInfo.yMax);
		}
		virtual bool getTitle(RenderInfo& renderInfo) const;
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

	class GraphPlot : public TwoDPlot{
		protected:
			std::vector<unsigned int> edges;
			std::vector<double> points;
			double pointSize = 3;
			virtual const char *getString(unsigned int index) const{
				return nullptr; 
			}
		public:
			void show() override;
	};

	class TreePlot : public GraphPlot{
		public:
			template<class T, class F>
				void calculate(T t, F f){
					edges.clear();
					points.clear();
					calculate(t, f, 0.0, 0.2, -0.1);
				}
		private:
			struct NodePosition{
				unsigned int l = 0, r = 0;
				double lB = 0.0, rB = 0.0;
				double pos() const{
					return (lB + rB) / 2;
				}
			};
			template<class T, class F>
			void calculate(T t, F f, double lB, double interval, double hInterval){
				std::vector<NodePosition> positions;
				positions.push_back(NodePosition());
				calculate(t, f, 0, lB, interval, positions);

				std::vector<unsigned int> nodes;
				nodes.push_back(0);
				points.push_back(positions[0].pos());
				points.push_back(0.0);

				for(size_t i = 0;i < nodes.size();++ i){
					NodePosition& node = positions[nodes[i]];
					if(node.l != 0){
						edges.push_back(i);
						edges.push_back(nodes.size());
						nodes.push_back(node.l);
						points.push_back(positions[node.l].pos());
						points.push_back(points[2 * i + 1] + hInterval);
					}
					if(node.r != 0){
						edges.push_back(i);
						edges.push_back(nodes.size());
						nodes.push_back(node.r);
						points.push_back(positions[node.r].pos());
						points.push_back(points[2 * i + 1] + hInterval);
					}
				}
			}
			template<class T, class F>
			void calculate(T t, F f, unsigned int nIndex, double lB, double interval, 
				std::vector<NodePosition>& positions){
				T c;

				NodePosition* n = &positions[nIndex];
				n->lB = lB;
				n->rB = lB + interval;
				if(f.left(t, c)){
					unsigned int l = positions.size(); 
					positions.push_back(NodePosition());
					calculate(c, f, l, lB, interval, positions);
					lB = positions[l].rB;
					n = &positions[nIndex];
					n->rB = lB;
					n->l = l;
				}
				if(f.right(t, c)){
					unsigned int r = positions.size();
					positions.push_back(NodePosition());
					calculate(c, f, r, lB, interval, positions);
					n = &positions[nIndex];
					n->rB = positions[r].rB;
					n->r = r;
				}
			}
	};

	class UndirGraphPlot : public GraphPlot {
	public:
		template<class Graph, class CoordinateTraits, class Stream>
		void init(const Graph& graph, Stream& error, const CoordinateTraits& traits){
			pointSize = 6;
			points.clear();
			edges.clear();
			std::unordered_map<unsigned int, unsigned int> vertexToPoint;
			for(auto& vertex : *graph.getCoordinates()){
				vertexToPoint[vertex.first] = points.size() / 2;
				points.push_back(traits.x(vertex.second));
				points.push_back(traits.y(vertex.second));
			}
			for(auto& vertex : graph.getEdges()){
				auto pos = vertexToPoint.find(vertex.first);
				if(pos == vertexToPoint.end()){
					error << "vertex " << vertex.first << " does not have coordinates\n";
					continue;
				}
				for(auto& v : vertex.second){
					auto pos1 = vertexToPoint.find(v);
					if(pos1 == vertexToPoint.end()){
						error << "vertex " << v << " does not have coordinates\n";
						error << "edge (" << vertex.first << ", " << v << ") ignored\n";
						continue;
					}
					edges.push_back(pos->second);
					edges.push_back(pos1->second);
				}
			}
		}
	};

}
#endif //MY_LIB_MY_TEST_BASE_H     
