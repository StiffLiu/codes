#ifndef _MYCPPCODE_COMMON_H
#define _MYCPPCODE_COMMON_H
#include <vector>
#include <sstream>
#include <GL/glut.h>
#include <ctime>
#include <random>
#include <algorithm>

class PerformancePloter {
public:
	static void openGLInit() {
		glClearColor(0.000, 0.110, 0.392, 0.0); // JMU Gold
		glColor3f(0.314, 0.314, 0.000); // JMU Purple
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glPointSize(2.0);
		gluOrtho2D(-1.0, 1.0, -1.0, 1.0);
	}
	static void drawString(double x, double y, const std::string& s) {
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glTranslatef(x, y, 0.0);
		glScaled(0.0003, 0.0003, 0.0003);
		for (size_t i = 0; i < s.size(); ++i)
			glutStrokeCharacter(GLUT_STROKE_ROMAN, s[i]);
		glPopMatrix();
	}
	static void drawString(double x, double y, const char *str,
			unsigned int maxValue) {
		if (str != NULL) {
			std::ostringstream os;
			os << str;
			os << ", ";
			os << "max value : " << maxValue;
			drawString(x, y, os.str());
		}
	}
	template<class T>
	static void drawBars(double xStart, double xEnd, double yStart, double yEnd,
			const std::vector<T>& values, T maxValue) {
		double xLen = xEnd - xStart;
		double yLen = yEnd - yStart;
		glBegin(GL_LINES);
		for (size_t i = 0; i < values.size(); ++i) {
			double x = xStart + i * xLen / values.size();
			double y = yStart + values[i] * yLen / maxValue;
			glVertex2d(x, yStart);
			glVertex2d(x, y);
		}
		glEnd();
	}
	template<class T>
	static void plotPoints(const char *title, const std::vector<T>& data,
			double maxValue, size_t count) {
		double xStart = -0.95;
		double yStart = -0.95;
		double xEnd = 0.95;
		double yEnd = 0.95;
		double xLen = xEnd - xStart;
		double yLen = yEnd - yStart;
		if (title != NULL)
			drawString(-0.8, 0.8, title, maxValue);
		glBegin(GL_POINTS);
		for (size_t i = 0; i < count; ++i) {
			double x = xStart + i * xLen / data.size();
			double y = yStart + yLen * data[i] / maxValue;
			//assert(data[i] <= maxValue);
			glVertex2f(x, y);
		}
		glEnd();
	}
	template<class T>
	static void display(const char *title, const std::vector<T>& data) {
		double xStart = -0.95;
		double yStart = -0.95;
		double xEnd = 0.95;
		double yEnd = 0.95;
		double xLen = xEnd - xStart;
		double yLen = yEnd - yStart;
		double maxValue = data[data.size() - 1];
		glClear(GL_COLOR_BUFFER_BIT);
		drawString(-0.8, 0.8, title, maxValue);
		// amortized cost.
		glColor3f(1.0, 0.0, 0.000);
		glBegin(GL_POINTS);
		for (size_t i = 0; i < data.size(); ++i) {
			double x = xStart + (i + 1) * xLen / data.size();
			double y = yStart + yLen * data[i] / ((i + 1) * maxValue);
			//assert(data[i] <= maxValue);
			glVertex2f(x, y);
		}
		glEnd();
		glColor3f(0.314, 0.314, 0.000); // JMU Purple
		glBegin(GL_LINES);
		glVertex2f(xStart, yStart);
		glVertex2f(xEnd, yStart);
		glVertex2f(xStart, yStart);
		glVertex2f(xStart, yEnd);
		glEnd();
		glBegin(GL_POINTS);
		for (size_t i = 0; i < data.size(); ++i) {
			double x = xStart + i * xLen / data.size();
			double y = yStart + yLen * data[i] / maxValue;
			//assert(data[i] <= maxValue);
			glVertex2f(x, y);
		}
		glEnd();
		glFlush();
	}
};
template<class U>
void doublingTest(unsigned int start, unsigned int end, U& u) {
	while (start < end) {
		u.problemSize(start);
		clock_t begin = clock();
		u();
		u.add(start, clock() - begin);
		start *= 2;
	}
}
template<class ForwardIterator>
void randUInts(ForwardIterator begin, ForwardIterator end,
		unsigned int maxValue) {
	std::random_device rd;
	std::uniform_int_distribution<decltype(maxValue)> generator(0, maxValue);
	std::generate(begin, end, [&] {return generator(rd);});
}
static unsigned int ceilLog2(unsigned long long x) {
	static const unsigned long long t[6] = { 0xFFFFFFFF00000000ull,
			0x00000000FFFF0000ull, 0x000000000000FF00ull, 0x00000000000000F0ull,
			0x000000000000000Cull, 0x0000000000000002ull };

	unsigned int y = (((x & (x - 1)) == 0) ? 0 : 1);
	unsigned int j = 32;
	unsigned int i;

	for (i = 0; i < 6; i++) {
		unsigned int k = (((x & t[i]) == 0) ? 0 : j);
		y += k;
		x >>= k;
		j >>= 1;
	}

	return y;
}
void randUInts(std::vector<unsigned int>& values, unsigned int count,
		unsigned int maxValue) {
	values.resize(count);
	randUInts(values.begin(), values.end(), maxValue);
}
#endif
