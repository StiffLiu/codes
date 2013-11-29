#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cassert>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
using namespace std;

/**
 * This interface comes from section 1.5,  "algorithms 4th Edition" book.
 * And many implementations in this file are copied from that section,
 * except that they are written in "C++".
 *
 *
 * A brief introduction of some concepts:
 *      site : An abstraction of entity, uniquely identified by an unsigned integer
 *              which is called an identifier or the name of the site.
 *      connected : A symmetric and transitive relation. That is
 *                   (1) if site a is connected to site b, then site b is connected to site a.
 *                   (2) if site a is connected to site b and site b is connected to site c, then site a is connected to site c
 *      component : Subset of all the sites, where every two sites in the subset are connected.
 *               uniquely identified by an integer which is called an identifier or the name of the component.
 *      It's trivial to show that two sites are connected if and only if they belong to the same component.
 *      We can partition the set of all sites into components, where every site must belong to one component and every two
 *      sites in a component are connected.
 */
class UF {
public:

	/**
	 * Make two sites connected.
	 * @param p the name of one site.
	 * @param q the name of another site.
	 *@return true if the two sets are not connected before the operation
	 */
	virtual bool connect(unsigned int p, unsigned int q) = 0;

	/**
	 * Find the component of a site
	 * @param p the integer identifier of the site
	 * @return the identifier of the component of the site
	 */
	virtual unsigned int find(unsigned int p) = 0;

	/**
	 * Determine whether two sites are connected.
	 * @param p the name of one site.
	 * @param q the name of another site.
	 * @return {@code true} if the two sites are connected, else {@code false}.
	 *
	 */
	virtual bool connected(int p, int q) {
		return find(p) == find(q);
	}

	/**
	 * @return the number of components.
	 */
	virtual int count() = 0;

	virtual ~UF() {
	}
};

/**
 * This class is based class.
 * Which represents n sites by integers from 0 to (n - 1).
 */
template<class SiteArray>
class UFIndexedT: public UF {
public:
	/**
	 * Representation of the connected relation between two sites.
	 * Where initially no two sites are connected.
	 * @param n number of sites.
	 */
	UFIndexedT(unsigned int n) :
			ids(n) {
		reset();
	}
	virtual void reset() {
		//since no two sites are connected, we have n components.
		unsigned int n = ids.length();
		c = n;
		for (unsigned int i = 0; i < n; ++i)
			ids[i] = i;

	}
	virtual int count() {
		return c;
	}
	SiteArray& getArray() {
		return ids;
	}
protected:
	//an array that stores some information about sites
	//the actual meaning of values stored in this array is determined by
	//the subclass.
	SiteArray ids;
	//number of components
	unsigned int c;
private:
	//Avoid assignment and copy construction.
	UFIndexedT& operator=(const UFIndexedT&);
	UFIndexedT(const UFIndexedT&);
};
template<class T>
class ArrayT {
	T *ar;
	unsigned int n;
	ArrayT(const ArrayT&);
	ArrayT& operator=(const ArrayT&);
public:
	ArrayT(unsigned int n) {
		ar = new T[n];
		this->n = n;
	}
	T& operator[](unsigned int index) {
		return ar[index];
	}
	const T& operator[](unsigned int index) const {
		return ar[index];
	}
	unsigned int length() {
		return n;
	}
	~ArrayT() {
		delete[] ar;
	}
};
typedef ArrayT<unsigned int> UIntArray;
class AccessCountedUIntArray: public UIntArray {
	mutable unsigned int accessCount;
	UIntArray& operator=(const UIntArray&);
public:
	AccessCountedUIntArray(unsigned int n) :
			UIntArray(n), accessCount(0) {
	}
	unsigned int getAccessCount() const {
		return accessCount;
	}
	void setAccessCount(unsigned int count) {
		accessCount = count;
	}
	void clearAccessCount() {
		setAccessCount(0);
	}
	unsigned int& operator[](unsigned int index) {
		++accessCount;
		return UIntArray::operator[](index);
	}
	unsigned int operator[](unsigned int index) const {
		++accessCount;
		return UIntArray::operator[](index);
	}
};
typedef UFIndexedT<AccessCountedUIntArray> UFIndexed;
/**
 * This implementation stores the component identifier of a site with name {@code i} at {@code ids[i]}
 * As the name suggests, retrieving the component of a site is fast.
 */
class QuickFindUF: public UFIndexed {
public:
	QuickFindUF(unsigned int n) :
			UFIndexed(n) {
	}
	//Obviously the cost of this function is O(1).
	virtual unsigned int find(unsigned int p) {
		return ids[p];
	}

	//the cost of this function is O(n)
	virtual bool connect(unsigned int p, unsigned int q) {
		//the component of site p.
		unsigned int cp = find(p);
		//the component of site q.
		unsigned int cq = find(q);
		unsigned int n = getArray().length();
		//if the two sites are not connected, then make them connected.
		if (cp != cq) {
			//If site p and q are connected, then every two sites in the component "cp"
			//and component "cq" are connected.
			for (unsigned int i = 0; i < n; ++i)
				if (ids[i] == cp)
					ids[i] = cq;
			--c;
			return true;
		}
		return false;
	}
};
/**
 * quick-union find.
 * Forest-of-trees representation.
 * 	if {@code ids[i] == i}, then {@code ids[i]} is the component name of site {@code i}.
 * 	if {@code ids[i] != i}, then {@code ids[i]} is the name another site(call the parent site)
 * 	which is in the same component with site {@code i}.
 *
 */
class QuickUnionUF: public UFIndexed {
public:
	QuickUnionUF(unsigned int n) :
			UFIndexed(n) {
	}
	//The cost is O(n) in the worst case.
	//The cost is O(lg(n)) in the average case.
	unsigned int find(unsigned int p) {
		while (ids[p] != p)
			p = ids[p];
		return p;
	}

	//The cost is O(n) in the worst case.
	//The cost is O(lg(n)) in the average case.
	bool connect(unsigned int p, unsigned int q) {
		unsigned int cp = find(p);
		unsigned int cq = find(q);
		if (cp == cq)
			return false;
		merge(cp, cq);
		--c;
		return true;
	}
	/**
	 * path compression union
	 */
	bool uniteWithCompress(unsigned int p, unsigned int q) {
		unsigned int cp = find(p);
		unsigned int cq = find(q);
		if (cp == cq)
			return false;
		if (merge(cp, cq)) {
			while (p != cp) {
				unsigned int next = ids[p];
				ids[p] = cq;
				p = next;
			}
		} else {
			while (q != cq) {
				unsigned int next = ids[q];
				ids[q] = cp;
				q = next;
			}
		}
		--c;
		return true;
	}
protected:
	/**
	 * merge two trees
	 * @param cp root of the first tree
	 * @param cq root of the second tree
	 * @return {@code true}, If the first tree become a child of the second tree.
	 *             {@code false}, otherwise.
	 *
	 */
	virtual bool merge(unsigned int cp, unsigned int cq) {
		ids[cp] = cq;
		return true;
	}
};
/**
 * weighted quick-union find.
 * We record the number of sites in each component.
 * when connect two sites,  the tree with less sites is added
 *      as a child of the root of the tree with more sites.
 *
 *the number of sites rooted at site {@code i} is recored in the array {@code sz}.
 */
class WeightedQuickUnionUF: public QuickUnionUF {
	void resetSizes() {
		unsigned int n = getArray().length();
		for (unsigned int i = 0; i < n; ++i)
			sz[i] = 1;
	}
public:
	WeightedQuickUnionUF(unsigned int n) :
			QuickUnionUF(n) {
		sz = new unsigned int[n];
		resetSizes();
	}
	~WeightedQuickUnionUF() {
		delete[] sz;
	}
	virtual void reset() {
		QuickUnionUF::reset();
		resetSizes();
	}
	//The cost is O(lg(n)).
	bool merge(unsigned int cp, unsigned int cq) {
		if (sz[cp] > sz[cq]) {
			ids[cq] = cp;
			sz[cp] += sz[cq];
			return false;
		} else {
			ids[cp] = cq;
			sz[cq] += sz[cp];
			return true;
		}
	}
private:
	unsigned int *sz;
};
class WeightedByHeightQuickUnionUF: public QuickUnionUF {
	void resetHeights() {
		unsigned int n = getArray().length();
		for (unsigned int i = 0; i < n; ++i)
			ht[i] = 0;
	}
public:
	WeightedByHeightQuickUnionUF(unsigned int n) :
			QuickUnionUF(n) {
		ht = new unsigned int[n];
		resetHeights();
	}
	~WeightedByHeightQuickUnionUF() {
		delete[] ht;
	}
	virtual void reset() {
		QuickUnionUF::reset();
		resetHeights();
	}
	//The cost is O(lg(n)).
	bool merge(unsigned int cp, unsigned int cq) {
		if (ht[cp] > ht[cq]) {
			ids[cq] = cp;
			return false;
		} else if (ht[cq] > ht[cp]) {
			ids[cp] = cq;
			return true;
		} else {
			//prefer the smaller identifier
			if (cp < cq) {
				ids[cq] = cp;
				ht[cq]++;
				return false;
			} else {
				ids[cp] = cq;
				ht[cp]++;
				return true;
			}
		}
	}
private:
	unsigned int *ht;
};

namespace {
struct Tester {
	vector<unsigned int> qfResult;
	vector<unsigned int> quResult, quResultPC;
	vector<unsigned int> wquResult, wquResultPC;
	vector<unsigned int> wbhquResult, wbhquResultPC;
public:
	void assertConnected(UF& uf, unsigned int *sites, unsigned int n) {
		for (unsigned int i = 0; i < n; ++i) {
			for (unsigned int j = i + 1; j < n; ++j)
				assert(uf.connected(sites[i], sites[j]));
		}
	}
	void assertNotConnected(UF& uf, unsigned int *sites1, unsigned int n1,
			unsigned int *sites2, unsigned int n2) {
		for (unsigned int i = 0; i < n1; ++i) {
			for (unsigned int j = i + 1; j < n2; ++j)
				assert(!uf.connected(sites1[i], sites2[j]));
		}
	}
	void validateAlgorithm(bool withPathCompression) {
#define ARSIZE(a) (sizeof (a) / sizeof (*(a)))
		unsigned int pairs[] = { 4, 3, 3, 8, 6, 5, 9, 4, 2, 1, 8, 9, 5, 0, 7, 2,
				6, 1, 1, 0, 6, 7 };
		unsigned int component1[] = { 1, 0, 2, 7, 5, 6 };
		unsigned int component2[] = { 3, 4, 8, 9 };
		unsigned int siteCount = 10;
		QuickFindUF qfUF(siteCount);
		QuickUnionUF quUF(siteCount);
		WeightedQuickUnionUF wquUF(siteCount);
		WeightedByHeightQuickUnionUF wbhquUF(siteCount);
		unsigned int n = ARSIZE(pairs);
		if (!withPathCompression)
			for (unsigned int i = 0; i < n; i += 2) {
				qfUF.connect(pairs[i], pairs[i + 1]);
				quUF.connect(pairs[i], pairs[i + 1]);
				wquUF.connect(pairs[i], pairs[i + 1]);
				wbhquUF.connect(pairs[i], pairs[i + 1]);
			}
		else
			for (unsigned int i = 0; i < n; i += 2) {
				qfUF.connect(pairs[i], pairs[i + 1]);
				quUF.uniteWithCompress(pairs[i], pairs[i + 1]);
				wquUF.uniteWithCompress(pairs[i], pairs[i + 1]);
				wbhquUF.uniteWithCompress(pairs[i], pairs[i + 1]);
			}
		assert(qfUF.count() == 2);
		assert(quUF.count() == 2);
		assert(wquUF.count() == 2);
		assert(wbhquUF.count() == 2);
		assertConnected(qfUF, component1, ARSIZE(component1));
		assertConnected(quUF, component1, ARSIZE(component1));
		assertConnected(wquUF, component1, ARSIZE(component1));
		assertConnected(wbhquUF, component1, ARSIZE(component1));
		assertConnected(qfUF, component2, ARSIZE(component2));
		assertConnected(quUF, component2, ARSIZE(component2));
		assertConnected(wquUF, component2, ARSIZE(component2));
		assertConnected(wbhquUF, component2, ARSIZE(component2));
		assertNotConnected(qfUF, component1, ARSIZE(component1), component2,
				ARSIZE(component2));
		assertNotConnected(quUF, component1, ARSIZE(component1), component2,
				ARSIZE(component2));
		assertNotConnected(wquUF, component1, ARSIZE(component1), component2,
				ARSIZE(component2));
		assertNotConnected(wbhquUF, component1, ARSIZE(component1), component2,
				ARSIZE(component2));
	}

	static void experimentOne(UFIndexed& uf, const vector<unsigned int>& pairs,
			vector<unsigned int>& arrayAccesses) {
		AccessCountedUIntArray& qfAr = uf.getArray();
		uf.reset();
		qfAr.clearAccessCount();
		arrayAccesses.resize(pairs.size() / 2);
		for (size_t i = 0; i < pairs.size(); i += 2) {
			uf.connect(pairs[i], pairs[i + 1]);
			arrayAccesses[i / 2] = qfAr.getAccessCount();
		}

	}
	static void experimentOne(QuickUnionUF& uf,
			const vector<unsigned int>& pairs,
			vector<unsigned int>& arrayAccesses,
			vector<unsigned int>& arrayAccessesPC) {
		experimentOne(uf, pairs, arrayAccesses);
		AccessCountedUIntArray& qfAr = uf.getArray();
		uf.reset();
		qfAr.clearAccessCount();
		arrayAccessesPC.resize(pairs.size() / 2);
		for (size_t i = 0; i < pairs.size(); i += 2) {
			uf.uniteWithCompress(pairs[i], pairs[i + 1]);
			arrayAccessesPC[i / 2] = qfAr.getAccessCount();
		}
	}
	void run(const char *dataFile) {
		vector<unsigned int> pairs;
		ifstream is(dataFile);
		while (is) {
			unsigned int val;
			is >> val;
			pairs.push_back(val);
		}
		run(pairs);
	}
	void run(const vector<unsigned int>& pairs) {
		qfResult.clear();
		quResult.clear();
		quResultPC.clear();
		wquResult.clear();
		wquResultPC.clear();
		wbhquResult.clear();
		wbhquResultPC.clear();

		if (pairs.size() < 2)
			return;
		unsigned int siteCount = *std::max_element(pairs.begin(), pairs.end())
				+ 1;
		QuickFindUF qfUF(siteCount);
		QuickUnionUF quUF(siteCount);
		WeightedQuickUnionUF wquUF(siteCount);
		WeightedByHeightQuickUnionUF wbhquUF(siteCount);
		experimentOne(qfUF, pairs, qfResult);
		experimentOne(quUF, pairs, quResult, quResultPC);
		experimentOne(wquUF, pairs, wquResult, wquResultPC);
		experimentOne(wbhquUF, pairs, wbhquResult, wbhquResultPC);
	}
};
#include <windows.h>
#include <GL/glut.h>

class PerformancePloter {
public:
	static void drawString(double x, double y, const char *str,
			unsigned int maxValue) {
		ostringstream os;
		if (str != NULL) {
			os << str;
			os << ", ";
		}
		os << "max value : " << maxValue;
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glTranslatef(x, y, 0.0);
		glScaled(0.0003, 0.0003, 0.0003);
		const string& s = os.str();
		/*while (*str) {
		 glutStrokeCharacter(GLUT_STROKE_ROMAN, *str);
		 ++str;
		 }*/
		for (size_t i = 0; i < s.size(); ++i)
			glutStrokeCharacter(GLUT_STROKE_ROMAN, s[i]);
		glPopMatrix();
	}
	template<class T>
	void display(const char *title, const vector<T>& data) {
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
	void init() {
		glClearColor(0.000, 0.110, 0.392, 0.0); // JMU Gold
		glColor3f(0.314, 0.314, 0.000); // JMU Purple
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glPointSize(2.0);
		gluOrtho2D(-1.0, 1.0, -1.0, 1.0);
	}
};
class UFPerformancePloter {
	Tester tester;
	PerformancePloter p;
	vector<const char *> titles;
	vector<vector<unsigned int>*> datas;
	unsigned int current;
	static UFPerformancePloter*instance;
	static void display() {
		instance->show();
	}
	static void keyboard(unsigned char ch, int x, int y) {
		instance->kbd(ch, x, y);
	}
	void randTest() {
		const unsigned int siteCount = 5000;
		unsigned int connectionCount = 10 + rand() % (siteCount * siteCount);
		vector<unsigned int> pairs;
		connectionCount = connectionCount / 2 * 2;
		for (unsigned int i = 0; i < connectionCount; ++i) {
			unsigned int val1 = rand() % siteCount;
			unsigned int val2 = 0;
			while ((val2 = rand() % siteCount) == val1)
				;
			pairs.push_back(val1);
			pairs.push_back(val2);
		}
		tester.run(pairs);
	}
public:
	void show() {
		if (current < titles.size())
			p.display(titles[current], *datas[current]);
		else
			cout << "won't show" << endl;
	}
	void kbd(unsigned char ch, int x, int y) {
		if (ch == 'n') {
			current = (current + 1) % titles.size();
			glutPostRedisplay();
		} else if (ch == 'r') {
			current = 0;
			randTest();
			glutPostRedisplay();
		}
	}
	int run(int argc, char *argv[]) {
		instance = this;
		current = 0;
		titles.clear();
		datas.clear();
		titles.push_back("quick find");
		datas.push_back(&tester.qfResult);
		titles.push_back("quick union");
		datas.push_back(&tester.quResult);
		titles.push_back("quick union with path compression");
		datas.push_back(&tester.quResultPC);
		titles.push_back("weighted quick union");
		datas.push_back(&tester.wquResult);
		titles.push_back("weighted quick union with path compression");
		datas.push_back(&tester.wquResultPC);
		titles.push_back("weighted quick union by height");
		datas.push_back(&tester.wbhquResult);
		titles.push_back(
				"weighted quick union by height with path compression");
		datas.push_back(&tester.wbhquResultPC);
		tester.validateAlgorithm(false);
		tester.validateAlgorithm(true);
		tester.run("F:\\source\\java\\algos4\\algs4-data\\mediumUF.txt");
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize(640, 480);
		glutInitWindowPosition(0, 0);
		glutCreateWindow("Test");
		glutDisplayFunc(display);
		glutKeyboardFunc(keyboard);
		p.init();
		glutMainLoop();
		return 0;
	}
};
UFPerformancePloter *UFPerformancePloter::instance = NULL;
/**
 * If we randomly choose two sites from a set of sites and connect them if they are not connected.
 * how many times should we choose before all sites are connected?
 * The hypothesis is that the number of choices  is ~ ½*N*ln N where N is the number of sites.
 * This class counts the number of random choices for a site set with at most n elements.
 * Then plot it on a graph.
 */
class ErdosRenyi {
	//"counts[i]" the number of random choices for a set with "i" sites.
	vector<unsigned int> counts;
	//measure the performance by number of array accesses;
	vector<unsigned int> arrayAccesses;
	//measure the performance by time;
	static ErdosRenyi *instance;
	static void display() {
		instance->show();
	}
	void show() {
		double xStart = -0.95;
		double yStart = -0.95;
		double xEnd = 0.95;
		double yEnd = 0.95;
		double xLen = xEnd - xStart;
		double yLen = yEnd - yStart;
		double maxValue = counts.size() * log(counts.size()) / 2 * 1.5;
		maxValue = arrayAccesses[arrayAccesses.size() - 1];
		glClear(GL_COLOR_BUFFER_BIT);
		glColor3f(0.314, 0.314, 0.000); // JMU Purple
		glBegin(GL_LINES);
		glVertex2f(xStart, yStart);
		glVertex2f(xEnd, yStart);
		glVertex2f(xStart, yStart);
		glVertex2f(xStart, yEnd);
		glEnd();
		glBegin(GL_LINE_STRIP);
		for (size_t i = 0; i < counts.size(); ++i) {
			double x = xStart + xLen * (i + 1) / counts.size();
			double y = yStart + yLen * counts[i] / maxValue;
			glVertex2f(x, y);
		}
		glEnd();
		glBegin(GL_LINE_STRIP);
		glColor3f(1.0, 0.0, 0.0); //Red
		for (size_t i = 0; i < counts.size(); ++i) {
			double x = xStart + xLen * (i + 1) / counts.size();
			double y = yStart + yLen * (i + 1) * log(i + 1) / (2 * maxValue);
			glVertex2f(x, y);
		}
		glEnd();
		glBegin(GL_LINE_STRIP);
		for (size_t i = 0; i < arrayAccesses.size(); ++i) {
			double x = xStart + xLen * (i + 1) / arrayAccesses.size();
			double y = yStart + yLen * arrayAccesses[i] / maxValue;
			glVertex2f(x, y);
		}
		glEnd();
		glColor3f(0.314, 0.314, 0.000); // JMU Purple
		glBegin(GL_LINE_STRIP);
		for (size_t i = 0; i < counts.size(); ++i) {
			double x = xStart + xLen * (i + 1) / counts.size();
			double y = yStart
					+ yLen * (i + 1) * log(i + 1) * log(i + 1) / (2 * maxValue);
			glVertex2f(x, y);
		}
		glEnd();
		glBegin(GL_LINE_STRIP);
		for (size_t i = 0; i < counts.size(); ++i) {
			double x = xStart + xLen * (i + 1) / counts.size();
			double y = yStart
					+ yLen * (i + 1)  * (i + 1) / (2 * maxValue);
			glVertex2f(x, y);
		}
		glEnd();
		glFlush();
	}
public:
	ErdosRenyi(unsigned int n) {
		//If there are zero or one site, no need to choose
		counts.push_back(0);
		counts.push_back(0);
		arrayAccesses.push_back(0);
		arrayAccesses.push_back(0);

		for (unsigned int i = 2; i < n; ++i) {
			WeightedQuickUnionUF wquUF(i);
			counts.push_back(0);
			unsigned int& count = counts.back();
			//random choose sites from a set of "i" sites,
			//utill all sites are connected.
			while (wquUF.count() > 1) {
				unsigned int val1 = rand() % i;
				unsigned int val2 = rand() % i;
				wquUF.connect(val1, val2);
				++count;
			}
			arrayAccesses.push_back(wquUF.getArray().getAccessCount());
		}
	}
	void init() {
		glClearColor(0.000, 0.110, 0.392, 0.0); // JMU Gold
		glColor3f(0.314, 0.314, 0.000); // JMU Purple
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glPointSize(2.0);
		gluOrtho2D(-1.0, 1.0, -1.0, 1.0);
	}
	int run(int argc, char *argv[]) {
		instance = this;
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize(640, 480);
		glutInitWindowPosition(0, 0);
		glutCreateWindow("Test");
		glutDisplayFunc(display);
		glutMainLoop();
		return 0;
	}
};
class RandGridGenerator {

};
ErdosRenyi* ErdosRenyi::instance = NULL;
}
int main(int argc, char *argv[]) {
	/*UFPerformancePloter ploter;
	 return ploter.run(argc, argv);*/
	ErdosRenyi erdosRenyi(10000);
	erdosRenyi.run(argc, argv);
}
