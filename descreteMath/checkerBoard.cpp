#include <cstdlib>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <vector>
namespace{
template<class T>
class Sliced2DArray{
	unsigned int startr, startc, lenc, lenr, r, c;
	T *vals;
public:
	Sliced2DArray(T * vals_, 
		unsigned int startr_, unsigned int startc_, 
		unsigned int lenr_, unsigned int lenc_, 
		unsigned int r_, unsigned int c_):
	vals(vals_), startr(startr_), startc(startc_),
		lenc(lenc_), lenr(lenr_), r(r_), c(c_){
			assert(startc + lenc <= c && startr + lenr <= r);
	}
	Sliced2DArray(const Sliced2DArray& ar,
		unsigned int startr_, unsigned int startc_, 
		unsigned int lenr_, unsigned int lenc_):
	vals(ar.vals), startr(ar.startr + startr_),
		startc(ar.startc + startc_), lenc(lenc_), 
		lenr(lenr_), r(ar.r), c(ar.c){
			assert(startr + lenr <= ar.startr + ar.lenr &&
				startc + lenc <= ar.startc + ar.lenc);
	}
	unsigned int cLen()const{
		return lenc;
	}
	unsigned int rLen()const{
		return lenr;
	}
	class Sliced1DArray{
		T *vals;
		unsigned int lenc;
		friend class Sliced2DArray<T>;
		Sliced1DArray(T* vals_, unsigned int lenc_)
			:vals(vals_), lenc(lenc_){
		}
	public:
		T& operator[](unsigned int c){
			assert(c < lenc);
			return vals[c];
		}
	};
	class ConstSliced1DArray{
		const T *vals;
		unsigned int lenc;
		friend class Sliced2DArray<T>;
		ConstSliced1DArray(const T* vals_, unsigned int lenc_)
			:vals(vals_), lenc(lenc_){
		}
	public:
		const T& operator[](unsigned int c) const{
			assert(c < lenc);
			return vals[c];
		}
	};
	Sliced1DArray operator[](unsigned int r){
		assert(r < lenr);
		return Sliced1DArray(vals + (r + startr) * c + startc, lenc);
	}
	ConstSliced1DArray operator[](unsigned int r)const{
		assert(r < lenr);
		return ConstSliced1DArray(vals + (r + startr) * c + startc, lenc);
	}
	friend std::ostream& operator<<(std::ostream&os, const Sliced2DArray& ar){
		for(int i = 0;i < ar.rLen();++ i){
			for(int j = 0;j < ar.cLen();++ j)
				os << std::setw(5) << ar[i][j];
			os << std::endl;
		}
		return os;
	}
};
void testSliced2DArray(int argc, char *argv[]){
	using namespace std;
	const int c = 30, r = 20;
	int tmp[r * c];
	int lenc = 10, lenr = 5, lenc1 = 8, lenr1 = 3;
	Sliced2DArray<int> ar(tmp, 1, 1, lenr, lenc, r, c);
	Sliced2DArray<int> ar1(ar, 2, 2, lenr1, lenc1);
	int current = 1;
	for(int i = 0;i < r;++ i)
		for(int j = 0;j < c;++ j)
			tmp[i * c + j] = (current ++);
	for(int i = 0;i < r;++ i){
		for(int j = 0;j < c;++ j)
			cout << setw(5) << tmp[i * c + j];
		cout << endl;
	}
	cout << "----------------" << endl << ar << "----------------" << endl <<  ar1 << endl;

}
template<class T>
class AlignedSliced2DArray{
	unsigned int startr, startc, pc, pt;
	T *vals;
public:
	AlignedSliced2DArray(T * vals_, 
		unsigned int startr_, unsigned int startc_, 
		unsigned int pc_, unsigned int pt_):
	vals(vals_), startr(startr_), startc(startc_),
		pc(pc_), pt(pt_){
			assert(startc + (1 << pc) <= (1 << pt) && startr + (1 << pc) <= (1 << pt));
	}
	AlignedSliced2DArray(const AlignedSliced2DArray& ar,
		unsigned int startr_, unsigned int startc_, 
		unsigned int pc_):
	vals(ar.vals), startr(ar.startr + startr_),
		startc(ar.startc + startc_), pc(pc_), pt(ar.pt){
			assert(startr + (1 << pc) <= ar.startr + (1 << ar.pc) &&
				startc + (1 << pc) <= ar.startc + (1 << ar.pc));
	}
	unsigned int len()const{
		return 1 << pc;
	}
	unsigned int count()const{
		return pc;
	}
	class Sliced1DArray{
		T *vals;
		unsigned int pc;
		friend class AlignedSliced2DArray<T>;
		Sliced1DArray(T* vals_, unsigned int pc_)
			:vals(vals_), pc(pc_){
		}
	public:
		T* operator&(){
			return vals;
		}
		T& operator[](unsigned int c){
			assert(c < (1 << pc));
			return vals[c];
		}
	};
	class ConstSliced1DArray{
		const T *vals;
		unsigned int pc;
		friend class AlignedSliced2DArray<T>;
		ConstSliced1DArray(const T* vals_, unsigned int pc_)
			:vals(vals_), pc(pc_){
		}
	public:
		const T* operator&(){
			return vals;
		}
		const T& operator[](unsigned int c) const{
			assert(c < (1 << pc));
			return vals[c];
		}
	};
	Sliced1DArray operator[](unsigned int r){
		assert(r < (1 << pc));
		return Sliced1DArray(vals + (r + startr) * (1 << pt) + startc, pc);
	}
	ConstSliced1DArray operator[](unsigned int r)const{
		assert(r < (1 << pc));
		return ConstSliced1DArray(vals + (r + startr) * (1 << pt) + startc, pc);
	}
	friend std::ostream& operator<<(std::ostream&os, const AlignedSliced2DArray& ar){
		for(int i = 0;i < ar.len();++ i){
			for(int j = 0;j < ar.len();++ j)
				os << std::setw(5) << ar[i][j];
			os << std::endl;
		}
		return os;
	}
};
class TileCheckerboard{
	mutable unsigned int r, c, n;
	mutable int *vals;
	std::vector<std::pair<int, int> > *seq;
	TileCheckerboard(const TileCheckerboard&);
	TileCheckerboard& operator=(const TileCheckerboard&);
	bool search(const int* val, int i) const{
		unsigned int index = (val - vals);
		unsigned int count = (1 << n);
		unsigned int row = index / count, col = index % count;
		if(row  > 0){
			if(vals[index - count] == i)
				return true;
			if(vals[index - count] == -1 && row > 1 && vals[index - count - count] == i)
				return true;
		}
		if(row < count - 1){
			if(vals[index + count] == i)
				return true;
			if(vals[index + count] == -1 && row < count - 2 && vals[index + count + count] == i)
				return true;
		}
		if(col > 0){
			if(vals[index - 1] == i)
				return true;
			if(vals[index - 1] == -1 && col > 1 && vals[index - 2] == i)
				return true;
		}
		if(col < count - 1){
			if(vals[index + 1] == i)
				return true;
			if(vals[index + 1] == -1 && col < count -2 && vals[index + 2] == i)
				return true;
		}
		return false;
	}
	void setIndex(int* val1, int* val2, int* val3) const{
		int i = 1;
		while(search(val1, i) || search(val2, i) || search(val3, i)) ++i;
		*val1 = *val2 = *val3 = i;
		//assert(i <= 3);
		if(seq != NULL){
			seq->push_back(std::make_pair(val1 - vals, i));
			seq->push_back(std::make_pair(val2 - vals, i));
			seq->push_back(std::make_pair(val3 - vals, i));
		}
	}
	void calculate(unsigned int r, unsigned int c, AlignedSliced2DArray<int>& ar) const{
		if(ar.len() == 2){
			if(r == 0){
				if(c == 0){
					setIndex(&ar[0][1], &ar[1][0], &ar[1][1]);
				}else if(c == 1){
					setIndex(&ar[0][0], &ar[1][0], &ar[1][1]);
				}else
					assert(false);
			}else if(r == 1){
				if(c == 0){
					setIndex(&ar[0][1], &ar[0][0], &ar[1][1]);
				}else if(c == 1){
					setIndex(&ar[0][1], &ar[0][0], &ar[1][0]);
				}else
					assert(false);
			}else
				assert(false);
			return;
		}
		unsigned int halfLen = (ar.len() / 2);
		AlignedSliced2DArray<int> lu(ar, 0, 0, ar.count() - 1);
		AlignedSliced2DArray<int> ru(ar, 0, halfLen, ar.count() - 1);
		AlignedSliced2DArray<int> ll(ar, halfLen, 0, ar.count() - 1);
		AlignedSliced2DArray<int> rl(ar, halfLen, halfLen, ar.count() - 1);
		if(r < halfLen){
			if(c < halfLen){
				setIndex(&ru[halfLen-1][0], &ll[0][halfLen-1], &rl[0][0]);
				calculate(r, c, lu);
				calculate(halfLen - 1, 0, ru);
				calculate(0 , halfLen - 1, ll);
				calculate(0, 0, rl);
			}else if(c < ar.len()){
				setIndex(&lu[halfLen-1][halfLen - 1], &ll[0][halfLen-1], &rl[0][0]);
				calculate(halfLen - 1, halfLen - 1, lu);
				calculate(r, c - halfLen, ru);
				calculate(0 , halfLen - 1, ll);
				calculate(0, 0, rl);
			}else
				assert(false);
		}else if(r < ar.len()){
			if(c < halfLen){
				setIndex(&lu[halfLen-1][halfLen - 1], &ru[halfLen - 1][0], &rl[0][0]);
				calculate(halfLen - 1, halfLen - 1, lu);
				calculate(halfLen - 1, 0, ru);
				calculate(r - halfLen , c, ll);
				calculate(0, 0, rl);
			}else if(c < ar.len()){
				setIndex(&lu[halfLen-1][halfLen - 1], &ru[halfLen - 1][0], &ll[0][halfLen-1]);
				calculate(halfLen - 1, halfLen - 1, lu);
				calculate(halfLen - 1, 0, ru);
				calculate(0 , halfLen - 1, ll);
				calculate(r - halfLen, c - halfLen, rl);
			}else
				assert(false);
		}else
			assert(false);
	}
	void begin() const{
		if(vals != NULL)
			return;
		if(vals == NULL)
			vals = new int[(1 << n) * (1 << n)];
		for(int i = 0;i < (1 << (2*n));++ i)
			vals[i] = 0;
		AlignedSliced2DArray<int> ar(vals, 0, 0, n, n);
		ar[r][c] = -1;
		calculate(r, c, ar);
	}
public:
	TileCheckerboard(unsigned int r_, unsigned int c_, 
		unsigned int n_):r(r_), c(c_), n(n_), seq(NULL){
			assert(r < (1 << n) && c < (1 << n));
			vals = NULL;	
	}
	void setSeq(std::vector<std::pair<int, int> >* seq){
		this->seq = seq;
		delete [] vals;
		vals = NULL;
		begin();
	}
	~TileCheckerboard(){
		delete []vals;
	}
	friend std::ostream& operator<<(std::ostream&os, TileCheckerboard& t){
		t.begin();
		for(int i = 0;i < (1 << t.n);++ i){
			for(int j = 0;j < (1 << t.n);++ j)
				os << std::setw(2) << t.vals[i * (1 << t.n) + j];
			os << std::endl;
		}
		return os;
	}

};
void testAlignedSliced2DArray(int argc, char *argv[]){
	using namespace std;
	const int p = 5;
	int tmp[(1 << p) * (1 << p)];
	int p1 = 4, p2 = 2;
	AlignedSliced2DArray<int> ar(tmp, 1, 1, p1, p);
	AlignedSliced2DArray<int> ar1(ar, 2, 2, p2);
	int current = 1;
	for(int i = 0;i < (1 << p);++ i)
		for(int j = 0;j < (1 << p);++ j)
			tmp[i * (1 << p) + j] = (current ++);
	for(int i = 0;i < (1 << p);++ i){
		for(int j = 0;j < (1 << p);++ j)
			cout << setw(5) << tmp[i * (1 << p) + j];
		cout << endl;
	}
	cout << "----------------" << endl << ar << "----------------" << endl <<  ar1 << endl;

}
void getResult(unsigned int r, unsigned int c, unsigned int n, std::vector<std::pair<int, int> >& ret){
	if(n == 0){
		assert(false);
		return;
	}
	assert(r < (1 << n) && c < (1 << n));
	r %= (1 << n);
	c %= (1 << n);
	TileCheckerboard t(r, c, n);
	t.setSeq(&ret);
}

std::vector<std::pair<int, int> > values;
int count, position, current = 0, width = 800, height = 800, division = 4;;
bool started = false;
clock_t lastTime;
#include <glut.h>
void init(){
	glClearColor(1.0, 1.0, 1.0, 1.0);
}
void initData(){
	srand(time(0));
	current = 0;
	count = 3 + rand() % 4;
	position = (rand() % (1 << count) )*(1 << count);// * (1 << count));
	if(rand () % 2 == 0)
		position += (1 << count) - 1;
	values.clear();
	getResult(position / (1 << count), position % (1 << count), count, values);
}
void display(){
	int len = std::min(width, height) - 10;
	int sLen = len / (1 << count);
	if(sLen < 10)
		sLen = 10;
	glClear(GL_COLOR_BUFFER_BIT);
	glColor4f(0.0, 0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	for(int i = 0;i <= (1 << count);++ i){
		glVertex2i(i * sLen, 0);
		glVertex2i(i * sLen, (1 << count) * sLen);
	}
	for(int i = 0;i <= (1 << count);++ i){
		glVertex2i(0, i * sLen);
		glVertex2i((1 << count) * sLen, i * sLen);
	}
	glEnd();
	glBegin(GL_QUADS);
	{
		int row = position / (1 << count);
		int col = position % (1 << count);
		glVertex2i(row * sLen, col * sLen);
		glVertex2i(row * sLen + sLen, col * sLen);
		glVertex2i(row * sLen + sLen, col * sLen + sLen);
		glVertex2i(row * sLen, col * sLen + sLen);
	}
	for(int i = 0;i < current;++ i){
		for(int j = 3 * i;j < 3 * i + 3;++ j){
			int row = values[j].first / (1 << count);
			int col = values[j].first % (1 << count);
			switch(values[j].second){
			case 1:
				glColor4f(1.0, 0.0, 0.0, 1.0);break;
			case 2:
				glColor4f(0.0, 1.0, 0.0, 1.0);break;
			case 3:
				glColor4f(0.0, 0.0, 1.0, 1.0);break;
			case 4:
				glColor4f(1.0, 1.0, 0.0, 1.0);break;
			default:
				assert(false);
			}
			glVertex2i(row * sLen, col * sLen);
			glVertex2i(row * sLen + sLen, col * sLen);
			glVertex2i(row * sLen + sLen, col * sLen + sLen);
			glVertex2i(row * sLen, col * sLen + sLen);
		}
	}
	glEnd();
	glutSwapBuffers();
}
void reshape(int w, int h){
	width = w;
	height = h;
	int len = std::min(width, height);
	glViewport(10, 10, len - 10, len - 10);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, len, 0, len);
}
void animate(){
	if(!started || (clock() - lastTime)  < CLOCKS_PER_SEC / division)
		return;
	lastTime = clock();
	assert(values.size() % 3 == 0);
	if(current < values.size() / 3)
		++current;
	glutPostRedisplay();
}
void mouseFunc(int mouse, int state, int x, int y){
	started = true;
	initData();
}
void keyboardFunc(unsigned char ch, int x, int y){
	switch(ch){
		case 'i':++division;break;
		case 'd':--division;if(division <= 0)division = 1;break;
	}
}
}
int testCheckerBoard(int argc, char *argv[]){
	using namespace std;
	initData();
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(width, height);
	glutCreateWindow("Checkerboard Tiling");
	glutDisplayFunc(display);
	glutIdleFunc(animate);
	glutReshapeFunc(reshape);
	glutMouseFunc(mouseFunc);
	glutKeyboardFunc(keyboardFunc);
	init();
	glutMainLoop();
	return 0;
}