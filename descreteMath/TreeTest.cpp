#include <vector>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <queue>
#include <algorithm>
#include <fstream>
#include <iterator>
#include <cmath>
#include <set>
#include <cassert>
using namespace std;

#define ARSIZE(ar) (sizeof (ar) / sizeof (*(ar)))
class BitArray {
	typedef unsigned char Byte;
	unsigned int bitCount;
	Byte *bits;
	unsigned int capacity;
	static const int ByteBytes = (sizeof(Byte));
	static const int ByteBits = ByteBytes * 8;
	void assure(int size) {
		unsigned int bytesCount = size / ByteBits;
		if (size % ByteBits != 0)
			++bytesCount;
		if (bytesCount > capacity) {
			capacity = bytesCount;
			if (bits != NULL)
				bits = (Byte*) realloc(bits, capacity * ByteBytes);
			else
				bits = (Byte*) malloc(capacity * ByteBytes);
		}
	}
public:
	~BitArray() {
		free(bits);
	}
	BitArray() :
		bitCount(0), bits(NULL), capacity(0) {
	}
	BitArray(unsigned int bitCount) :
		bitCount(bitCount), bits(NULL) {
		capacity = 0;
		assure(bitCount);
	}
	BitArray(const BitArray& bitArray) :
		bitCount(0), bits(NULL), capacity(0) {
		append(bitArray);
	}
	BitArray& operator=(const BitArray& bitArray) {
		if (this != &bitArray) {
			bitCount = 0;
			append(bitArray);
		}
		return *this;
	}
	void set(unsigned int index, bool value) {
		unsigned int i = index / ByteBits;
		unsigned int j = index % ByteBits;
		if (value)
			bits[i] |= (1 << j);
		else
			bits[i] &= ~(1 << j);
	}
	bool operator[](unsigned int index) const {
		unsigned int i = index / ByteBits;
		unsigned int j = index % ByteBits;
		return (bits[i] & (1 << j)) != 0;
	}
	bool at(unsigned int index) const {
		if (index >= bitCount)
			return false;
		return this->operator[](index);
	}
	unsigned int getCapacity() const {
		return capacity;
	}
	unsigned int getBitCount() const {
		return bitCount;
	}
	void remove(unsigned int startPos, unsigned int count) {
		if (count == 0)
			return;
		unsigned int endPos = startPos + count;
		if (count % ByteBits == 0 && endPos < bitCount) {
			while (startPos % ByteBits != 0) {
				set(startPos, this->operator[](endPos));
				++startPos;
				++endPos;
				--count;
			}
			unsigned int bytesCount = count / ByteBits;
			if (bytesCount > 0) {
				unsigned int bitsCount = bytesCount * ByteBits;
				memmove(&bits[startPos / ByteBits], &bits[endPos / ByteBits],
						bytesCount * ByteBytes);
				startPos += bitsCount;
				endPos += bitsCount;
			}
		}
		while (endPos < bitCount) {
			set(startPos, this->operator[](endPos));
			++startPos;
			++endPos;
		}
		bitCount = startPos;
	}
	void append(void *bs, unsigned int bCount) {
		assure(bitCount + bCount);
		unsigned int startPos = bitCount;
		unsigned int startPos0 = 0;
		bitCount += bCount;
		if (startPos % ByteBits == 0) {
			unsigned int bytesCount = bCount / ByteBits;
			if (bytesCount > 0) {
				unsigned int bitsCount = bytesCount * ByteBits;
				memmove(&bits[startPos / ByteBits], bs, bytesCount * ByteBytes);
				startPos += bitsCount;
				startPos0 += bitsCount;
			}
		}
		while (startPos0 != bCount) {
			set(startPos, get(startPos0, bs));
			++startPos;
			++startPos0;
		}
	}
	void append(const BitArray& val) {
		append(val.bits, val.bitCount);
	}
	void fill(bool val) {
		if (bitCount != 0 && bits == NULL)
			assure(bitCount);
		if (bits != NULL)
			memset(bits, (val ? -1 : 0), capacity);
	}
	unsigned int countBitOne() const {
		const int intBits = (sizeof(unsigned int) * 8);
		unsigned int intsCount = bitCount / intBits;
		unsigned int c = 0;
		unsigned int *vals = (unsigned int*) bits;
		for (unsigned int i = 0; i < intsCount; ++i) {
			unsigned int val = vals[i];
			while (val != 0) {
				val &= (val - 1);
				++c;
			}
		}
		unsigned int startPos = intsCount * intBits;
		while (startPos < bitCount) {
			if (this->operator[](startPos))
				++c;
			++startPos;
		}
		return c;
	}
	unsigned int countBitZero() const {
		return bitCount - countBitOne();
	}
	bool operator ==(const BitArray& val) const {
		if (bitCount != val.bitCount)
			return false;

		const int intBits = (sizeof(unsigned int) * 8);
		unsigned int intsCount = bitCount / intBits;
		unsigned int *vals0 = (unsigned int*) bits, *vals1 =
				(unsigned int*) val.bits;
		for (unsigned int i = 0; i < intsCount; ++i)
			if (vals0[i] != vals1[i])
				return false;
		unsigned int startPos = intsCount * intBits;
		while (startPos < bitCount) {
			if (this->operator[](startPos) != val[startPos])
				return false;
			++startPos;
		}
		return true;
	}
	bool operator<(const BitArray& val) const {
		unsigned int bCount = min(bitCount, val.bitCount);
		if (bCount == 0) {
			return bitCount < val.bitCount;
		}
		const int intBits = (sizeof(unsigned int) * 8);
		unsigned int intsCount = bCount / intBits;
		unsigned int *vals0 = (unsigned int*) bits, *vals1 =
				(unsigned int*) val.bits;

		unsigned int i = 0;
		for (; i < intsCount; ++i)
			if (vals0[i] != vals1[i])
				break;
		unsigned int startPos = i * intBits;
		while (startPos < bitCount) {
			bool b1 = this->operator[](startPos), b2 = val[startPos];
			if (b1 != b2)
				return !b1;
			++startPos;
		}
		return bitCount < val.bitCount;
	}
	void swap(BitArray& val) {
		if (this == &val)
			return;
		std::swap(bitCount, val.bitCount);
		std::swap(bits, val.bits);
	}
	void random() {
		for (unsigned int i = 0; i < bitCount; ++i)
			set(i, rand() % 2 == 0);
	}
	static bool get(unsigned int index, void *bits) {
		unsigned int i = index / ByteBits;
		unsigned int j = index % ByteBits;
		return (((Byte*) bits)[i] & (1 << j)) != 0;
	}
	friend ostream& operator<<(ostream& os, const BitArray& val) {
		unsigned int bitCount = val.getBitCount();
		for (unsigned int i = 0; i < bitCount; ++i) {
			os << (val[i] ? '1' : '0');
		}
		return os;
	}
};
class BitArray2D {
	unsigned int row;
	unsigned int col;
	BitArray bitArray;
	friend class Row;
	friend class Col;
	void outputRow(unsigned int rowNo, ostream& os) const {
		unsigned int startIndex = rowNo * col;
		for (unsigned int j = 0; j < col; ++j)
			os << (bitArray[startIndex + j] ? '1' : '0');
	}
public:
	class Row {
		const BitArray2D& array;
		unsigned int rowNo;
		friend class BitArray2D;
		Row(const BitArray2D& array, unsigned int rowNo) :
			array(array), rowNo(rowNo) {
		}
	public:
		bool operator[](unsigned int i) const {
			return array.at(rowNo, i);
		}
		bool at(unsigned int i) const {
			return array.at(rowNo, i);
		}
		inline ostream& output(ostream& os) const {
			array.outputRow(rowNo, os);
			return os;
		}
		inline friend ostream& operator<<(ostream& os, const Row& val) {
			return val.output(os);
		}
	};
	class Col {
		const BitArray2D& array;
		unsigned int colNo;
		friend class BitArray2D;
		Col(const BitArray2D& array, unsigned int colNo) :
			array(array), colNo(colNo) {
		}
	public:
		bool operator[](unsigned int i) const {
			return array.get(i, colNo);
		}
		bool at(unsigned int i) const {
			return array.at(i, colNo);
		}
	};
	BitArray2D(unsigned int row, unsigned int col) :
		row(row), col(col), bitArray(row * col) {
	}
	unsigned int getRow() const {
		return row;
	}
	unsigned int getCol() const {
		return col;
	}
	Row operator[](unsigned int i) const {
		return Row(*this, i);
	}
	bool get(unsigned int i, unsigned int j) const {
		return bitArray[i * col + j];
	}
	bool at(unsigned int i, unsigned int j) const {
		if (i >= row)
			return false;
		if (j >= col)
			return false;
		return this->get(i, j);
	}
	unsigned int countBitOne() const {
		return bitArray.countBitOne();
	}
	unsigned int countBitZero() const {
		return bitArray.countBitZero();
	}
	inline void set(unsigned int i, unsigned int j, bool val) {
		bitArray.set(i * col + j, val);
	}
	inline void fill(bool val) {
		bitArray.fill(val);
	}
	bool operator==(const BitArray2D& val) const {
		return row == val.row && col == val.col && bitArray == val.bitArray;
	}
	void swap(BitArray2D& val) {
		std::swap(row, val.row);
		std::swap(col, val.col);
		bitArray.swap(val.bitArray);
	}
	friend ostream& operator<<(ostream& os, const BitArray2D& val) {
		for (unsigned int i = 0; i < val.row; ++i) {
			val.outputRow(i, os);
			os << endl;
		}
		return os;
	}
};
template<class T>
struct BinTreeNode{
	T p, l, r;
	BinTreeNode(const T& p, const T& l, const T& r) : p(p), l(l), r(r){
	}
};
template<class T, class Comparator=less<T> >
class BinSearchTree{
	typedef BinTreeNode<unsigned int> TreeNode;
	vector<TreeNode> tree;
	vector<T> datas;
	Comparator comparator;
	void output(unsigned int level, unsigned int node, ostream& os)const{
		for(unsigned int i = 0;i < level;++ i)
			os << ' ';
		os << datas[node] << '\n';
		if(tree[node].l != static_cast<unsigned int>(-1))
			output(level + 1, tree[node].l, os);
		if(tree[node].r != static_cast<unsigned int>(-1))
			output(level + 1, tree[node].r, os);

	}
	int compare(const T& item1, const T& item2){
		if(comparator(item1, item2))
			return -1;
		if(comparator(item2, item1))
			return 1;
		return 0;
	}
public:
	bool add(const T& item){
		if(tree.empty()){
			tree.push_back(TreeNode(0, -1, -1));
			datas.push_back(item);
			return true;
		}
		T node = 0;
		while(true){
			int ret = compare(datas[node], item);
			if(ret == 0)
				return false;
			else if(ret < 0){
				if(tree[node].r == static_cast<unsigned int>(-1)){
					tree[node].r = tree.size();
					tree.push_back(TreeNode(node, -1, -1));
					datas.push_back(item);
					return true;
				}
				node = tree[node].r;
			}else{
				if(tree[node].l == static_cast<unsigned int>(-1)){
					tree[node].l = tree.size();
					tree.push_back(TreeNode(node, -1, -1));
					datas.push_back(item);
					return true;
				}
				node = tree[node].l;
			}
		}
		return false;
	}
	friend ostream& operator<<(ostream& os, const BinSearchTree tree){
		if(!tree.tree.empty())
			tree.output(0, 0, os);
		return os;
	}
	BinSearchTree(const T*vals, unsigned int n){
		for(unsigned int i = 0;i < n;++ i)
			add(vals[i]);
	}
};
class AdjacencyMatrix {
	BitArray2D matrix;

public:
	AdjacencyMatrix(unsigned int vertices) :
		matrix(vertices, vertices) {
		matrix.fill(false);
	}
	bool isAdajacent(unsigned int v1, unsigned int v2) const {
		return matrix.get(v1, v2);
	}
	bool at(unsigned int v1, unsigned int v2) const {
		return matrix.at(v1, v2);
	}
	void set(unsigned int v1, unsigned int v2, bool val) {
		matrix.set(v1, v2, val);
	}
	unsigned int getDegree() const {
		return matrix.countBitOne();
	}
	bool getNb(unsigned int v, unsigned int startIndex, unsigned int& nb) const {
		unsigned int col = matrix.getCol();
		for (unsigned int j = startIndex; j < col; ++j)
			if (matrix.get(v, j)) {
				nb = j;
				return true;
			}
		return false;
	}
	bool getFirstNeighbour(unsigned int v, unsigned int& nb) const {
		return getNb(v, 0, nb);
	}
	bool getNextNeighbour(unsigned int v, unsigned int n, unsigned int& nb) const {
		return getNb(v, n + 1, nb);
	}
	unsigned int getDegreeOfVertex(unsigned int v) const {
		unsigned int col = matrix.getCol();
		unsigned int c = 0;
		for (unsigned int j = 0; j < col; ++j)
			if (matrix.get(v, j))
				++c;
		return c;
	}
	unsigned int getVertexCount() const {
		return matrix.getRow();
	}
	friend ostream& operator<<(ostream& os, const AdjacencyMatrix& val) {
		return os << val.matrix;
	}
};
class AdjacencyList {
protected:
	unsigned int vertexCount;
	unsigned int *adjacentVertex;
	unsigned int *vertexDegrees;
	void set(const AdjacencyList& val) {
		if (this != &val) {
			delete[] adjacentVertex;
			delete[] vertexDegrees;

			unsigned int degree = val.getDegree();
			vertexCount = val.vertexCount;
			adjacentVertex = new unsigned int[degree];
			vertexDegrees = new unsigned int[vertexCount + 1];
			memcpy(adjacentVertex, val.adjacentVertex, degree
					* sizeof(unsigned int));
			memcpy(vertexDegrees, val.vertexDegrees, (vertexCount + 1)
					* sizeof(unsigned int));
		}
	}
public:
	AdjacencyList() :
		vertexCount(0), adjacentVertex(NULL), vertexDegrees(NULL) {

	}
	AdjacencyList(const AdjacencyMatrix& matrix) {
		adjacentVertex = vertexDegrees = NULL;
		set(matrix);
	}
	~AdjacencyList() {
		delete[] adjacentVertex;
		delete[] vertexDegrees;
	}
	void set(const AdjacencyMatrix& matrix) {
		unsigned int sum = 0;
		vertexCount = matrix.getVertexCount();
		vertexDegrees = new unsigned int[vertexCount + 1];
		vertexDegrees[0] = 0;
		for (unsigned int i = 0; i < vertexCount; ++i) {
			unsigned int degree = matrix.getDegreeOfVertex(i);
			vertexDegrees[i + 1] = vertexDegrees[i] + degree;
			sum += degree;
		}
		adjacentVertex = new unsigned int[sum];

		unsigned int start = 0;
		for (unsigned int i = 0; i < vertexCount; ++i) {
			unsigned int nb = i;
			if (matrix.getFirstNeighbour(i, nb)) {
				adjacentVertex[start] = nb;
				++start;
				while (matrix.getNextNeighbour(i, nb, nb)) {
					adjacentVertex[start] = nb;
					++start;
				}
			}
		}
	}
	unsigned int getDegree() const {
		return vertexDegrees[vertexCount];
	}
	unsigned int getDegreeOfVertex(unsigned int v) const {
		return vertexDegrees[v + 1] - vertexDegrees[v];
	}
	unsigned int* getNb(unsigned int v) {
		if (v > vertexCount)
			return NULL;
		return &adjacentVertex[vertexDegrees[v]];
	}
	bool isAdjacent(unsigned int v1, unsigned int v2) {
		int end = vertexDegrees[v1] + vertexDegrees[v1 + 1];
		for (int i = vertexDegrees[v1]; i < end; ++i)
			if (adjacentVertex[i] == v2)
				return true;
		return false;
	}
	void set(unsigned int vertexCount, unsigned int *edges,
			unsigned int edgeCount, bool isEdgeOrdered = false) {
		delete[] vertexDegrees;
		delete[] adjacentVertex;

		this->vertexCount = vertexCount;
		vector<vector<unsigned int> > adjVerts;
		adjVerts.resize(vertexCount);
		for (unsigned int i = 0; i < edgeCount; ++i) {
			unsigned int index = i << 1;
			unsigned int v1 = edges[index];
			unsigned int v2 = edges[index + 1];
			if (v1 >= vertexCount || v2 >= vertexCount)
				continue;
			vector<unsigned int>& adj1 = adjVerts[v1];
			if (find(adj1.begin(), adj1.end(), v2) == adj1.end())
				adj1.push_back(v2);
			vector<unsigned int>& adj2 = adjVerts[v2];
			if (find(adj2.begin(), adj2.end(), v1) == adj2.end())
				adj2.push_back(v1);
		}

		unsigned int sum = 0;
		vertexDegrees = new unsigned int[vertexCount + 1];
		vertexDegrees[0] = 0;
		for (unsigned int i = 0; i < vertexCount; ++i) {
			unsigned int degree = adjVerts[i].size();
			vertexDegrees[i + 1] = vertexDegrees[i] + degree;
			sum += degree;
		}
		adjacentVertex = new unsigned int[sum];
		unsigned int start = 0;
		for (unsigned int i = 0; i < vertexCount; ++i) {
			vector<unsigned int>& adj = adjVerts[i];
			unsigned int degree = adj.size();
			for (unsigned int j = 0; j < degree; ++j, ++start)
				adjacentVertex[start] = adj[j];
		}
	}
	AdjacencyList(unsigned int vertexCount, unsigned int *edges,
			unsigned int edgeCount, bool isEdgeOrdered = false) {
		adjacentVertex = NULL;
		vertexDegrees = NULL;
		set(vertexCount, edges, edgeCount, isEdgeOrdered);
	}
	AdjacencyList(const AdjacencyList& val) {
		vertexCount = 0;
		adjacentVertex = vertexDegrees = NULL;
		set(val);
	}
	AdjacencyList& operator=(const AdjacencyList& val) {
		set(val);
		return *this;
	}
	unsigned int outputOneVertex(ostream& os, unsigned int v) const {
		unsigned int end = vertexDegrees[v + 1];
		unsigned int count = 0;
		for (unsigned int i = vertexDegrees[v]; i < end; ++i, ++count)
			os << adjacentVertex[i] << ' ';
		return count;
	}
	friend ostream& operator<<(ostream& os, const AdjacencyList& adjList) {
		for (unsigned int i = 0; i < adjList.vertexCount; ++i) {
			if (adjList.getDegreeOfVertex(i) > 0) {
				os << i << " : ";
				if (adjList.outputOneVertex(os, i) > 0)
					os << endl;
			}
		}
		return os;
	}
};
class UndirectedGraph {
	AdjacencyMatrix matrix;
	bool colorGraph(unsigned int colors, unsigned int v, unsigned int* coloring) {
		unsigned int i = 0;
		for (; i < colors; ++i) {
			unsigned int nb = 0;
			while (matrix.getNb(v, nb, nb)) {
				if (coloring[nb] == i)
					goto nextLoop;
				++nb;
			}
			coloring[v] = i;
			nb = 0;
			while (matrix.getNb(v, nb, nb)) {
				if (coloring[nb] == colors && !colorGraph(colors, nb, coloring))
					goto nextLoop;
				++nb;
			}
			return true;
			nextLoop: ;
		}
		return false;
	}
public:
	UndirectedGraph(unsigned int vertexCount = 0) :
		matrix(vertexCount) {
	}
	bool isAdjacent(unsigned int v1, unsigned int v2) const {
		return matrix.isAdajacent(v1, v2);
	}
	unsigned int getVertexCount() const {
		return matrix.getVertexCount();
	}
	unsigned int getEdgeCount() const {
		return matrix.getDegree() / 2;
	}
	void set(unsigned int i, unsigned int j, bool val) {
		if (i == j)
			return;
		matrix.set(i, j, val);
		matrix.set(j, i, val);
	}
	bool isConnected() const {
		unsigned int vertexCount = getVertexCount();
		vector<bool> isUsed(vertexCount, false);
		vector<unsigned int> vertices;
		if (isUsed.empty())
			return true;
		unsigned int count = 1;
		vertices.push_back(0);
		isUsed[0] = true;
		for (unsigned int i = 0; i < count; ++i) {
			unsigned int v = vertices[i];
			unsigned int nb = 0;
			while (matrix.getNb(v, nb, nb)) {
				if (!isUsed[nb]) {
					++count;
					vertices.push_back(nb);
					isUsed[nb] = true;
				}
				++nb;
			}
		}
		return vertices.size() == vertexCount;

	}
	bool isTree() const {
		return isConnected() && (getEdgeCount() == getVertexCount() - 1);
	}
	void dfsTree(unsigned int root) {
		vector<unsigned int> vertices;
		dfsTree(root, vertices);
	}
	void dfsTree(unsigned int root, vector<unsigned int>& vertices) {
		vector<bool> isUsed;
		dfsTree(root, isUsed, vertices);
	}
	void bfsTree(unsigned int root, vector<bool>& isUsed,
			vector<unsigned int>& vertices) {
		unsigned int i = vertices.size();
		unsigned int count = i + 1;
		unsigned int j = 0;
		vector<unsigned int> parent;
		vertices.push_back(root);
		parent.push_back(root);
		isUsed[root] = true;
		for (; i < count; ++i, ++j) {
			unsigned int v = vertices[i];
			unsigned int nb = 0;
			while (matrix.getNb(v, nb, nb)) {
				if (!isUsed[nb]) {
					++count;
					vertices.push_back(nb);
					parent.push_back(v);
					isUsed[nb] = true;
				} else if (nb != parent[j]) {
					set(v, nb, false);
				}
				++nb;
			}
		}
	}
	void bfsTree(unsigned int root) {
		vector<unsigned int> vertices;
		bfsTree(root, vertices);
	}
	void bfsTree(unsigned int root, vector<unsigned int>& vertices) {
		vector<bool> isUsed;
		bfsTree(root, isUsed, vertices);
	}
	void dfsTree(unsigned int root, vector<bool>& isUsed,
			vector<unsigned int>& vertices) {
		int i = vertices.size();
		int startIndex = i;
		unsigned int count = i + 1;
		unsigned int lastNb = 0;
		vector<unsigned int> parent;
		vertices.push_back(root);
		parent.push_back(-1);
		isUsed[root] = true;
		while (i >= startIndex) {
			unsigned int v = vertices[i];
			bool b = true;
			while (matrix.getNb(v, lastNb, lastNb)) {
				if (!isUsed[lastNb]) {
					parent.push_back(i);
					i = count;
					++count;
					vertices.push_back(lastNb);
					isUsed[lastNb] = true;
					lastNb = 0;
					b = false;
					break;
				}
				if (lastNb != vertices[parent[i - startIndex]])
					set(v, lastNb, false);
				lastNb += 1;
			}
			if (b) {
				lastNb = v + 1;
				i = parent[i - startIndex];
			}
		}
	}
	unsigned int dfsForest(unsigned int root = 0) {
		unsigned int vertexCount = getVertexCount();
		vector<bool> isUsed(vertexCount, false);
		vector<unsigned int> vertices;
		if (isUsed.empty())
			return 0;
		if (root > vertexCount)
			root = 0;
		unsigned int count = 1;
		dfsTree(root, isUsed, vertices);
		for (unsigned int i = 0; i < vertexCount; ++i)
			if (!isUsed[i]) {
				dfsTree(i, isUsed, vertices);
				++count;
			}
		return count;
	}
	unsigned int bfsForest(unsigned int root = 0) {
		unsigned int vertexCount = getVertexCount();
		vector<bool> isUsed(vertexCount, false);
		vector<unsigned int> vertices;
		if (isUsed.empty())
			return 0;
		unsigned int count = 1;
		if (root > vertexCount)
			root = 0;
		bfsTree(root, isUsed, vertices);
		for (unsigned int i = 0; i < vertexCount; ++i)
			if (!isUsed[i]) {
				bfsTree(i, isUsed, vertices);
				++count;
			}
		return count;
	}
	bool colorGraph(unsigned int colors, unsigned int* coloring = NULL) {
		unsigned int vertexCount = getVertexCount();
		if (colors == 0 && vertexCount != 0)
			return false;
		unsigned int *tmpColoring = NULL;
		if (coloring == NULL)
			tmpColoring = coloring = new unsigned int[vertexCount];
		for (unsigned int i = 0; i < vertexCount; ++i)
			coloring[i] = colors;
		for (unsigned int i = 0; i < vertexCount; ++i)
			if (coloring[i] == colors && !colorGraph(colors, i, coloring))
				return false;
		delete[] tmpColoring;
		return true;
	}
	unsigned int connectedComponents() {
		unsigned int vertexCount = getVertexCount();
		vector<bool> isUsed(vertexCount, false);
		vector<unsigned int> vertices;
		if (isUsed.empty())
			return 0;
		unsigned int count = 0;
		unsigned int componentCount = 0;
		while (count != vertexCount) {
			for (unsigned int i = 0; i < vertexCount; ++i)
				if (!isUsed[i]) {
					vertices.push_back(i);
					isUsed[i] = true;
					++count;
					break;
				}
			for (unsigned int i = count - 1; i < count; ++i) {
				unsigned int v = vertices[i];
				unsigned int nb = 0;
				while (matrix.getNb(v, nb, nb)) {
					if (!isUsed[nb]) {
						++count;
						vertices.push_back(nb);
						isUsed[nb] = true;
					}
					++nb;
				}
			}
			++componentCount;
		}
		return componentCount;
	}
	friend ostream& operator<<(ostream& os, const UndirectedGraph& graph) {
		return os << graph.matrix << endl;
	}
	AdjacencyList getAdjList() {
		return AdjacencyList(matrix);
	}
};
struct Node {
	unsigned int frequency;
	unsigned int height;
	unsigned int count;
	int l, r;
	Node(unsigned int f = 0, unsigned int h = 0, unsigned int cnt = 0, int l =
			-1, int r = -1) :
		frequency(f), height(h), count(cnt), l(l), r(r) {
	}
	bool operator<(const Node& node) const {
		if (frequency < node.frequency)
			return true;
		if (frequency == node.frequency) {
			if (height < node.height)
				return true;
			if (height == node.height)
				return count < node.count;
		}
		return false;
	}
	bool operator>(const Node& node) const {
		if (frequency > node.frequency)
			return true;
		if (frequency == node.frequency) {
			if (height > node.height)
				return true;
			if (height == node.height)
				return count > node.count;
		}
		return false;
	}
};
class Maze {
	BitArray2D data;
	void build(unsigned int start, double ratio = 0.6){
		unsigned int row = data.getRow();
		unsigned int col = data.getCol();
		if(row <= 2 || col <= 2)
			return;
		vector<unsigned int> usedCells;

		usedCells.push_back(start);
		unsigned int rowDirect[4];
		unsigned int colDirect[4];
		unsigned int count = 0;
		unsigned int total = round(row * col * ratio);
		data.set(start / col, start % col, true);
		while(count < total && !usedCells.empty()){
			unsigned int r = usedCells.back() / col;
			unsigned int c = usedCells.back() % col;
			unsigned int dc = 0;
			if(r >= 2 && !data.get(r - 2, c) && !data.get(r - 1, c)
					&& (c <= 0 || !data.get(r - 1, c - 1)) && (c >= col - 1 || !data.get(r - 1, c + 1))){
				rowDirect[dc] = r - 1;
				colDirect[dc] = c;
				++ dc;
			}
			if(r <= row - 3 && !data.get(r + 2, c) && !data.get(r + 1, c)
					&& (c <= 0 || !data.get(r + 1, c - 1)) && (c >= col - 1 || !data.get(r + 1, c + 1))){
				rowDirect[dc] = r + 1;
				colDirect[dc] = c;
				++ dc;
			}
			if(c >= 2 && !data.get(r, c - 2) && !data.get(r, c - 1)
					&& (r <= 0 || !data.get(r - 1, c - 1)) && (r >= row - 1 || !data.get(r + 1, c - 1))){
				rowDirect[dc] = r;
				colDirect[dc] = c - 1;
				++ dc;
			}
			if(c <= col - 3 && !data.get(r, c + 2) && !data.get(r, c + 1)
					&& (r <= 0 || !data.get(r - 1, c + 1)) && (r >= row - 1 || !data.get(r + 1, c + 1))){
				rowDirect[dc] = r;
				colDirect[dc] = c + 1;
				++ dc;
			}
			if(dc > 0){
				unsigned int index = rand() % dc;
				data.set(rowDirect[index], colDirect[index], true);
				usedCells.push_back(rowDirect[index] * col + colDirect[index]);
				++ count;
			}else{
				usedCells.pop_back();
			}
		}

	}
public:
	Maze(unsigned int row = 0, unsigned int col = 0) :
		data(row, col) {
		data.fill(false);
	}
	void init(unsigned int row, unsigned int col){
		data = BitArray2D(row, col);
		data.fill(false);
	}
	bool isWall(unsigned int i, unsigned int j){
		return !data.get(i, j);
	}
	unsigned getWallCount(){
		return data.countBitZero();
	}
	unsigned int getRow(){
		return data.getRow();
	}
	unsigned int getCol(){
		return data.getCol();
	}
	friend ostream& operator<<(ostream& os, const Maze& maze) {
		unsigned int row = maze.data.getRow();
		unsigned int col = maze.data.getCol();
		for (unsigned int i = 0; i < row; ++i) {
			for (unsigned int j = 0; j < col; ++j)
				os << (maze.data.get(i, j) ? ' ' : '@');
			if (i != row - 1)
				os << endl;
		}
		return os;
	}
	//not the shortest path, but there are algorithms to find shortest path.
	bool findPath(set<unsigned int>& path, unsigned int start, unsigned int end){
		vector<unsigned int> usedCells;
		vector<char> lastDirection;
		unsigned int row = getRow();
		unsigned int col = getCol();
		assert(data.get(start / col ,start % col));
		assert(data.get(end / col, end % col));
		usedCells.push_back(start);
		path.insert(start);
		lastDirection.push_back(-1);
		while(!usedCells.empty() && usedCells.back() != end){
			unsigned int r = usedCells.back() / col;
			unsigned int c = usedCells.back() % col;
			for(int i = lastDirection.back() + 1;i < 4;++ i){
				unsigned int tmp = -1;
				switch(i){
				case 0:
					if(r > 0 && data.get(r - 1, c)){
						tmp = usedCells.back() - col;
					}
					break;
				case 1:
					if(c < col - 1 && data.get(r, c + 1)){
						tmp = usedCells.back() + 1;
					}
					break;
				case 2:
					if(c > 0 && data.get(r, c - 1)){
						tmp = usedCells.back() - 1;
					}
					break;
				case 3:
					if(r < row - 1 && data.get(r + 1, c)){
						tmp = usedCells.back() + col;
					}
				}
				if(tmp != static_cast<unsigned int>(-1)){
					unsigned int cnt = usedCells.size();
					if(cnt <= 1 || path.find(tmp) == path.end()){
						lastDirection.back() = i;
						lastDirection.push_back(-1);
						usedCells.push_back(tmp);
						path.insert(tmp);
						goto endLoop;
					}
				}
			}
			path.erase(usedCells.back());
			usedCells.pop_back();
			lastDirection.pop_back();
			endLoop:;
		}
		if(!usedCells.empty()){
			set<unsigned int > tmp(usedCells.begin(), usedCells.end());
			path.swap(tmp);
			return true;
		}
		return false;
	}
	friend void outputWithPath(ostream& os, const Maze& maze, set<unsigned int >& path){
		unsigned int row = maze.data.getRow();
		unsigned int col = maze.data.getCol();
		for (unsigned int i = 0; i < row; ++i) {
			for (unsigned int j = 0; j < col; ++j){
				if(path.find(i * col + j) != path.end())
					cout << "@";
				else
					os << (maze.data.get(i, j) ? ' ' : '*');
			}
				if (i != row - 1)
					os << endl;
		}
	}
	void random() {
		unsigned int row = data.getRow();
		unsigned int col = data.getCol();
		srand(time(0));
		data.fill(false);
		build(0, 0.58);
		for(unsigned int i = 1;i < row - 1;++ i)
			for(unsigned int j = 1;j < col - 1;++ j)
				if(isWall(i,j)){
					unsigned hWall = 0, vWall = 0;
					if(!isWall(i - 1, j))
						++vWall;
					if(!isWall(i + 1, j))
						++vWall;
					if(!isWall(i, j - 1))
						++hWall;
					if(!isWall(i, j + 1))
						++hWall;
					if(vWall * hWall == 0 && vWall + hWall != 0 && rand () % 3 == 0)
						data.set(i, j, true);
				}
	}

};
struct NodeComparator {
	vector<Node>& nodes;
	NodeComparator(vector<Node>& nodes) :
		nodes(nodes) {
	}
	bool operator()(int i, int j) const {
		return nodes[i] > nodes[j];
	}
};
struct Token {
	int type;
	int data;
	Token(int type, int data) :
		type(type), data(data) {
	}
};
struct Lexer {
	enum TokenType {
		OPERATOR, NUMBER, END
	};
	vector<Token> tokens;
	vector<double> value;
	vector<Token>::iterator current;
	Token getToken() {
		if (current == tokens.end())
			return Token(END, -1);
		Token tmp = *current;
		++current;
		return tmp;
	}
	bool toTree(vector<BinTreeNode<unsigned int> >& tree){
		vector<Token>::iterator begin = tokens.begin(), end = tokens.end();
		unsigned int current = static_cast<unsigned int>(-1);
		/*while(begin != end && tree[current]){
			if(tree.empty()){

			}else{
			}
			++begin;
		}*/
		return begin == end;
	}
	bool getValue(Token token, double& val) const {
		if (token.type == NUMBER) {
			val = value[token.data];
			return true;
		}
		return false;
	}
	friend ostream& operator<<(ostream& os, const Lexer& lexer) {
		vector<Token>::const_iterator begin = lexer.tokens.begin(), end =
				lexer.tokens.end();
		while (begin != end) {
			Token token = *begin;
			if (token.type == Lexer::OPERATOR)
				os << (char) token.data << ' ';
			else if (token.type == Lexer::NUMBER)
				os << lexer.value[token.data] << ' ';
			++begin;
		}
		return os;
	}
};
class ThresholdGate {
	vector<double> weights;
	double threshold;
public:
	ThresholdGate(const vector<double>& weights, double threshold) :
		weights(weights), threshold(threshold) {
	}
	ThresholdGate(const double *weights, unsigned int n, double threshold) :
		weights(weights, weights + n), threshold(threshold) {
	}
	bool getValue(const BitArray& input) const {
		double sum = 0;
		unsigned int count = min(input.getBitCount(), weights.size());
		for (unsigned int i = 0; i < count; ++i)
			if (input[i])
				sum += weights[i];
		return sum > threshold;
	}
	friend ostream& operator<<(ostream& os, const ThresholdGate& gate) {
		unsigned int count = gate.weights.size();
		os << "weights : ";
		for (unsigned int i = 0; i < count; ++i)
			cout << gate.weights[i] << ' ';
		os << endl;
		return os << "threshold : " << gate.threshold;
	}
};
struct Expr {
public:
	static void randomPrefix(Lexer& lexer) {
		lexer.tokens.clear();
		lexer.value.clear();
		randomPrefix0(lexer);
		lexer.current = lexer.tokens.begin();
	}
	static void randomPostfix(Lexer& lexer) {
		lexer.tokens.clear();
		lexer.value.clear();
		randomPostfix0(lexer);
		lexer.current = lexer.tokens.begin();

	}
	static bool getValue(int op, double val1, double val2, double& val) {
		switch (op) {
		case '+':
			val = val1 + val2;
			return true;
		case '-':
			val = val1 - val2;
			return true;
		case '*':
			val = val1 * val2;
			return true;
		case '/':
			if (val2 == 0.0)
				return false;
			val = val1 / val2;
			return true;
		}
		return false;
	}
	static bool prefix(Lexer& lexer, double& ret) {
		vector<Token> s;
		Token token = lexer.getToken();
		while (token.type != Lexer::END) {
			s.push_back(token);
			size_t count = s.size();
			while (count >= 3 && s[count - 1].type == Lexer::NUMBER && s[count
					- 2].type == Lexer::NUMBER) {
				if (s[count - 3].type != Lexer::OPERATOR)
					return false;
				double val;
				double val1;
				if (!lexer.getValue(s[count - 2], val) || !lexer.getValue(
						s[count - 1], val1) || !getValue(s[count - 3].data,
						val, val1, val))
					return false;
				s.pop_back();
				s.pop_back();
				s.back() = Token(Lexer::NUMBER, lexer.value.size());
				lexer.value.push_back(val);
				count -= 2;
			}
			token = lexer.getToken();
		}
		if (s.size() != 1 || s[0].type != Lexer::NUMBER || !lexer.getValue(
				s[0], ret))
			return false;
		return true;
	}
	static bool postfix(Lexer& lexer, double& ret) {
		vector<Token> s;
		Token token = lexer.getToken();
		while (token.type != Lexer::END) {
			if (token.type == Lexer::OPERATOR) {
				size_t count = s.size();
				if (count < 2)
					return false;
				if (s[count - 1].type != Lexer::NUMBER || s[count - 2].type
						!= Lexer::NUMBER)
					return false;

				double val;
				double val1;
				if (!lexer.getValue(s[count - 2], val) || !lexer.getValue(
						s[count - 1], val1) || !getValue(token.data, val, val1,
						val))
					return false;

				s.pop_back();
				s.back() = Token(Lexer::NUMBER, lexer.value.size());
				lexer.value.push_back(val);

			} else
				s.push_back(token);
			token = lexer.getToken();
		}
		if (s.size() != 1 || s[0].type != Lexer::NUMBER || !lexer.getValue(
				s[0], ret))
			return false;
		return true;
	}
private:
	static int randOp() {
		static char op[4] = { '+', '-', '*', '/' };
		return op[rand() % 4];
	}
	static double randDouble() {
		return rand() % 10000 / 1000.0;
	}
	static void randomPrefix0(Lexer& lexer, unsigned int maxSize = 10) {
		if (lexer.tokens.empty()) {
			lexer.tokens.push_back(Token(Lexer::NUMBER, lexer.value.size()));
			lexer.value.push_back(randDouble());
		}
		while (lexer.tokens.size() < maxSize) {
			Token op = Token(Lexer::OPERATOR, randOp());
			vector<Token> right;
			vector<Token> newTokens;
			unsigned int tmpSize = maxSize - lexer.tokens.size();
			lexer.tokens.swap(right);
			randomPrefix0(lexer, tmpSize);
			newTokens.push_back(op);
			if (rand() % 2 == 0)
				lexer.tokens.swap(right);
			newTokens.insert(newTokens.end(), lexer.tokens.begin(),
					lexer.tokens.end());
			newTokens.insert(newTokens.end(), right.begin(), right.end());
			lexer.tokens.swap(newTokens);
		}
	}
	static void randomPostfix0(Lexer& lexer, unsigned int maxSize = 10) {
		if (lexer.tokens.empty()) {
			lexer.tokens.push_back(Token(Lexer::NUMBER, lexer.value.size()));
			lexer.value.push_back(randDouble());
		}
		while (lexer.tokens.size() < maxSize) {
			Token op = Token(Lexer::OPERATOR, randOp());
			vector<Token> right;
			vector<Token> newTokens;
			unsigned int tmpSize = maxSize - lexer.tokens.size();
			lexer.tokens.swap(right);
			randomPostfix0(lexer, tmpSize);
			if (rand() % 2 == 0)
				lexer.tokens.swap(right);
			newTokens.insert(newTokens.end(), lexer.tokens.begin(),
					lexer.tokens.end());
			newTokens.insert(newTokens.end(), right.begin(), right.end());
			newTokens.push_back(op);
			lexer.tokens.swap(newTokens);
		}
	}
};
class BooleanFunction {
	vector<BitArray> function;
	unsigned int varNum;
	void nandRep(const BitArray& bitArray, unsigned int startIndex, unsigned int count, string& ret){
		if(count <= 1){
			if(count == 1){
				char ch = (char)('a' + startIndex);
				ret.push_back(ch);
				if(!bitArray[startIndex]){
					ret.push_back('|');
					ret.push_back((char)('a' + startIndex));
				}
			}
			return;
		}
		unsigned int c1 = count / 2;
		unsigned int c2 = count - c1;
		string ret1;
		string ret2;
		nandRep(bitArray, startIndex, c1, ret1);
		nandRep(bitArray, startIndex + c1, c2, ret2);
		ret.push_back('(');	ret.push_back('(');
		ret += ret1;
		ret.push_back(')');	ret.push_back('|');
		ret.push_back('(');
		ret += ret2;
		ret.push_back(')');	ret.push_back(')');
		ret.push_back('|');
		ret.push_back('(');	ret.push_back('(');
		ret += ret1;
		ret.push_back(')');	ret.push_back('|');
		ret.push_back('(');
		ret += ret2;
		ret.push_back(')');	ret.push_back(')');
	}
	void norRep(const BitArray& bitArray, unsigned int startIndex, unsigned int count, string& ret){
		if(count <= 1){
			if(count == 1){
				char ch = (char)('a' + startIndex);
				ret.push_back(ch);
				if(!bitArray[startIndex]){
					ret.push_back('!');
					ret.push_back((char)('a' + startIndex));
				}
			}
			return;
		}
		unsigned int c1 = count / 2;
		unsigned int c2 = count - c1;
		string ret1;
		string ret2;
		norRep(bitArray, startIndex, c1, ret1);
		norRep(bitArray, startIndex + c1, c2, ret2);
		ret.push_back('(');	ret.push_back('(');
		ret += ret1;
		ret.push_back(')');	ret.push_back('!');
		ret.push_back('(');
		ret += ret1;
		ret.push_back(')');	ret.push_back(')');
		ret.push_back('!');
		ret.push_back('(');	ret.push_back('(');
		ret += ret2;
		ret.push_back(')');	ret.push_back('!');
		ret.push_back('(');
		ret += ret2;
		ret.push_back(')');	ret.push_back(')');
	}
	void nandRep(unsigned int startIndex, unsigned int count, string& ret) {
		if (count <= 1) {
			if(count == 1)
				nandRep(function[startIndex], 0, varNum, ret);
			return;
		}
		unsigned int c1 = count / 2;
		unsigned int c2 = count - c1;
		string ret1;
		string ret2;
		nandRep(startIndex, c1, ret1);
		nandRep(startIndex + c1, c2, ret2);
		ret.push_back('(');	ret.push_back('(');
		ret += ret1;
		ret.push_back(')');	ret.push_back('|');
		ret.push_back('(');
		ret += ret1;
		ret.push_back(')');	ret.push_back(')');
		ret.push_back('|');
		ret.push_back('(');	ret.push_back('(');
		ret += ret2;
		ret.push_back(')');	ret.push_back('|');
		ret.push_back('(');
		ret += ret2;
		ret.push_back(')');	ret.push_back(')');

	}
	void norRep(unsigned int startIndex, unsigned int count, string& ret) {
		if (count <= 1) {
			if(count == 1)
				norRep(function[startIndex], 0, varNum, ret);
			return;
		}
		unsigned int c1 = count / 2;
		unsigned int c2 = count - c1;
		string ret1;
		string ret2;
		norRep(startIndex, c1, ret1);
		norRep(startIndex + c1, c2, ret2);
		ret.push_back('(');	ret.push_back('(');
		ret += ret1;
		ret.push_back(')');	ret.push_back('!');
		ret.push_back('(');
		ret += ret2;
		ret.push_back(')');	ret.push_back(')');
		ret.push_back('!');
		ret.push_back('(');	ret.push_back('(');
		ret += ret1;
		ret.push_back(')');	ret.push_back('!');
		ret.push_back('(');
		ret += ret2;
		ret.push_back(')');	ret.push_back(')');
	}
public:
	BooleanFunction(unsigned int varNum) :
		varNum(varNum) {
	}
	void random(unsigned int count) {
		for (unsigned int i = 0; i < count; ++i) {
			BitArray tmp(varNum);
			tmp.random();
			while (find(function.begin(), function.end(), tmp)
					!= function.end())
				tmp.random();
			function.push_back(BitArray());
			function.back().swap(tmp);
		}
	}
	void sumOfProductRep(ostream& os) {
		for (unsigned int i = 0; i < function.size(); ++i) {
			if (i != 0)
				os << " + ";
			BitArray& tmp = function[i];
			for (unsigned int j = 0; j < varNum; ++j) {
				if (!tmp[j])
					os << '~';
				os << (char) ('a' + j);
			}
		}
	}
	void andNegRep(ostream& os) {
		os << "~(";
		for (unsigned int i = 0; i < function.size(); ++i) {
			os << "~(";
			BitArray& tmp = function[i];
			for (unsigned int j = 0; j < varNum; ++j) {
				if (!tmp[j])
					os << '~';
				os << (char) ('a' + j);
			}
			os << ")";
		}
		os << ")";
	}
	void orNegRep(ostream& os) {
		for (unsigned int i = 0; i < function.size(); ++i) {
			if (i != 0)
				os << " + ";
			BitArray& tmp = function[i];
			os << "~(";
			for (unsigned int j = 0; j < varNum; ++j) {
				if (j != 0)
					os << '+';
				if (tmp[j])
					os << '~';
				os << (char) ('a' + j);
			}
			os << ")";
			if ((i + 1) % 6 == 0)
				os << endl;
		}
	}
	void nandRep(ostream& os) {
		string tmp;
		this->nandRep(0, function.size(), tmp);
		os << tmp;
	}
	void addFunction(const BitArray& tmp){
		if(tmp.getBitCount() == varNum){
			if(find(function.begin(), function.end(), tmp) == function.end())
				function.push_back(tmp);
			return;
		}
		unsigned int bCount = min(tmp.getBitCount(), varNum);
		BitArray b(varNum);
		b.fill(false);
		for(unsigned int i = 0;i < bCount;++ i)
			b.set(i, tmp[i]);
		if(find(function.begin(), function.end(), b) != function.end())
			return;
		function.push_back(BitArray());
		function.back().swap(b);
	}
	void norRep(ostream& os) {
		string tmp;
		this->norRep(0, function.size(), tmp);
		os << tmp;
	}
};
void huffamEncoding(unsigned int *fs, unsigned int n, vector<BitArray>& result,
		vector<Node>& nodes) {
	vector<unsigned int> allIndices;
	for (unsigned int i = 0; i < n; ++i) {
		nodes.push_back(Node(fs[i], 0, 1));
		allIndices.push_back(i);
	}

	NodeComparator tmp(nodes);
	std::priority_queue<int, vector<int> , NodeComparator> all(
			allIndices.begin(), allIndices.end(), tmp);
	while (all.size() > 1) {
		int l = all.top();
		all.pop();
		int r = all.top();
		all.pop();
		Node& n1 = nodes[l];
		Node& n2 = nodes[r];
		int index = nodes.size();
		nodes.push_back(Node(n1.frequency + n2.frequency, max(n1.height,
				n2.height) + 1, n1.count + n2.count, l, r));
		all.push(index);
	}
	unsigned int root = all.top();

	BitArray bitOne(1);
	BitArray bitZero(1);
	BitArray tmpTmp;
	bitOne.set(0, true);
	bitZero.set(0, false);

	cout << tmpTmp << endl;
	queue<unsigned int> nodeQueue;
	result.resize(nodes.size());
	nodeQueue.push(root);
	while (!nodeQueue.empty()) {
		unsigned int n = nodeQueue.front();
		Node& node = nodes[n];
		nodeQueue.pop();
		if (node.l != -1 && node.r != -1) {
			result[node.l] = result[n];
			result[node.r] = result[n];
			result[node.l].append(bitOne);
			result[node.r].append(bitZero);
			nodeQueue.push(node.l);
			nodeQueue.push(node.r);
		}
	}
}
bool nQueen(int n, vector<int>* ns = NULL) {
	if (n <= 0)
		return false;
	vector<int> nums;
	int count = 0;
	int startValue = 0;
	nums.resize(n);
	if(ns != NULL){
		for(size_t i = 0;i < ns->size();++ i){
			if((*ns)[i] >= n)
				break;
			nums[i] = (*ns)[i];
			++ count;
		}
		if(count > 0){
			--count;
			startValue = nums[count] + 1;
		}
	}
	while (count < n && count >= 0) {
		int i = startValue;
		for (; i < n; ++i) {
			int j = 0;
			for (; j < count; ++j)
				if (nums[j] == i || abs(count - j) == abs(nums[j] - i))
					break;
			if (j == count)
				break;
		}
		if (i != n) {
			nums[count] = i;
			startValue = 0;
			++count;
		} else {
			--count;
			if (count < 0)
				break;
			startValue = nums[count] + 1;
		}
	}
	if(count == n){
		ns->swap(nums);
		return true;
	}
	return false;
}
template<class T>
void randomPermute(T *vals, unsigned int n){
	while(n > 1){
		unsigned int i = rand() % n;
		if(i > 0){
			T tmp = vals[0];
			vals[0] = vals[i];
			vals[i] = tmp;
		}
		++vals;
		--n;
	}
}
template<class T>
void badRandomPermute(T *vals, unsigned int n){
	if(n <= 1)
		return;
	for(unsigned int i = 0;i < n;++ i){
		unsigned int j = rand() % n;
		if(i != j){
			T tmp = vals[j];
			vals[j] = vals[i];
			vals[i] = tmp;
		}
	}
}
template<class T>
bool findSum(const T *vals, unsigned int n, T sum, vector<const T*>& ret) {
	for (unsigned int i = 0; i < n; ++i)
		if (vals[i] == sum || findSum(vals + i + 1, n - i - 1, sum - vals[i],
				ret)) {
			ret.push_back(vals + i);
			return true;
		}
	return false;
}
int test0(int argc, char *argv[]) {
	vector<Node> nodes;
	vector<BitArray> result;
	ifstream in("/home/tmlg/src/descreteMath/TreeTest.cpp", ios::binary);
	string content((istreambuf_iterator<char> (in)),
			(istreambuf_iterator<char> ()));
	const int size = 256;
	unsigned int fs[size] = { 0 };
	for (unsigned int i = 0; i < content.size(); ++i)
		++fs[static_cast<unsigned int>(content[i])];
	vector<unsigned int> realFs;
	vector<char> chars;
	for (int i = 0; i < size; ++i)
		if (fs[i] != 0) {
			realFs.push_back(fs[i]);
			chars.push_back((char) i);
		}
	huffamEncoding(&realFs[0], realFs.size(), result, nodes);

	unsigned int totalBits = 0;
	for (unsigned int i = 0; i < chars.size(); ++i) {
		switch (chars[i]) {
		case '\n':
			cout << "\\n";
			break;
		case '\r':
			cout << "\\r";
			break;
		case '\t':
			cout << "\\t";
			break;
		case ' ':
			cout << "\\b";
			break;
		default:
			cout << chars[i];
			break;
		}
		cout << '\t' << realFs[i] << '\t' << result[i] << endl;
		totalBits += (realFs[i] * result[i].getBitCount());
	}
	int totalBytes = totalBits / 8;
	if (totalBits % 8 != 0)
		++totalBytes;
	cout << "bytes after compression : " << totalBytes << endl;
	cout << "bytes before compression : " << content.size() << endl;
	return 0;
}
int test1(int argc, char *argv[]) {
	for(int i = 0;i < 12;++ i){
		int count = 0;
		vector<int> tmp;
		while(nQueen(i, &tmp)){
			for(size_t i = 0;i < tmp.size();++ i)
				cout << tmp[i] << ' ';
			cout << endl;
			++ count;
		}
		cout << "num : " << i << ", count : " << count << endl << endl;
	}
	return 0;
}
int test2(int argc, char *argv[]) {
	BitArray bitOne(1);
	BitArray bits(100);
	for (unsigned int i = 0; i < bits.getBitCount(); ++i)
		bits.set(i, i % 2 == 0);
	cout << bits << endl;
	bits.set(59, true);
	bits.remove(80, 30);
	cout << bits << endl;
	bits.remove(60, 20);
	cout << bits << endl;
	bits.remove(1, 1);
	cout << bits << endl;
	bits.remove(10, 20);
	cout << bits << endl;
	cout << "number of bit one : " << bits.countBitOne() << endl;
	BitArray tmp = bits;
	assert(tmp == bits);
	tmp.set(0, false);
	cout << "clear the first bit : " << tmp << endl;
	assert(tmp < bits);
	tmp.set(0, true);
	tmp.set(tmp.getBitCount() - 1, false);
	cout << "clear the last bit : " << tmp << endl;
	assert(tmp < bits);
	cout << "remove the last 10 bits : " << tmp << endl;
	tmp.remove(tmp.getBitCount() - 10, 10);
	assert(tmp < bits);
	cout << "remove the last 10 bits : " << tmp << endl;

	return 0;
}
int test3(int argc, char *argv[]) {
	int vals[] = { 27, 24, 11, 19, 14, 8 };
	int count = ARSIZE(vals);
	int sum = 19;
	vector<const int*> ret;
	copy(vals, vals + count, (ostream_iterator<int> (cout, " ")));
	cout << endl;
	if (findSum(vals, count, sum, ret)) {
		for (int i = ret.size() - 1; i >= 0; --i)
			cout << *ret[i] << ' ';
		cout << endl;
	} else {
		cout << "does not exist!" << endl;
	}
	return 0;
}
int test4(int argc, char *argv[]) {
	const unsigned int verticesCount = 15;
	unsigned int edges[] = { 0, 1, 0, 2, 0, 3, 0, 4, 1, 5, 1, 6, 2, 7, 3, 8, 8,
			9, 4, 10, 4, 11, 4, 12, 4, 13, 4, 14, };
	const unsigned int edgeCount = ARSIZE(edges) / 2;
	UndirectedGraph graph(verticesCount);
	AdjacencyList adjList0(verticesCount, edges, edgeCount);
	cout << graph;
	cout << "adjacency list with constructor call:\n" << adjList0;
	for (unsigned int i = 0; i < edgeCount; ++i)
		graph.set(edges[2 * i], edges[2 * i + 1], true);
	assert(graph.getVertexCount() == verticesCount);
	for (unsigned int i = 0; i < edgeCount; ++i) {
		assert(graph.isAdjacent(edges[2 * i], edges[2 * i + 1]));
		assert(adjList0.isAdjacent(edges[2 * i], edges[2 * i + 1]));
	}
	cout << "adjacency list with conversion from adjacency matrix:\n"
			<< graph.getAdjList();
	cout << graph;
	cout << graph.getEdgeCount() << endl;
	assert(graph.getEdgeCount() == edgeCount);
	graph.set(2, 3, true);
	assert(graph.isConnected());
	assert(!graph.isTree());
	graph.set(2, 3, false);
	assert(graph.isConnected());
	assert(graph.connectedComponents() == 1);
	assert(graph.isTree());
	assert(graph.colorGraph(2));
	assert(graph.getEdgeCount() == graph.getVertexCount() - 1);
	graph.set(0, 1, false);
	assert(!graph.isConnected());
	assert(graph.connectedComponents() == 2);
	assert(graph.getEdgeCount() == edgeCount - 1);
	assert(!graph.isConnected());
	graph.set(0, 1, true);
	return 0;
}
int test5(int argc, char *argv[]) {
	const unsigned int verticesCount = 13;
	unsigned int edges[] = { 0, 2, 0, 4, 0, 5, 0, 7, 1, 2, 1, 4, 3, 2, 3, 5, 6,
			4, 6, 7, 6, 9, 8, 7, 8, 5, 8, 10, 10, 5, 10, 12, 9, 7, 9, 11, };
	const unsigned int edgeCount = ARSIZE(edges) / 2;
	UndirectedGraph graph(verticesCount);
	UndirectedGraph backup, backup0, backup1, backup2;
	for (unsigned int i = 0; i < edgeCount; ++i)
		graph.set(edges[2 * i], edges[2 * i + 1], true);
	cout << "graph : \n" << graph.getAdjList();
	backup = backup0 = backup1 = graph;
	assert(graph.colorGraph(3));
	assert(!graph.colorGraph(2));
	assert(backup0.bfsForest() == 1);
	cout << "bfs tree : \n" << backup0.getAdjList();
	assert(graph.dfsForest() == 1);
	cout << "dfs tree : \n" << graph.getAdjList();
	backup1.set(9, 7, false);
	backup1.set(9, 6, false);
	backup1.set(10, 5, false);
	backup1.set(10, 8, false);
	assert(backup1.connectedComponents() == 3);
	backup2 = backup1;
	assert(backup2.connectedComponents() == 3);
	assert(backup1.bfsForest() == 3);
	cout << "bfs tree : \n" << backup1.getAdjList();
	assert(backup2.dfsForest() == 3);
	cout << "dfs tree : \n" << backup2.getAdjList();
	return 0;
}

#include <GL/glut.h>
class MazeDemo{
	static int width;
	static int height;
	static unsigned int target;
	static const unsigned int size = 50;
	static set<unsigned int> path;
	static Maze maze;
public:
	static void reset(){
		unsigned int dest = 0;
		unsigned int count = 0;
		maze.init(size, size);
		maze.random();
		dest = rand() % (size * size - maze.getWallCount());
		for(unsigned int i = 0;i < size;++ i)
			for(unsigned int j = 0;j < size;++ j)
				if(!maze.isWall(i, j)){
					++count;
					if(count == dest){
						target = i * size + j;
						break;
					}
				}
		path.clear();
		maze.findPath(path, 0, target);
		cout << path.size() << endl;
	}
	static void init(){
		glClearColor(1.0, 1.0, 1.0, 1.0);
		reset();
	}
	static 	void reshape(int w, int h){
		width = w;
		height = h;
		glViewport(10, 10, width - 20, height - 20);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluOrtho2D(0, width, 0, height);
		glutPostRedisplay();
	}
	static void keyboardFunc(unsigned char ch, int x, int y){
		switch(ch){
			case 'n':reset();break;
		}
		glutPostRedisplay();
	}
	static void display(){
		unsigned int row = maze.getRow();
		unsigned int col = maze.getCol();
		double rLen = height / (double)row;
		double cLen = width / (double)col;
		glClear(GL_COLOR_BUFFER_BIT);
		glColor4f(0.0, 0.0, 0.0, 1.0);
		glBegin(GL_LINES);

		for(unsigned int i = 0;i <= col;++ i){
			glVertex2i(i * cLen, 0);
			glVertex2i(i * cLen, row * rLen);
		}
		for(unsigned int i = 0;i <= row;++ i){
			glVertex2i(0, i * rLen);
			glVertex2i(col * cLen, i * rLen);
		}
		glEnd();
		glBegin(GL_QUADS);
		for(unsigned int j = 0;j < col;++ j){
			unsigned int k = col - 1 - j;
			for(unsigned int i = 0;i < col;++ i){
				if(maze.isWall(i, k))
					glColor4f(0.0, 1.0, 0.0, 1.0);
				else
					glColor4f(0.0, 0.0, 1.0, 1.0);
				glVertex2i(i * cLen, j * rLen);
				glVertex2i(i * cLen + cLen, j * rLen);
				glVertex2i(i * cLen + cLen, j * rLen + rLen);
				glVertex2i(i * cLen, j * rLen + rLen);
			}
		}
		set<unsigned int>::iterator begin = path.begin(), end = path.end();
		double tol = min(cLen, rLen) / 4;
		while(begin != end){
			int i = *begin / col;
			int j = col - 1 - *begin % col;
			glColor4f(1.0, 1.0, 0.0, 1.0);
			glVertex2i(i * cLen + tol, j * rLen + tol);
			glVertex2i(i * cLen + cLen - tol, j * rLen + tol);
			glVertex2i(i * cLen + cLen - tol, j * rLen + rLen - tol);
			glVertex2i(i * cLen + tol, j * rLen + rLen - tol);
			++begin;
		}
		{
			int i = target / col;
			int j = col - 1 - target % col;
			glColor4f(1.0, 0.0, 0.0, 1.0);
			glVertex2i(i * cLen - tol, j * rLen - tol);
			glVertex2i(i * cLen + cLen + tol, j * rLen - tol);
			glVertex2i(i * cLen + cLen + tol, j * rLen + rLen + tol);
			glVertex2i(i * cLen - tol, j * rLen + rLen + tol);
		}
		glEnd();
		glutSwapBuffers();
	}
	static int run(int argc, char *argv[]){
		using namespace std;
		init();
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
		glutInitWindowPosition(0, 0);
		glutInitWindowSize(width, height);
		glutCreateWindow("maze");
		glutDisplayFunc(display);
		glutKeyboardFunc(keyboardFunc);
		glutReshapeFunc(reshape);
		init();
		glutMainLoop();
		return 0;
	}
};

int MazeDemo::width = 700;
int MazeDemo::height = 700;
Maze MazeDemo::maze;
unsigned int MazeDemo::target = 0;
set<unsigned int> MazeDemo::path;
int test6(int argc, char *argv[]) {
	MazeDemo::run(argc, argv);
	return 0;
}
int test7(int argc, char *argv[]) {
	Lexer lexer;
	double val = 0.0;
	for (unsigned int i = 0; i < 10; ++i) {
		Expr::randomPrefix(lexer);
		assert(Expr::prefix(lexer, val));
		cout << round(val * 1000) / 1000.0 << "\t:\t" << lexer << endl;
	}
	cout << endl << endl;
	for (unsigned int i = 0; i < 10; ++i) {
		Expr::randomPostfix(lexer);
		assert(Expr::postfix(lexer, val));
		cout << round(val * 1000) / 1000.0 << "\t:\t" << lexer << endl;
	}
	return 0;
}
int test8(int argc, char *argv[]) {
	const char *reps[] = { "000", "001", "010", "011", "100", "101", "110",
			"111" };
	const unsigned int repsCount = ARSIZE(reps);
	for (unsigned int i = 0; i < repsCount; ++i)
		cout << ' ' << reps[i];
	cout << endl;
	for (unsigned int i = 0; i < (1 << 8); ++i) {
		for (unsigned int j = 0; j < repsCount; ++j)
			cout << "   " << ((i & (1 << j)) == 0 ? '0' : '1');
		cout << endl;
	}
	return 0;
}
int test9(int argc, char *argv[]) {
	BooleanFunction function(5);
	srand(time(0));
	function.random(12);
	function.sumOfProductRep(cout);
	cout << endl << endl;
	function.andNegRep(cout);
	cout << endl << endl;
	function.orNegRep(cout);
	cout << endl << endl;
	function.nandRep(cout);
	cout << endl << endl;
	function.norRep(cout);
	cout << endl;

	BooleanFunction func(3);
	BitArray tmp(3);
	tmp.fill(true);
	func.addFunction(tmp);
	//tmp.fill(false);
	//func.addFunction(tmp);
	func.sumOfProductRep(cout);
	cout << endl << endl;
	func.andNegRep(cout);
	cout << endl << endl;
	func.orNegRep(cout);
	cout << endl << endl;
	func.nandRep(cout);
	cout << endl << endl;
	func.norRep(cout);
	cout << endl;

	const unsigned int weightsCount = 10;
	double weights[weightsCount];
	for (unsigned int i = 0; i < weightsCount; ++i)
		weights[i] = (rand() % 1000 - 500) / 100.0;
	ThresholdGate gate(weights, weightsCount, (rand() % 1000 - 500) / 100.0);
	cout << gate << endl;
	for (unsigned int i = 0; i < 10; ++i) {
		BitArray tmp(weightsCount);
		tmp.random();
		cout << "input : " << tmp << endl;
		cout << "output : " << gate.getValue(tmp) << endl;
	}
	return 0;
}
int test10(int argc, char *argv[]){
	int vals[] = {3, 43, 9, 7, -4, 1, 57, -89};
	BinSearchTree<int> binSearchTree(vals, ARSIZE(vals));
	cout << binSearchTree << endl;
	assert(!binSearchTree.add(7));
	assert(binSearchTree.add(45));
	cout << binSearchTree << endl;
	return 0;
}
int test11(int argc, char *argv[]){
	const unsigned int count = 12;
	const unsigned int iteration = 1000000;
	int vals[count], backup[count];
	int counts[count][count];

	memset(counts, 0, sizeof (counts));
	for(unsigned int i = 0;i < count;++ i)
		backup[i] = i;
	for(unsigned int i = 0;i < iteration;++ i){
		memcpy(vals, backup, (sizeof (int)) * count);
		badRandomPermute(vals, count);
		for(unsigned int j = 0;j < count;++ j){
			++ counts[vals[j]][j];
		}
	}
	for(unsigned int i = 0;i < count;++ i){
		for(unsigned int j = 0;j < count;++ j)
			cout << counts[i][j] * count / (double)iteration << '\t' << "  ";
		cout << endl;
	}
	return 0;
}
/*int main(int argc, char *argv[]) {
	return test6(argc, argv);
}*/
