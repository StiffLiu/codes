/*
 * An illustration of swapable stack, described as following:
 * A region of disk space(external memory) is mapped to a limited amount 
 * of RAM(internal memory) in unit of page, where each page contains a specificed number 
 * of bytes.
 *
 * When an element(called a word) is pushed on to the stack it first checks the RAM,
 * if the RAM is full, a page will be swapped out to a disk page.
 * When an element is poped from the stack and the RAM for the stack is empty,
 * it will swap in a page from disk to RAM.
 *
 * Assuming that a page can store at most m elements, then the amortized disk access for
 * each push and pop operation takes 1/m disk access.
 *
 * This code is just for illustration purpose.
 * */
template<size_t N, class PageSetTraits>
class PageSet{
	char _buffer[N * PageSetTraits::pageSize()];
	size_t _current = 0;
	size_t _count = 0;
public:
	bool empty() const {return 0 == _count;}
	bool full() const {return N == _count;}
	char* front() {return &_buffer[_current * N];}
	char* back() {return &_buffer[(_current + _count - 1) % N *N];}
	void pop(){_current = (_current + 1) % N; assert(_count > 0);--_count;}
	void push(){assert(_count <= N);++_count;}
	size_t size(){return _count;}
	template<class PageNum>
	void write_front(PageNum pageNum_){
		assert(!empty());
		write_disk(pageNum_ - size(), front());
	}
	template<class PageNum>
	void read_back(PageNum pageNum_){
		assert(!empty());
		read_disk(pageNum_ - 1, back());
	}
};

template<class PageSetTraits, class Address = size_t>
class Stack{
	auto _pages = PageSetTraits::create();
	Address _current = Address();
	static_assert(PageSetTraits::pageSize() > 1);
public:

	void push(const Word& value_){
		auto offset = _current % PageSetTraits::pageSize();
		if(offset == 0){
			auto pageNum = _current / PageSetTraits::pageSize();
			if(_pages.full()){
				assert(pageNum > _pages.size());
				_pages.write_front(pageNum);
				_pages.pop();
			}
			_pages.push();
		}
		_pages.back()[offset] = value_;
		++_current;
	}

	bool pop(Word& word_){
		if(_current == 0) return false;
		auto offset = _current % PageSetTraits::pageSize();
		if(offset == 1){
			assert(!_pages.empty());
			word_ = _pages.back()[0];
			_pages.pop();
			--_current;
			return true;
		}
		
		if(offset == 0 && _pages.empty()){
			auto pageNum = _current / PageSetTraits::pageSize();
			assert(pageNum > 1);
			_pages.push();
			_pages.read_back(pageNum - 1);
		}
		--_current;
		word_ = _pages.back()[_current % PageSetTraits::pageSize()];
		return true;
	}
};
