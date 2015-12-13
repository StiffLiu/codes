#include <limits>
#include <set>
namespace my_lib{

template<typename T>
class MyAllocator {
public : 
	//    typedefs
	typedef T value_type;
	typedef value_type* pointer;
	typedef const value_type* const_pointer;
	typedef value_type& reference;
	typedef const value_type& const_reference;
	typedef std::size_t size_type;
	typedef std::ptrdiff_t difference_type;

public : 
	//    convert an allocator<T> to allocator<U>
	template<typename U>
	struct rebind {
		typedef MyAllocator<U> other;
	};

public : 
	inline explicit MyAllocator() {}
	inline ~MyAllocator() {}
	inline explicit MyAllocator(MyAllocator const&) {}
	template<typename U>
	inline explicit MyAllocator(MyAllocator<U> const&) {}

	//    address
	inline pointer address(reference r) { return &r; }
	inline const_pointer address(const_reference r) { return &r; }

	//    memory allocation
	inline pointer allocate(size_type cnt, 
			typename std::allocator<void>::const_pointer = 0) { 
		return reinterpret_cast<pointer>(malloc(cnt * sizeof (T))); 
	}
	inline void deallocate(pointer p, size_type) { 
		free(p); 
	}

	//    size
	inline size_type max_size() const { 
		return std::numeric_limits<size_type>::max() / sizeof(T);
	}

	//    construction/destruction
	inline void construct(pointer p, const T& t) {
	       	new(p) T(t);
       	}
	inline void destroy(pointer p) { p->~T(); }

	inline bool operator==(MyAllocator const&) { return true; }
	inline bool operator!=(MyAllocator const& a) { return !operator==(a); }
};    //    end of class Allocator 

struct MemAllocator{
	std::set<void*, std::less<void*>, MyAllocator<void*> > allocated;
	~MemAllocator(){
		assert(allocated.empty());
	}
	static MemAllocator& get(){
		static MemAllocator instance;
		return instance;
	}
	static void* allocate(unsigned int size){
		auto p = malloc(size);
		//printf("allocate memory : %p\n", p);
		get().allocated.insert(p);
		return p;
	}
	
	static void deallocate(void *p){
		if(p == nullptr)
			return;
		MemAllocator& m = get();
		assert(m.allocated.find(p) != m.allocated.end());
		m.allocated.erase(p);
		//printf("deallocate memory : %p\n", p);
		free(p);
	}
};
}
void* operator new(size_t size){
	return my_lib::MemAllocator::allocate(size);
}
void operator delete(void *p){
	return my_lib::MemAllocator::deallocate(p);
}
