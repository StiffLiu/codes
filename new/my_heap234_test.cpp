#include <iostream>
#include <cstdlib>
#include <cassert>
#include <ctime>
#include <cstdlib>
#include "my_mem_allocator.h"
#include "my_heap234.h"
int main(int argc, char *argv[]){
	{
		using namespace my_lib;
		using namespace std;
		srand(time(0));
		Heap234<unsigned int> h1, h2;
		for(unsigned int i = 0;i < 1000;++ i){
			int key = rand() % 100;
			//cout << " push " << key << std::endl;
			h1.push(key);
			if(rand() % 4 == 0) h2.push(key);
			//cout << "h1 : \n" << h1 << std::endl;
			//cout << "h2 : \n" << h2 << std::endl;
			try{
				assert(h1.isValid());
				assert(h2.isValid());
				if(h2.size() > 10 && rand() % 4 == 0){
					//std::cout << "merge two 2-3-4 node heap" << std::endl;
					if(rand() % 2 == 0) h2.merge(h1);
					else h1.merge(h2);
					h1.merge(h2);
					//cout << h1 << std::endl;
					assert(h1.isValid());
				}
			}catch(int lineNo){
				std::cout << "invalid at line : " << lineNo << std::endl;
				break;
			}
		}
		std::cout << "end of test" << std::endl;
	}
	std::cout << "unfreed : " << my_lib::MemAllocator::get().allocated.size() << std::endl;
	return 0;
}
