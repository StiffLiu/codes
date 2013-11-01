#include <iostream>
#include <map>
using namespace std;

void generateUlam(int n, int *numbers){
  if(n < 1)
    return;
  numbers[0] = 1;
  if(n < 2)
    return;
  numbers[1] = 2;
  map<int, int> counts;
  map<int, int>::iterator pos = counts.begin();
  for(int i = 2;i < n;++ i){
    for(int k = 0;k < i - 1;++ k){
      ++ counts[numbers[k] + numbers[i - 1]];
    }
    if(pos == counts.end()){
      pos = counts.begin();
      while(pos != counts.end()){
	if(pos->first > numbers[i - 1])
	  break;
	++pos;
      }
    }
    while(pos != counts.end()){
      if(pos->second == 1){
	numbers[i] = pos->first;
	++pos;
	break;
      }
      ++pos;
    }
  }
  
}

int main(int argc, char *argv[]){
  int tmp[1000];
  generateUlam(1000, tmp);
  for(int i = 0;i < 1000;++ i){
    cout << tmp[i] << ' ';
    if((i + 1) % 8 == 0)
      cout << endl;
  }
}
