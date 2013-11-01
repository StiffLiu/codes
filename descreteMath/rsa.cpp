#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
using namespace std;

template<class T>
T modularExponentiation(T e, T b, T n){
  T r = 1;
  T p = n % b;
  T o = 1;
  while(e != 0){
    if((e & o) != 0) r = (r * p) % b;
    p = (p*p) % b;
    e >>= 1;
  }
  return r;
}

template<class T>
void rsa(T *m, T *o, int c, T n, T e){
  for(int i = 0;i < c;++ i)
    o[i] = modularExponentiation(e, n, m[i]);
}
template<class T>
T gcd(T a, T b){
  while(b != 0){
    T tmp = b;
    b = a % b;
    a = tmp;
  }
  if(a < 0)
    return -a;
  return a;
}
template<class T>
T extGcd(T a, T b, T& s, T& t){
  T sn[3], tn[3];
  int i = 0;
  T q;
  sn[0] = tn[1] = 1;
  sn[1] = tn[0] = 0;
  sn[2] = tn[2] = 0;
  while(b != 0){
    T tmp = b;
    ++i;
    if(i > 1){
      sn[i%3] = sn[(i-2)%3] - q * sn[(i-1)%3];
      tn[i%3] = tn[(i-2)%3] - q * tn[(i-1)%3];
    }
    q = a / b;
    b = a % b;
    a = tmp;
  }
  s = sn[i%3], t = tn[i%3];
  if(a < 0){
    s = -s;
    t = -t;
    return -a;
  }
  return a;
}
void func(int a, int b){
  int s, t, cd;
  cd = extGcd(a, b, s, t);
  cout << cd << '=' << a << '*';
  if(s < 0)
    cout << '(' << s << ')';
  else
    cout << s;
  cout << '+' << b << '*';
  if(t < 0)
    cout << '(' << t << ')';
  else
    cout << t;
}
template<class T>
bool chineseRemainder(T *modolees, T * moduluses, int n, T& ret, T& tot){
  if(n < 1)
    return false;
  T total = moduluses[0];
  T current = modolees[0] % total;
  cout << current << endl;
  for(int i = 1;i < n;++ i){
    T s, t;
    T tmp = total * moduluses[i];
    if(extGcd(total, moduluses[i], s, t) != 1)
      return false;
    func(total, moduluses[i]);cout << endl;
    cout << current << '*' << moduluses[i] << '*' << t << '+' << modolees[i] << '*' << total << '*' << s;
    current = ((((current * moduluses[i]) % tmp) * t)%tmp + (((modolees[i] * total) % tmp) * s)%tmp) % tmp;
    cout << '=' << current << endl;
    total = tmp;
  }
  ret = current;
  tot = total;
  return true;
}
class ShiftCypher{
	int shift;
	int modulee;
	int moduli;
	int inverseModulee;
public:
	ShiftCypher(int modulee, int moduli, int shift){
		int tmp;
		this->shift = shift % moduli;
		this->modulee = modulee % moduli;
		this->moduli = moduli;
		assert(extGcd(this->modulee, moduli, inverseModulee, tmp)==1);
		this->inverseModulee = inverseModulee % moduli;
	}
	
	int encrypt(int input){
		int ret = (modulee*input+shift)%moduli;
		if(ret < 0)
			ret += moduli;
		return ret;
	}
	int decrypt(int input){
		int ret = (inverseModulee*(input-shift))%moduli;
		if(ret < 0)
			ret += moduli;
		return ret;
	}
	void encrypt(int *input, int *output, int n){
		for(int i = 0;i < n;++ i){
			output[i] = (modulee*input[i]+shift)%moduli;
		}
	}
	void decrypt(int *input, int *output, int n){
		for(int i = 0;i < n;++ i){
			output[i] = inverseModulee*(input[i]-shift)%moduli;
		}
	}
};
void outMsg(const char *input, bool dprpt = true){
  ShiftCypher cypher(7, 26, 10);
  while(*input){
    char ch = *input;
    if(ch >= 'a' && ch <= 'z')
      ch -= 32;
      int diff = ch - 'A';
      if(diff >= 0 && diff < 26){
        if(dprpt)
          cout << (char)(cypher.decrypt(diff) + 'A');
        else
          cout << (char)(cypher.encrypt(diff) + 'A');
      }else
        cout << ch;
      ++ input;
  }
}
void findPrime(vector<bool>& ret){
	int n = ret.size();
	int i = 0;
	if(n < 1)
		return;
	ret[0] = false;
	if(n < 2)
		return;
	ret[1] = true;
	while(i < n){
		while(i < n && !ret[i]) ++i;
		for(int j = 2;j * (i + 1) <= n;++ j)
			ret[j * (i + 1) - 1] = false;
		++i;
	}
}
template<class T>
bool isPrime(T p){
	T primes[]={2, 3, 5, 7, 11};
	for(int i = 0;i < (sizeof primes / sizeof *primes);++ i)
		if((p % primes[i] == 0) || modularExponentiation((p - 1), p, primes[i]) != 1)
			return false;
	return true;
}
ostream& outValue(int i){
	if(i >= 0)
		cout << i;
	else
		cout << '(' << i << ')';
	return cout;
}
void out(int *values, int *cf, int n){
  for(int i = 0;i < n;++ i){
	outValue(values[i]) << '*';
	outValue(cf[i]);
	if(i != n -1)
		cout << '+';
  }
}
int gcd(int *values, unsigned int n){
	assert(n > 0);
	int tmp = abs(values[0]);
	for(int i = 1;i < n;++ i){
		tmp = gcd(tmp, values[i]);
	}
	return tmp;
}
int extGcd(int *values, int *ret, unsigned int n){
	assert(n > 0);
	int tmp = abs(values[0]);
	ret[0] = (values[0] > 0 ? 1 : -1);
	for(int i = 1;i < n;++ i){
		int s, t;
		tmp = extGcd(tmp, values[i], s, t);
		ret[i] = t;
		for(int j = 0;j < i;++ j)
			ret[j] *= s;
	}
	return tmp;
}
int test1(int argc, char *argv[]){
	int values[] = {-64, 32, 48, 40, -44, 79};
	const int count = (sizeof values / sizeof *values);
	int ret[count];
	cout << "gcd(";
	for(int i = 0;i < count;++ i)
		cout << values[i] << (i < count - 1 ? ',':')');
	cout << "=" <<  extGcd(values, ret, count) << "=";
	out(values, ret, count);
	cout << endl;
	return 0;
}
int test0(int argc, char *argv[]){
  unsigned int p = 43, q = 59;
  unsigned int n = p*q;
  unsigned int e = 13;
  unsigned int c = 937;
  unsigned int temp = 0;
  if(extGcd(e, (p-1)*(q-1), c, temp) == 1){
    unsigned int message[]={667,1947,671};
    const int count = (sizeof message / sizeof *message);
    unsigned int out[count];
    if(c < 0)
      c += (p-1)*(q-1);
    cout << c << endl;
    rsa(message, out, count, n, c);
    for(int i = 0;i < count;++ i){
      short a1 = out[i] / 100, a2 = out[i] % 100;
      cout << (char)('A'+a1) << (char)('A'+a2);
    }
    cout << endl;
  }
  {
    cout << gcd(135, 96) << endl;
    cout << gcd(-135, 96) << endl;
    cout << gcd(-135, -96) << endl;
    cout << gcd(135, -96) << endl;
  }
  cout << "------------------" << endl;
  {
    int s, t, a = 135, b = 96;
    int cd = extGcd(a, b, s, t);
    cout << cd << ' ' << s << ' ' << t << endl;
    cd = extGcd(-a, b, s, t);
    cout << cd << ' ' << s << ' ' << t << endl;
    cd = extGcd(-a, -b, s, t);
    cout << cd << ' ' << s << ' ' << t << endl;
    cd = extGcd(a, -b, s, t);
    cout << cd << ' ' << s << ' ' << t << endl;
    cd = extGcd(-b, a, s, t);
    cout << cd << ' ' << s << ' ' << t << endl;
    func(144, 89);cout << endl;
    func(1001, 100001);cout << endl;
  }
  cout << "------------------" << endl;
  {
    long long modus[]= {99, 98, 97, 95};
    const int count = (sizeof modus/ sizeof *modus);
    long long modes[count]={65, 2, 51, 10};
    long long ret, total;
    if(chineseRemainder(modes, modus, count, ret, total)){
      cout << total << ":" << ret << endl;
    }
  }
  cout << "------------------" << endl;
  {
    outMsg("LJMKGM GMXF QEXMW");
    cout << endl;
    outMsg("TUYC WEU");
    cout << endl;
    bool test = false;
    if(test){
      string msg;
      do{
        getline(cin, msg);
	outMsg(msg.c_str(), false);
        cout << endl;
      }while(!msg.empty());
    }
  }
  cout << "------------------" << endl;
  cout << isPrime((unsigned long long)172947529) << endl;
  {
    unsigned int total = 100000000;
    vector<bool> tmp(total, true);
    int count = 0, count1 = 0;
    findPrime(tmp);
    unsigned int i = 1;
    unsigned int segments[] = {1000000,10000000-1000000,total-10000000};
    unsigned int current = 0;
    for(int j = 0;j < (sizeof segments / sizeof *segments);++ j){
      current += segments[j];
      for(;i < current;++ i){
        if(tmp[i]) ++ count1;
        if(!tmp[i] && isPrime((unsigned long long)(i + 1))){
          //cout << (i + 1) << ' ';
          ++ count;
        //if(count % 8 == 0)
        //  cout << endl;
        }
      }
      cout << "\n---------counts prime:" << count1 << " pesudoprimes:" << count << endl;
    }
    count = 0;
    for(unsigned int i = 1;i * i + 1<= total;++ i){
      if(tmp[i*i]){
        cout << (i*i+1) << ' ';
        ++count;
        if(count % 8 == 0)
          cout << endl;       
      }
    }
  }
  return 0;
}

int main(int argc, char *argv[]){
	test1(argc, argv);
}
