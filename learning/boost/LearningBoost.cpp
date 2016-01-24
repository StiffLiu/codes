// LearningBoost.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/deque.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/date_time.hpp>
#include <boost/locale.hpp>
#include <vector>
#include <deque>

using namespace std;
namespace test_archive{
  class gps_position{
   public:
     gps_position(int d = 0, int m = 0, float s = 0.0) : degrees(d), minutes(m), seconds(s){
     }

     template<class Archive>
     void serialize(Archive &ar, const unsigned int version){
       ar & degrees;
       ar & minutes;
       ar & seconds;
     }
     int degrees;
     int minutes;
     float seconds;
  };
  const int cnt = 100;
  void save(const char *filename){
    std::ofstream ofs(filename, ios::binary);
    std::vector<gps_position> positions(cnt);
    for(int i = 0;i < cnt;++ i){
      gps_position& g = positions[i];
      g.degrees = rand() % 1000;
      g.minutes = rand() % 60;
      g.seconds = rand() % 1000;
      std::cout << "save : d=>" << g.degrees << ", m=>" << g.minutes <<", s=>" << g.seconds << std::endl;
    }

    boost::archive::binary_oarchive toa(ofs);
    toa << positions;
  }

  void read(const char *filename){
    std::ifstream ifs(filename, ios::binary);
    std::deque<gps_position> positions;
    boost::archive::binary_iarchive toa(ifs);

    toa >> positions;
    for(int i = 0;i < cnt;++ i){
      gps_position& g = positions[i];
      std::cout << "read : d=>" << g.degrees << ", m=>" << g.minutes <<", s=>" << g.seconds << std::endl;
    }
  }
}

namespace test_datetime{
	void test(int argc, char *argv[]){
		// extra words in format
		using namespace boost::gregorian;
		using namespace boost::posix_time;
		using namespace boost::local_time;
		std::stringstream ss;
		local_date_time ldt(not_a_date_time);
		string my_format("%Y-%m-%d %H:%M:%s");
		local_time_facet output_facet;// = new local_time_facet();
		local_time_input_facet* input_facet = new local_time_input_facet();
		ptime pt = boost::posix_time::time_from_string("2014-01-08 21:00:22.530");
		ss.imbue(locale(locale::classic(), &output_facet));
		ss.imbue(locale(ss.getloc(), input_facet));
		output_facet.format(my_format.c_str());
		input_facet->format(my_format.c_str());
		ss.str("");
		ss << ldt;
		cout << ss.str() << endl;
		ss.str("");
		ss << pt;
		cout << ss.str() << endl;
		ss.str("2014-01-08 21:00:22.530");
		ss >> ldt;
		ss.str("");
		ss << ldt;
		ptime tmp(date(1970, 01, 01), milliseconds(1388653201000));
		cout << ss.str() << endl;
		cout << tmp << endl;
		// matching extra words in input 
		ss.str("The extended ordinal time 2005-128T12:15 can also be \
			     represented as Sunday May 08, 2005");
		ss >> ldt;
		cout << ldt << endl;

	}
}

namespace test_network{
struct Info{
  std::wstring name;
  std::wstring tel;
  std::wstring addr;
  template<class Stream>
  friend Stream& operator<<(Stream& os, const Info& info){
    os << info.name << '\t' << info.addr << '\t' << info.tel << '\n';
    return os;
  }
};
 void extract_info(boost::property_tree::wptree& tree, std::vector<Info>& infos){
   typedef boost::property_tree::wptree PTree;
   for(PTree::iterator begin = tree.begin(); begin != tree.end();++ begin){
     if(begin->first != _T("content"))
       continue;
     for(PTree::iterator begin1 = begin->second.begin(); begin1 != begin->second.end();++ begin1){
       infos.push_back(Info());

       Info& info = infos.back();
       for(PTree::iterator begin2 = begin1->second.begin(); begin2 != begin1->second.end();++ begin2){
         if(begin2->first  == _T("name")){
           info.name = begin2->second.data();
         }else if(begin2->first == _T("tel")){
           info.tel = begin2->second.data();
         }else if(begin2->first == _T("addr")){
           info.addr = begin2->second.data();
         }
       }
     }
   }
 }
 int readFromFile(const char *filePath, std::vector<Info>& infos)
 {
   try{
     std::wifstream in(filePath);
     boost::property_tree::wptree pt;
     boost::property_tree::read_json(in, pt);
     extract_info(pt, infos);
   }catch (std::exception& e){
     std::cout << "Exception when reading file " << filePath << ", reason" << e.what() << "\n";
     return -1;
   }
   return 0;
 }
 int test(int argc, const char* argv[])
{
  try
  {
    if (argc < 3)
    {
      std::cout << "Usage: http_client <server> <path> <dest>\n";
      std::cout << "Example:\n";
      std::cout << "  http_client www.boost.org /LICENSE_1_0.txt\n";
      return 1;
    }

    boost::asio::ip::tcp::iostream s;

    // The entire sequence of I/O operations must complete within 60 seconds.
    // If an expiry occurs, the socket is automatically closed and the stream
    // becomes bad.
    s.expires_from_now(boost::posix_time::seconds(80));

    // Establish a connection to the server.
    s.connect(argv[1], "http");
    if (!s)
    {
      std::cout << "Unable to connect: " << s.error().message() << "\n";
      return 1;
    }

    // Send the request. We specify the "Connection: close" header so that the
    // server will close the socket after transmitting the response. This will
    // allow us to treat all data up until the EOF as the content.
    s << "GET " << argv[2] << " HTTP/1.0\r\n";
    //s << "Host: " << argv[1] << "\r\n";
    s << "Accept: */*\r\n";
    s << "Connection: close\r\n\r\n";

    // By default, the stream is tied with itself. This means that the stream
    // automatically flush the buffered output before attempting a read. It is
    // not necessary not explicitly flush the stream at this point.

    // Check that response is OK.
    std::string http_version;
    s >> http_version;
    unsigned int status_code;
    s >> status_code;
    std::string status_message;
    std::getline(s, status_message);
    if (!s || http_version.substr(0, 5) != "HTTP/")
    {
      std::cout << "Invalid response\n";
      return 1;
    }
    if (status_code != 200)
    {
      std::cout << "Response returned with status code " << status_code << "\n";
      return 1;
    }
    

    // Process the response headers, which are terminated by a blank line.
    std::string header;
    while (std::getline(s, header) && header != "\r")
      std::cout << header << "\n";
    std::cout << "\n";
    
    std::ostringstream os;
    os << s.rdbuf();
    
    std::string content = os.str();
    boost::algorithm::trim(content);
    std::ofstream ofs(argv[3]);
    ofs << content;
  }
  catch (std::exception& e)
  {
    std::cout << "Exception when downloading from url " << argv[2] << "  reason : " << e.what() << "\n";
    return -1;
  }

  return 0;
}
int downloaDatas(int argc, char *argv[]){
  std::vector<test_network::Info> infos;
  const char *args[]={"test_network", "ditu.baidu.com", NULL, NULL};
  
  for(int i = 76;i < 100;++ i){
     std::string url;
     std::string filePath;
    {
      std::ostringstream os;
      os << "http://ditu.baidu.com/?"
              "newmap=1&reqflag=pcmap&biz=1&pcevaname=pc2&da_par=direct&from=webmap&qt=s&da_src=pcmappg.searchBox.button&"
              "wd=%E5%9B%BE%E4%B9%A6%E9%A6%86&c=289&src=0&wd2=&sug=0&l=6&b=(13482217.92,3641546.71;13537897.92,3660234.71)&"
              "from=webmap&tn=B_NORMAL_MAP&ie=utf-8&t=1426811220222&nn="
         << (i - 76) * 10;
      url = os.str();      
      args[2] = url.c_str();
    }
    {
      std::ostringstream os;
      os << "D:\\dat\\" << i << ".json";
      filePath = os.str();
      args[3] = filePath.c_str();
    }
    for(int j = 0;j < 4 && test_network::test(3, args) != 0;++ j);
    readFromFile(args[3], infos);
  }
  std::ofstream ofs("D:\\dat\\result.txt");
  std::string encoding("gb2312");
  for(int i = 0;i < infos.size();++ i){
    wostringstream os;
    os << infos[i];
    ofs << boost::locale::conv::from_utf(os.str(), encoding);
  }
  return 0;
}
}
int _tmain(int argc, _TCHAR* argv[]){
  /*test_archive::save("filename");
  test_archive::read("filename");
  cin.get();
  const char *args[]={"test_network", "news.163.com", "http://news.163.com", "D:\\news.txt"};
  test_network::test(4, args);*/
  test_datetime::test(0, nullptr);
  return 0;
}

