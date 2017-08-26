#include <iostream>
#include <sstream>
#include <string>


int main(int argc, char** argv) {


	std::string s = "blubb!!";

	std::size_t start = 2;
	std::size_t end = 4;
	

	long level = 0;
	long N = 0;

	std::string options("2 200 -");

	std::string::size_type loc = options.find("-",0);

	std::cout << "loc: " << loc << std::endl;

	std::istringstream str_stream ( std::string (options, 0, loc) );

	str_stream >> level;

	std::cout << "level: " << level << std::endl;

	str_stream >> N;

	std::cout << "N: " << N << std::endl;
	
	return 0;

}