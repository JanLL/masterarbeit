#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <functional>


int main(int argc, char** argv) {

	std::function<double (double)> func = [](double x) {return x*x;};

	std::cout << func(2) << std::endl;

	return 0;

}