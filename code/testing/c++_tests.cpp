#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <vector>
#include <algorithm>


void func() {

	static std::vector<double> vec(3, 0.);

	for (int i=0; i<3; ++i) {
		std::cout << vec[i] << std::endl;
	}

	for (int i=0; i<3; ++i) {
		vec[i] = i; 
	}


	for (int i=0; i<3; ++i) {
		std::cout << vec[i] << std::endl;
	}




}



int main(int argc, char** argv) {

	func();
	func();

	return 0;

}