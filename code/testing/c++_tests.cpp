#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <sstream>


#include <time.h>



int main(int argc, char** argv) {


	std::vector<double> v1(3, 0.);
	std::vector<double> v2(3, 6.);
	
	v1[0] = 1.;
	v1[1] = 1.;
	v1[2] = 1.;
	
	v2 = v1;

	v1[0] = 999.;

	for (int i=0; i<3; ++i) {
		std::cout << v1[i] << "\t" << v2[i] << std::endl;
	}


	return 0;

}