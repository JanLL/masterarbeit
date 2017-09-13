#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <cmath>

#include <time.h>

#include <adolc/adouble.h>



template<typename T>
T fraser_suzuki(T x, T h, T r, T wr, T sr, T z) {

	T value;

	T log_value = log(1 + (x-z)*(sr*sr-1)/(wr*sr));

	return h*exp(-log(r)/(log(sr)*log(sr)) * log_value*log_value);

} 


namespace TEST {

	enum Component {

		Component_X = 0,
		Component_Y = 1,
		
	};

}


//using namespace TEST;


int main(int argc, char** argv) {


	std::cout << Component_X << "\t" << Component_Y << std::endl;


	return 0;

}