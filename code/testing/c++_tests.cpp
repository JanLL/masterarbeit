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


int main(int argc, char** argv) {


	adouble x;

	x = 5.;


	//std::cout << x.val << std::endl;

	return 0;

}