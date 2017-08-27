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


	double a_const = 1.32313123213123; 
	std::ostringstream oss;



	double time_start, time_duration;

	time_start = clock();  // start time measurement
	for (size_t i=0; i<200; ++i) {
		oss << a_const;
		std::string s = oss.str();
		double a_const_double = std::stod(s, NULL);
	}
	time_duration = (clock() - time_start) / CLOCKS_PER_SEC;


	std::cout << time_duration << " s\n";


	return 0;

}