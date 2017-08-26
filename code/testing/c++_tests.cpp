#include <iostream>
#include <fstream>
#include <string>
#include <set>




int main(int argc, char** argv) {

	std::set<double> s1;
	std::set<double>::iterator it1, it2;

	s1.insert(3.);
	s1.insert(3.5);
	s1.insert(10.);
	s1.insert(-2.);
	s1.insert(6.);
	s1.insert(7.);

	/*
	it1 = s1.lower_bound(3.6);
	it2 = it1;
	it2--;

	std::cout << *it2 << "\t" << *it1 << std::endl;
	*/

	for(it1=s1.begin(); it1!=s1.end(); ++it1) {
		std::cout << *it1 << "\t" << *(++it1) << std::endl;
	}

	return 0;

}