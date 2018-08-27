#include "VCFstream.hpp"
#include <iostream>
#include <vector>
#include <string>
int main()
{
	cstatgen::VCFstream vs("test.vcf.gz");

	for (auto i : vs.GetSampleNames()) std::cout << i << ' ';
	std::cout << std::endl;
	// vs.Extract("1", 69269, 801942);
	vs.Extract("1", 32768, 801942);
	while (vs.Next()) {
		std::cout << vs.GetChrom() << " " << vs.GetPosition() << " " << vs.IsBiAllelic() << " " << vs.GetInfo("AF") << std::endl;
		for (auto i : vs.GetGenotypes()) std::cout << i << ' ';
		std::cout << std::endl;
		vector<int> idx { 1, 2 };
		for (auto i : vs.GetGenotypes(idx)) std::cout << i << ' ';
		std::cout << std::endl;
	}
	return 0;
}


