// $File: HaplotypingUtils.hpp $
// $LastChangedDate:  $
// $Rev:  $
// Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)

#ifndef _HTPU_HPP_
#define _HTPU_HPP_

#include "Pedigree.h"
#include "MerlinFamily.h"
#include "MerlinHaplotype.h"
#include "MerlinSort.h"

#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <iterator>

#include "Exception.hpp"

namespace cstatgen {

typedef std::vector<int> VecInt;
typedef std::vector<std::vector<int> > VecVecInt;
typedef std::vector<std::vector<std::vector<int> > > VecVecVecInt;
typedef std::vector<std::string> VecString;
typedef std::vector<std::vector<std::string> > VecVecString;
typedef std::vector<std::vector<std::vector<std::string> > > VecVecVecString;
typedef std::vector<std::vector<double> > VecVecDouble;
typedef std::map<std::string, std::vector<double> > VecDoubleDict;
typedef std::map<std::string, std::vector<std::vector<double> > > VecVecDoubleDict;

class DataLoader
{
public:
	DataLoader() {};
	~DataLoader() {};
	DataLoader * clone() const { return new DataLoader(*this); }
	void LoadVariants(Pedigree * & ped,
		const VecString & names,
		const VecInt & positions,
		const std::string & chrom,
		double positionAdjustment = 0.01);

	void LoadSamples(Pedigree * & ped,
		const VecVecString & samples,
		const VecString & names);

private:
	void __AddPerson(Pedigree * & ped,
		const VecString & fam_info,
		const VecString & genotypes,
		const VecString & names);

};

class MendelianErrorChecker
{
public:
	MendelianErrorChecker() : __errorCount(0) {};
	~MendelianErrorChecker() {};
	MendelianErrorChecker * clone() const { return new MendelianErrorChecker(*this); }
	void Apply(Pedigree * & ped);

	int CountMendelianErrors() { return __errorCount; }

private:
	int __errorCount;

};

class GeneticHaplotyper
{
public:
	GeneticHaplotyper(const std::string & chrom) : data(0), __chrom(chrom) {}
	~GeneticHaplotyper() {};
	GeneticHaplotyper * clone() const { return new GeneticHaplotyper(*this); }
	// [family][sample][haplotypes]
	VecVecVecString data;
	// Apply haplotyping. Missing data are imputed as possible
	void Apply(Pedigree * & ped, double Rsq, const char * logname, bool reorder=true);

	void Print();

private:
	std::string __chrom;
};


class HaplotypeCoder
{
public:
	HaplotypeCoder(const double size) : __data(0), __freqs(), __recombCount(0), __size(size) {}
	~HaplotypeCoder() {};
	HaplotypeCoder * clone() const { return new HaplotypeCoder(*this); }

	// each element of haploVecs is a family's data
	// each element of haploVecs[i] is a haplotype with the first 2 items being fid and sid
	// each element of maf is a family's data
	// each element of maf[i] is founder population MAF of the corresponding variant
	void Execute(const VecVecVecString & haploVecsConst, const VecVecDouble & mafVecsConst,
		const VecVecVecInt & markerIdxClusters);

	void Print();

	VecVecString GetHaplotypes() { return __data; }
	VecVecDouble GetAlleleFrequencies(const std::string & family) { return __freqs[family]; }
	int CountRecombinations() { return __recombCount; }

private:
	VecVecString __data;
	VecVecDoubleDict __freqs;
	int __recombCount;
	double __size;
};
}
#endif
