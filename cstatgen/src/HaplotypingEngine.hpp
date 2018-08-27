// $File: HaplotypingEngine.hpp $
// $LastChangedDate:  $
// $Rev:  $
// Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)
#ifndef _HPE_HPP_
#define _HPE_HPP_

#include <string>
#include <stdexcept>
#include "Exception.hpp"
#include "HaplotypingUtils.hpp"

namespace cstatgen {
inline void reset_ped(Pedigree & ped)
{

	// FIXME: It does not make sense to me but I have to reset ped data object manually,
	// otherwise in a Python program that calls CHP::Apply multiple times in a loop,
	// ped.GetMarkerInfo(i)->freq.dim will not equal 0 after a few runs
	// complaining the same marker name has been previously used.
	// cannot yet figure out why as this is suppose to be a brand new ped object here!
	// UPDATE:
	// It might due to using it from Python. The swig generated wrapper did not delete it after use
	// anyways let me just manually clean up everything instead of wrestling with swig
	// UPDATE 2:
	// seems it's just the markerInfo part gives problems; ped.count and ped.familyCount are always 0

	// delete old pointer
	for (int i = 0; i < ped.markerInfoCount; i++)
		if (ped.markerInfo[i]) delete ped.markerInfo[i];
	if (ped.markerInfo) delete [] ped.markerInfo;
	if (ped.markerInfoByInteger) delete [] ped.markerInfoByInteger;
	// FIXME: Only Clear() is not going to delete the char * buffer pointer which unfortunately is protected ...
	// so there will still be memory leak here
	ped.markerNames.Clear();
	ped.markerLookup.Clear();
	ped.markerInfoByName.Clear();
	ped.markerCount = ped.markerInfoCount = ped.markerInfoSize = 0;
	// reset pointer
	ped.GrowMarkerInfo();
	// delete old pointer
	for (int i = 0; i < ped.count; i++)
		if (ped.persons[i]) delete ped.persons[i];

	for (int i = 0; i < ped.familyCount; i++)
		if (ped.families[i]) delete ped.families[i];

	if (ped.families) delete [] ped.families;
	if (ped.persons) delete [] ped.persons;
	ped.size = 10000;
	ped.count = ped.familyCount = ped.haveTwins = 0;
	// reset pointer
	ped.persons = new Person *[ped.size];
	ped.families = new Family * [1];
}


class HaplotypingEngine
{
public:
	HaplotypingEngine(const int verbose = 0) :
		__verbose(verbose), __mendelianErrorCount(0) { __ped = new Pedigree(); };
	~HaplotypingEngine() { delete __ped; };
	HaplotypingEngine * clone() const { return new HaplotypingEngine(*this); }
	// input chrom must be 1 .. 22 and X
	// input samples follows PED convention, e.g., { "1", "1", "0", "0", "1", "21", "21", "21" }
	// first 5 cols are fid, sid, pid, mid and sex (0/1/2 coding), followed by genotypes
	// (1/2 coding, 0 for missing)
	// positionAdjustment: adjust physical distance to map distance, 1 / 100 million
	VecVecVecString Execute(const std::string & chrom, const VecString & marker_names,
	                        const VecInt & marker_positions, const VecVecString & samples, 
							double Rsq, const char* logname, bool reorder=true,
	                        double positionAdjustment = 1E-8)
	{
		reset_ped(*__ped);
		try {
			DataLoader dl;
			dl.LoadVariants(__ped, marker_names, marker_positions, chrom, positionAdjustment);
			dl.LoadSamples(__ped, samples, marker_names);
			MendelianErrorChecker mc;
			mc.Apply(__ped);
			__mendelianErrorCount += mc.CountMendelianErrors();
			GeneticHaplotyper gh(chrom);
			gh.Apply(__ped,Rsq,logname,reorder);
			if (__verbose) gh.Print();
			return gh.data;
			// } catch (...) {
		} catch (std::exception e) {
			// std::clog << e.what() << std::endl;
			const VecVecVecString nulldata(0);
			return nulldata;
		}
	}


	int CountMendelianErrors() { return __mendelianErrorCount; }

private:
	int __verbose;
	int __mendelianErrorCount;
	Pedigree * __ped;
};

}
#endif
