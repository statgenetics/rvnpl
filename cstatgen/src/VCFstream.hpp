// $File: VCFstream.hpp $
// $LastChangedDate:  $
// $Rev:  $
// Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)

#ifndef _VCFSTM_HPP_
#define _VCFSTM_HPP_
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include "Exception.hpp"
#include "VcfFileReader.h"

namespace cstatgen {

typedef std::vector<int> VecInt;
typedef std::vector<std::string> VecString;
typedef std::vector<std::vector<std::string> > VecVecString;

class VCFstream
{
public:
	/// Initialize VCFstream with VCF file
	/// \param vcf file
	/// \param vcf index
	VCFstream(const std::string & vcf) : __vcf(vcf)
	{
		std::string msg;

		try {
			__reader.open(vcf.c_str(), __header);
		} catch (std::exception & e) {
			msg.assign(e.what());
			throw RuntimeError("[BAD VCF file] " + msg);
		}
		try {
			__reader.readVcfIndex();
		} catch (std::exception & e) {
			msg.assign(e.what());
			throw RuntimeError("[BAD index file] " + msg);
		}

		const Tabix * tabixPtr = __reader.getVcfIndex();
		if (tabixPtr == NULL || tabixPtr->getFormat() != Tabix::FORMAT_VCF)
			throw RuntimeError("Failed to load a proper VCF index");
		sampleCount = __header.getNumSamples();
	}


	~VCFstream()
	{
		__reader.close();
	}


	/// Get list of sample names
	/// \return list of sample names
	VecString GetSampleNames()
	{
		VecString samples(sampleCount);

		for (unsigned i = 0; i < samples.size(); ++i) {
			samples[i].assign(__header.getSampleName(i));
		}
		return samples;
	}


	/// Get all variants in VCF file
	/// \return list of (chr, pos)
	VecVecString GetGenomeCoordinates()
	{
		VecVecString res(0);
		VcfRecord line;
		VcfFileReader reader;
		VcfHeader header;

		reader.open(__vcf.c_str(), header);
		reader.readVcfIndex();
		while (reader.readRecord(line)) {
			VecString variant(2);
			variant[0] = line.getChromStr();
			variant[1] = std::to_string(line.get1BasedPosition());
			res.push_back(variant);
		}
		return res;
	}


	/// Extract VCF region
	/// \param chrom
	/// \param start pos
	/// \param end pos
	void Extract(const std::string & chrom, int start, int end)
	{
		try {
			__reader.set1BasedReadSection(chrom.c_str(), start, end);
		} catch (std::exception & e) {
			throw RuntimeError("Failed to extract VCF region");
		}
	}


	/// Point to next record
	/// \return true if __line is valid otherwise false
	bool Next()
	{
		return __reader.readRecord(__line);
	}


	std::string GetChrom()
	{
		std::string chrom = __line.getChromStr();

		return chrom;
	}


	int GetPosition()
	{
		return __line.get1BasedPosition();
	}


	std::string GetInfo(const std::string & key)
	{
		VcfRecordInfo & info = __line.getInfo();
		const std::string * s = info.getString(key.c_str());

		if (s) return *s;
		else throw ValueError("VCF line does not have info field " + key);
	}


	int CountSampleGenotypes()
	{
		return __line.getNumSamples();
	}


	bool IsBiAllelic()
	{
		return __line.getNumAlts() == 1;
	}


	/// Get sample genotype
	/// \param variant ID
	/// \param sample IDs
	/// \return list of sample genotypes
	VecString GetGenotypes(const VecInt & sid)
	{
		VecString genotypes(sid.size());

		for (unsigned i = 0; i < sid.size(); ++i) {
			int allele1 = __line.getGT(sid[i], 0);
			int allele2 = __line.getGT(sid[i], 1);
			allele1 = (allele1 == 0 || allele1 == 1) ? allele1 + 1 : 0;
			allele2 = (allele2 == 0 || allele2 == 1) ? allele2 + 1 : 0;
			genotypes[i] = std::to_string(allele1) + std::to_string(allele2);
		}
		return genotypes;
	}


	VecString GetGenotypes()
	{
		VecString genotypes(sampleCount);

		for (unsigned i = 0; i < sampleCount; ++i) {
			int allele1 = __line.getGT(i, 0);
			int allele2 = __line.getGT(i, 1);
			allele1 = (allele1 == 0 || allele1 == 1) ? allele1 + 1 : 0;
			allele2 = (allele2 == 0 || allele2 == 1) ? allele2 + 1 : 0;
			genotypes[i] = std::to_string(allele1) + std::to_string(allele2);
		}
		return genotypes;
	}


	unsigned sampleCount;

private:
	VCFstream(const VCFstream &);
	VCFstream & operator=(const VCFstream &);

	std::string __vcf;
	VcfRecord __line;
	VcfFileReader __reader;
	VcfHeader __header;
};
}
#endif
