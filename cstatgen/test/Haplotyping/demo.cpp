#include <iostream>
#include "HaplotypingUtils.hpp"
using namespace cstatgen;

void showPed(Pedigree & ped)
{
	printf("Loaded %d individuals\n", ped.count);
	for (int i = 0; i < std::min(ped.count, 10); i++) {
		printf("[%s]: %s, %s, %s, %d\t|\t",
			(const char *)ped[i].pid, (const char *)ped[i].famid,
			(const char *)ped[i].fatid, (const char *)ped[i].motid,
			ped[i].sex);
		for (int j = 0; j < ped.markerCount; ++j) {
			printf("%d%d\t", ped[i].markers[j].one, ped[i].markers[j].two);
		}
		printf("\n");
	}
	//
	printf("Loaded %d markers\n", ped.markerCount);
	// Estimate allele frequencies for all markers, verbose mode
	ped.EstimateFrequencies(1, false);
	// Get genotype statistics for markers
	for (int i = 0; i < ped.markerNames.Length(); i++) {
		printf("Statistics for marker [%s]\n", (const char *)ped.markerNames[i]);
		// Allele index starts with 1 not 0
		for (int j = 1; j <= ped.GetMarkerInfo(i)->CountAlleles(); j++) {
			printf("\tFrequency for allele %d: %f\n", j, ped.GetMarkerInfo(i)->freq[j]);
			printf("\tName for allele %d: %s\n", j, (const char *)ped.GetMarkerInfo(i)->GetAlleleLabel(j));
		}
	}
	return;
}


void readData(Pedigree & ped,
              const char * datfile, const char * pedfile, const char * mapfile)
{
	// The data file contains a description of the contents of the
	// pedigree file, including for example, a list of marker and
	// trait names
	ped.Prepare(datfile);
	// The pedigree file contains a list of individuals, stored one
	// per row, with specific information about each individual as
	// detailed in the data file.
	ped.Load(pedfile);
	SortFamilies(ped);
	ped.LoadMarkerMap(mapfile);
	return;
}


int main(int argc, char ** argv)
{
	if (argc != 3 && argc != 4) {
		printf("usage: %s <data source code: 1, 2, 3> <task: 1 or 2> -v\n", argv[0]);
		return 0;
	}

	Pedigree ped;
	std::string chrom = "1";
	if (atoi(argv[1]) == 1) readData(ped, "haplo.dat", "haplo.ped", "haplo.map");
	else if (atoi(argv[1]) == 2) readData(ped, "gene.dat", "gene.ped", "gene.map");
	else if (atoi(argv[1]) == 3) {
		std::vector<std::string> marker_ids { "V1", "V2", "V3" };
		std::vector<int> marker_positions { 1, 2, 3 };
		std::vector< std::vector<std::string> > samples;
		std::vector<std::string> s0 { "1", "1", "0", "0", "1", "21", "21", "21" };
		samples.push_back(s0);
		std::vector<std::string> s1 { "1", "2", "0", "0", "2", "11", "11", "11" };
		samples.push_back(s1);
		std::vector<std::string> s2 { "1", "3", "1", "2", "1", "21", "21", "21" };
		samples.push_back(s2);
		std::vector<std::string> s3 { "2", "1", "0", "0", "1", "22", "21", "00" };
		samples.push_back(s3);
		std::vector<std::string> s4 { "2", "2", "0", "0", "2", "11", "11", "11" };
		samples.push_back(s4);
		std::vector<std::string> s5 { "2", "3", "1", "2", "1", "21", "21", "21" };
		samples.push_back(s5);
		std::vector<std::string> s6 { "3", "1", "0", "0", "1", "22", "21", "21" };
		samples.push_back(s6);
		std::vector<std::string> s7 { "3", "2", "0", "0", "2", "11", "11", "21" };
		samples.push_back(s7);
		std::vector<std::string> s8 { "3", "3", "1", "2", "1", "21", "21", "21" };
		samples.push_back(s8);
		//
		DataLoader dl;
		dl.LoadVariants(ped, marker_ids, marker_positions, chrom);
		dl.LoadSamples(ped, samples);
	} else ;
	if (atoi(argv[2]) == 1) showPed(ped);
	else if (atoi(argv[2]) == 2) {
		MendelianErrorChecker mc;
		mc.Apply(ped);
		std::cout << "Mendelian Errors " << mc.errorCount << std::endl;
		GeneticHaplotyper gh(chrom);
		gh.Apply(ped);
		if (argc == 4) gh.Print();
		HaplotypeCoder hc(1);
		hc.Apply(gh.data);
		if (argc == 4) hc.Print();
	}   else ;
}


