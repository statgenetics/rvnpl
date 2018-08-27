/*
    Copyright 2009 Stéphane De Mita, Mathieu Siol

    This file is part of the EggLib library.

    EggLib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    EggLib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with EggLib.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "HFStatistics.hpp"
#include "EggException.hpp"
#include <cstdlib>


namespace egglib {


    void HFStatistics::d_init() {
        d_reserved = 0;
        d_flag = false;
        d_numberOfGenotypes = 0;
        d_genotypes = NULL;
        d_populationLabels = NULL;
    }
    
    
    void HFStatistics::d_clear() {
        if (d_genotypes) free(d_genotypes);
        if (d_populationLabels) free(d_populationLabels);
        d_init();
    }


    void HFStatistics::s_init() {
        s_flag = false;
        s_numberOfAlleles = 0;
        s_alleleValueMapping = NULL;
        s_numberOfPopulations = 0;
        s_populationLabelMapping = NULL;
        
        s_populationFrequencies = NULL;
        s_alleleFrequenciesTotal = NULL;
        s_alleleFrequenciesPerPopulation = NULL;
    }


    void HFStatistics::s_clear() {

        if (s_alleleValueMapping) free(s_alleleValueMapping);
        if (s_populationLabelMapping) free(s_populationLabelMapping);
        if (s_populationFrequencies) free(s_populationFrequencies);
        if (s_alleleFrequenciesTotal) free(s_alleleFrequenciesTotal);
        for (unsigned int i=0; i<s_numberOfPopulations; i++) {
            if (s_alleleFrequenciesPerPopulation[i]) free(s_alleleFrequenciesPerPopulation[i]);
        }
        if (s_alleleFrequenciesPerPopulation) free(s_alleleFrequenciesPerPopulation);
        s_init();
    }


    void HFStatistics::w_init() {
        w_flag = false;
        w_T = 0.;
        w_T1 = NULL;
        w_T2 = NULL;
        w_nbar = 0.;
        w_nc = 0.;
        w_pbar = NULL;
        w_ssquare = NULL;
        w_sum_T1 = 0.;
        w_sum_T2 = 0.;
    }
    
    
    void HFStatistics::w_clear() {
        if (w_T1) free(w_T1);
        if (w_T2) free(w_T2);
        if (w_pbar) free(w_pbar);
        if (w_ssquare) free(w_ssquare);
        w_init();
    }


/* ****************************************************************** */

    HFStatistics::HFStatistics() {
        d_init();
        s_init();
        w_init();
    }


    HFStatistics::~HFStatistics() {
        d_clear();
        s_clear();
        w_clear();
    }

    
/* ****************************************************************** */

    void HFStatistics::reserve(unsigned int numberOfIndividuals) {
        if (numberOfIndividuals==0) {
            throw EggArgumentValueError("HFStatistics: cannot reserve memory for zero sites for F-statistics analysis");
        }

        if (numberOfIndividuals<d_numberOfGenotypes) {
            throw EggArgumentValueError("HFStatistics: cannot reserve less memory than what is currently used");
        }
        
        d_genotypes = (unsigned int*) realloc(d_genotypes, numberOfIndividuals*sizeof(unsigned int));
        if (!d_genotypes) throw EggMemoryError();

        d_populationLabels = (unsigned int*) realloc(d_populationLabels, numberOfIndividuals*sizeof(unsigned int));
        if (!d_populationLabels) throw EggMemoryError();
            
        d_reserved = numberOfIndividuals;
    }
    

    void HFStatistics::loadIndividual(unsigned int genotype, unsigned int populationLabel) {
    
        // clears statistics if already computed
        if (s_flag) s_clear();
        d_flag = true;

        // allocates memory
        d_numberOfGenotypes++;

        if (d_numberOfGenotypes>d_reserved) {
            d_genotypes = (unsigned int*) realloc(d_genotypes, d_numberOfGenotypes*sizeof(unsigned int));
            if (!d_genotypes) throw EggMemoryError();

            d_populationLabels = (unsigned int*) realloc(d_populationLabels, d_numberOfGenotypes*sizeof(unsigned int));
            if (!d_populationLabels) throw EggMemoryError();
        }

        // loads data
        d_genotypes[d_numberOfGenotypes-1] = genotype;
        d_populationLabels[d_numberOfGenotypes-1] = populationLabel;

    }


/* ****************************************************************** */

    void HFStatistics::s_compute() {
        if (!d_flag) throw EggRuntimeError("HFStatistics: cannot compute frequency statistics: data not loaded");
        if (s_flag) throw EggRuntimeError("HFStatistics: inconsistency: frequency statistics already computed");

        // gets population data (number, mapping and frequencies)
        processPopulations();

        // makes memory allocation accordingly
        s_alleleFrequenciesPerPopulation = (unsigned int**) malloc(s_numberOfPopulations*sizeof(unsigned int*));
        if (!s_alleleFrequenciesPerPopulation) {
            s_clear();
            throw EggMemoryError();
        }
        for (unsigned int i=0; i<s_numberOfPopulations; i++) {
            s_alleleFrequenciesPerPopulation[i] = NULL;
        }
        if (s_numberOfPopulations<2) throw EggRuntimeError("HFStatistics: cannot compute frequency statistics: not enough different populations");
        
        // gets alleles, mapping and frequencies
        processAlleles();
        
        s_flag = true;
    }
    
    
    void HFStatistics::processPopulations() {
        for (unsigned int genotype=0; genotype<d_numberOfGenotypes; genotype++) {
            unsigned int populationLabel = d_populationLabels[genotype];
            unsigned int populationIndex = getPopulationIndex(populationLabel);
            
            // new population
            if (populationIndex==s_numberOfPopulations) {
                
                s_numberOfPopulations++;

                // memory
                s_populationLabelMapping = (unsigned int*) realloc(s_populationLabelMapping, s_numberOfPopulations*sizeof(unsigned int));
                if (!s_populationLabelMapping) {
                    s_clear();
                    throw EggMemoryError();
                }
                s_populationFrequencies = (unsigned int*) realloc(s_populationFrequencies, s_numberOfPopulations*sizeof(unsigned int));
                if (!s_populationFrequencies) {
                    s_clear();
                    throw EggMemoryError();
                }
                
                // data
                s_populationLabelMapping[populationIndex] = populationLabel;
                s_populationFrequencies[populationIndex] = 1;
            }
            
            // one of previous populations
            else {
                s_populationFrequencies[populationIndex]++;
            }
        }
    }
    
    
    void HFStatistics::processAlleles() {
        for (unsigned int genotype=0; genotype<d_numberOfGenotypes; genotype++) {

            unsigned int populationIndex = getPopulationIndex(d_populationLabels[genotype]);

            unsigned int alleleValue = d_genotypes[genotype];
            unsigned int alleleIndex = getAlleleIndex(alleleValue);

            // new allele

            if (alleleIndex == s_numberOfAlleles) {

                s_numberOfAlleles++;

                // adds the allele
                s_alleleValueMapping = (unsigned int*) realloc(s_alleleValueMapping, s_numberOfAlleles*sizeof(unsigned int));
                if (!s_alleleValueMapping) throw EggMemoryError();
                s_alleleValueMapping[alleleIndex] = alleleValue;
                    
                // adds the frequency

                // memory and initialization of new cells
                    
                // ... for the allele frequency tables
                s_alleleFrequenciesTotal = (unsigned int*) realloc(s_alleleFrequenciesTotal, s_numberOfAlleles*sizeof(unsigned int));
                if (!s_alleleFrequenciesTotal) {
                    s_clear();
                    throw EggMemoryError();
                }
                
                s_alleleFrequenciesTotal[alleleIndex] = 0;
                    
                for (unsigned int i=0; i<s_numberOfPopulations; i++) {
                    s_alleleFrequenciesPerPopulation[i] = (unsigned int*) realloc(s_alleleFrequenciesPerPopulation[i], s_numberOfAlleles*sizeof(unsigned int));
                    if (!s_alleleFrequenciesPerPopulation[i]) {
                        s_clear();
                        throw EggMemoryError();
                    }
                    s_alleleFrequenciesPerPopulation[i][alleleIndex] = 0;
                }
                    
            }

            // increments allele frequencies (even if new)
            s_alleleFrequenciesTotal[alleleIndex]++;
            s_alleleFrequenciesPerPopulation[populationIndex][alleleIndex]++;
            
        }
    }


    unsigned int HFStatistics::getAlleleIndex(unsigned int alleleValue) const {
        for (unsigned int allele=0; allele<s_numberOfAlleles; allele++) {
            if (s_alleleValueMapping[allele] == alleleValue) {
                return allele;
            }
        }
        return s_numberOfAlleles;
    }
    
    
    unsigned int HFStatistics::getPopulationIndex(unsigned int populationLabel) const {
        for (unsigned int population=0; population<s_numberOfPopulations; population++) {
            if (s_populationLabelMapping[population] == populationLabel) {
                return population;
            }
        }
        return s_numberOfPopulations;
    }


/* ****************************************************************** */

    void HFStatistics::w_compute() {

        /* implementation of Weir & Cockerham 1984 */
        
        if (!s_flag) s_compute();
        if (w_flag) throw EggRuntimeError("HFStatistics: inconsistency: F-statistics already computed");
        
        // average sample size (per population)
        
        w_nbar = 0.;
        unsigned int sum_nisquare = 0;
        for (unsigned int i=0; i<s_numberOfPopulations; i++) {
            sum_nisquare += s_populationFrequencies[i] * s_populationFrequencies[i];
            w_nbar += s_populationFrequencies[i];
        }
        w_nbar /= s_numberOfPopulations;
        
        // sample size variation
        
        w_nc = (s_numberOfPopulations * w_nbar - sum_nisquare / (s_numberOfPopulations * w_nbar))
                                        / (s_numberOfPopulations - 1.);


        // average allele frequency

        w_pbar = (double*) malloc(s_numberOfAlleles * sizeof(double));
        if (!w_pbar) {
            w_clear();
            throw EggMemoryError();
        }

        for (unsigned int allele=0; allele<s_numberOfAlleles; allele++) {
            w_pbar[allele] = 0.;
            for (unsigned int population=0; population<s_numberOfPopulations; population++) {
                w_pbar[allele] += s_alleleFrequenciesPerPopulation[population][allele];       // absolule frquency

            }
            w_pbar[allele] /= ( s_numberOfPopulations * w_nbar );
        }

        // allele frequency variation

        w_ssquare = (double*) malloc(s_numberOfAlleles * sizeof(double));
        if (!w_ssquare) {
            w_clear();
            throw EggMemoryError();
        }
        
        for (unsigned int allele=0; allele<s_numberOfAlleles; allele++) {
            w_ssquare[allele] = 0.;
            for (unsigned int population=0; population<s_numberOfPopulations; population++) {
                w_ssquare[allele] += ( s_populationFrequencies[population] *
        (1.*s_alleleFrequenciesPerPopulation[population][allele] / s_populationFrequencies[population] - w_pbar[allele]) *
        (1.*s_alleleFrequenciesPerPopulation[population][allele] / s_populationFrequencies[population] - w_pbar[allele]) );
            }
            w_ssquare[allele] /= ( (s_numberOfPopulations - 1.) * w_nbar );
        }
        
        // between-population component of variance
        
        w_T1 = (double*) malloc(s_numberOfAlleles * sizeof(double));
        if (!w_T1) {
            w_clear();
            throw EggMemoryError();
        }

        for (unsigned int allele=0; allele<s_numberOfAlleles; allele++) {

            w_T1[allele] = //(w_nbar/w_nc)* () ???
                w_ssquare[allele] - (1./(w_nbar-1)) * (
                w_pbar[allele]*(1-w_pbar[allele])
                - w_ssquare[allele]*(s_numberOfPopulations-1.)/s_numberOfPopulations ) ;
        }
        
        // within-population, between individual component of variance
        
        w_T2 = (double*) malloc(s_numberOfAlleles * sizeof(double));
        if (!w_T2) {
            w_clear();
            throw EggMemoryError();
        }

        // total variance
        
        for (unsigned int allele=0; allele<s_numberOfAlleles; allele++) {
            w_T2[allele] = ((w_nc-1)/(w_nbar-1)) * ( w_pbar[allele]*(1-w_pbar[allele])
                + (w_ssquare[allele]/s_numberOfPopulations) * 
                    (1+((s_numberOfPopulations-1)*(w_nbar-w_nc)/(w_nbar-1))) );
        }

        
        // F-statistics  equation modified for multi allele

        for (unsigned int allele=0; allele<s_numberOfAlleles; allele++) {
            w_sum_T1 += w_T1[allele];
            w_sum_T2 += w_T2[allele];
        }
        
        w_T = w_sum_T1 / w_sum_T2;
        
        // conclusion
        
        w_flag = true;
        
    }

    unsigned int HFStatistics::allele(unsigned int individualIndex) const {
        return d_genotypes[individualIndex];
    }

    unsigned int HFStatistics::individualLabel(unsigned int individualIndex) const {
        return d_populationLabels[individualIndex];
    }


/* ****************************************************************** */

    unsigned int HFStatistics::numberOfGenotypes() const {
        return d_numberOfGenotypes;
    }

    unsigned int HFStatistics::numberOfPopulations() {
        if (!s_flag) s_compute();
        return s_numberOfPopulations;
    }
    
    
    unsigned int HFStatistics::numberOfAlleles() {
        if (!s_flag) s_compute();
        return s_numberOfAlleles;
    }
    
    
    unsigned int HFStatistics::populationLabel(unsigned int i) {
        if (!s_flag) s_compute();
        return s_populationLabelMapping[i];
    }


    unsigned int HFStatistics::alleleValue(unsigned int i) {
        if (!s_flag) s_compute();
        return s_alleleValueMapping[i];
    }


    unsigned int HFStatistics::alleleFrequencyTotal(unsigned int alleleIndex) {
        if (!s_flag) s_compute();
        return s_alleleFrequenciesTotal[alleleIndex];
    }

    
    unsigned int HFStatistics::alleleFrequencyPerPopulation(unsigned int populationIndex, unsigned int alleleIndex) {
        if (!s_flag) s_compute();
        return s_alleleFrequenciesPerPopulation[populationIndex][alleleIndex];
    }

    
    unsigned int HFStatistics::populationFrequency(unsigned int populationIndex) {
        if (!s_flag) s_compute();
        return s_populationFrequencies[populationIndex];
    }


/* ****************************************************************** */

    double HFStatistics::theta() {
        if (!w_flag) w_compute();
        return w_T;
    }
    
    double HFStatistics::T1() {
        if (!w_flag) w_compute();
        return w_sum_T1;
    }
    
    double HFStatistics::T2() {
        if (!w_flag) w_compute();
        return w_sum_T2;
    }
       
}

