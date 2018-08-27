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


#include "FStatistics.hpp"
#include "EggException.hpp"
#include <cstdlib>


namespace egglib {


    void FStatistics::d_init() {
        d_reserved = 0;
        d_flag = false;
        d_numberOfGenotypes = 0;
        d_genotypes = NULL;
        d_populationLabels = NULL;
    }
    
    
    void FStatistics::d_clear() {
        if (d_genotypes) free(d_genotypes);
        if (d_populationLabels) free(d_populationLabels);
        d_init();
    }


    void FStatistics::s_init() {
        s_flag = false;
        s_numberOfAlleles = 0;
        s_alleleValueMapping = NULL;
        s_numberOfPopulations = 0;
        s_populationLabelMapping = NULL;
        
        s_populationFrequencies = NULL;
        s_alleleFrequenciesTotal = NULL;
        s_alleleFrequenciesPerPopulation = NULL;
        s_genotypeFrequenciesTotal = NULL;
        s_genotypeFrequenciesPerPopulation = NULL;
    }


    void FStatistics::s_clear() {

        if (s_alleleValueMapping) free(s_alleleValueMapping);
        if (s_populationLabelMapping) free(s_populationLabelMapping);
        if (s_populationFrequencies) free(s_populationFrequencies);
        if (s_alleleFrequenciesTotal) free(s_alleleFrequenciesTotal);
        for (unsigned int i=0; i<s_numberOfPopulations; i++) {
            if (s_alleleFrequenciesPerPopulation[i]) free(s_alleleFrequenciesPerPopulation[i]);
        }
        if (s_alleleFrequenciesPerPopulation) free(s_alleleFrequenciesPerPopulation);
        for (unsigned int i=0; i<s_numberOfAlleles; i++) {
            if (s_genotypeFrequenciesTotal[i]) free(s_genotypeFrequenciesTotal[i]);
        }
        if (s_genotypeFrequenciesTotal) free(s_genotypeFrequenciesTotal);
        for (unsigned int i=0; i<s_numberOfPopulations; i++) {
            for (unsigned int j=0; j<s_numberOfAlleles; j++) {
                if (s_genotypeFrequenciesPerPopulation[i][j]) {
                    free(s_genotypeFrequenciesPerPopulation[i][j]);
                }
            }
            if (s_genotypeFrequenciesPerPopulation[i]) {
                free(s_genotypeFrequenciesPerPopulation[i]);
            }
        }
        if (s_genotypeFrequenciesPerPopulation) {
            free(s_genotypeFrequenciesPerPopulation);
        }
        s_init();
    }


    void FStatistics::w_init() {
        w_flag = false;
        w_F = 0.;
        w_T = 0.;
        w_f = 0.;
        w_a = NULL;
        w_b = NULL;
        w_c = NULL;
        w_nbar = 0.;
        w_nc = 0.;
        w_pbar = NULL;
        w_ssquare = NULL;
        w_hbar = NULL;
        w_sum_a = 0.;
        w_sum_b = 0.;
        w_sum_c = 0.;
        w_sum_abc = 0.;
        w_sum_bc = 0.;
    }
    
    
    void FStatistics::w_clear() {
        if (w_a) free(w_a);
        if (w_b) free(w_b);
        if (w_c) free(w_c);
        if (w_pbar) free(w_pbar);
        if (w_ssquare) free(w_ssquare);
        if (w_hbar) free(w_hbar);
        w_init();
    }


/* ****************************************************************** */

    FStatistics::FStatistics() {
        d_init();
        s_init();
        w_init();
    }


    FStatistics::~FStatistics() {
        d_clear();
        s_clear();
        w_clear();
    }

    
/* ****************************************************************** */

    void FStatistics::reserve(unsigned int numberOfIndividuals) {
        if (numberOfIndividuals==0) {
            throw EggArgumentValueError("FStatistics: cannot reserve memory for zero sites for F-statistics analysis");
        }

        if (numberOfIndividuals<d_numberOfGenotypes) {
            throw EggArgumentValueError("FStatistics: cannot reserve less memory than what is currently used");
        }
        
        d_genotypes = (unsigned int*) realloc(d_genotypes, 2*numberOfIndividuals*sizeof(unsigned int));
        if (!d_genotypes) throw EggMemoryError();

        d_populationLabels = (unsigned int*) realloc(d_populationLabels, numberOfIndividuals*sizeof(unsigned int));
        if (!d_populationLabels) throw EggMemoryError();
            
        d_reserved = numberOfIndividuals;
    }
    

    void FStatistics::loadIndividual(unsigned int genotype1,
                    unsigned int genotype2, unsigned int populationLabel) {
    
    
        // clears statistics if already computed
        if (s_flag) s_clear();
        d_flag = true;

        // allocates memory
        d_numberOfGenotypes++;

        if (d_numberOfGenotypes>d_reserved) {
            d_genotypes = (unsigned int*) realloc(d_genotypes, 2*d_numberOfGenotypes*sizeof(unsigned int));
            if (!d_genotypes) throw EggMemoryError();

            d_populationLabels = (unsigned int*) realloc(d_populationLabels, d_numberOfGenotypes*sizeof(unsigned int));
            if (!d_populationLabels) throw EggMemoryError();
        }

        // loads data
        d_genotypes[2*(d_numberOfGenotypes-1)] = genotype1;
        d_genotypes[2*(d_numberOfGenotypes-1) + 1] = genotype2;
        d_populationLabels[d_numberOfGenotypes-1] = populationLabel;

    }


/* ****************************************************************** */

    void FStatistics::s_compute() {
        if (!d_flag) throw EggRuntimeError("FStatistics: cannot compute frequency statistics: data not loaded");
        if (s_flag) throw EggRuntimeError("FStatistics: inconsistency: frequency statistics already computed");

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
        s_genotypeFrequenciesPerPopulation = (unsigned int***) malloc(s_numberOfPopulations*sizeof(unsigned int**));
        if (!s_genotypeFrequenciesPerPopulation) {
            s_clear();
            throw EggMemoryError();
        }
        for (unsigned int i=0; i<s_numberOfPopulations; i++) {
            s_genotypeFrequenciesPerPopulation[i] = NULL;
        }

        if (s_numberOfPopulations<2) throw EggRuntimeError("FStatistics: cannot compute frequency statistics: not enough different populations");
        
        // gets alleles, mapping and frequencies
        processAlleles();
        
        s_flag = true;
    }
    
    
    void FStatistics::processPopulations() {
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
    
    
    void FStatistics::processAlleles() {
        for (unsigned int genotype=0; genotype<d_numberOfGenotypes; genotype++) {

            unsigned int populationIndex = getPopulationIndex(d_populationLabels[genotype]);

            // both allele (for genotype indexing)
            unsigned int alleleIndex1 = 0;
            unsigned int alleleIndex2 = 0;
            
            for (unsigned int offset=0; offset<2; offset++) {

                unsigned int alleleValue = d_genotypes[2*genotype+offset];
                unsigned int alleleIndex = getAlleleIndex(alleleValue);

                // stores both alleles (works also if alleleIndex==s_numberOfAlleles)
                if (offset==0) alleleIndex1 = alleleIndex;
                if (offset==1) alleleIndex2 = alleleIndex;
                
                // new alleles

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
                    
                    // ... for the genotype frequency tables
                    
                    // total table
                    
                    // new row
                    
                    s_genotypeFrequenciesTotal = (unsigned int**) realloc(s_genotypeFrequenciesTotal, s_numberOfAlleles*sizeof(unsigned int*));
                    if (!s_genotypeFrequenciesTotal) {
                        s_clear();
                        throw EggMemoryError();
                    }
                    s_genotypeFrequenciesTotal[alleleIndex] = NULL;
                    
                    // new column
                    for (unsigned int i=0; i<s_numberOfAlleles; i++) {
                        s_genotypeFrequenciesTotal[i] = (unsigned int*) realloc(s_genotypeFrequenciesTotal[i], s_numberOfAlleles*sizeof(unsigned int));
                        if (!s_genotypeFrequenciesTotal[i]) {
                            s_clear();
                            throw EggMemoryError();
                        }
                        s_genotypeFrequenciesTotal[i][alleleIndex] = 0;
                    }
                    
                    // initialization of the new line (except new column)
                    for (unsigned int i=0; i<(s_numberOfAlleles-1); i++) {
                        s_genotypeFrequenciesTotal[alleleIndex][i] = 0;
                    }
                    
                    // the meta-table
                    
                    for (unsigned int i=0; i<s_numberOfPopulations; i++) {

                        // adds the new row
                        s_genotypeFrequenciesPerPopulation[i] = (unsigned int**) realloc(s_genotypeFrequenciesPerPopulation[i], s_numberOfAlleles*sizeof(unsigned int*));
                        if (!s_genotypeFrequenciesPerPopulation[i]) {
                            s_clear();
                            throw EggMemoryError();
                        }
                        s_genotypeFrequenciesPerPopulation[i][alleleIndex] = NULL;

                        // adds the new column (one cell per row, incl. new)
                        for (unsigned int j=0; j<s_numberOfAlleles; j++) {
                            s_genotypeFrequenciesPerPopulation[i][j] = (unsigned int*) realloc(s_genotypeFrequenciesPerPopulation[i][j], s_numberOfAlleles*sizeof(unsigned int));
                            if (!s_genotypeFrequenciesPerPopulation[i][j]) {
                                s_clear();
                                throw EggMemoryError();
                            }
                            // initializes the new column
                            s_genotypeFrequenciesPerPopulation[i][j][alleleIndex] = 0;
                        }

                        // initializes the new row (but not the new column)
                        for (unsigned int j=0; j<(s_numberOfAlleles-1); j++) {
                            s_genotypeFrequenciesPerPopulation[i][alleleIndex][j] = 0;
                        }
                    }
                }
                
                // increments allele frequencies (even if new)
                s_alleleFrequenciesTotal[alleleIndex]++;
                s_alleleFrequenciesPerPopulation[populationIndex][alleleIndex]++;
            }
            
            // increments frequencies of genotype
            s_genotypeFrequenciesTotal[alleleIndex1][alleleIndex2]++;
            s_genotypeFrequenciesPerPopulation[populationIndex][alleleIndex1][alleleIndex2]++;
        }
    }


    unsigned int FStatistics::getAlleleIndex(unsigned int alleleValue) const {
        for (unsigned int allele=0; allele<s_numberOfAlleles; allele++) {
            if (s_alleleValueMapping[allele] == alleleValue) {
                return allele;
            }
        }
        return s_numberOfAlleles;
    }
    
    
    unsigned int FStatistics::getPopulationIndex(unsigned int populationLabel) const {
        for (unsigned int population=0; population<s_numberOfPopulations; population++) {
            if (s_populationLabelMapping[population] == populationLabel) {
                return population;
            }
        }
        return s_numberOfPopulations;
    }


/* ****************************************************************** */

    void FStatistics::w_compute() {

        /* implementation of Weir & Cockerham 1984 */
        
        if (!s_flag) s_compute();
        if (w_flag) throw EggRuntimeError("FStatistics: inconsistency: F-statistics already computed");
        
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
                w_pbar[allele] += 0.5 * s_alleleFrequenciesPerPopulation[population][allele];       // absolule frquency

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
        (0.5* s_alleleFrequenciesPerPopulation[population][allele] / s_populationFrequencies[population] - w_pbar[allele]) *
        (0.5* s_alleleFrequenciesPerPopulation[population][allele] / s_populationFrequencies[population] - w_pbar[allele]) );
            }
            w_ssquare[allele] /= ( (s_numberOfPopulations - 1.) * w_nbar );
        }
        
        // average heterozygosity
        
        w_hbar = (double*) malloc(s_numberOfAlleles * sizeof(double));
        if (!w_hbar) {
            w_clear();
            throw EggMemoryError();
        }

        for (unsigned int allele=0; allele<s_numberOfAlleles; allele++) {
            w_hbar[allele] = 0.;
            for (unsigned int population=0; population<s_numberOfPopulations; population++) {
                
                // computes the frequency of heterozygotes
                double H = 0.; // equivalent to ni*hi
                for (unsigned int otherAllele = 0; otherAllele<s_numberOfAlleles; otherAllele++) {
                    if (allele!=otherAllele) {
                        H += ((
                        s_genotypeFrequenciesPerPopulation[population][allele][otherAllele] +
                        s_genotypeFrequenciesPerPopulation[population][otherAllele][allele]));
                    }
                }
                
                w_hbar[allele] += H;
            }
            w_hbar[allele] /= ( s_numberOfPopulations * w_nbar );
        }
            

        // between-population component of variance
        
        w_a = (double*) malloc(s_numberOfAlleles * sizeof(double));
        if (!w_a) {
            w_clear();
            throw EggMemoryError();
        }

        /* equation (2) of WC84 */

        for (unsigned int allele=0; allele<s_numberOfAlleles; allele++) {
            w_a[allele] = (w_nbar/w_nc)*( 
            w_ssquare[allele] - (1./(w_nbar-1)) * (
                w_pbar[allele]*(1-w_pbar[allele])
                - w_ssquare[allele]*(s_numberOfPopulations-1.)/s_numberOfPopulations
                - w_hbar[allele] / 4.) );
        }
        
        // within-population, between individual component of variance
        
        w_b = (double*) malloc(s_numberOfAlleles * sizeof(double));
        if (!w_b) {
            w_clear();
            throw EggMemoryError();
        }

        /* equation (3) of WC84) */
        
        for (unsigned int allele=0; allele<s_numberOfAlleles; allele++) {
            w_b[allele] = (w_nbar/(w_nbar-1)) * ( w_pbar[allele]*(1-w_pbar[allele])
                - w_ssquare[allele] * ((s_numberOfPopulations-1.)/s_numberOfPopulations)
                - w_hbar[allele] * (2*w_nbar-1)/(4*w_nbar) );
        }

        // within individual component of variance
        
        w_c = (double*) malloc(s_numberOfAlleles * sizeof(double));
        
        if (!w_c) {
            w_clear();
            throw EggMemoryError();
        }
        
        for (unsigned int allele=0; allele<s_numberOfAlleles; allele++) {
            w_c[allele] = 0.5 * w_hbar[allele];
        }

        // sums (for multi-allele computation, also available for multi-locus)
        
        for (unsigned int allele=0; allele<s_numberOfAlleles; allele++) {
            w_sum_a += w_a[allele];
            w_sum_b += w_b[allele];
            w_sum_c += w_c[allele];
            w_sum_bc += (w_b[allele] + w_c[allele]);
            w_sum_abc += (w_a[allele] + w_b[allele] + w_c[allele]);
        }
        
        // F-statistics  WC equation (1) modified for multi allele
        
        w_F = 1 - w_sum_c / w_sum_abc;
        w_T = w_sum_a / w_sum_abc;
        w_f = 1 - w_sum_c / w_sum_bc;
        
        // conclusion
        
        w_flag = true;
        
    }

    unsigned int FStatistics::firstAllele(unsigned int individualIndex) const {
        return d_genotypes[2*individualIndex];
    }

    unsigned int FStatistics::secondAllele(unsigned int individualIndex) const {
        return d_genotypes[2*individualIndex+1];
    }

    unsigned int FStatistics::individualLabel(unsigned int individualIndex) const {
        return d_populationLabels[individualIndex];
    }


/* ****************************************************************** */

    unsigned int FStatistics::numberOfGenotypes() const {
        return d_numberOfGenotypes;
    }

    unsigned int FStatistics::numberOfPopulations() {
        if (!s_flag) s_compute();
        return s_numberOfPopulations;
    }
    
    
    unsigned int FStatistics::numberOfAlleles() {
        if (!s_flag) s_compute();
        return s_numberOfAlleles;
    }
    
    
    unsigned int FStatistics::populationLabel(unsigned int i) {
        if (!s_flag) s_compute();
        return s_populationLabelMapping[i];
    }


    unsigned int FStatistics::alleleValue(unsigned int i) {
        if (!s_flag) s_compute();
        return s_alleleValueMapping[i];
    }


    unsigned int FStatistics::alleleFrequencyTotal(unsigned int alleleIndex) {
        if (!s_flag) s_compute();
        return s_alleleFrequenciesTotal[alleleIndex];
    }

    
    unsigned int FStatistics::alleleFrequencyPerPopulation(unsigned int populationIndex, unsigned int alleleIndex) {
        if (!s_flag) s_compute();
        return s_alleleFrequenciesPerPopulation[populationIndex][alleleIndex];
    }

    
    unsigned int FStatistics::genotypeFrequencyTotal(unsigned int alleleIndex1, unsigned int alleleIndex2) {
        if (!s_flag) s_compute();
        return s_genotypeFrequenciesTotal[alleleIndex1][alleleIndex2];
    }

    
    unsigned int FStatistics::genotypeFrequencyPerPopulation(unsigned int populationIndex, unsigned int alleleIndex1, unsigned int alleleIndex2) {
        if (!s_flag) s_compute();
        return s_genotypeFrequenciesPerPopulation[populationIndex][alleleIndex1][alleleIndex2];
    }

    unsigned int FStatistics::populationFrequency(unsigned int populationIndex) {
        if (!s_flag) s_compute();
        return s_populationFrequencies[populationIndex];
    }


/* ****************************************************************** */

    double FStatistics::F() {
        if (!w_flag) w_compute();
        return w_F;
    }
    
    double FStatistics::theta() {
        if (!w_flag) w_compute();
        return w_T;
    }
    
    double FStatistics::f() {
        if (!w_flag) w_compute();
        return w_f;
    }
    
    double FStatistics::Vpopulation() {
        if (!w_flag) w_compute();
        return w_sum_a;
    }
    
    double FStatistics::Vindividual() {
        if (!w_flag) w_compute();
        return w_sum_b;
    }

    double FStatistics::Vallele() {
        if (!w_flag) w_compute();
        return w_sum_c;
    }
        
}

