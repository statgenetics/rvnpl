/*
    Copyright 2008-2010 Stéphane De Mita, Mathieu Siol

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

#include "MicrosatelliteDiversity.hpp"
#include "EggException.hpp"


namespace egglib {


    MicrosatelliteDiversity::MicrosatelliteDiversity() {
        init();
    }
    
    
    MicrosatelliteDiversity::~MicrosatelliteDiversity() {
        clear();
    }

    
/* ****************************************************************** */


    void MicrosatelliteDiversity::init() {
        v_numberOfSites = 0;
        v_He = NULL;
        v_numberOfAlleles = NULL;
        v_sizeVariance = NULL;
        v_thetaAssumingIAM = NULL;
        v_thetaAssumingSMMfromHe = NULL;
        v_thetaAssumingSMMfromSizeVariance = NULL;
    }


    void MicrosatelliteDiversity::clear() {
        free( v_He );
        free( v_numberOfAlleles );
        free( v_sizeVariance );
        free( v_thetaAssumingIAM );
        free( v_thetaAssumingSMMfromHe );
        free( v_thetaAssumingSMMfromSizeVariance );
    }


/* ****************************************************************** */

    void MicrosatelliteDiversity::load(const DataMatrix& dataMatrix,
            int missingData, bool noMissingData) {

        clear();
        init();
        if (dataMatrix.numberOfSites()==0) return;
        
        // allocates memory
        v_numberOfSites = dataMatrix.numberOfSites();

        v_He = (double*) malloc(v_numberOfSites * sizeof(double));
        if (!v_He) throw EggMemoryError();
        
        v_numberOfAlleles = (unsigned int*) malloc(v_numberOfSites * sizeof(unsigned int));
        if (!v_numberOfAlleles) throw EggMemoryError();

        v_sizeVariance = (double*) malloc(v_numberOfSites * sizeof(double));
        if (!v_sizeVariance) throw EggMemoryError();

        v_thetaAssumingIAM = (double*) malloc(v_numberOfSites * sizeof(double));
        if (!v_thetaAssumingIAM) throw EggMemoryError();

        v_thetaAssumingSMMfromHe = (double*) malloc(v_numberOfSites * sizeof(double));
        if (!v_thetaAssumingSMMfromHe) throw EggMemoryError();

        v_thetaAssumingSMMfromSizeVariance = (double*) malloc(v_numberOfSites * sizeof(double));
        if (!v_thetaAssumingSMMfromSizeVariance ) throw EggMemoryError();

        // process all sites
        for (unsigned int site=0; site<v_numberOfSites; site++) {

            v_numberOfAlleles[site] = 0;
            int* alleles = NULL;
            unsigned int* frequencies = NULL;
            unsigned int numberOfValidSamples = 0;

            // process all alleles
            for (unsigned int ind=0; ind<dataMatrix.numberOfSequences(); ind++) {

                // gets the characters
                int allele = dataMatrix.get(ind, site);
                if (!noMissingData && allele == missingData) continue;
                numberOfValidSamples++;

                unsigned int all;
                for (all=0; all<v_numberOfAlleles[site]; all++) {

                    // already known allele
                    if (allele==alleles[all]) {
                        frequencies[all]+=1;
                        break;
                    }
                }
                
                // new allele
                if (all==v_numberOfAlleles[site]) {
                    
                    v_numberOfAlleles[site]++;
                    
                    // reallocations
                    alleles = (int*) realloc(alleles, v_numberOfAlleles[site]*sizeof(int));
                    if (!alleles) {
                        if (frequencies) free(frequencies);
                        throw EggMemoryError();
                    }
                    
                    frequencies = (unsigned int*) realloc(frequencies, v_numberOfAlleles[site]*sizeof(unsigned int));
                    if (!frequencies) {
                        if (alleles) free(alleles);
                        throw EggMemoryError();
                    }       

                    // assignment
                    alleles[v_numberOfAlleles[site]-1] = allele;
                    frequencies[v_numberOfAlleles[site]-1] = 1;  
                }
            }
            
            // computes heterozygoty
            double sumfreqsq = 0.;
            for (unsigned int all=0; all<v_numberOfAlleles[site]; all++) {
                sumfreqsq+= 
                    ((1.*frequencies[all])/numberOfValidSamples) 
                *   ((1.*frequencies[all])/numberOfValidSamples);
            }
         
            v_He[site] = 1. - sumfreqsq;
            v_He[site]*= numberOfValidSamples/(numberOfValidSamples-1.);
            
            // variance nb repet
            double sum = 0.;
            double sumsq = 0.;
            for (unsigned int all=0; all<v_numberOfAlleles[site]; all++) {
                sum+= alleles[all]*(int)frequencies[all];
                sumsq+= (alleles[all]*alleles[all])*frequencies[all];
            }        
            sum /= numberOfValidSamples;
            sumsq /= numberOfValidSamples;

            // computes final statistics
            v_sizeVariance[site] = sumsq-sum*sum;
            v_thetaAssumingIAM[site] = v_He[site]/(1.-v_He[site]);
            v_thetaAssumingSMMfromHe[site] = 0.5*(1/((1.-v_He[site])*(1.-v_He[site]))-1.);
            v_thetaAssumingSMMfromSizeVariance[site] = 2.*v_sizeVariance[site];
            
            // final freeing
            free(frequencies);
            free(alleles);
        }
    }


/* ****************************************************************** */

    unsigned int MicrosatelliteDiversity::numberOfSites() const {
        return v_numberOfSites;
    }
            
    double MicrosatelliteDiversity::He(unsigned int siteIndex) const {
        return v_He[siteIndex];
    }
            
    unsigned int MicrosatelliteDiversity::numberOfAlleles(unsigned int siteIndex) const {
        return v_numberOfAlleles[siteIndex];
    }
            
    double MicrosatelliteDiversity::sizeVariance(unsigned int siteIndex) const {
        return v_sizeVariance[siteIndex];
    }
            
    double MicrosatelliteDiversity::thetaAssumingIAM(unsigned int siteIndex) const {
        return v_thetaAssumingIAM[siteIndex];
    }
            
    double MicrosatelliteDiversity::thetaAssumingSMMfromHe(unsigned int siteIndex) const {
        return v_thetaAssumingSMMfromHe[siteIndex];
    }

    double MicrosatelliteDiversity::thetaAssumingSMMfromSizeVariance(unsigned int siteIndex) const {
        return v_thetaAssumingSMMfromSizeVariance[siteIndex];
    }

}
