/*
    Copyright 2008,2009,2011 Stéphane De Mita, Mathieu Siol

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


#include "HaplotypeDiversity.hpp"
#include "EggException.hpp"
#include <cstdlib>

namespace egglib {

    // helps

    void HaplotypeDiversity::init() {

        BaseDiversity::init();

        m_loaded = false;
        m_K = 0;
        m_He = 0.;
        m_Kst = 0.;
        m_Fst = 0.;
        m_Gst = 0.;
        m_Hst = 0.;
        m_Snn = 0.;
        m_haplotypeIndex = NULL;
    }


    void HaplotypeDiversity::clear() {
        BaseDiversity::clear();
        if (m_haplotypeIndex) {
            free(m_haplotypeIndex);
        }
    }


    // constructor

    HaplotypeDiversity::HaplotypeDiversity()  {
        init();
    }


    HaplotypeDiversity::~HaplotypeDiversity() {
        clear();
    }


    // main method

    void HaplotypeDiversity::load(CharMatrix& data, bool allowMultipleMutations,
           unsigned int ignoreFrequency, std::string characterMapping) {

        // initialization
        clear();
        init();
        m_loaded = true;
        // create the site vector
        importSites(data, allowMultipleMutations, 1, ignoreFrequency,
                                        characterMapping, false, true);

        m_haplotypeIndex = (unsigned int*) malloc(v_ns * sizeof(unsigned int));
        if (!m_haplotypeIndex) throw EggMemoryError();
        for (unsigned int i=0; i<v_ns; i++) {
            m_haplotypeIndex[i] = 0;
        }

        if (!v_S) {
            m_K = 1.;
            return;
        }
        if (!v_npop) throw EggRuntimeError("unexpected situation: polymorphic sites but not populations");


        // initializes tables

        // size of each population
        unsigned int* ni = (unsigned int*) malloc(v_npop*sizeof(unsigned int));
        if (!ni) throw EggMemoryError();
        for (unsigned int i=0; i<v_npop; i++) {
            ni[i] = 0;
        }

        // marginal frequencies of alleles
        unsigned int* nj = NULL;

        // allele frequencies in populations
        unsigned int** nij = (unsigned int**) malloc(v_npop*sizeof(unsigned int*));
        if (!nij) throw EggMemoryError();
        for (unsigned int i=0; i<v_npop; i++) {
            nij[i] = NULL;
        }

        // index of sequence representing each allele
        unsigned int* haplotypeRef = NULL;

        // population heterozygosy
        double* Hi = NULL;

        // number of pairwise differences in population
        double* Ki = NULL;

        // population indices
        unsigned int* groups = (unsigned int*) malloc(v_ns * sizeof(unsigned int));
        if (!groups) throw EggMemoryError();

        // runs through the sequences skipping outgroups
        unsigned int c=0;
        for (unsigned int individual=0; individual<data.numberOfSequences(); individual++) {

            unsigned int popLabel = data.populationLabel(individual);
            if (popLabel==999) {
                continue;
            }

            // identify pop, count 1 individual in this pop
            unsigned int popIndex = getPopIndex(popLabel);
            if (popIndex==v_npop) {
                throw EggRuntimeError("inconsistency haplotype analysis");
            }
            ni[popIndex]++;
            groups[c] = popIndex;

            // finds allele in the haplotype list
            // checks if already in allele list
            unsigned int h;
            for (h=0; h<m_K; h++) {

                if (!diff(data, individual, haplotypeRef[h])) {
                    // already known allele
                    nj[h]++;
                    nij[popIndex][h]++;
                    break;
                }
            }

            // new allele
            if (h==m_K) {
                m_K++;

                // adds the new haplotype
                haplotypeRef = (unsigned int*) realloc(haplotypeRef, m_K*sizeof(unsigned int));
                if (!haplotypeRef) throw EggMemoryError();
                haplotypeRef[h] = individual;

                // adds the new marginal haplotype frequency
                nj = (unsigned int*) realloc(nj, m_K*sizeof(unsigned int));
                if (!nj) throw EggMemoryError();
                nj[h] = 1;

                // adds the new frequency
                for (unsigned int i=0; i<v_npop; i++) {
                    nij[i] = (unsigned int*) realloc(nij[i], m_K*sizeof(unsigned int));
                    if (!nij[i]) throw EggMemoryError();
                    nij[i][h]= 0;
                }
                nij[popIndex][h] = 1;
            }

            m_haplotypeIndex[c]= h;
            c++;
    
        }

        if (haplotypeRef) free(haplotypeRef);

        // computations of various stats

        // computes He (alias Ht total heterozygosy)
        for (unsigned int i=0; i<m_K; i++) {
            m_He+= ((double) nj[i]/v_ns)*((double) nj[i]/v_ns);
        }
        double Httilde = (1- m_He); // stored for use by Gst
        if (m_K==1) m_He= 0.;
        else m_He = Httilde * v_ns/(v_ns-1.);

        // computes Fst estimator, the Hudson et al.'s way (Hst)
        Hi = (double*) malloc(v_npop * sizeof(double));
        if (!Hi) throw EggMemoryError();
        for (unsigned int pop=0; pop<v_npop; pop++) {
            Hi[pop] = 0.;
            for (unsigned int haplotype=0; haplotype<m_K; haplotype++) {
                Hi[pop]+= ((double)nij[pop][haplotype]/ni[pop])*
                            ((double)nij[pop][haplotype]/ni[pop]);
// beware, small sized population generate NAN around here
            }
            Hi[pop] = (ni[pop]/(ni[pop]-1.))*(1-Hi[pop]);
        }

        // computes Hst itself
        double Hs = 0.;
        for (unsigned int pop=0; pop<v_npop; pop++) {
            Hs += (double)Hi[pop] * (double)ni[pop]/v_ns;
        }
        if (m_He) m_Hst = 1-Hs/m_He;

        // the Nei's way (Gst)
        double ntilde = 0.;
        for (unsigned int pop=0; pop<v_npop; pop++) {
            if (ni[pop]) ntilde+=1./ni[pop];
        }
        if (ntilde) ntilde= v_npop*v_npop/ntilde;
        Httilde += Hs/ntilde;
        if (Httilde && m_K!=1) m_Gst = 1-Hs/Httilde;

        // Snn, Fst and Kst (needs iteration over pairs)

        // initialization
        // number of pairwise differences in each population
        Ki = (double*) malloc(v_npop * sizeof(double));
        if (!Ki) throw EggMemoryError();
        for (unsigned int i=0; i<v_npop; i++) Ki[i] = 0.;
        double Kt = 0.;
        double Hw= 0.;  // will be: mean number of pairwise differences, within subpopulations
        double Hb= 0;   // will be: mean number of pairwise differences, between subpopulations
        unsigned int nw=0, nb=0;  // counters (for convenience)

        // iteration over all possible pairs
        if (v_nseff>=2) {

            unsigned int c1=0;
            for (unsigned int ind1=0; ind1<data.numberOfSequences(); ind1++) {
                if (data.populationLabel(ind1)==999) continue;    
                unsigned int min=0;
                unsigned int nni=0; // nearest neighbors from the same group
                unsigned int nnt=0; // nearest neighbors (total)
                unsigned int dist;
                unsigned int c2=0;
                for (unsigned int ind2=0; ind2<data.numberOfSequences(); ind2++) {
                    if (ind1==ind2) {
                        c2++;
                        continue;
                    }
                    if (data.populationLabel(ind2)==999) continue;
                    dist= diff(data, ind1, ind2);

                    // increments number of pairwise differences counters
                    if (ind2>ind1) {
                        unsigned int pop1= groups[c1];
                        unsigned int pop2= groups[c2];
                        if (pop1==pop2) {
                            Ki[pop1]+= dist;
                            Hw+= dist;
                            nw++;
                        }
                        else {
                            Hb+= dist;
                            nb++;
                        }
                        Kt+= dist;
                    }

                    if (nnt==0 || dist<min) {
                        min= dist;
                        nnt=1;
                        nni= (groups[c1]==groups[c2]);
                        c2++;
                        continue;
                    }
                    if (dist==min) {
                        nnt++;
                        nni+= (groups[c1]==groups[c2]);
                    }
                    c2++;
                }
                if (nnt) m_Snn+= (double)nni/nnt;
                c1++;
            }
            if (v_nseff) m_Snn/=v_nseff;

            // normalisation of all K's
            for (unsigned int pop=0; pop<v_npop; pop++) {
                if (ni[pop]) {
                    Ki[pop]= 2.*Ki[pop]/(ni[pop]*(ni[pop]-1));
                }
            }
            Kt = 2.*Kt/(v_ns*(v_ns-1));

            // Ks computation
            double Ks= 0.;
            for (unsigned int i=0; i<v_npop; i++) {
                Ks+= ni[i]*Ki[i];
            }
            if (v_ns) Ks/=v_ns;

            // and Kst itself
            if (Kt) m_Kst = 1-Ks/Kt;

            // and finally Fst
            if (nw && nb && Hb) m_Fst = 1- (Hw/nw)/(Hb/nb);

        } // end computation of Snn, Fst, Kst
                    
        // cleans memory
        free(ni);
        free(nj);
        for (unsigned int i=0; i<v_npop; i++) {
            free(nij[i]);
        }
        free(nij);
        free(Hi);
        free(Ki);
        free(groups);
    
    }
            
            
          
    
            
    // accessors
    
    unsigned int HaplotypeDiversity::K() const {
        if (!m_loaded) throw EggRuntimeError("cannot access HaplotypeDiversity statistics: data not loaded");
        return m_K;
    }
    
    double HaplotypeDiversity::He() const {
        if (!m_loaded) throw EggRuntimeError("cannot access HaplotypeDiversity statistics: data not loaded");
        return m_He;        
    }
    
    unsigned int HaplotypeDiversity::haplotypeIndex(unsigned int index) const {
        if (!m_loaded) throw EggRuntimeError("cannot access HaplotypeDiversity statistics: data not loaded");
        if (index>=v_ns) throw EggRuntimeError("cannot access HaplotypeDiversity statistics - invalid haplotype index");
        return m_haplotypeIndex[index];
    }
    
    double HaplotypeDiversity::Kst() const {
        if (!m_loaded) throw EggRuntimeError("cannot access HaplotypeDiversity statistics: data not loaded");
        return m_Kst;
    }
    
    double HaplotypeDiversity::Fst() const {
        if (!m_loaded) throw EggRuntimeError("cannot access HaplotypeDiversity statistics: data not loaded");
        return m_Fst;        
    }
    
    double HaplotypeDiversity::Gst() const {
        if (!m_loaded) throw EggRuntimeError("cannot access HaplotypeDiversity statistics: data not loaded");
        return m_Gst;
    }
    
    double HaplotypeDiversity::Hst() const {
        if (!m_loaded) throw EggRuntimeError("cannot access HaplotypeDiversity statistics: data not loaded");
        return m_Hst;
    }
    
    double HaplotypeDiversity::Snn() const {
        if (!m_loaded) throw EggRuntimeError("cannot access HaplotypeDiversity statistics: data not loaded");
        return m_Snn;
    }

    inline unsigned int HaplotypeDiversity::diff(CharMatrix& data, unsigned int ind1, unsigned int ind2) const {

        unsigned int numberOfDifferences = 0;
        for (unsigned int site=0; site<v_S; site++) {
            char c= data.character(ind1, v_sitePositions[site]);
            char d= data.character(ind2, v_sitePositions[site]);
            if (p_characterMapping.find(c) < p_pos_sep_mapping &&
                p_characterMapping.find(d) < p_pos_sep_mapping && c!=d)
            {
                numberOfDifferences+=1;
            }
        }
        return numberOfDifferences;
    }

}

