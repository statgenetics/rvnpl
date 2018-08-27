/*
    Copyright 2008,2009,2012 Stéphane De Mita, Mathieu Siol

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


#include "SitePolymorphism.hpp"
#include "EggException.hpp"
#include <cstdlib>



namespace egglib {


    SitePolymorphism::SitePolymorphism() {
        init();
    }

    
    SitePolymorphism::~SitePolymorphism() {
        clear();
    }
    

    SitePolymorphism::SitePolymorphism(unsigned int npop) {
        init();
        numberOfPopulations(npop);
    }


    SitePolymorphism::SitePolymorphism(const SitePolymorphism& source) {
        init();
        copy(source);
    }


    SitePolymorphism& SitePolymorphism::operator=(const SitePolymorphism& source) {
        init();
        copy(source);
        return *this;
    }

    
    void SitePolymorphism::init() {
        m_numberOfPopulations = 0;
        m_numberOfStates = 0;
        m_states = NULL;
        m_frequencies = NULL;
        m_ns = 0;
        m_pop_ns = NULL;
        m_numberOfOutgroups = 0;
        m_outgroups = NULL;
    }
    
    
    void SitePolymorphism::copy(const SitePolymorphism& source) {
        m_numberOfPopulations = source.m_numberOfPopulations;
        m_numberOfStates = source.m_numberOfStates;
        
        m_frequencies = (unsigned int**) malloc(m_numberOfPopulations*sizeof(unsigned int*));
        if (!m_frequencies) throw EggMemoryError();
        for (unsigned int i=0; i<m_numberOfPopulations; i++) {
            m_frequencies[i] = (unsigned int*) malloc(m_numberOfStates*sizeof(unsigned int));
            if (!m_frequencies[i]) throw EggMemoryError();
            for (unsigned int j=0; j<m_numberOfStates; j++) {
                m_frequencies[i][j] = source.m_frequencies[i][j];
            }
        }
        m_states = (char*) malloc(m_numberOfStates*sizeof(char));
        if (!m_states) throw EggMemoryError();
        for (unsigned int i=0; i<m_numberOfStates; i++) {
            m_states[i] = source.m_states[i];
        }
        
        m_numberOfOutgroups = source.m_numberOfOutgroups;
        if (m_numberOfOutgroups) {
            m_outgroups = (char*) malloc(m_numberOfOutgroups*sizeof(char));
            if (!m_outgroups) throw EggMemoryError();
                for (unsigned int i=0; i<m_numberOfOutgroups; i++) {
                m_outgroups[i] = source.m_outgroups[i];
            }
        }

        m_ns = source.m_ns;
        
        m_pop_ns = (unsigned int*) malloc(m_numberOfPopulations*sizeof(unsigned int));
        if (!m_pop_ns) throw EggMemoryError();
        for (unsigned int i=0; i<m_numberOfPopulations; i++) {
            m_pop_ns[i] = source.m_pop_ns[i];
        }
    }


    void SitePolymorphism::clear() {
        if ((m_numberOfPopulations && !m_frequencies) ||
            (!m_numberOfPopulations && m_frequencies)) {
                throw EggRuntimeError("memory management error in SitePolymorphism");
        }
        if ((m_numberOfStates && !m_states) ||
            (!m_numberOfStates && m_states)) {
                throw EggRuntimeError("memory management error in SitePolymorphism");
        }
        if (m_numberOfPopulations) {
            if (m_numberOfStates) {
                for (unsigned int i=0; i<m_numberOfPopulations; i++) {
                    if (!m_frequencies[i]) {
                        throw EggRuntimeError("memory management error in SitePolymorphism");
                    }
                    free(m_frequencies[i]);
                }
            }
            free(m_frequencies);
        }
        if (m_numberOfStates) {
            free(m_states);
        }
        if (m_numberOfOutgroups) {
            free(m_outgroups);
        }
        if (m_pop_ns) {
            free(m_pop_ns);
        }
    }
   


    void SitePolymorphism::numberOfPopulations(unsigned int npop) {
        clear();
        init();
        if (npop==0) throw EggArgumentValueError("a site cannot have zero populations");
        m_numberOfPopulations = npop;
        m_frequencies = (unsigned int**) malloc(m_numberOfPopulations * sizeof(unsigned int*));
        if (!m_frequencies) throw EggMemoryError();
        for(unsigned int i=0; i<m_numberOfPopulations; i++) {
            m_frequencies[i] = NULL;
        }
        m_pop_ns = (unsigned int*) malloc(m_numberOfPopulations * sizeof(unsigned int));
        if (!m_pop_ns) throw EggMemoryError();
        for (unsigned int i=0; i<m_numberOfPopulations; i++) {
            m_pop_ns[i] = 0;
        }
    }


    void SitePolymorphism::load(unsigned int populationIndex, char character) {
        m_ns++;
        m_pop_ns[populationIndex]++;

        // checks whether allele in list
        for (unsigned int i=0; i<m_numberOfStates; i++) {
            if (m_states[i]==character) {
                m_frequencies[populationIndex][i]++;
                return;
            }
        }
        
        // adds a new allele
        m_numberOfStates++;
        m_states = (char*) realloc(m_states, m_numberOfStates*sizeof(char));
        if (!m_states) throw EggMemoryError();
        m_states[m_numberOfStates-1] = character;
        for (unsigned int i=0; i<m_numberOfPopulations; i++) {
            m_frequencies[i] = (unsigned int*) realloc(m_frequencies[i], m_numberOfStates*sizeof(unsigned int));
            if (!m_frequencies[i]) throw EggMemoryError();
            m_frequencies[i][m_numberOfStates-1] = (i==populationIndex)?1:0;
        }



    }
    
    
    void SitePolymorphism::outgroup(char character) {
        m_numberOfOutgroups++;
        m_outgroups = (char*) realloc(m_outgroups, m_numberOfOutgroups*sizeof(char));
        if (!m_outgroups) throw EggMemoryError();
        m_outgroups[m_numberOfOutgroups-1] = character;
    }
    
    
    char SitePolymorphism::allele(unsigned int index) const {
        return m_states[index];
    }

    
    unsigned int SitePolymorphism::numberOfAlleles() const {
        return m_numberOfStates;
    }
    
    
    unsigned int SitePolymorphism::alleleFrequency(unsigned int alleleIndex) const {
        unsigned int p = 0;
        for (unsigned int i=0; i<m_numberOfPopulations; i++) {
            p+=m_frequencies[i][alleleIndex];
        }
        return p;
    }


    unsigned int SitePolymorphism::alleleFrequency(unsigned int popIndex, unsigned int alleleIndex) const {
        return m_frequencies[popIndex][alleleIndex];
    }


    unsigned int SitePolymorphism::derivedAlleleFrequency() const {
        unsigned int p = 0;
        for (unsigned i=0; i<m_numberOfStates; i++) {
            if (m_states[i]!=m_outgroups[0]) {
                p+=alleleFrequency(i);
            }
        }
        return p;
    }


    unsigned int SitePolymorphism::ns() const {
        return m_ns;
    }


    unsigned int SitePolymorphism::ns(unsigned int popIndex) const {
        return m_pop_ns[popIndex];
    }


    bool SitePolymorphism::isOrientable() const {
        // at least one outgroup
        if (m_numberOfOutgroups==0) return false;
        
        // all outgroups the same
        for (unsigned int i=1; i<m_numberOfOutgroups; i++) {
            if (m_outgroups[i]!=m_outgroups[0]) return false;
        }

        // GOOD IF EITHER
        //     outgroup is one the allele present
        for (unsigned int i=0; i<m_numberOfStates; i++) {
            if (m_outgroups[0] == m_states[i]) return true;
        }

        // OR
        //     there is no polymorphism (2 alleles overall)
        if (numberOfAlleles()==1) return true;

        return false;
    }


    bool SitePolymorphism::isPolymorphic(unsigned int popIndex) const {
        bool flag = false;
        for (unsigned int state=0; state<m_numberOfStates; state++) {
            if (m_frequencies[popIndex][state]>0) {
                if (!flag) {
                    flag=true;
                }
                else {
                    return true;  // true if at least two states have freq >=1
                }
            }
        }
        return false;  // false otherwise
    }


    bool SitePolymorphism::hasSpecificAllele(unsigned int popIndex, bool restrictToDerived) const {
        for (unsigned int state=0; state<m_numberOfStates; state++) {
            if (restrictToDerived && m_states[state]==m_outgroups[0]) continue;
            unsigned int pop=0;
            for (pop=0; pop<m_numberOfPopulations; pop++) {
                    if ((pop==popIndex && m_frequencies[pop][state]==0) 
                     || (pop!=popIndex && m_frequencies[pop][state]!=0))  {
                    break;
                }
            }
            if (pop==m_numberOfPopulations) return true;
        }
        return false;
    }


    bool SitePolymorphism::haveFixedDifference(unsigned int pop1, unsigned int pop2) const {
        if (isPolymorphic(pop1)) return false;
        if (isPolymorphic(pop2)) return false;
        unsigned int allele1=0, allele2=0;
        for (unsigned int allele=0; allele<m_numberOfStates; allele++) {
            if (m_frequencies[pop1][allele]>0) allele1=allele;
            if (m_frequencies[pop2][allele]>0) allele2=allele;
            if (m_frequencies[pop1][allele]>0 && m_frequencies[pop2][allele]>0) break;
        }
        if (allele1!=allele2) return true;
        return false;
    }


    bool SitePolymorphism::haveCommonAllele(unsigned int pop1, unsigned int pop2) const {
        for (unsigned int allele=0; allele<m_numberOfStates; allele++) {
            if (m_frequencies[pop1][allele]>0 && m_frequencies[pop2][allele]>0) {
                return true;
            }
        }
        return false;
    }


    bool SitePolymorphism::haveSharedAllele(unsigned int pop1, unsigned int pop2) const {
        for (unsigned int allele=0; allele<m_numberOfStates; allele++) {
            if (m_frequencies[pop1][allele]>0 && m_frequencies[pop2][allele]>0
             && m_frequencies[pop1][allele]<m_pop_ns[pop1]
             && m_frequencies[pop2][allele]<m_pop_ns[pop2]) {
                    return true;
            }
        }
        return false;
    }


}

