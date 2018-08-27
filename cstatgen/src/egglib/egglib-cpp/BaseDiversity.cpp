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


#include "BaseDiversity.hpp"
#include "EggException.hpp"
#include <cstdlib>


namespace egglib {

    const std::string BaseDiversity::dnaMapping = "ACGT MRWSYKBDHVN?-";
    const std::string BaseDiversity::rnaMapping = "ACGU MRWSYKBDHVN?-";
    const std::string BaseDiversity::aaMapping = "GALMFWKQESPVICYHRNDT X-";

    void BaseDiversity::init() {
        v_reserved = 0;
        v_S = 0;
        v_sites = NULL;
        v_orientables = NULL;
        v_sitePositions = NULL;
        v_ns = 0;
        v_So = 0;
        v_eta = 0;
        v_nseff = 0.;
        v_nseffo = 0.;
        v_lseff = 0;
        v_lseffo = 0;
        v_npop = 0;
        v_popLabel = NULL;

        p_allowMultipleMutations = false;
        p_minimumExploitableData = 0.;
        p_characterMapping.clear();
        p_pos_sep_mapping = 0;
        p_useZeroAsAncestral = false;
        p_ignoreFrequency = 0;
    }
    
    
    void BaseDiversity::clear() {
        if (v_sites) {
            for (unsigned int i=0; i<v_S; i++) {
                delete v_sites[i];
            }
            free(v_sites);
            v_sites = NULL;
        }
        if (v_orientables) {
            free(v_orientables);
            v_orientables = NULL;
        }
        if (v_sitePositions) {
            free(v_sitePositions);
            v_sitePositions = NULL;
        }
        if (v_popLabel) {
            free(v_popLabel);
            v_popLabel = NULL;
        }

    }


    void BaseDiversity::reset() {
        clear();
        init();
    }


    BaseDiversity::BaseDiversity() {
        init();
    }


    BaseDiversity::~BaseDiversity() {
        clear();
    }

    const SitePolymorphism* BaseDiversity::get_site(unsigned int index) const {
        if (index>=v_S) {
            throw EggArgumentValueError("invalid site index");
        }
        return v_sites[index];
    }

    unsigned int BaseDiversity::get_position(unsigned int index) const {
        if (index>=v_S) {
            throw EggArgumentValueError("invalid site index");
        }
        return v_sitePositions[index];
    }
    
    void BaseDiversity::reserve(unsigned int numberOfSites) {
        if (numberOfSites==0) {
            return;
        }

        if (numberOfSites<v_S) {
            throw EggArgumentValueError("cannot reserve less memory than what is currently used for diversity analysis");
        }

        v_sites = (SitePolymorphism**) realloc(v_sites, numberOfSites*sizeof(SitePolymorphism*));
        if (!v_sites) throw EggMemoryError();
            
        v_orientables = (bool*) realloc(v_orientables, numberOfSites*sizeof(bool));
        if (!v_orientables) throw EggMemoryError();
            
        v_sitePositions = (unsigned int*) realloc(v_sitePositions, numberOfSites*sizeof(unsigned int));
        if (!v_sitePositions) throw EggMemoryError();
        
        v_reserved = numberOfSites;
    }
    
    
    void BaseDiversity::importSites(CharMatrix& data,
        bool allowMultipleMutations,  double minimumExploitableData,
        unsigned int ignoreFrequency,  std::string characterMapping,
        bool useZeroAsAncestral, bool ignoreOutgroups) {

            // copies parameters
            p_allowMultipleMutations = allowMultipleMutations;
            if (minimumExploitableData<0. || minimumExploitableData>1.) {
                throw EggArgumentValueError("invalid value for argument minimumExploitableData of Diversity::load: only values from the range [0.,1.] are accepted");
            }
            p_minimumExploitableData = minimumExploitableData;
            p_characterMapping = characterMapping;
            p_pos_sep_mapping = characterMapping.find(' ');
            p_useZeroAsAncestral = useZeroAsAncestral;
            p_ignoreFrequency = ignoreFrequency;
            
            // presets threshold of missing data
            v_ns=0;
            for (unsigned int i=0; i<data.numberOfSequences(); i++) {
                if (data.populationLabel(i)!=999) {
                    v_ns++;
                }
            }
            double maxMissingData = (1-p_minimumExploitableData) * v_ns;
    
            // number of populations
            for (unsigned int i=0; i<data.numberOfSequences(); i++) {
                unsigned int label = data.populationLabel(i);
                if (label==999 && !useZeroAsAncestral) continue;
                unsigned int index = getPopIndex(label);
                if (index==v_npop) {
                    v_npop++;
                    v_popLabel = (unsigned int*) realloc(v_popLabel, v_npop*sizeof(unsigned int));
                    if (!v_popLabel) throw EggMemoryError();
                    v_popLabel[index] = label;
                }
            }
                    
            // iterates over sites
            for (unsigned int i=0; i<data.numberOfSites(); i++) {
                analyzeSite(data, i, maxMissingData, ignoreOutgroups);
            }
                
            // computes averages
            if (v_lseff) {
                v_nseff/=v_lseff;
                v_nseffo/=v_lseffo;
            }

    }

    unsigned int BaseDiversity::getPopIndex(unsigned int label) const {
        for (unsigned int i=0; i<v_npop; i++) {
            if (v_popLabel[i]==label) return i;
        }
        return v_npop;
    }

    void BaseDiversity::analyzeSite(CharMatrix& data, unsigned int index, double maxMissingData, bool ignoreOutgroup) {

        // create a site
        SitePolymorphism* site = new(std::nothrow) SitePolymorphism(v_npop);
        if (!site) throw EggMemoryError();
        if (p_useZeroAsAncestral) site->outgroup('0');
        
        // loads the sequences (only valid characters, and separating outgroup(s))
        // in case missingData is reach in the loop, exits immediately
        unsigned int missingData = 0;
        for (unsigned int i=0; i<data.numberOfSequences(); i++) {
            char c = data.character(i, index);
            std::size_t pos = p_characterMapping.find(c);
            if (pos<p_pos_sep_mapping) {
                if (p_useZeroAsAncestral || data.populationLabel(i)!=999) {
                    site->load(getPopIndex( data.populationLabel(i) ), c);
                }
                else {
                    if (!ignoreOutgroup) {
                        site->outgroup(c);
                    }
                }
            }
            else {
                if (pos>p_pos_sep_mapping && pos<std::string::npos) {
                    // the if statement ensure that doesnt eliminate the site if only the outgroup has missing data
                    if (p_useZeroAsAncestral || data.populationLabel(i)!=999) {
                        missingData+=1;
                        if (missingData>maxMissingData) {
                            delete site;
                            return;
                        }
                    }
                }
                else {
                    delete site;
                    throw EggInvalidCharacterError(c, i+1, index+1);
                }
            }
        }
        
        // ignore site without data
        if (site->ns()==0) {
            delete site;
            return;
        }
        
        // computes number of mutations anyway
        v_eta+= site->numberOfAlleles()-1;
        bool orientable = site->isOrientable();
        
        // safety check
        if (site->numberOfAlleles()==0) {
            delete site;
            throw EggRuntimeError("site has data but no alleles");
        }
        
        // if too many alleles, don't count the site
        if (site->numberOfAlleles()>2 && (!p_allowMultipleMutations)) {
            delete site;
            return;
        }

        // from here the site is "analyzed"
        v_nseff+=site->ns();
        v_lseff++;
        if (orientable) {
            v_nseffo+=site->ns();
            v_lseffo++;
        }

        // if all alleles (but one) is a singleton (according to parameter)
        // don't count it
        
        bool frequentEnough = true;
        
        if (p_ignoreFrequency>0) {
            if (p_ignoreFrequency>0) {
                unsigned int numberOfCopyFrequentEnough=0;
                for (unsigned int i=0; i<site->numberOfAlleles(); i++) {
                    if (site->alleleFrequency(i)>p_ignoreFrequency) {
                        numberOfCopyFrequentEnough++;
                        if (numberOfCopyFrequentEnough>=2) break;
                    }
                }
                if (numberOfCopyFrequentEnough<2) {
                    frequentEnough = false;
                }
            }
        }
        
        // if the site is polymorph, keeps it
        if (frequentEnough && site->numberOfAlleles()>1) {

            // counts the mutation
            v_S++;
            if (orientable) v_So++;
            
            // adds to the site vector
            if (v_S>v_reserved) {
                v_sites = (SitePolymorphism**) realloc(v_sites, v_S*sizeof(SitePolymorphism*));
                if (!v_sites) throw EggMemoryError();
            }
            v_sites[v_S-1] = site;
            
            // adds the orientable flag
            if (v_S>v_reserved) {
                v_orientables = (bool*) realloc(v_orientables, v_S*sizeof(bool));
                if (!v_orientables) throw EggMemoryError();
            }
            v_orientables[v_S-1] = orientable;
            
            // adds the position
            if (v_S>v_reserved) {
                v_sitePositions = (unsigned int*) realloc(v_sitePositions, v_S*sizeof(unsigned int));
                if (!v_sitePositions) throw EggMemoryError();
            }
            v_sitePositions[v_S-1] = index;
        }
        
        // otherwise, trashs it
        else {
            delete site;
        }
        
    }

}

