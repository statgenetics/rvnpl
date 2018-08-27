/*
    Copyright 2008-2009 Stéphane De Mita, Mathieu Siol

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


#ifndef EGGLIB_HAPLOTYPEDIVERSITY_HPP
#define EGGLIB_HAPLOTYPEDIVERSITY_HPP

#include "BaseDiversity.hpp"

namespace egglib {


   /** \brief Computes diversity based on haplotype analysis
    *
    * \ingroup polymorphism
    * 
    * This class relies on detection of polymorphic sites, as does
    * NucleotideDiversity, with the exception that sites with missing
    * data cannot be processed (minimumExploitableData is enforced to
    * 1.).
    * 
    * Like NucleotideDiversity, the same object can be used to analyze
    * different data sets. Only the call to load() is required before
    * accessing the data.
    * 
    * Hst, Gst and Kst are between population differenciation indices.
    * They are respectively defined in equations 2, 5-6 and 9 of Hudson
    * et al. 1992a (Molecular Biology and Evolution 9:138-151). Also,
    * Fst is defined in equation 3 of Hudson et al. 1992b (Genetics
    * 132:583-589). Finally, Snn is from Hudson 2000 Genetics. It is
    * computed as the average of Xi for all sequences. Where Xi is the
    * ratio of nearest neighbours from the same group to the number of
    * nearest neighbours. Nearest neigbours are all the sequences with
    * the lowest number of differences to the focal sequence. NOTE: 
    * Gst/Hst are quite similar, but Fst and Kst are more different. Snn
    * is a different statistic. Gst and Hst are two ways to estimate the
    * between-population fraction of haplotypic diversity.
    * 
    */
    class HaplotypeDiversity : public BaseDiversity {

        public:

           /** \brief Constructor
            * 
            */
            HaplotypeDiversity();
            
           /** \brief Destructor
            * 
            */ 
            virtual ~HaplotypeDiversity();

           /** \brief Identifies polymorphic sites and computes basis
            * statistics
            * 
            * \param data an alignment object (subclass of CharMatrix).
            * The presence of outgroup or of different populations will
            * be detected based on the populationLabel members of the
            * passed object. The populationLabel 999 will be interpreted
            * as outgroups. If several outgroups are passed, sites were
            * the outgroups are not consistent will be treated as "non-
            * orientable".
            * 
            * \param allowMultipleMutations if true, sites with more
            * than two alleles will not be ignored. The sum of the
            * frequencies of all alleles not matching the outgroup will
            * treated as the derived allele frequency (for orientable
            * sites).
            * 
            * \param ignoreFrequency removes sites that are polymorph
            * because of an allele at absolute frequency smaller than or
            * equal to this value. If ignoreFrequency=1, no sites are
            * removed, if ignoreFrequency=1, singleton sites are
            * ignored. Such sites are completely removed from the
            * analysis (not counted in lseff). Note that if more than
            * one mutation is allowed, the site is removed only if all
            * the alleles but one are smaller than or equal to this
            * value. For example, an alignment column AAAAAAGAAT is
            * ignored with an ignoreFrequency of 1, but AAAAAAGGAT is
            * conserved (including the third allele T which is a
            * singleton).
            * 
            * \param characterMapping a string giving the list of
            * characters that should be considered as valid data. If a
            * space is present in the string, the characters left of the
            * space will be treated as valid data and the characters
            * right of the space will be treated as missing data, that
            * is tolerated but ignored. All characters not in the string
            * will cause an EggInvalidCharacterError to be raised.
            * 
            */
            void load(CharMatrix& data,
                bool allowMultipleMutations=false,
                unsigned int ignoreFrequency=0,
                std::string characterMapping=dnaMapping
            );
            
            /// Number of distinct haplotypes
            unsigned int K() const;
            
            /// Haplotype diversity (unbiased)
            double He() const;
            
            /** \brief Returns the allele number of a given sequence
             * 
             * The passed index must be given ignoring any outgroup
             * sequence.
             * 
             */
            unsigned int haplotypeIndex(unsigned int) const;
            
            /// Population differenciation, based on nucleotides (Hudson 1992a)
            double Kst() const;

            /// Population differenciation, based on nucleotides (Hudson 1992b)
            double Fst() const;

            /// Population differenciation, based on haplotypes (Nei version)
            double Gst() const;

            /// Population differenciation, based on haplotypes (Hudson et al. version)
            double Hst() const;
            
            /// Hudson's Snn (nearest neighbor statistics)
            double Snn() const;


        protected:
        
            void init();
            void clear();
                        
            inline unsigned int diff(CharMatrix& data, unsigned int ind1, unsigned int ind2) const;

            bool m_loaded;
            unsigned int m_K;
            double m_He;
            double m_Kst;
            double m_Fst;
            double m_Gst;
            double m_Hst;
            double m_Snn;
            unsigned int *m_haplotypeIndex;


        private:
        
            HaplotypeDiversity(const HaplotypeDiversity& source) {
                
            }
            
            HaplotypeDiversity& operator=(const HaplotypeDiversity& source) {
                return *this;
            }

    };
}

#endif
