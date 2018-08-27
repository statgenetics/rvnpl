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


#ifndef EGGLIB_NUCLEOTIDEDIVERSITY_HPP
#define EGGLIB_NUCLEOTIDEDIVERSITY_HPP


#include "BaseDiversity.hpp"
#include <string>
#include <vector>



namespace egglib {


   /** \brief Performs analyzes of population genetics
    *
    * \ingroup polymorphism
    * 
    * This class computes several summary statistics based on
    * nucleotide analysis. Note that it is possible to use the same
    * object to analyze different data set. Calling the load() method
    * erases all data preivously computed (if any). Calling the load()
    * method is absolutely required to compute any statistics. Some
    * statistics are not computed by default, but are if the
    * corresponding accessor is used (only load() is required).
    * 
    * Note that "unsecure" accessors don't perform out-of-bound checks.
    * 
    * S is the number of varying sites (only in sites that were not
    * rejected).
    * 
    * eta is the minimum number of mutations, that is the sum of the
    * number of alleles minus 1 for each varying site. eta = S if all
    * sites have no variant or 2 alleles. eta is computed independently
    * of the option multiple and IS NOT computed over lseff sites.
    *
    * Pi is the average number of pairwise differences between sequences
    * (expressed here per site) or (as computed here) the mean per site
    * (unbiased) heterozygosity. Pi is zero if no polymorphic sites.
    *
    * D is the Tajima's test of neutrality
    * Ref. Tajima F.: Statistical method for testing the neutral
    * mutation hypothesis by DNA polymorphism. Genetics 1989, 123:585-595.
    * It is arbitrary set to 0 if no polymorphic sites.
    *
    * tW: thetaW: estimator of theta based on polymorphic sites (ref.
    * e.g. Watterson 1975 Theor. Pop. Biol.).
    * Both D and thetaW are computed assuming that rounded nseff samples
    * have been sampled.
    * The variance of D is computed using rounded nseff instead of ns.
    *
    * H is the Fay and Wu's test of neutrality.
    * Z is the standardized version and E a similar test.
    * Ref. Fay J. C., Wu C.-I.: Hitchhiking under positive Darwinian
    * selection. Genetics 2000, 155:1405-1413. and Zeng K., Fu Y. X.,
    * Shi S., Wu C.-I.: Statistical tests for detecting positive
    * selection by utilizing high-frequency variants. Genetics 2006,
    * 174:1431-9. Both are arbitrary set to 0 if no polymorphic or
    * orientable sites.
    *
    * tH and tL: theta H: estimators of theta based on derived
    * polymorphic sites (ref in Fay and Wu and Zeng al.). The variance
    * of H/Z are computed assuming that rounded nseff samples have
    * been sampled.
    * 
    */
    class NucleotideDiversity : public BaseDiversity {

        public:

           /** \brief Builds an object
            * 
            */
            NucleotideDiversity();


           /** \brief Destroys an object
            * 
            */
            virtual ~NucleotideDiversity();


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
            * \param minimumExploitableData sites where the non-missing
            * data (as defined by characterMapping) are at a frequency
            * larger than this value will be removed from the analysis.
            * Use 1. to take only 'complete' sites into account and 0.
            * to use all sites. (The outgroup is not considered in this
            * computation.)
            * 
            * \param ignoreFrequency removes sites that are polymorph
            * because of an allele at absolute frequency smaller than or
            * equal to this value. If ignoreFrequency=1, no sites are
            * removed, if ignoreFrequency=0, singleton sites are
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
            * \param useZeroAsAncestral if true, all outgroups (if
            * present) will be ignored and the character "0" will be
            * considered as ancestral for all sites, whatever the
            * character mapping.
            * 
            */
            virtual void load(
                CharMatrix& data,
                bool allowMultipleMutations=false,
                double minimumExploitableData=1.,
                unsigned int ignoreFrequency=0,
                std::string characterMapping=dnaMapping,
                bool useZeroAsAncestral=false
            );


        // accessors for the "site analysis" section

            /// Number of polymorphic sites
            unsigned int S() const;
            
            /// Number of polymorphic orientable sites
            unsigned int So() const;
            
            /// Minimum number of mutations
            unsigned int eta() const;
            
            /// Average of per-site number of sequences effectively used
            double nseff() const;
            
            /// Number of sites effectively used
            unsigned int lseff() const;
            
            /// Average of number of sequences effectively used at orientable sites
            double nseffo() const;
            
            /// Number of orientable sites
            unsigned int lseffo() const;
            
            /// Number of detected populations
            unsigned int  npop() const;
            
            /// Label of the population with given index (unsecure)
            unsigned int popLabel(unsigned int popIndex) const; // no check!


        // accessors for the "diversity" section

            /// Nucleotide diversity
            double Pi();
            
            /// Watterson estimator of theta
            double thetaW();
            
            /// Average of Pi over populations
            double average_Pi();
            
            /// Pi of a given population (unsecure)
            double pop_Pi(unsigned int popIndex);  // no check!
            
            /// Tajima's D
            double D();

        // accessors for the "outgroup diversity" section
        
            /// Fay and Wu estimator of theta
            double thetaH();
            
            /// Zeng et al. estimator of theta
            double thetaL();
            
            /// Fay and Wu's H
            double H();
            
            /// Standardized H
            double Z();
            
            /// Zeng et al.'s E
            double E();

         // accessors for the "differentiation" section
         
            /// Number of sites with at least one fixed difference
            unsigned int FixedDifferences();
            
            /// Number of sites with at least one allele shared among at least two populations
            unsigned int CommonAlleles();
            
            /// Number of sites with at least one non-fixed allele shared among at least two populations
            unsigned int SharedAlleles();
            
            /// Number of sites with at least one allele specific to one population
            unsigned int SpecificAlleles();
            
            /// Number of sites with at least one derived allele specific to one population
            unsigned int SpecificDerivedAlleles();
            
            /// Number of polymorphisms in a given population (unsecure)
            unsigned int Polymorphisms(unsigned int pop);
            
            /// Number of specific alleles for a given population (unsecure)
            unsigned int SpecificAlleles(unsigned int pop);
            
            /// Number of specific derived allele for a given population (unsecure)
            unsigned int SpecificDerivedAlleles(unsigned int pop);
            
            /// Number of fixed differences between a given pair of populations (unsecure; pop2 must be larger than pop1)
            unsigned int FixedDifferences(unsigned int pop1, unsigned int pop2);

            /// Number of common alleles between a given pair of populations (unsecure; pop2 must be larger than pop1)
            unsigned int CommonAlleles(unsigned int pop1, unsigned int pop2);

            /// Number of shared non-fixed alleles between a given pair of populations (unsecure; pop2 must be larger than pop1)
            unsigned int SharedAlleles(unsigned int pop1, unsigned int pop2);


        // accessor for the "triConfigurations" section

           /** \brief Number falling into one of the possible site configurations
            *
            * The statistics are limited to three populations.
            * Assuming an unrooted A/G polymorphism (A and G can be
            * substitued), the site configurations are:
            *     -  0: A&G  A   A  specific 1
            *     -  1: A&G  A   G  specific 1 + fixed 2-3
            *     -  2:  A  A&G  A  specific 2
            *     -  3:  A  A&G  G  specific 2 + fixed 1-3
            *     -  4:  A   A  A&G specific 3
            *     -  5:  A   G  A&G specific 3 + fixed 1-2
            *     -  6: A&G A&G  A  shared 1-2
            *     -  7: A&G  A  A&G shared 1-3
            *     -  8:  A  A&G A&G shared 2-3
            *     -  9: A&G A&G A&G shared 1-2-3
            *     - 10:  A   G   G  fixed 1
            *     - 11:  A   G   A  fixed 2
            *     - 12:  A   A   G  fixed 3
            *
            * \param index must be an index from 0 to 12.
            * 
            */
            unsigned int triConfiguration(unsigned int index);


        /// Builds and returns the vector of positions of all polymorphic sites
        std::vector<unsigned int> polymorphic_positions() const;


        /** \brief Builds and returns the vector of positions of all singleton sites
         * 
         * A site singleton when it is polymorphic according to
         * parameter of the diversity analysis, when it has exactly two
         * alleles and one of them is at absolute frequency 1 (one
         * copy) disregarding the outgroup.
         * 
         */
        std::vector<unsigned int> singleton_positions() const;


        protected:

           /** \brief This class cannot be copied
            * 
            */
            NucleotideDiversity(const NucleotideDiversity& source) { }


           /** \brief This class cannot be copied
            * 
            */
            NucleotideDiversity& operator=(const NucleotideDiversity& source) { return *this; }


            void init();  // initializes values
            void clear();  // free memory but doesn't initializes
            
            // diversity (without outgroup)
            void diversity();
            
            // diversity with outgroup
            void outgroupDiversity();
            
            // site patterns
            void differentiation();
            
            // triconfigurations
            void triConfigurations();
            

            // holders for statistics, with booleans flagging groups of stats
            
            bool b_analysisSites;
            
            bool b_diversity;
            
            double  v_Pi;             // nucleotide diversity
            double  v_thetaW;         // theta (Watterson estimator)
            double  v_average_Pi;     // average diversity across populations
            double *v_pop_Pi;         // diversity per population
            double  v_D;              // Tajima's D
            
            bool b_outgroupDiversity;
            
            double v_thetaH;        // theta (Fay and Wu estimator)
            double v_thetaL;        // theta (Zeng estimator)
            double v_H;             // Fay and Wu's H
            double v_Z;             // normalized Fay and Wu's H
            double v_E;             // Zeng et al.'s E
            
            bool b_differentiation;
            
            unsigned int  *v_pairwiseFixedDifferences;
            unsigned int  *v_pairwiseCommonAlleles;
            unsigned int  *v_pairwiseSharedAlleles;
            unsigned int  *v_popPolymorphic;
            unsigned int  *v_popSpecific;
            unsigned int  *v_popSpecificDerived;
            unsigned int   v_countFixedDifferences;
            unsigned int   v_countCommonAlleles;
            unsigned int   v_countSharedAlleles;
            unsigned int   v_countSpecificAlleles;
            unsigned int   v_countSpecificDerivedAlleles;
        
            
            bool b_triConfigurations;
            
            unsigned int  *v_triConfigurations;

    };
}

#endif
