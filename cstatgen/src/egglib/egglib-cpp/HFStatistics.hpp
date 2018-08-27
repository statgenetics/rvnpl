/*
    Copyright 2010 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_HFSTATISTICS_HPP
#define EGGLIB_HFSTATISTICS_HPP



namespace egglib {


    /** \brief Computes Fst and Fit from haploid data
    *
    * The class requires loading data. Data are loaded by haploid
    * (one genotype per individual). The analyses are cached: they are
    * performed upon the first call to statistics accessors. The cache
    * is emptied whenever a datum is loaded.
    * 
    * The computations are performed after Weir and Cockerham. The
    * statistic theta is generalized for multiple alleles. To allow
    * computation of multi-locus statistics, variance components are
    * also available. The two components of the variance are T1 and T2
    * and theta is T1/T2 (from Weir 1996 "Genetic Data Analysis II",
    * Sinauer associates, Sunderland MA).
    * 
    * \ingroup polymorphism
    *
    */
    class HFStatistics {
    
        public:
    
           /** \brief Constructor
            * 
            */ 
            HFStatistics();

            
           /** \brief Destructor
            * 
            */ 
            virtual ~HFStatistics();

            
           /** \brief Reserve sufficient memory for a given number of
            * individuals.
            * 
            * This method makes the load function faster by allocating
            * all required memory at once.
            * 
            * \param numberOfIndividuals a strictly positive integer.
            * 
            */
            void reserve(unsigned int numberOfIndividuals);


           /** \brief Loads the data for one individual
            * 
            * \param genotype an integer giving the allele.
            * \param populationLabel an integer indication belonging to
            * a population.
            * 
            * Genotypes and population labels are not required to be
            * consecutive (both are labels, not indices). They are
            * internally mapped to indices (the mapping can be obtained
            * by accessors populationLabel and allele).
            * 
            * All genotypes are considered to be valid (no missing data).
            * If statistics were computed previous to call to this
            * function, all data will be erased.
            * 
            */
            void loadIndividual(unsigned int genotype, unsigned int populationLabel);


           /** \brief Label of a population
            * 
            * The index corresponds to the local mapping of populations
            * regardless of the ranking of population labels. (No out
            * of bound checking.)
            * 
            */
            unsigned int populationLabel(unsigned int populationIndex);


           /** \brief Value of an allele
            * 
            * The index corresponds to the local mapping of alleles
            * regardless of the ranking of allele values. (No out of
            * bound checking.)
            * 
            */
            unsigned int alleleValue(unsigned int alleleIndex);


            /// Allele of a given individual (no checking)
            unsigned int allele(unsigned int individualIndex) const;

            /// Population label of a given individual (no checking)
            unsigned int individualLabel(unsigned int individualIndex) const;


           /** \brief Number of alleles
            * 
            */
            unsigned int numberOfAlleles();


           /** \brief Number of populations
            * 
            */
            unsigned int numberOfPopulations();


           /** \brief Number of loaded genotypes
            * 
            */
            unsigned int numberOfGenotypes() const;


           /** \brief Absolute total allele frequency
            * 
            */
            unsigned int alleleFrequencyTotal(unsigned int alleleIndex);
            

           /** \brief Absolute allele frequency in a population
            * 
            */
            unsigned int alleleFrequencyPerPopulation(unsigned int populationIndex, unsigned int alleleIndex);


           /** \brief Sample size of a population
            * 
            */
            unsigned int populationFrequency(unsigned int populationIndex);


           /** \brief Weir-Cockerham theta-statistic
            * 
            * Note: equivalent to Fst.
            * 
            */
            double theta();


           /** \brief Between-population component of variance
            * 
            */
            double T1();


           /** \brief Total variance
            * 
            */
            double T2();
            

        protected:
    
            bool d_flag;
            void d_init();
            void d_clear();
            unsigned int  d_reserved;
            unsigned int  d_numberOfGenotypes;
            unsigned int *d_genotypes;
            unsigned int *d_populationLabels;

            bool s_flag;
            void s_init();
            void s_clear();
            void s_compute();
            void processPopulations();
            void processAlleles();
            unsigned int getPopulationIndex(unsigned int) const;
            unsigned int getAlleleIndex(unsigned int) const;
            unsigned int    s_numberOfAlleles;
            unsigned int   *s_alleleValueMapping;
            unsigned int    s_numberOfPopulations;
            unsigned int   *s_populationLabelMapping;
            unsigned int   *s_populationFrequencies;
            unsigned int   *s_alleleFrequenciesTotal;
            unsigned int  **s_alleleFrequenciesPerPopulation;

            bool w_flag;
            void w_init();
            void w_clear();
            void w_compute();
            double  w_T;
            double *w_T1;
            double *w_T2;
            double  w_nbar;
            double  w_nc;
            double *w_pbar;
            double *w_ssquare;
            double  w_sum_T1;
            double  w_sum_T2;

    
        private:
            
            HFStatistics(const HFStatistics& source) { }
            
            HFStatistics& operator=(const HFStatistics& source) {
                return *this;
            }

    };
}

#endif
