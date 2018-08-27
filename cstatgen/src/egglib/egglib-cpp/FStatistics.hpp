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

#ifndef EGGLIB_FSTATISTICS_HPP
#define EGGLIB_FSTATISTICS_HPP



namespace egglib {


    /** \brief Computes Fis, Fst and Fit from diploid data
    *
    * The class requires loading data. Data are loaded by individual
    * (two genotypes per individual). The analyses are cached: they are
    * performed upon the first call to statistics accessors. The cache
    * is emptied whenever a datum is loaded.
    * 
    * The computations are performed after Weir and Cockerham. The
    * statistics F, theta and f are generalized for multiple alleles.
    * To allow computation of multi-locus statistics, variance
    * components are also available. The three components of the
    * variance are Vpopulation (between-population), Vindividual
    * (within-population, between-individual) and Vallele (within-
    * individual). The formulas to compute the F-statistics are as
    * follows:
    *       - 1-F = Vallele/(Vpopulation+Vindividual+Vallele)
    *       - theta = Vpopulation/(Vpopulation+Vindividual+Vallele)
    *       - 1-f = Vallele/(Vindividual+Vallele).
    * 
    * \ingroup polymorphism
    *
    */
    class FStatistics {
    
        public:
    
           /** \brief Constructor
            * 
            */ 
            FStatistics();

            
           /** \brief Destructor
            * 
            */ 
            virtual ~FStatistics();

            
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
            * \param genotype1 an integer giving the first allele.
            * \param genotype2 an integer giving the second allele.
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
            * function, all data will be erase.
            * 
            */
            void loadIndividual(unsigned int genotype1,
                    unsigned int genotype2, unsigned int populationLabel);


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


            /// First allele of a given individual (no checking)
            unsigned int firstAllele(unsigned int individualIndex) const;

            /// Second allele of a given individual (no checking)
            unsigned int secondAllele(unsigned int individualIndex) const;

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


           /** \brief Absolute genotype frequency
            * 
            * Note that allele AB is considered different to BA (this
            * means that values can be accessed both sides of the
            * diagonal.
            * 
            */
            unsigned int genotypeFrequencyTotal(unsigned int alleleIndex1, unsigned int alleleIndex2);


           /** \brief Absolute genotype frequency in a population
            * 
            * Note that allele AB is considered different to BA (this
            * means that values can be accessed both sides of the
            * diagonal.
            * 
            */
            unsigned int genotypeFrequencyPerPopulation(unsigned int populationIndex, unsigned int alleleIndex1, unsigned int alleleIndex2);

            
           /** \brief Sample size of a population
            * 
            */
            unsigned int populationFrequency(unsigned int populationIndex);


           /** \brief Weir-Cockerham F-statistic
            * 
            * Note: equivalent to Fit.
            * 
            */
            double F();


           /** \brief Weir-Cockerham theta-statistic
            * 
            * Note: equivalent to Fst.
            * 
            */
            double theta();


           /** \brief Weir-Cockerham f-statistic
            * 
            * Note: equivalent to Fis.
            * 
            */
            double f();
            

           /** \brief Between-population component of variance
            * 
            */
            double Vpopulation();


           /** \brief Within-population, between-individual component of variance
            * 
            */
            double Vindividual();
            
            
           /** \brief Within-individual component of variance
            * 
            */
            double Vallele();


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
            unsigned int  **s_genotypeFrequenciesTotal;
            unsigned int ***s_genotypeFrequenciesPerPopulation;

            bool w_flag;
            void w_init();
            void w_clear();
            void w_compute();
            double  w_F;
            double  w_T;
            double  w_f;
            double *w_a;
            double *w_b;
            double *w_c;
            double  w_nbar;
            double  w_nc;
            double *w_pbar;
            double *w_ssquare;
            double *w_hbar;
            double  w_sum_a;
            double  w_sum_b;
            double  w_sum_c;
            double  w_sum_abc;
            double  w_sum_bc;

    
        private:
            
            FStatistics(const FStatistics& source) { }
            
            FStatistics& operator=(const FStatistics& source) {
                return *this;
            }

    };
}

#endif
