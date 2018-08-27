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


#ifndef EGGLIB_SITEPOLYMORPHISM_HPP
#define EGGLIB_SITEPOLYMORPHISM_HPP



namespace egglib {


   /** \brief Implements diversity analysis at the site level
    *
    * \ingroup polymorphism
    * 
    * Data are loaded along with a population index. It is necessary to
    * set the number of populations prior to use.
    * 
    * Outgroup sequence must be loaded separetedly. There can be any
    * number of outgroups, but they must be all consistent otherwise the
    * site will be considered as not orientable.
    * 
    */
    class SitePolymorphism {

        public:

           /** \brief Builds an object
            * 
            */
            SitePolymorphism();


           /** \brief Builds an object
            * 
            * \param npop number of populations
            * 
            */
            SitePolymorphism(unsigned int npop);


           /** \brief Destroys an object
            * 
            */
            virtual ~SitePolymorphism();


           /** \brief Copy constructor
            * 
            */
            SitePolymorphism(const SitePolymorphism& source);


           /** \brief Assignment operator
            * 
            */
            SitePolymorphism& operator=(const SitePolymorphism& source);


           /** \brief Sets the number of populations
            * 
            * NOTE THAT all previous data is lost.
            * 
            */
            void numberOfPopulations(unsigned int npop);


           /** \brief Adds a character
            * 
            * \param populationIndex the index of the population from
            * which is sampled this character (do not use "population
            * label").
            * 
            * \param character the character value (it is assumed it
            * represents a valid character.
            * 
            */
            void load(unsigned int populationIndex, char character);


           /** \brief Loads outgroup state
            * 
            * There can be any number of outgroup states. Only
            * characters that are considered as valid (whatever the list
            * is) should be loaded.
            * 
            */
            void outgroup(char state);


           /** \brief Number of different alleles
            * 
            */
            unsigned int numberOfAlleles() const;
            
            
           /** \brief Gets an allele (unsecure)
            * 
            * Assumes that the index provided lies in the valid range
            * 
            */
            char allele(unsigned int index) const;


           /** \brief Gets a frequency (unsecure)
            * 
            * The sum of of frequencies of the allele over populations
            * is computed. Not out-of-bounds check is performed.
            * 
            */
            unsigned int alleleFrequency(unsigned int alleleIndex) const;


           /** \brief Gets the frequency of an allele in one pop (unsecure)
            * 
            * The frequency of the allele in the given population is
            * returned. Not out-of-bounds check is performed.
            * 
            */
            unsigned int alleleFrequency(unsigned int popIndex, unsigned int alleleIndex) const;


           /** \brief Sums the frequency of derived allele(s)
            * 
            * This method assumes that the site is orientable. It will
            * use as outgroup the first outgroup character entered,
            * assuming at least one was entered and that all (if more
            * than one) were identical.
            * 
            */
            unsigned int derivedAlleleFrequency() const;


           /** \brief Number of sequences that were analyzed
            * 
            */
            unsigned int ns() const;


           /** \brief Gets the number of analyzed sequences for a population
            * 
            * No out-of-bound check is performed
            * 
            */
            unsigned int ns(unsigned int popIndex) const;


           /** \brief Checks if the site can be oriented
            * 
            * Returns true if at least one outgroup datum has been
            * loaded, if all outgroup data are identical (regardless of
            * their value) and if the outgroup allele is one of the
            * allele in the sample.
            * 
            */
            bool isOrientable() const;

            bool isPolymorphic(unsigned int popIndex) const;
            bool hasSpecificAllele(unsigned int popIndex, bool restrictToDerived) const;
            bool haveFixedDifference(unsigned int pop1, unsigned int pop2) const;
            bool haveCommonAllele(unsigned int pop1, unsigned int pop2) const;
            bool haveSharedAllele(unsigned int pop1, unsigned int pop2) const;




        protected:

            // helpers
            void init();
            void clear();
            void copy(const SitePolymorphism& site);


            // data
            unsigned int m_numberOfPopulations;
            unsigned int m_numberOfStates;
            char * m_states;
            unsigned int ** m_frequencies;
            unsigned int m_numberOfOutgroups;
            char * m_outgroups;
            unsigned int m_ns;
            unsigned int * m_pop_ns;
            
            bool m_cache_orientable;

    };
}

#endif
