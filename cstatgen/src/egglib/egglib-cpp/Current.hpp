/*
    Copyright 2009-2010 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_CURRENT_HPP
#define EGGLIB_CURRENT_HPP


namespace egglib {

    class Population;
    class ParamSet;

   /** \brief Represents the current set of populations
    *
    * \ingroup coalesce
    * 
    */
    class Current {

        public:

           /** \brief Default constructor
            * 
            */
            Current();
            
           /** \brief Standard constructor
            *
            * \param paramSet allows to initiate the correct structure
            * of populations.
            * 
            */
            Current(ParamSet* paramSet);

           /** \brief Rebuilds the object
            *
            * \param paramSet allows to initiate the correct structure
            * of populations.
            * 
            */
            void reset(ParamSet* paramSet);
            
           /** \brief Destructor
            * 
            */
            virtual ~Current();

           /** \brief Copy constructor
            * 
            */
            Current(const Current&);

           /** \brief Assignment operator
            * 
            */
            Current& operator=(const Current&);

           /** \brief Gets the current number of populations
            * 
            */
            unsigned int numberOfPopulations() const;

            
           /** \brief Adds an empty population to the system
            *
            */
            void addPopulation();

            
           /** \brief Gets the number of lineages contained by a given
            * population
            * 
            */
            unsigned int populationNumberOfLineages(unsigned int populationIndex) const;

            
           /** \brief Provides access to a given population
            * 
            * The returned pointer can be used to modify the object.
            * 
            */
            Population* population(unsigned int populationIndex);
            
            
           /** \brief Total number of lineages
            * 
            */
            unsigned int totalNumberOfLineages() const;


           /** \brief Efficient number of lineages
            * 
            * This sums the number of covered segments of each lineage.
            * 
            */
            unsigned int efficientNumberOfLineages() const;


        private:

            void setPopulationArray();
            void copy(const Current&);
            void clear();

            unsigned int _numberOfPopulations;
            unsigned int _numberOfSegments;
            Population** populations;
    };

}

#endif
