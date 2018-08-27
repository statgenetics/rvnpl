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

#ifndef EGGLIB_POPULATION_HPP
#define EGGLIB_POPULATION_HPP


#include "Edge.hpp"

namespace egglib {

    class Random;

   /** \brief Handles a single population
    *
    * \ingroup coalesce
    *
    */
    class Population {

        public:
    
           /** \brief Default constructor
            *
            * Generates an empty population.
            * 
            */
            Population();
            
           /** \brief Copy constructor
            * 
            */
            Population(const Population& source);
            
           /** \brief Assignment operator
            * 
            */
            Population& operator=(const Population& source);
            
           /** \brief Destructor
            * 
            * The object only cleans Edge objects currently stored in it.
            *
            */
            ~Population();

           /** \brief Standard constructor
            * 
            * The Edge instances will be handled by address and they
            * MUST be passed using the method set().
            * 
            * \param numberOfSegments number of recombining segments.
            * 
            * \param numberOfLineages the number of lineages contained
            * in this population.
            * 
            * \param firstIndex the absolute index (or ID) of the first
            * lineage (the other will have consecutive incremented
            * ID's).
            *
            */
            Population(unsigned int numberOfSegments,
                    unsigned int numberOfLineages, unsigned firstIndex);

           /** \brief Gets the number of lineages
            * 
            */
            unsigned int numberOfLineages() const;
            
           /** \brief Gets the efficient number of lineages
            * 
            * The number of lineages is multiplied by the number of
            * covered segments of each lineages.
            * 
            */
            unsigned int efficientNumberOfLineages() const;            
            
           /** \brief Sets the Edge of a lineage
            * 
            * \param index the index of the lineage within the
            * population.
            * \param edge the address of the Edge instance representing
            * the lineage.
            *
            */
            void set(unsigned int index, Edge* edge);

           /** \brief Removes and returns a random lineage.
            * 
            * \param random pointer to simulator's random generator
            * instance.
            * 
            */
            Edge* extractRandomly(Random* random);

           /** \brief Removes and returns a given lineage.
            * 
            * \param index the relative index of the lineage.
            * 
            */
            Edge* extractByIndex(unsigned int index);

           /** \brief Appends a lineage to the object
            * 
            */
            void push(Edge* edge);

           /** \brief Gets coverage
            * 
            */
            unsigned int coverage(unsigned int edgeIndex) const;


        private:

            void copy(const Population& source);
            void clear();
            Edge* pick(unsigned int index);
            void init();
            unsigned int _numberOfLineages;
            unsigned int _efficientNumberOfLineages;
            Edge** lineages;
    };

}

#endif
