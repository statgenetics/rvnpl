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

#ifndef EGGLIB_CONTROLLER_HPP
#define EGGLIB_CONTROLLER_HPP

#include "Current.hpp"
#include "Arg.hpp"
#include "ParamSet.hpp"

namespace egglib {
    
    class Random;

   /** \brief Controls a coalescent simulation
    *
    * \ingroup coalesce
    * 
    * This class generates the gene genealogy, based on the parameters
    * stocked in a ParamSet object.
    *
    */
    class Controller {

        public:
    
           /** \brief Default constructor
            *
            * Uses a default ParamSet object that will not allow
            * performing simulations.
            *
            */
            Controller();

           /** \brief Destructor
            * 
            */
            ~Controller();
            
           /** \brief Standard constructor
            * 
            * \param paramSet a ParamSet object containing run
            * parameters (it is taken as a reference and stored as this
            * so it must not be modified during simulations).
            * 
            * \param random the address of the random number generator.
            *
            */
            Controller(const ParamSet* paramSet, Random* random);
            
          /** \brief Reset for a new simulation
           * 
           * Object is reinitiliazed for a new simulation.
           * 
           */
           void reset();
            
           /** \brief Increments the coalescent model
            * 
            * \return The number of lineages.
            *
            */
            unsigned int step();
            
           /** \brief Gets the Ancestral Recombination Graph
            * 
            * \return The address of the ARG contained in the object.
            * 
            */
            Arg* getArg();

           /** \brief Applies a bottleneck to a given population
            * 
            * The bottleneck is applied following Galtier, Depaulis and
            * Barton (Genetics, 2000): the general time counter is
            * stopped, and coalescence events are performed during a
            * time (of normal coalescent process) given by the parameter
            * strength. All coalescent events are instantaneous.
            * 
            * \param populationIndex index of the population concerned
            * by the event.
            * 
            * \param strength strength of the bottleneck given by a
            * number of time units (2N generations times the size of
            * the population).
            * 
            */
            void bottleneck(unsigned int populationIndex, double strength);

           /** \brief Migrate a complete population
            * 
            * Takes all the lineages currently present in the population
            * source to the population dest.
            * 
            */
            void moveAllLineages(unsigned int source, unsigned int dest);
            
           /** \brief Migrate a complete population
            * 
            * Takes all the lineages currently present in the population
            * source to the population dest.
            * 
            * \param source source population.
            * \param dest destination population.
            * \param probability the probability that a lineage of
            * source migrates to dest.
            * 
            */
            void moveSomeLineages(unsigned int source, unsigned int dest, double probability);

            /// Adds an empty population
            void addPopulation();

        private:
        
            /// The copy constructor is disabled
            Controller(const Controller& source) {}
            
            /// The assignment operator is disabled
            Controller& operator=(const Controller& source) {return *this;}

            void diploids();
            double getMigrationTime(double& migrationParameterDestination);
            void getCoalescenceTime(double& destTime, unsigned int& destPopIndex);
            double getCoalescenceTimeForPopulation(unsigned int populationIndex);
            double getRecombinationTime() const;
            void migrate(double migrationParameter);

            const ParamSet* initialParamSet;
            ParamSet paramSet;
            Current current;
            Arg arg;
            
            Random* random;

    };

}

#endif
