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

#ifndef EGGLIB_PARAMSET_HPP
#define EGGLIB_PARAMSET_HPP


#include "DataMatrix.hpp"


namespace egglib {

    class Change;
    class Controller;


   /** \brief Set of parameters
    *
    * \ingroup coalesce
    *
    */
    class ParamSet {

        public:
    
           /** \brief Default constructor
            *
            * Initializes all parameters to reasonnable values (except
            * that the sample size is null: 1 population, 0 samples,
            * selfing rate of 0, recombination rate of 0, growth rate of
            * 0, population size of 1 and no changes.
            *
            */
            ParamSet();

           /** \brief Destructor
            * 
            */
            ~ParamSet();
            
           /** \brief Copy constructor
            * 
            */
            ParamSet(const ParamSet&);
            
           /** \brief Assignment operator
            * 
            */
            ParamSet& operator=(const ParamSet&);

           /** \brief Restores default value of all parameters
            * 
            */
            void reset();

           /** \brief Gets the number of populations
            * 
            */
            unsigned int numberOfPopulations() const;
            
           /** \brief Gets a pairwise migration rate
            * 
            * It is allowed to access a diagonal value. Diagonal
            * values contain the sum of values of the corresponding
            * line (diagonal cell excepted, of course).
            * 
            */
            double pairwiseMigrationRate(unsigned int source, unsigned int dest) const;
            
           /** \brief Sets a pairwise migration rate
            * 
            * It is not allowed to set a value on the diagonal (this
            * would raise an exception). The method takes care of
            * modifying the diagonal accordingly (still this is not
            * relevant for the client);
            * 
            */
            void pairwiseMigrationRate(unsigned int source, unsigned int dest, double value);
            
           /** \brief Sets the migration rate for all matrix
            * 
            */
            void migrationRate(double value);
            
           /** \brief Gets a population size
            * 
            */
            double populationSize(unsigned int populationIndex) const;
            
           /** \brief Sets a population size
            * 
            * The size must be strictly positive.
            * 
            */
            void populationSize(unsigned int populationIndex, double value);
            
           /** \brief Gets a growth rate
            * 
            */
            double growthRate(unsigned int populationIndex) const;
            
           /** \brief Sets a growth rate
            * 
            */
            void growthRate(unsigned int populationIndex, double value);
            
           /** \brief Gets the recombination rate
            * 
            */
            double recombinationRate() const;
            
           /** \brief Sets the recombination rate
            * 
            */
            void recombinationRate(double value);

           /** \brief Gets the number of recombining segments
            * 
            */
            unsigned int numberOfSegments() const;
            
           /** \brief Sets the number of recombining segments
            * 
            */
            void numberOfSegments(unsigned int value);

           /** \brief Gets the selfing rate
            * 
            */
            double selfingRate() const;
            
           /** \brief Sets the selfing rate
            * 
            */
            void selfingRate(double value);
            
           /** \brief Adds a population to the model
            * 
            * \param migr pairwise migration rate between other
            * population and the new one.
            *
            * The parameters for the new population are equivalent to
            * the default parameters.
            * 
            */
            void addPopulation(double migr);
            
           /** \brief Adds a change
            * 
            * The change can be of any type derived from Change. A
            * const Change* must be passed, so ParamSet will neither
            * modify or delete the object itself. All passed Change
            * object must be kept available for use by ParamSet.
            *
            */
            void addChange(const Change* change);

           /** \brief Gets the date of the next change
            * 
            * Returns -1 if no change is planned.
            * 
            */
            double nextChangeDate() const;
            
           /** \brief Applies the next change event
            * 
            * \param controller the Change event might need to have
            * access to simulation controller (to trigger coalescent
            * events, for example).
            * 
            */
            void nextChangeDo(Controller* controller);
            
           /** \brief Gets the number of single sample from a population
            * 
            */
            unsigned int singles(unsigned int populationIndex) const;

           /** \brief Sets the number of single sample from a population
            * 
            */
            void singles(unsigned int populationIndex, unsigned int value);

           /** \brief Gets the number of double sample from a population
            * 
            */
            unsigned int doubles(unsigned int populationIndex) const;

           /** \brief Sets the number of double sample from a population
            * 
            */
            void doubles(unsigned int populationIndex, unsigned int value);
            
           /** \brief Computes the total number of samples
            * 
            */
            unsigned int numberOfSamples() const;
            
           /** \brief Gives the date of the last size change
            * 
            * \param populationIndex the index of the population.
            * \return The date where the last change occurred, or 0. if
            * no change occurred during the simulation.
            *
            */
            double dateOfLastChange(unsigned int populationIndex) const;


           /** \brief Sets the date of the last size change
            * 
            * \param populationIndex the index of the population.
            * \param date the date where the last change occurred, or 0.
            * if no change occurred during the simulation.
            *
            */
            void dateOfLastChange(unsigned int populationIndex, double date) const;

            
           /** \brief Set groups labels
            * 
            * Sets the group labels of the DataMatrix, according to the
            * current state of population structure, and assuming that
            * the DataMatrix was generated by the class Arg.
            * 
            * \param dataMatrix the DataMatrix object to modify. The
            * number of sequences must match the total number of samples
            * defined by the ParamSet object this method is called on.
            * 
            * \param labelIndividuals by default, labels the different
            * samples depending on the population they come from (0
            * being the label of the first population). If this flag is
            * set to true, then the samples are labelled depending on
            * the individual they come from, regardless of populations.
            * In that case there can be only one or two genes for a
            * given group label.
            * 
            */
            void setGroups(DataMatrix& dataMatrix, bool labelIndividuals=false);

        private:

            void clear();
            void init();
            void copy(const ParamSet&);
        
            double _selfingRate;
            double _recombinationRate;
            unsigned int _numberOfSegments;
            unsigned int _numberOfPopulations;
            unsigned int* _singles;
            unsigned int* _doubles;
            double* _growthRates;
            double* _populationSize;
            double* _dateOfLastChange;
            double** migrationMatrix;
            unsigned int _numberOfChanges;
            unsigned int nextChangeIndex;
            Change const** changes;
    };

}

#endif
