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

#ifndef EGGLIB_CHANGETYPES_HPP
#define EGGLIB_CHANGETYPES_HPP

#include "ParamSet.hpp"
#include "Controller.hpp"

namespace egglib {

/**********************************************************************/

   /** \brief Pure virtual base class for parameter changes
    *
    * \ingroup coalesce
    *
    */
    class Change {
        public:
    
           /** \brief Default constructor
            *
            * The default date is 0.
            *
            */
            Change();
            
           /** \brief Standard constructor
            *
            * \param date the event date.
            *
            */
            Change(double date);

            /// Gets the event date value
            double date() const;
            
            /// Sets the event date value
            void date(double value);
            
           /** \brief Applies the event
            *
            * \param paramSet the ParamSet instance to which the Change
            * instance is attached.
            * \param controller the Controller instance of the
            * simulation.
            *
            */
            virtual void apply(ParamSet* paramSet, Controller* controller) const = 0;
        
        protected:
            double _date;
    
    };

/**********************************************************************/

   /** \brief Pure virtual base class for single parameter changes
    *
    * \ingroup coalesce
    *
    */
    class SingleParamChange : public Change {
        public:
    
           /** \brief Default constructor
            *
            * The default date is 0., the default parameter value is 0.
            *
            */
            SingleParamChange();
            
           /** \brief Standard constructor
            *
            * \param date the event date.
            * \param value the parameter value.
            *
            */
            SingleParamChange(double date, double value);

            /// Gets the parameter value
            double value() const;
            
            /// Sets the parameter value
            void value(double value);

        protected:
            double _value;

    };

/**********************************************************************/

   /** \brief Single parameter changes applied to a single population
    *
    * \ingroup coalesce
    *
    */
    class PopulationParamChange : public SingleParamChange {
        public:

           /** \brief Default constructor
            *
            * The default date is 0., the default parameter value is 0.,
            * the  default population is 0
            *
            */
            PopulationParamChange();
            
           /** \brief Standard constructor
            *
            * \param date the event date.
            * \param population the population index.
            * \param value the parameter value.
            *
            */
            PopulationParamChange(double date, unsigned int population, double value);

            /// Gets the population index
            unsigned int population() const;
            
            /// Sets the population index
            void population(unsigned int value);
    
        protected:
        
            unsigned int _population;

    };

            
/**********************************************************************/
    
   /** \brief Bottleneck event
    *
    * \ingroup coalesce
    * 
    * The bottleneck parameter is its strength, corresponding to an
    * amount of time where time is locked and only coalescences are
    * allowed (resulting in a given - and random - number of
    * instantaneous coalescence with branches).
    *
    */
    class Bottleneck : public SingleParamChange {
        public:
            Bottleneck(double date, double param) : SingleParamChange(date, param) {}
            void apply(ParamSet* paramSet, Controller* controller) const;
    };

/**********************************************************************/
    
   /** \brief Population-specific bottleneck event
    *
    * \ingroup coalesce
    *
    */
    class PopulationBottleneck : public PopulationParamChange {
        public:
            PopulationBottleneck(double date, unsigned int population, double value) : PopulationParamChange(date, population, value) {}
            void apply(ParamSet* paramSet, Controller* controller) const;
    };


/**********************************************************************/
    
   /** \brief Change of the size of all populations
    *
    * The parameter is the new size (applied to all populations)
    *
    * \ingroup coalesce
    *
    */
    class AllPopulationSizeChange : public SingleParamChange {
        public:
            AllPopulationSizeChange(double date, double value) : SingleParamChange(date, value) {}
            void apply(ParamSet* paramSet, Controller* controller) const; 
    };
    
    
/**********************************************************************/
    
   /** \brief Change of a single population size
    *
    * \ingroup coalesce
    *
    */
    class SinglePopulationSizeChange : public PopulationParamChange {
        public:
            SinglePopulationSizeChange(double date, unsigned int population, double value) : PopulationParamChange(date, population, value) {}
            void apply(ParamSet* paramSet, Controller* controller) const;    
    };
    

    
/**********************************************************************/
    
   /** \brief Change of the growth rate of all populations
    *
    * The parameter is the new rate (applied to all populations)
    *
    * \ingroup coalesce
    *
    */
    class GrowthRateChange : public SingleParamChange {
        public:
            GrowthRateChange(double date, double value) : SingleParamChange(date, value) {}
            void apply(ParamSet* paramSet, Controller* controller) const;
    };
    
    
/**********************************************************************/
    
   /** \brief Change of a single population's growth rate
    *
    * \ingroup coalesce
    *
    */
    class PopulationGrowthRateChange : public PopulationParamChange {
        public:
            PopulationGrowthRateChange(double date, unsigned int population, double value) : PopulationParamChange(date, population, value) {}
            void apply(ParamSet* paramSet, Controller* controller) const;
    };
    
/**********************************************************************/
    
   /** \brief Change of the selfing rate
    *
    * The parameter is the new rate
    *
    * \ingroup coalesce
    *
    */
    class SelfingRateChange : public SingleParamChange {
        public:
            SelfingRateChange(double date, double value) : SingleParamChange(date, value) {}
            void apply(ParamSet* paramSet, Controller* controller) const;
    };
    
/**********************************************************************/
    
   /** \brief Fusion of two populations
    *
    * \ingroup coalesce
    *
    */
    class PopulationFusion : public Change {
        public:
           /** \brief Default constructor
            *
            * The default date is 0., the default mother is 0, the
            * default daughter is 0.
            *
            */
            PopulationFusion();

           /** \brief Standard constructor
            *
            * \param date the date of the event.
            * \param mother first population to merge.
            * \param daughter second population to merge.
            *
            * A time date, all the lineages from the daughter population
            * are moved to the mother population and all mutation rates
            * to the daughter population are cancelled. This functions
            * emulates a population split (forward in time).
            *
            */
            PopulationFusion(double date, unsigned int mother, unsigned int daughter);
            
            void apply(ParamSet* paramSet, Controller* controller) const;
        
            /// Sets the daughter population
            void daughter(unsigned int);
            
            /// Gets the daughter population
            unsigned int daughter() const;
            
            /// Sets the mother population
            void mother(unsigned int);
            
            /// Gets the mother population
            unsigned int mother() const;
        
        protected:
            unsigned int _mother;
            unsigned int _daughter;
    };    
    
    
/**********************************************************************/
    
   /** \brief Split of a population
    *
    * \ingroup coalesce
    *
    */
    class PopulationSplit : public Change {
        public:
           /** \brief Default constructor
            *
            * The default date is 0., the default population is 0, the
            * default probability is 0.5.
            *
            */
            PopulationSplit();

           /** \brief Standard constructor
            *
            * A the time given by date, the specified population is
            * split in two. An additional population (whose index is
            * incremented from the current total number of population)
            * is created and lineages are randomly picked and moved to
            * the new population. The parameter proba gives the
            * probability that a lineage from the population number pop
            * moves instantly to the new population. If proba is 0,
            * the program emulates the creation of an empty population
            * (thinking forward in time, this is a population
            * extinction). In general, forward in time, this is a
            * population fusion.
            * 
            * \param date the date of the event.
            * \param pop population index.
            * \param proba the probability that lineages move to the
            * new population.
            *
            */
            PopulationSplit(double date, unsigned int pop, double proba);
            
            void apply(ParamSet* paramSet, Controller* controller) const;
            
            /// Gets the population index
            unsigned int population() const;
            
            /// Sets the population index
            void population(unsigned int);
            
            /// Gets the probability of instant migration
            double probability() const;

            /// Sets the probability of instant migration
            void probability(double);
            
        protected:
            unsigned int _population;
            double _probability;
    };

/**********************************************************************/

   /** \brief Change of the migration rate of all population pairs
    *
    * The parameter is the new rate (applied to all population pairs)
    *
    * \ingroup coalesce
    *
    */
    class AllMigrationRateChange : public SingleParamChange {
        public:
        AllMigrationRateChange(double date, double value) : SingleParamChange(date, value) {}
        void apply(ParamSet* paramSet, Controller* controller) const;
    };

/**********************************************************************/
    
   /** \brief Change of a single migration rate
    *
    * \ingroup coalesce
    *
    */
    class SingleMigrationRateChange : public SingleParamChange {
        public:
           /** \brief Default constructor
            *
            * The default date is 0., the default parameter value is 0.,
            * the  default source population is 0, the default
            * destination population 1.
            *
            */
            SingleMigrationRateChange();

           /** \brief Standard constructor
            *
            * \param date the date of the event.
            * \param source index of the source population.
            * \param dest index of the destination population.
            * \param migr new value of the pairwise migration rate.
            *
            */
            SingleMigrationRateChange(double date, unsigned int source, unsigned int dest, double migr);
            
            /// Gets the source population index
            unsigned source() const;
            
            /// Sets the source population index
            void source(unsigned int);

            /// Gets the dest population index
            unsigned dest() const;
            
            /// Sets the dest population index
            void dest(unsigned int);

            void apply(ParamSet* paramSet, Controller* controller) const;
        
        protected:
            unsigned int _source;
            unsigned int _dest;
    };
}

#endif
