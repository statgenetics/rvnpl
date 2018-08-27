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


#include "EggException.hpp"
#include "ChangeTypes.hpp"
#include <cmath>


namespace egglib {


    /**********************************************************************/

    Change::Change() {
        _date = 0.;
    }

    Change::Change(double date) {
        if (date<0.) throw EggArgumentValueError("change date must not be negative");
        _date = date;
    }

    double Change::date() const {
        return _date;
    }

    void Change::date(double value) {
        if (value<0.) throw EggArgumentValueError("change date must not be negative");
        _date =  value;
    }


    /**********************************************************************/

    SingleParamChange::SingleParamChange() : Change() {
        _value = 0.;
    }
                
    SingleParamChange::SingleParamChange(double date, double value) : Change(date) {
        _value = value;
    }

    double SingleParamChange::value() const {
        return _value;
    }

    void SingleParamChange::value(double value) {
        _value = value;
    }

    /**********************************************************************/

    PopulationParamChange::PopulationParamChange() : SingleParamChange() {
        _population = 0;
    }
                
    PopulationParamChange::PopulationParamChange(double date, unsigned int population, double value) : SingleParamChange(date, value) {
        _population = population;
    }

    unsigned int PopulationParamChange::population() const {
        return _population;
    }

    void PopulationParamChange::population(unsigned int value) {
        _population = value;
    }


    /**********************************************************************/

    void Bottleneck::apply(ParamSet* paramSet, Controller* controller) const {
        for (unsigned int i=0; i<paramSet->numberOfPopulations(); i++) {
            controller->bottleneck(i, _value);
        }
    }


    /**********************************************************************/
        
    void PopulationBottleneck::apply(ParamSet* paramSet, Controller* controller) const {
        if (_population>=paramSet->numberOfPopulations()) throw EggArgumentValueError("invalid population index for bottleneck event");
        controller->bottleneck(_population, _value);
    }


    /**********************************************************************/

    void AllPopulationSizeChange::apply(ParamSet* paramSet, Controller* controller) const {
        for (unsigned int i=0; i<paramSet->numberOfPopulations(); i++) {
            paramSet->populationSize(i, _value);
            paramSet->dateOfLastChange(i, _date);
        }
    }
        
        
    /**********************************************************************/
        
    void  SinglePopulationSizeChange::apply(ParamSet* paramSet, Controller* controller) const {
        if (_population>=paramSet->numberOfPopulations()) throw EggArgumentValueError("invalid population index for population size change");
        paramSet->populationSize(_population, _value);
        paramSet->dateOfLastChange(_population, _date);
    }
        
        
    /**********************************************************************/

    void GrowthRateChange::apply(ParamSet* paramSet, Controller* controller) const {
        for (unsigned int i=0; i<paramSet->numberOfPopulations(); i++) {
            paramSet->populationSize(i,  paramSet->populationSize(i)
                    * exp(-paramSet->growthRate(i)*
                            (_date - paramSet->dateOfLastChange(i))) );
            paramSet->growthRate(i, _value);
            paramSet->dateOfLastChange(i, _date);
        }
    }
        
        
    /**********************************************************************/
        
    void  PopulationGrowthRateChange::apply(ParamSet* paramSet, Controller* controller) const {
        if (_population>=paramSet->numberOfPopulations()) throw EggArgumentValueError("invalid population index for growth rate change");
        paramSet->populationSize(_population,  paramSet->populationSize(_population)
                * exp(-paramSet->growthRate(_population)*
                        (_date - paramSet->dateOfLastChange(_population))) );
        paramSet->growthRate(_population, _value);
        paramSet->dateOfLastChange(_population, _date);
    }
        
        
    /**********************************************************************/
        
    void SelfingRateChange::apply(ParamSet* paramSet, Controller* controller) const {
        if (_value<0. || _value>1.) throw EggArgumentValueError("invalid value for selfing rate");
        paramSet->selfingRate(_value);
    }
        
        
    /**********************************************************************/
        
    PopulationFusion::PopulationFusion() : Change() {
        _mother = 0;
        _daughter = 0;
    }

    PopulationFusion::PopulationFusion(double date, unsigned int mother, unsigned int daughter) : Change(date) {
        _mother = mother;
        _daughter = daughter;
    }

    void PopulationFusion::daughter(unsigned int value) {
        _daughter = value;
    }

    unsigned int PopulationFusion::daughter() const {
        return _daughter;
    }

    void PopulationFusion::mother(unsigned int value) {
        _mother = value;
    }

    unsigned int PopulationFusion::mother() const {
        return _mother;
    }

    void PopulationFusion::apply(ParamSet* paramSet, Controller* controller) const {

        // checks that the population indices are correct
        if (_mother==_daughter) throw EggArgumentValueError("invalid population fusion parameters (population indices are identical)");
        if (_mother>=paramSet->numberOfPopulations()) throw EggArgumentValueError("invalid population fusion parameters (mother population out of range)");
        if (_daughter>=paramSet->numberOfPopulations()) throw EggArgumentValueError("invalid population fusion parameters (daughter population out of range)");
        
        // cancels all migration rates to daughter (from wherever)
        for (unsigned int i=0; i<paramSet->numberOfPopulations(); i++) {
            if (i==_daughter) continue;
            paramSet->pairwiseMigrationRate(i, _daughter, 0.);
        }
        
        // migrates all lineages
        controller->moveAllLineages(_daughter, _mother);
    }

        
    /**********************************************************************/
        
    PopulationSplit::PopulationSplit() : Change() {
        _population = 0;
        _probability = 0.5;
    }

    PopulationSplit::PopulationSplit(double date, unsigned int pop, double proba) : Change(date) {
        _population = pop;
        _probability = proba;
    }
                
    unsigned int PopulationSplit::population() const {
        return _population;
    }
                
    void PopulationSplit::population(unsigned int pop) {
        _population = pop;
    }
                
    double PopulationSplit::probability() const {
        return _probability;
    }

    void PopulationSplit::probability(double proba) {
        _probability = proba;
    }
                
    void PopulationSplit::apply(ParamSet* paramSet, Controller* controller) const {
        if (_population>=paramSet->numberOfPopulations()) throw EggArgumentValueError("invalid population fusion parameters (mother population out of range)");
        double migr=0.;
        if (paramSet->numberOfPopulations()==1) migr = paramSet->pairwiseMigrationRate(0, 0);
        else migr = paramSet->pairwiseMigrationRate(0, 1);
        paramSet->addPopulation(migr);
        controller->addPopulation();
        controller->moveSomeLineages(_population, paramSet->numberOfPopulations()-1, _probability);
        paramSet->dateOfLastChange(paramSet->numberOfPopulations()-1, _date);
    }

        
    /**********************************************************************/

    void AllMigrationRateChange::apply(ParamSet* paramSet, Controller* controller) const {
        if (_value<0.) throw EggArgumentValueError("invalid value for migration rate");
        for (unsigned int i=0; i<paramSet->numberOfPopulations(); i++) {
            paramSet->migrationRate(_value);
        }
    }


    /**********************************************************************/
        
    SingleMigrationRateChange::SingleMigrationRateChange() : SingleParamChange() {
        _source = 0;
        _dest = 1;
    }

    SingleMigrationRateChange::SingleMigrationRateChange(double date, unsigned int source, unsigned int dest, double migr) : SingleParamChange(date, migr) {
        _source = source;
        _dest = dest;
    }
                
    unsigned SingleMigrationRateChange::source() const {
        return _source;
    }
                
    void SingleMigrationRateChange::source(unsigned int value) {
        _source = value;
    }

    unsigned SingleMigrationRateChange::dest() const {
        return _dest;
    }
                
    void SingleMigrationRateChange::dest(unsigned int value) {
        _dest = value;
    }

    void SingleMigrationRateChange::apply(ParamSet* paramSet, Controller* controller) const {
        if (_value<0.) throw EggArgumentValueError("invalid value for migration rate");
        if (_source>=paramSet->numberOfPopulations()) throw EggArgumentValueError("cannot set migration rate (source population out of bounds)");
        if (_dest>=paramSet->numberOfPopulations()) throw EggArgumentValueError("cannot set migration rate (dest population out of bounds)");
        if (_source==_dest) throw EggArgumentValueError("cannot set migration rate (source and dest populations are the same!)");
        paramSet->pairwiseMigrationRate(_source, _dest, _value);
    }

}
