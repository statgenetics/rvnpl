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


#include "Current.hpp"
#include "ParamSet.hpp"
#include "Population.hpp"
#include "EggException.hpp"
#include <cstdlib>


namespace egglib {

    void Current::setPopulationArray() {
        if (!_numberOfPopulations) {
            throw EggRuntimeError("Current cannot process zero populations");
        }
        populations = (Population**) malloc(sizeof(Population*) * _numberOfPopulations);
        if (!populations) throw EggMemoryError();
    }


    Current::Current(ParamSet* paramSet) {
        reset(paramSet);
    }


    Current::Current() {
        _numberOfPopulations = 0;
        populations = NULL;
        _numberOfSegments = 1;
    }

    Current::Current(const Current& source) {
        _numberOfPopulations = 0;
        _numberOfSegments = 1;
        populations = NULL;
        copy(source);
    }


    void Current::reset(ParamSet* paramSet) {
        clear();
        _numberOfSegments = paramSet->numberOfSegments();
        _numberOfPopulations = paramSet->numberOfPopulations();
        setPopulationArray();
        unsigned int c=0;
        for (unsigned int i=0; i<_numberOfPopulations; i++) {
            unsigned int n = paramSet->singles(i) + 2*paramSet->doubles(i);
            populations[i] = new Population(
                paramSet->numberOfSegments(), n, c);
            c+=n;
        }
    }


    Current& Current::operator=(const Current& source) {
        clear();
        copy(source);
        return *this;
    }


    void Current::copy(const Current& source) {
        _numberOfPopulations = source._numberOfPopulations;
        _numberOfSegments  = source._numberOfSegments;
        setPopulationArray();
        for (unsigned int i=0; i<_numberOfPopulations; i++) {
            populations[i] = new Population(*source.populations[i]);
        }
    }


    void Current::clear() {
        if (_numberOfPopulations) {
            if (populations) {
                for (unsigned int i=0; i<_numberOfPopulations; i++) delete populations[i];
                free(populations);
            }
        }
    }

                
    Current::~Current() {
        clear();
    }


    unsigned int Current::numberOfPopulations() const {
        return _numberOfPopulations;
    }


    void Current::addPopulation() {
        if (!_numberOfPopulations) throw EggRuntimeError("Current cannot process zero populations");
        _numberOfPopulations++;
        populations = (Population**) realloc(populations, sizeof(Population*) * _numberOfPopulations);
        if (!populations) throw EggMemoryError();
        populations[_numberOfPopulations-1] =
                new Population(_numberOfSegments, 0, 0);
    }

                
    unsigned int Current::populationNumberOfLineages(unsigned int populationIndex) const {
        if (populationIndex>=_numberOfPopulations) throw EggArgumentValueError("invalid population index");
        return populations[populationIndex]->numberOfLineages();
    }


    unsigned int Current::totalNumberOfLineages() const {
        unsigned int total=0;
        for (unsigned int i=0; i<_numberOfPopulations; i++) {
            total+= populations[i]->numberOfLineages();
        }
        return total;
    }


    unsigned int Current::efficientNumberOfLineages() const {
        unsigned int total=0;
        for (unsigned int i=0; i<_numberOfPopulations; i++) {
            total+= populations[i]->efficientNumberOfLineages();
        }
        return total;
    }


    Population* Current::population(unsigned int populationIndex) {
        if (populationIndex>=_numberOfPopulations) throw EggArgumentValueError("invalid population index");
        return populations[populationIndex];
    }

}
