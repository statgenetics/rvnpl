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
#include "ParamSet.hpp"
#include "ChangeTypes.hpp"
#include "Controller.hpp"
#include <cstdlib>



namespace egglib {

    ParamSet::ParamSet() {
        init();
        addPopulation(0.);
    }


    ParamSet::~ParamSet() {
        clear();
    }


    void ParamSet::init() {
        // sets the default values / assuming the object is new or was cleared
        _selfingRate = 0.;
        _recombinationRate = 0.;
        _numberOfSegments = 1;
        _numberOfPopulations = 0;
        _singles = NULL;
        _doubles = NULL;
        _growthRates = NULL;
        _populationSize = NULL;
        _dateOfLastChange = NULL;
        migrationMatrix = NULL;
        _numberOfChanges = 0;
        nextChangeIndex = 0;
        changes = NULL;
    }


    void ParamSet::clear() {
        // free everything
        if (_numberOfPopulations) {
            if (migrationMatrix) {
                for (unsigned int i=0; i<_numberOfPopulations; i++) if (migrationMatrix[i]) free(migrationMatrix[i]);
                free(migrationMatrix);
            }
            if (_singles) free(_singles);
            if (_doubles) free(_doubles);
            if (_growthRates) free(_growthRates);
            if (_populationSize) free(_populationSize);
            if (_dateOfLastChange) free(_dateOfLastChange);
        }
        if (!_numberOfChanges && changes) throw EggRuntimeError("possible memory leak in ParamSet");
        if (_numberOfChanges) free(changes);
    }


    ParamSet::ParamSet(const ParamSet& source) {
        init();
        copy(source);
    }

                
    ParamSet& ParamSet::operator=(const ParamSet& source) {
        clear();
        init();
        copy(source);
        return *this;
    }


    void ParamSet::reset() {
        clear();
        init();
    }

                        
    void ParamSet::copy(const ParamSet& source) {
        // everything  must have been cleared (or be new)
        // copies the data
        _numberOfPopulations = source._numberOfPopulations;
        _selfingRate = source._selfingRate;
        _recombinationRate = source._recombinationRate;
        _numberOfSegments = source._numberOfSegments;
            
        // allocates all the population-based arrays
        _singles = (unsigned int*) malloc(_numberOfPopulations * sizeof(unsigned int));
        if (!_singles) throw EggMemoryError();
            
        _doubles = (unsigned int*) malloc(_numberOfPopulations * sizeof(unsigned int));
        if (!_doubles) throw EggMemoryError();

        _growthRates = (double*) malloc(_numberOfPopulations * sizeof(double));
        if (!_growthRates) throw EggMemoryError();

        _populationSize = (double*) malloc(_numberOfPopulations * sizeof(double));
        if (!_populationSize) throw EggMemoryError();

        _dateOfLastChange = (double*) malloc(_numberOfPopulations * sizeof(double));
        if (!_dateOfLastChange) throw EggMemoryError();

        migrationMatrix = (double**) malloc(_numberOfPopulations * sizeof(double*));
        if (!migrationMatrix) throw EggMemoryError();
        for (unsigned int i=0; i<_numberOfPopulations; i++) {
            migrationMatrix[i] = (double*) malloc(_numberOfPopulations * sizeof(double));
            if (!migrationMatrix[i]) throw EggMemoryError();        
        }
        
        // copies the data
        for (unsigned int i=0; i<_numberOfPopulations; i++) {
            _singles[i] = source._singles[i];
            _doubles[i] = source._doubles[i];
            _growthRates[i] = source._growthRates[i];
            _populationSize[i] = source._populationSize[i];
            _dateOfLastChange[i] = source._dateOfLastChange[i];
            for (unsigned int j=0; j<_numberOfPopulations; j++) {
                migrationMatrix[i][j] = source.migrationMatrix[i][j];
            }
        }
        
        // copies Change data
        _numberOfChanges = source._numberOfChanges;
        nextChangeIndex = source.nextChangeIndex;
        if (_numberOfChanges) {
            changes = (Change const**) malloc(_numberOfChanges * sizeof(Change const*));
            if (!changes) throw EggMemoryError();
            for (unsigned int i=0; i<_numberOfChanges; i++) changes[i] = source.changes[i];
        }
        else changes=NULL;
    }


    void ParamSet::addPopulation(double migrationRate) {
        
        // if the current matrix is 1x1, the value in the migration matrix is irrelevant
        if (_numberOfPopulations==1) migrationMatrix[0][0]=0;
        
        _numberOfPopulations++;

        // reallocs population-specific arrays and sets reasonnable default value for new cell
        _singles = (unsigned int*) realloc(_singles, _numberOfPopulations * sizeof(unsigned int));
        if (!_singles) throw EggMemoryError();
        _singles[_numberOfPopulations-1] = 0;

        _doubles = (unsigned int*) realloc(_doubles, _numberOfPopulations * sizeof(unsigned int));
        if (!_doubles) throw EggMemoryError();
        _doubles[_numberOfPopulations-1] = 0;

        _populationSize = (double*) realloc(_populationSize, _numberOfPopulations * sizeof(double));
        if (!_populationSize) throw EggMemoryError();
        _populationSize[_numberOfPopulations-1] = 1.;
        
        _growthRates = (double*) realloc(_growthRates, _numberOfPopulations * sizeof(double));
        if (!_growthRates) throw EggMemoryError();
        _growthRates[_numberOfPopulations-1] = 0.;

        _dateOfLastChange = (double*) realloc(_dateOfLastChange, _numberOfPopulations * sizeof(double));
        if (!_dateOfLastChange) throw EggMemoryError();
        _dateOfLastChange[_numberOfPopulations-1] = 0.;
        
        // enlarges the migration matrix. but the passed rate in new cells and increments all diagonals
        migrationMatrix = (double**) realloc(migrationMatrix, _numberOfPopulations * sizeof(double*));
        if (!migrationMatrix) throw EggMemoryError();
        migrationMatrix[_numberOfPopulations-1] = (double*) malloc(_numberOfPopulations * sizeof(double));
        if (!migrationMatrix[_numberOfPopulations-1]) throw EggMemoryError();
        for (unsigned int i=0; i<(_numberOfPopulations-1); i++) {
            migrationMatrix[i] = (double*) realloc(migrationMatrix[i], _numberOfPopulations * sizeof(double));
            if (!migrationMatrix[i]) throw EggMemoryError();
            migrationMatrix[i][_numberOfPopulations-1] = migrationRate;
            migrationMatrix[_numberOfPopulations-1][i] = migrationRate;
            migrationMatrix[i][i]+= migrationRate;
        }
        migrationMatrix[_numberOfPopulations-1][_numberOfPopulations-1] = migrationRate * (_numberOfPopulations-1);    

    }


    unsigned int ParamSet::numberOfPopulations() const {
        return _numberOfPopulations;
    }

            
    double ParamSet::pairwiseMigrationRate(unsigned int source, unsigned int dest) const {
        if (source>=_numberOfPopulations || dest>=_numberOfPopulations) throw EggArgumentValueError("invalid population index");
        return migrationMatrix[source][dest];
    }

     
    void ParamSet::pairwiseMigrationRate(unsigned int source, unsigned int dest, double value) {
        if (source>=_numberOfPopulations || dest>=_numberOfPopulations || source==dest) {
            throw EggArgumentValueError("invalid population index");
        }
        migrationMatrix[source][source] += (value - migrationMatrix[source][dest]);
        migrationMatrix[source][dest] = value;
    }

                
    void ParamSet::migrationRate(double value) {
        for (unsigned int i=0; i<_numberOfPopulations; i++) {
            for (unsigned int j=0; j<_numberOfPopulations; j++) {
                if (i==j) migrationMatrix[i][j] = value;
                else migrationMatrix[i][j] = value / (_numberOfPopulations-1.);
            }
        }
    }

                
    double ParamSet::populationSize(unsigned int populationIndex) const {
        if (populationIndex>=_numberOfPopulations) throw EggArgumentValueError("invalid population index");
        return _populationSize[populationIndex];
    }

                
    void ParamSet::populationSize(unsigned int populationIndex, double value) {
        if (populationIndex>=_numberOfPopulations) throw EggArgumentValueError("invalid population index");
        if (value==0) throw EggArgumentValueError("population size must be non-null");
        _populationSize[populationIndex] = value;
    }

                
    double ParamSet::growthRate(unsigned int populationIndex) const {
        if (populationIndex>=_numberOfPopulations) throw EggArgumentValueError("invalid population index");
        return _growthRates[populationIndex];
    }


    void ParamSet::growthRate(unsigned int populationIndex, double value) {
        if (populationIndex>=_numberOfPopulations) throw EggArgumentValueError("invalid population index");
        _growthRates[populationIndex] = value;
    }


    double ParamSet::recombinationRate() const {
        return _recombinationRate;
    }
                

    void ParamSet::recombinationRate(double value) {
        if (value<0.) throw EggArgumentValueError("invalid recombination rate");
        if (_numberOfSegments<2) throw EggArgumentValueError("it is not allowed to set the recombination rate when there is <2 segments");
        _recombinationRate = value;
    }


    unsigned int ParamSet::numberOfSegments() const {
        return _numberOfSegments;
    }
                

    void ParamSet::numberOfSegments(unsigned int value) {
        if (value<1) throw EggArgumentValueError("invalid number of recombining segmentss");
        if (value<2 && _recombinationRate>0) throw EggArgumentValueError("it is not allowed to set <2 segments when the recombination rate is >0");
        _numberOfSegments = value;
    }


    double ParamSet::selfingRate() const {
        return _selfingRate;
    }


    void ParamSet::selfingRate(double value) {
        if (value<0. || value>1.) throw EggArgumentValueError("invalid selfing rate");
        _selfingRate = value;
    }
                

    void ParamSet::addChange(const Change* change) {
        _numberOfChanges++;
        // sets the pointer to this change, if required
        if (_numberOfChanges==1 || change->date()<changes[nextChangeIndex]->date()) {
            nextChangeIndex = _numberOfChanges - 1;
        }
        // enlarges the array
        changes = (Change const**) realloc(changes, _numberOfChanges * sizeof(Change const*));
        if (!changes) throw EggMemoryError();
        // copy the address
        changes[_numberOfChanges-1] = change;
    }


    double ParamSet::nextChangeDate() const {
        if (_numberOfChanges==0) return -1;
        return changes[nextChangeIndex]->date();
    }
                

    void ParamSet::nextChangeDo(Controller* controller) {

        if (_numberOfChanges==0) throw EggRuntimeError("try to do a demographic change when there is none");
        
        // applies the change
        changes[nextChangeIndex]->apply(this, controller);

        // removes the change (only disconnecting)
        changes[nextChangeIndex] = NULL;

        // searches for the next change
        unsigned int acc=0;
        unsigned int c=nextChangeIndex; // makes an iteration with shift
        for (unsigned int i=0; i<_numberOfChanges; i++) {
            c++;
            if (c>=_numberOfChanges) c=0;
            if (changes[c]==NULL) continue;
            acc++;
            if (acc==1 || changes[c]->date()<changes[nextChangeIndex]->date()) nextChangeIndex = c;
        }

        // if no next change, cleans the array
        if (acc==0) {
            _numberOfChanges = 0;
            free(changes);
            changes = NULL;
            nextChangeIndex = 0;
        }

    }


    unsigned int ParamSet::singles(unsigned int populationIndex) const {
        if (populationIndex>=_numberOfPopulations) throw EggArgumentValueError("invalid population index");
        return _singles[populationIndex];
    }
        
        
    void ParamSet::singles(unsigned int populationIndex, unsigned int value) {
        if (populationIndex>=_numberOfPopulations) throw EggArgumentValueError("invalid population index");
        _singles[populationIndex] = value;
    }
        

    unsigned int ParamSet::doubles(unsigned int populationIndex) const {
        if (populationIndex>=_numberOfPopulations) throw EggArgumentValueError("invalid population index");
        return _doubles[populationIndex];
    }


    void ParamSet::doubles(unsigned int populationIndex, unsigned int value) {
        if (populationIndex>=_numberOfPopulations) throw EggArgumentValueError("invalid population index");
        _doubles[populationIndex] = value;
    }

    unsigned int ParamSet::numberOfSamples() const {
        unsigned int c=0;
        for (unsigned int i=0; i<_numberOfPopulations; i++) {
            c+= _singles[i] + 2*_doubles[i];
        }
        return c;
    }
                
    double ParamSet::dateOfLastChange(unsigned int populationIndex) const {
        if (populationIndex>=_numberOfPopulations) throw EggArgumentValueError("invalid population index");
        return _dateOfLastChange[populationIndex];
    }


    void ParamSet::dateOfLastChange(unsigned int populationIndex, double date) const {
        if (populationIndex>=_numberOfPopulations) throw EggArgumentValueError("invalid population index");
        _dateOfLastChange[populationIndex] = date;
    }


    void ParamSet::setGroups(DataMatrix& dataMatrix, bool labelIndividuals) {

        if (dataMatrix.numberOfSequences()!=numberOfSamples()) {
            throw EggRuntimeError("cannot set groups of DataMatrix object: the number of samples is incorrect");
        }
        
        // case 1: label populations
        unsigned int sample=0;
        if (!labelIndividuals) {
            for (unsigned int pop=0; pop<_numberOfPopulations; pop++) {
                for (unsigned int i=0; i<(_singles[pop]+2*_doubles[pop]); i++) {
                    dataMatrix.populationLabel(sample++, pop);
                }
            }
        }
        
        // case 2: label individuals
        else {
            unsigned int group=0;
            for (unsigned int pop=0; pop<_numberOfPopulations; pop++) {
                for (unsigned int i=0; i<_doubles[pop]; i++) {
                    dataMatrix.populationLabel(sample++, group);
                    dataMatrix.populationLabel(sample++, group++);
                }
                for (unsigned int i=0; i<_singles[pop]; i++) {
                    dataMatrix.populationLabel(sample++, group++);
                }
            }
        }
    }

}
