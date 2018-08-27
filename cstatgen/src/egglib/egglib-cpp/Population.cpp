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


#include "Population.hpp"
#include "EggException.hpp"
#include "Random.hpp"
#include <cstdlib>


namespace egglib {

    Population::Population() {
        init();
    }

    Population::Population(const Population& source) {
        init();
        copy(source);
    }

    Population& Population::operator=(const Population& source) {
        clear();
        init();
        copy(source);
        return *this;
    }
                
    Population::~Population() {
        clear();
    }


    void Population::copy(const Population& source) {
        _numberOfLineages = source._numberOfLineages;
        _efficientNumberOfLineages = source._efficientNumberOfLineages;
        if (!_numberOfLineages) {
            return;
        }
        lineages = (Edge**) malloc(_numberOfLineages * sizeof(Edge*));
        if (!lineages) {
            throw EggMemoryError();
        }
        for (unsigned int i=0; i<_numberOfLineages; i++) {
            lineages[i] = source.lineages[i];
        }
    }


    void Population::init() {
        _numberOfLineages = 0;
        _efficientNumberOfLineages = 0;
        lineages = NULL;
    }


    void Population::clear() {
        if (_numberOfLineages) {
           if (!lineages) throw EggRuntimeError("inconsistent content of Population");
            free(lineages);
            lineages=NULL;
        }
        else if (lineages) throw EggRuntimeError("possible memory leak in Population");
    }



    Population::Population(unsigned int numberOfSegments,
                    unsigned int numberOfLineages, unsigned firstIndex) {

        init();

        _numberOfLineages = numberOfLineages;

        if (!_numberOfLineages) {
            return;
        }

        lineages = (Edge**) malloc(_numberOfLineages * sizeof(Edge*));
        if (!lineages) {
            throw EggMemoryError();
        }

        for (unsigned int i=0; i<_numberOfLineages; i++) {
            lineages[i] = NULL;
        }
    }


    unsigned int Population::numberOfLineages() const {
        return _numberOfLineages;
    }

    unsigned int Population::efficientNumberOfLineages() const {
        return _efficientNumberOfLineages;
    }

    void Population::set(unsigned int index, Edge* edge) {
        if (index>=_numberOfLineages) throw EggArgumentValueError("invalid lineage index");
        lineages[index] = edge;
        if (edge->coverage>1) _efficientNumberOfLineages += (edge->coverage-1);
    }


    Edge* Population::extractRandomly(Random* random) {
        if (_numberOfLineages==0) throw EggRuntimeError("cannot pick lineage: empty population");
        if (_numberOfLineages==1) return pick(0);
        return pick(random->irand(_numberOfLineages));
    }


    Edge* Population::extractByIndex(unsigned int index) {
        if (index>=_numberOfLineages) throw EggArgumentValueError("invalid lineage index");
        return pick(index);
    }


    Edge* Population::pick(unsigned int index) {
        Edge* cache = lineages[index];
        for (unsigned int i=index; i<(_numberOfLineages-1); i++) {
            lineages[i] = lineages[i+1];
        }
        _numberOfLineages--;
        if (_numberOfLineages==0) {
            free(lineages);
            lineages=NULL;
        }
        else {
            lineages = (Edge**) realloc(lineages, _numberOfLineages * sizeof(Edge*));
            if (!lineages) {
                throw EggMemoryError();
            }
        }
        if (cache->coverage>1) _efficientNumberOfLineages -= (cache->coverage-1);
        return cache;
    }


    void Population::push(Edge* edge) {
        if (!edge) throw EggRuntimeError("received a null address as edge");
        _numberOfLineages++;
        lineages = (Edge**) realloc(lineages, _numberOfLineages * sizeof(Edge*));
        if (!lineages) {
            throw EggMemoryError();
        }
        if (edge->coverage>1) _efficientNumberOfLineages += (edge->coverage-1);
        lineages[_numberOfLineages-1] = edge;
    }

    unsigned int Population::coverage(unsigned int edgeIndex) const {
        if (edgeIndex>=_numberOfLineages) throw EggArgumentValueError("invalid lineage index");
        return lineages[edgeIndex]->coverage;
    }

}

