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

#include <sstream>
#include <cstdlib>
#include "Random.hpp"
#include "EggException.hpp"
#include "Arg.hpp"
#include "Population.hpp"

namespace egglib {

    void Arg::init_stable_parameters() {
        numberOfSegments = 0;
        _MRCA = NULL;
        segmentLengths = NULL;
        numberOfEdgesPerSegment = NULL;
    }
    
    
    void Arg::init_variable_parameters() {
        numberOfEdges = 0;
        numberOfSamples= 0;
        edges = NULL;
        _time = 0.;
        totalLength = 0;
        numberOfRecombinationEvents = 0;
        current = NULL;
    }


    Arg::Arg() {
        init_stable_parameters();
        init_variable_parameters();
    }


    Arg::Arg(Current* current, unsigned int numberOfSegments) {
        init_stable_parameters();
        init_variable_parameters();
        set(current, numberOfSegments);
    }


    Arg::~Arg() {
        clear();
    }


    void Arg::clear() {
        if (edges) free(edges);
        if (segmentLengths) free(segmentLengths);
        if (_MRCA) free(_MRCA);
        if (numberOfEdgesPerSegment) free(numberOfEdgesPerSegment);
        edgePool.clear();
        init_stable_parameters();
        init_variable_parameters();
    }


    void Arg::set(Current* current, unsigned int numberOfSegments) {

        edgePool.clear();

        edgePool.set(numberOfSegments, 2*current->totalNumberOfLineages());
        this->numberOfSegments = numberOfSegments;
        
        // initializes the list of segment lengths (their sum will be L)
        
        _MRCA = (Edge**) malloc(numberOfSegments*sizeof(Edge*));
        if (!_MRCA) throw EggMemoryError();
        
        segmentLengths = (double*) malloc(numberOfSegments*sizeof(double));
        if (!segmentLengths) throw EggMemoryError();

        // sets the list of ancestral nodes per segments
        
        numberOfEdgesPerSegment = (unsigned int*)
                        malloc(numberOfSegments* sizeof(unsigned int));
        if (!numberOfEdgesPerSegment) throw EggMemoryError();

        // prepares the ARG for simulation
        
        reset(current);
    }


    void Arg::reset(Current* current) {

        if (edges) free(edges);
        edgePool.releaseAll();
        init_variable_parameters();

        this->current = current;
        
        // collects the Edge addresses

        unsigned int c=0;
        for (unsigned int i=0; i<current->numberOfPopulations(); i++) {
            for (unsigned int j=0; j<current->populationNumberOfLineages(i); j++) {
                Edge* edge = edgePool.deliver();
                edge->set_terminal(c++);
                addEdge(edge);                // increments numberOfEdges
                current->population(i)->set(j, edge);
            }
        }
        
        numberOfSamples = numberOfEdges;

        for (unsigned int i=0; i<numberOfSegments; i++) {
            numberOfEdgesPerSegment[i] = numberOfEdges;
            _MRCA[i] = NULL;
            segmentLengths[i]  = 0;
        }

    }


    void Arg::addEdge(Edge* edge) {
        numberOfEdges++;
        edges = (Edge**) realloc(edges, numberOfEdges * sizeof(Edge*));
        if (!edges) throw EggMemoryError();
        edges[numberOfEdges-1] = edge;
    }


    double Arg::time() const {
        return _time;
    }
                

    void Arg::addTime(double increment) {
        _time+=increment;
    }


    void Arg::coalescence(double incr, unsigned int pop, unsigned int index1, unsigned int index2) {
        _time += incr;

        
        // argument check
        if (index1==index2) throw EggArgumentValueError("invalid lineage index (identical indices)");

        // picks coalescing lineages
        
        Edge* son1 = current->population(pop)->extractByIndex(index1);
        if (index2>index1) index2--;
        Edge* son2 = current->population(pop)->extractByIndex(index2);
        
        // coalesces
        Edge* father = edgePool.deliver();
        father->coalescence(_time, son1, son2, numberOfEdgesPerSegment,
                            _MRCA,  totalLength, segmentLengths);
        addEdge(father);
        current->population(pop)->push(father);
    }


    void Arg::coalescence(double incr, unsigned int pop, Random* random) {

        _time += incr;

        // argument check

        if (current->population(pop)->numberOfLineages()<2) throw EggRuntimeError("invalid population for coalescence: not enough lineages");

        // picks coalescing lineages

        Edge* son1 = current->population(pop)->extractRandomly(random);
        Edge* son2 = current->population(pop)->extractRandomly(random);

        // coalesces

        Edge* father = edgePool.deliver();
        father->coalescence(_time, son1, son2, numberOfEdgesPerSegment,
                                _MRCA, totalLength, segmentLengths);

        addEdge(father);
        current->population(pop)->push(father);

    }


    void Arg::recombination(double incr, Random* random) {

        _time += incr;
        numberOfRecombinationEvents += 1;
        
        // draws a random lineage
        unsigned int N = current->totalNumberOfLineages();
        unsigned int EN = current->efficientNumberOfLineages();
        unsigned int populationIndex= 9999;
        unsigned int lineageIndex = N; // this should generates an exception if Population performs a check
        unsigned int X = random->uniform()*EN;
        unsigned int c = 0, d;
        for (populationIndex=0; populationIndex<current->numberOfPopulations(); populationIndex++) {
            unsigned int n=current->populationNumberOfLineages(populationIndex);
            for (unsigned int individualIndex=0; individualIndex<n; individualIndex++) {
                d = current->population(populationIndex)->coverage(individualIndex);
                if (d>1) c += (d-1);
                if (X<c) {
                    lineageIndex = individualIndex;
                    break;
                }
            }
            if (lineageIndex != N) break;
        }
        
        // extracts the lineage
        Edge* son = current->population(populationIndex)->extractByIndex(lineageIndex);
        
        // performs the recombination
        Edge* father1 = edgePool.deliver();
        Edge* father2 = edgePool.deliver();
        son->recombination(_time, father1, father2, random,  totalLength, segmentLengths);

        // puts the things at the appropriate places
        addEdge(father1);
        addEdge(father2);
        current->population(populationIndex)->push(father1);
        current->population(populationIndex)->push(father2);

    }



    Edge* Arg::mute(unsigned int segment, double treePosition) {

        //Mutation mutation;
        // localizes the position of the mutation
        double acc = 0.;

        for (unsigned int i=0; i<numberOfEdges; i++) {

            if (edges[i]->segment(segment)) {
                    acc += edges[i]->length;
            }
            
            if (acc>treePosition) {
            
                //edges[i]->numberOfMutationsPerActualSite[segment]++;
            
//                mutation.age = edges[i]->bottom;
//                mutation.edge = edges[i];
//                mutation.segment = segment;
        
                // returns
                return edges[i];
            }
        }

        throw EggRuntimeError("hole in Arg::mute()");
        return NULL;
    }

        
    std::string Arg::newick(unsigned int segment) {
        const Edge* root = MRCA(segment);
        std::ostringstream sstream;
        sstream << "(";
        sstream << rnewick(root->son1, segment, 0.);
        sstream << ",";
        sstream << rnewick(root->son2, segment, 0.);
        sstream << ");";
        return sstream.str();
    }


    std::string Arg::rnewick(Edge* edge, unsigned int segment, double cache) {

        // increments cached branch length
        if ( edge->segment(segment)) {
                cache+= edge->length;
        }
        
        // predeclaration
        std::ostringstream sstream;
        Edge* son = NULL;
        bool seg1,seg2;
        
        switch (edge->numberOfSons) {

            // case of a terminal taxon
            case 0:

                // formats the string
                sstream << (edge->label())+1;
                sstream << ":";
                sstream << cache;

                return sstream.str();
        
            // case of an ancestral node
            case 2:
                
                // determines which descendant(s) have a segment with the specified position
                seg1 = edge->son1->segment(segment);
                seg2 = edge->son2->segment(segment);
                
                // if both sons are descending for this position, with make a new node
                if (seg1 && seg2) {
                    sstream << "(";
                    sstream << rnewick(edge->son1, segment, 0.);
                    sstream << ",";
                    sstream << rnewick(edge->son2, segment, 0.);
                    sstream << ")";
                    sstream << ":";
                    sstream << cache;
                    return sstream.str();
                }

                // if only one of sons is descending, then we just prolongate the current branch
                if (seg1) son = edge->son1;
                if (seg2) son = edge->son2;
                if (!son) throw EggRuntimeError("malformed ARG: interupted genealogy");
                return rnewick(son, segment, cache);
                
            // case of a recombined node, we prolongate the current branch
            case 1:
                return rnewick(edge->son1, segment, cache);
        }
        
        throw EggRuntimeError("maformed ARG: more than two descendants for a node");
        return "";
    }


    void Arg::set_actualNumberOfSites(unsigned int actualNumberOfSites) {
        for (unsigned int i=0; i<numberOfEdges; i++) {
            edges[i]->set_actualNumberOfSites(actualNumberOfSites);
        }
    }

}
