/*
    Copyright 2009-2011 Stéphane De Mita, Mathieu Siol

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


#include "Random.hpp"
#include "Edge.hpp"
#include <cstdlib>

namespace egglib {


    void Edge::init(unsigned int numberOfSegments) {

        this->numberOfSegments = numberOfSegments;

        if (numberOfSegments==0) {
            segments = NULL;
            numberOfMutationsPerActualSite = NULL;
        }
        else {
            segments = (unsigned int*) malloc(numberOfSegments*sizeof(unsigned int));

            segments[0] = numberOfSegments;
            if (!segments) throw EggMemoryError();
            for (unsigned int i=1; i<numberOfSegments; i++) {
                segments[i] = 0;
            }
            segbools = (bool*) malloc(numberOfSegments*sizeof(int));
            if (!segbools) {
                free (segments);
                throw EggMemoryError();
            }
            segbools[0] = false;
            
        }

        top = 0.;
        bottom = 0.;
        length = 0.;
        son1 = NULL;
        son2 = NULL;
        coverage = 0;
        _label = 0;
        numberOfSons = 0;
        numberOfMutationsPerActualSite = NULL;

    }
    
    Edge::Edge() {
        init(0);
    }


    Edge::Edge(unsigned int numberOfSegments) {
        init(numberOfSegments);
    }


    Edge::~Edge() {
        if (segments) free(segments);
        if (segbools) free(segbools);
        if (numberOfMutationsPerActualSite) free(numberOfMutationsPerActualSite);
    }
    

    void Edge::set_terminal(unsigned int leaf_index) {
        _label  = leaf_index;
        segments[0] = numberOfSegments;
        for (unsigned int i=1; i<numberOfSegments; i++) {
            segments[i] = -1;
        }
        segbools[0] = true;
        coverage = numberOfSegments;
    }


    void Edge::reset() {

        unsigned int cur = 0;
        unsigned int cache;
        while (true) {
            cache = segments[cur];
            if (segbools[cur]==true) {
                for (unsigned int i=1; i<cache; i++) {
                    segments[cur+i] = 0;
                }
                segbools[cur]=false;
            }
            segments[cur] = 0;
            cur += cache;
            if (cur>=numberOfSegments) break;
        } 

        segments[0] = numberOfSegments;
        segbools[0] = false;

        //segments[0] = numberOfSegments;
        //segbools[0] = true;

        top = 0.;
        bottom = 0.;
        length = 0.;
        son1 = NULL;
        son2 = NULL;
        coverage = 0;
        _label = 0;
        numberOfSons = 0;

        if (numberOfMutationsPerActualSite) free(numberOfMutationsPerActualSite);
        numberOfMutationsPerActualSite = NULL;
        
    }
    
        

    void Edge::coalescence(double date, Edge* son1, Edge* son2,
                        unsigned int* edgesPerSegments, Edge** MRCA,
                        double& totalLength, double* segmentLengths) {

        top = date;
        bottom = date;
        numberOfSons = 2;
        
        // connection with nodes below

        this->son1 = son1;
        this->son2 = son2;
        son1->top = bottom;
        son2->top = bottom;
        son1->length = son1->top - son1->bottom;
        son2->length = son2->top - son2->bottom;
        son1->update_lengths(totalLength, segmentLengths);
        son2->update_lengths(totalLength, segmentLengths);

        // merges sons segments

        coverage = 0;
        
        unsigned int i1, i2, incr, c1, c2;
        bool cur1, cur2, b1, b2;
        
        // initializes counters
        
        i1 = 0;
        i2 = 0;
        c1 = son1->segments[0];
        c2 = son2->segments[0];
        b1 = son1->segbools[0];
        b2 = son2->segbools[0];
        cur1 = b1|b2;

        // special case: lMRCA found
        if ((b1&b2) && edgesPerSegments[0]==2) {   /* will be 1 */
            cur1 = false; // count as non-ancestral
        }

        cur2 = cur1;

        while (i1<numberOfSegments) {

            // goes on while the status is consistent

            while (cur1==cur2) { // look ahead until change of status

                // step until end of shortest segment/domain

                incr = c1<c2?c1:c2;
                
                // if both segments are covered,
                //  decrement ancestral-lineage counter and check for MRCA
                
                if (b1&b2) {

                    // decrements the nb of ancestral lineages
                    
                    bool bMRCA = (edgesPerSegments[i2]==2);
                    
                    for (unsigned int j=0; j<incr; j++) {
                        
                        // possibility of change of status because of ancestrality
                        if (bMRCA != (edgesPerSegments[i2+j]==2)) {
                            incr = j;
                            break;   // don't process now
                        }

                        // decrement counter and detect lMRCA
                        edgesPerSegments[i2+j]--;
                        if (edgesPerSegments[i2+j]==1) MRCA[i2+j] = this;
                    }
                    
                    // if status regarding MRCA has changed, j will be smaller than incr

                }
             
                // moves the cursor (it will not necessarily fall on a domain)
                
                i2 += incr;
                c1 -= incr;
                c2 -= incr;
                
                if (i2==numberOfSegments) { // touches the end
                    break;
                }

                // if either one segment is finished, update info
                // if cx is 0, i2 must be on a new range

                if (c1==0) {
                    c1 = son1->segments[i2];
                    b1 = son1->segbools[i2];
                }

                if (c2==0) {
                    c2 = son2->segments[i2];
                    b2 = son2->segbools[i2];
                }
                
                // 'new' status (might not change) (used stored booleans)
                
                cur2 = b1|b2;
                if ((b1&b2) && edgesPerSegments[i2]==2) cur2 = false;
            }
            
            // enter cell (domain, range) data
            
            segments[i1] = i2-i1;
            if (cur1) segbools[i1] = true;
            else segbools[i1] = false;

            // forward main cursor

            i1 = i2;
            cur1 = cur2;

        }

        // VERY IMPORTANT

        fill();

    }




    void Edge::recombination(double date, Edge* dest1, Edge* dest2, Random* random,
                            double& totalLength, double* segmentLengths) {


        // draws a random position within the interval
        
        unsigned int X = random->irand(coverage-1);
        
        unsigned int i=0;
        unsigned int acc=0;
        unsigned int point=numberOfSegments;
        while(i<numberOfSegments) {

            if (segbools[i]) {
                acc+=segments[i];

                if (acc>X) {
                    point = i + X - (acc-segments[i]);
                    break;
                }
            }
            i += segments[i];

        }

       
        if (point>=numberOfSegments) {
            throw EggRuntimeError("hole in Edge::recombination()!");
        }

        // sets destinations (ancestors) information

        dest1->numberOfSons = 1;
        dest2->numberOfSons = 1;
        dest1->coverage = 0;
        dest2->coverage = 0;

        // copies all segments up to the recombination point
        
        unsigned int length_curr;
        bool sign;
        unsigned int pos = 0;
        
        while (pos<=point) {
            length_curr = segments[pos];
            sign = segbools[pos];
            dest1->segments[pos] = length_curr;
            dest1->segbools[pos] = sign;
            if (sign) dest1->coverage += length_curr;
            pos += length_curr;
        }
        
        // if a segment is covering the point, moves back to it
        
        if (pos>point) pos -= length_curr;

        // if last segment covered, crop it

        if (sign) {
            unsigned int overhang = pos+length_curr - (point+1);
            dest1->segments[pos] -= overhang;

            // sets the last segment
            dest1->segments[point+1] = numberOfSegments-point-1;
            dest1->segbools[point+1] = false;
        }
        
        // if last segment not covered, extends it to the end
        
        else {
            dest1->segments[pos] += (numberOfSegments - (pos+length_curr));
        }

        // no need to check that the recomb point is left-hand covered

        // now, pos should be at point or earlier

        // makes a pseudo-segment with only what is over than point
        
        length_curr = length_curr - (point-pos+1);
        pos = point+1;

        // if the pseudo segment is null (previous segment stops exactly at point)
        
        if (length_curr==0) {
            sign = !sign;
            length_curr = segments[pos];
        }

        // sets the left-hand empy segment
        
        dest2->segments[0] = point+1;
        dest2->segbools[0] = false;
        
        if (!sign) {
            dest2->segments[0] += length_curr;
            pos += length_curr;
            length_curr = segments[pos];
            sign = true;
        }
        
        // add the pseudo segment
        
        if (sign) {
            dest2->segments[pos] = length_curr;
            dest2->segbools[pos] = true;
        }
        else dest2->segments[pos] = -length_curr;
        pos += length_curr;

        // copies all other segments
        
        while (pos<numberOfSegments) {
            dest2->segments[pos] = segments[pos];
            dest2->segbools[pos] = segbools[pos];
            length_curr = segments[pos];
            if (segbools[pos]) dest2->coverage += length_curr;
            pos += length_curr;
        }
        
        // creates nodes and connects
        
        top = date;
        length = top - bottom;
        update_lengths(totalLength, segmentLengths);

        dest1->bottom = date;
        dest1->top = date;
        dest2->bottom = date;
        dest2->top = date;
        dest1->son1 = this;
        dest2->son1 = this;

        // VERY IMPORTANT
        dest1->fill();
        dest2->fill();

    }


    void Edge::set_actualNumberOfSites(unsigned int actualNumberOfSites) {
        numberOfMutationsPerActualSite = (unsigned int*)
                    malloc(actualNumberOfSites*sizeof(unsigned int));
        if (!numberOfMutationsPerActualSite) {
            throw EggMemoryError();
        }
        for (unsigned int i=0; i<actualNumberOfSites; i++) {
            numberOfMutationsPerActualSite[i] = 0;
        }
    }














    EdgePool::EdgePool() {
        used = 0;
        released = 0;
        ready = 0;
        cache = NULL;
    }


    void EdgePool::set(unsigned int numberOfSegments, unsigned int numberOfPreAllocated) {

        used = 0;
        released = 0;
        ready = numberOfPreAllocated;

        cache = (Edge**) malloc(numberOfPreAllocated * sizeof(Edge*));
        if (numberOfPreAllocated && !cache) throw EggMemoryError();
        for (unsigned int i=0; i<numberOfPreAllocated; i++) {
            try {
                cache[i] = new Edge(numberOfSegments);
            }
            catch (std::bad_alloc) {
                clear();
                throw EggMemoryError();
            }
        }

        this->numberOfSegments = numberOfSegments;
    }


    EdgePool::~EdgePool() {
        clear();
    }


    void EdgePool::clear() {
        if (used+released+ready) {
            for (unsigned int i=0; i<(used+released+ready); i++) {
                if (cache[i]) delete cache[i];
            }
            free(cache);
        }
        used=0;
        released=0;
        ready=0;
    }
            
            
    Edge* EdgePool::deliver() {

        if (released>0) {
            cache[used]->reset();
            released--;
        }
        else {
            if (ready==0) {
                cache = (Edge**) realloc(cache, (used+1)*sizeof(Edge*));
                if (!cache) throw EggMemoryError();
                try {
                    cache[used] = new Edge(numberOfSegments);
                }
                catch (std::bad_alloc) {
                    throw EggMemoryError();
                }
            }
            else {
                ready--;
            }
        }
        used++;

        return cache[used-1];
    }
            
            
    void EdgePool::releaseLast() {

        if (used==0) return;
        used--;
        released++;

    }

  
            
    void EdgePool::releaseAll() {
        released += used;
        used = 0;
    }

}
