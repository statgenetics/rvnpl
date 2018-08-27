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

#ifndef EGGLIB_EDGE_HPP
#define EGGLIB_EDGE_HPP

#include <vector>
#include <climits>
#include "EggException.hpp"

namespace egglib {

    class Random;

    /** \brief Edge of the ancestral recombination graph
    *
    * \ingroup coalesce
    *
    * Each Edge instance provides access to its 0, 1 or 2 descendants
    * (the former holds for a terminal node, the middle for the parent
    * of a recombined node and the latter for the parent of a coalesced
    * node (most classical node in the coalescent).The Edge also
    * provides to the edge length. Note that the Edge instance must be
    * understood as an ARG node and the branch above it (latter in the
    * coalescence process). Edge instances also keep track of the list
    * of descendants descending from this node (which may differ along
    * recombining segment). Edge instances *must* be created through one
    * of the "default" and "coalescence" constructors or through the 
    * recombination method. Edge instances should never be copied but
    * manipulated by references.
    * 
    */
    class Edge {
    
        public:
        
            /// Destructor
            virtual ~Edge();
        
           /** \brief Constructor
            * 
            * \param numberOfSegments the number of recombining segments
            * (one for a non-recombining region).
            * 
            * Use the Pool, instead. Objects are delivered with a
            * complete coverage.
            * 
            */ 
            Edge(unsigned int numberOfSegments);


            /// Restore object to `factory` state
            void reset();


           /** \brief Builds for internal node
            * 
            * \param date the date of creation of the edge.
            * \param son1 first edge descending from this edge.
            * \param son2 second edge descending from this edge.
            * \param edgesPerSegments counts the current number of
            * (non-coalesced lineages for each lineages); must have the
            * appropriate size and will be updated.
            * \param MRCA the list where to place the address of segment
            * MRCA, if it occurs.
            * \param totalLength the total length of the tree.
            * \param segmentLengths the table of tree lengths per
            * segment.
            *
            * Assumes the current object has the correct number of
            * segments.
            * 
            */
            void coalescence(double date, Edge* son1, Edge* son2,
                        unsigned int* edgesPerSegments, Edge** MRCA,
                        double& totalLength, double* segmentLengths);


           /** \brief Generates a recombination event
            * 
            * \param date the date of the event.
            * \param dest1 destination for the first resulting edge.
            * \param dest2 destination for the second resulting edge.
            * \param random pointer to the Random instance used by the 
            * simulator.
            * \param totalLength the total length of the tree.
            * \param segmentLengths the table of tree lengths per
            * segment.
            * 
            * dest1 and dest2 must be Edge address initialized with the
            * appropriate number of segments.
            * 
            */
            void recombination(double date,
                        Edge* dest1, Edge* dest2, Random* random,
                        double& totalLength, double* segmentLengths);


           /** \brief Define this edge to be terminal
            * 
            * The edge have only non-covered segments
            * 
            */
            void set_terminal(unsigned int leaf_index);
            
            /// Branch's raw length (doesn't take account segment coverage)
            double length;

            /// Number of covered segments
            unsigned int coverage;
            
            /// Time position of the branch's bottom
            double bottom;
            
            /// Time position of the branch's top
            double top;
            
            /// Address of the  first son
            Edge* son1;

            /// Address of the second son
            Edge* son2;

            /// Number of sons (0, 1 or 2)
            unsigned int numberOfSons;

            /// Checks if a given segment is covered
            inline bool segment(unsigned int segmentIndex) const {
                //while(segments[segmentIndex]==0) segmentIndex--;
                if (segments[segmentIndex]==0)  return false;
                if (segments[segmentIndex]==UINT_MAX) return true;
                return segbools[segmentIndex];
            }

            /// leaf index (0 for internal nodes)
            inline unsigned int label() const {
                return _label;
            }

            /// Number of mutations per segment
            unsigned int* numberOfMutationsPerActualSite;

            /// Sets the actual number of sites
            void set_actualNumberOfSites(unsigned int actualNumberOfSites);


        private:

            unsigned int _label;

            /// Default constructor is not available
            Edge();

            /// Copy constructor is not available
            Edge(const Edge& edge) {}

            /// Assignment operator is not available
            Edge& operator=(const Edge& edge) { return *this; }

            void init(unsigned int numberOfSegments);

            unsigned int numberOfSegments;
            unsigned int* segments;
            bool* segbools;

            /// complete the covered ranges with UINT_MAX's and the non-covered by 0's
            inline void fill() {
                coverage=0.;
                unsigned int i=0,j;
                while (i<numberOfSegments) {
                    if (segbools[i]==true) {
                        coverage += segments[i];
                        for (j=1; j<segments[i]; j++) {
                            segments[i+j] = UINT_MAX;
                        }
                    }
                    i+=segments[i];
                }
            }



            /// update containing Arg's branch lengths
            inline void update_lengths(double& totalLength, double* segmentLengths) {
                unsigned int i=0,j;
                while (i<numberOfSegments) {
                    if (segbools[i]==true) {
                        totalLength += segments[i]*length;
                        for (j=0; j<segments[i]; j++) {
                            segmentLengths[i+j] += length;
                        }
                    }
                    i+=segments[i];
                }
            }


    };












    /** \brief Pool of Edge objects
    *
    * \ingroup coalesce
    *
    * Holds a pool of Edge objects that can be recycled to spare the
    * building burden. A construction time, a number of Edge objects
    * equals to the predicted number of needed instances should be
    * requested. The Edge's will be prebuilt immediately and delivered
    * upon request. After use, the Edge's should be released. It is only
    * possible to release the last issued Edge instance or all of them
    * at once.
    * 
    */
    class EdgePool {

        public:
        
            /// Default constructor (nothing allocated)
            EdgePool();


            /// Destructor
            virtual ~EdgePool();


           /** \brief Configure pool
            * 
            * Pre-allocates a given number of Edge objects. The objects
            * will be immediately available.
            * 
            * Data previously allocated (by a previous call of this
            * function or by the deliver() method) will be lost so it
            * can be required to use clear() before.
            * 
            * \param numberOfSegments the number of segments of the
            * simulation; all Edge instances will use this value.
            * 
            * \param numberOfPreAllocated the number of Edge that should
            * be kept ready for immediate use.
            * 
            */
            void set(unsigned int numberOfSegments, unsigned numberOfPreAllocated);
            
            
           /** \brief Frees internally stored memory
            * 
            * This invalidate all points that have been delivered
            * previously. However, any previously set number of segments
            * (0, by default) is retained.
            * 
            */
            void clear();

            
           /** \brief Deliver an Edge
            * 
            * The object must not be freed by the client! This object is
            * allocated on the heap if the cache is not large enough,
            * only reset if it was previously released, or just delivered
            * if it is one of the initially allocated instances.
            * 
            */
            Edge* deliver();
            
            
           /** \brief Release an Edge
            * 
            * Release the last delivered Edge. The instance is only
            * cached for a potential future use; it is not freed nor
            * reset immediately.  If no Edge's are in use, nothing is
            * done.
            * 
            */
            void releaseLast();

            
           /** \brief Release all Edge's
            * 
            * Release all delivered Edges. The instances are only
            * cached for a potential future use; they are not freed nor
            * reset immediately. If no Edge's are in use, nothing is
            * done.
            * 
            */
            void releaseAll();

        private:
        
            /// Not available
            EdgePool(const EdgePool& ep) {}
            
            /// Not available
            EdgePool& operator=(const EdgePool& ep) { return *this; }

            unsigned int numberOfSegments;
            unsigned int used;
            unsigned int released;
            unsigned int ready;
            Edge** cache;
            
    };

}

#endif
