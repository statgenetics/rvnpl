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


#ifndef EGGLIB_ARG_HPP
#define EGGLIB_ARG_HPP


#include "Current.hpp"
#include "Edge.hpp"
#include <string>


/** \defgroup coalesce coalesce
 *
 * \brief Coalescent simulator
 *
 * The set of classes implements a three-scale coalescent simulator with
 * recombination, and a flexible mutation model. The main classes are
 * Controller (the starting point for generating genealogies), ParamSet
 * (that centralizes parameter specification), the Change hierarchy
 * (that implements demographic change specifications), Arg (ancestral
 * recombination graph; the result of generation a genealogy) and
 * Mutator (that generates genotype data from an ARG).
 * 
 */


namespace egglib {
    
    class Random;

   /** \brief Ancestral recombination graph
    *
    * \ingroup coalesce
    * 
    * This class stores the ARG (genealogical information). It is
    * progressively built by appropriate (especially regarding to the
    * timing) calls to coal() and recomb() methods. Then it can be
    * used by a mutator class to generates data, or it can also
    * generate newick trees (one tree by non-recombining segment).
    *
    */
    class Arg {

        public:
        
           /** \brief Default constructor
            *
            * Creates a null, useless, object.
            *
            */
            Arg();

            
           /** \brief Object initialization
            * 
            * \param current address of the Current instance used by
            * the simulator.
            * 
            * \param numberOfSegments number of recombining segments.
            * 
            */
            void set(Current* current, unsigned int numberOfSegments);

            
           /** \brief Object reset method
            * 
            * This method doesn't reset all parameters (the number of
            * segments and associated tables are retained, as well as
            * the Edge object pool).
            * 
            * \param current address of the Current instance used by
            * the simulator.
            * 
            */
            void reset(Current* current);


           /** \brief Standard constructor
            * 
            * \param current address of the Current instance used by
            * the simulator.
            * 
            * \param numberOfSegments number of recombining segments
            *
            */
            Arg(Current* current, unsigned int numberOfSegments);

            
           /** \brief Destructor
            * 
            * Clears all Edge instances referenced in the object.
            * 
            */
            virtual ~Arg();
            

           /** \brief Gets the current value of the time counter
            * 
            */ 
            double time() const;
            

           /** \brief Increments the time counter
            * 
            */
            void addTime(double increment);
            

           /** \brief Performs a coalescence event
            * 
            * For this version, the two lineages to coalesce are
            * predefined.
            * 
            * \param incr increment of the time counter.
            * \param pop index of the population.
            * \param index1 first lineage to coalesce.
            * \param index2 second lineage to coalesce.
            * 
            */ 
            void coalescence(double incr, unsigned int pop, unsigned int index1, unsigned int index2);


           /** \brief Performs a coalescence event
            * 
            * For this version, the two lineages to coalesce are
            * randomly picked in the given population
            * 
            * \param incr increment of the time counter.
            * \param pop index of the population.
            * \param random pointer to simulator's random generator
            * instance.
            * 
            */ 
            void coalescence(double incr, unsigned int pop, Random* random);


           /** \brief Performs a recombination event
            * 
            * \param incr increment of the time counter.
            * \param random pointer to simulator's random generator
            * instance.
            * 
            */
            void recombination(double incr, Random* random);
            

           /** \brief Places a mutation
            * 
            * \param segment index of the segment affected.
            * 
            * \param treePosition a random number placed on the
            * interval defined by the tree length at this position.
            * 
            * \return the concerned Edge's address.
            * 
            * \todo why this is not encapsulated?
            * 
            * Another nerve-taking point: calling this method assume 
            * that all Edge of have previously undergone a call of
            * branchLength(position) with intervalPosition - what 
            * should be done by the called (that is, Mutator) through
            * my (Arg's) treeLength of something. BEWARE WHEN MODIFYING
            * (enhancements should be directed to Edge in my view) 
            * 
            */
            Edge* mute(unsigned int segment, double treePosition);

            
           /** \brief Age of the uMRCA
            * 
            * The uMRCA is the ultimate Most Recent Common Ancestor,
            * that is the point where the last segment finds its most
            * recent common ancestor. This member will have a meaningful
            * value only if the coalescent process is completed.
            * 
            */
            inline double ageUltimateMRCA() const {
                return _time;
            }
            

           /** \brief Age of the MRCA for a given segment
            * 
            * The MRCA is the Most Recent Common Ancestor, that is the
            * point where the coalescent process is over (all lineages
            * have coalesced). This member will have a meaningful
            * value only if the coalescent process is completed.
            * 
            * Note that the value is cached; it is computed only one
            * upon first call and no again, even if the Arg is modified<
            * 
            */
            inline double ageMRCA(unsigned int segmentIndex) {
                return _MRCA[segmentIndex]->bottom;
            }

           /** \brief  MRCA for each segment
            * 
            * The MRCA is the Most Recent Common Ancestor, that is the
            * point where the coalescent process is over (all lineages
            * have coalesced). This member will have a meaningful
            * value only if the coalescent process is completed.
            * 
            * Note that the value is cached; it is computed only one
            * upon first call and no again, even if the Arg is modified
            * 
            */
            inline const Edge* MRCA(unsigned int segmentIndex) {
                return _MRCA[segmentIndex];
            }

            /// Ultimate MRCA
            
            inline const Edge* uMRCA() {
                return edges[numberOfEdges-1];
            }
            
            
            /// the number of recombining segments
            unsigned int numberOfSegments;

           /** \brief Formats the newick-formatted tree for a segment
            * 
            */
            std::string newick(unsigned int segment);
            

            /// Number of initial lineages
            unsigned int numberOfSamples;


           /** \brief Total tree length (summed over all segments)
            * 
            */
            double totalLength;

           /** \brief Segment-specific tree length
            * 
            */
            double* segmentLengths;

            /// Current number of Edges in the tree (including the MRCA node)
            unsigned int numberOfEdges;

            /// Total number of recombination events that occurred
            unsigned int numberOfRecombinationEvents;
           
            /// Set the number of actual sites in all branches
            void set_actualNumberOfSites(unsigned int actualNumberOfSites);
           
            
        private:
        
            /// Copy constructor not available
            Arg(const Arg&) { }
            
            /// Assignment operator not available
            Arg& operator=(const Arg&) { return *this; }

            void init_stable_parameters();
            void init_variable_parameters();
            void clear();
            void addEdge(Edge*);
            std::string rnewick(Edge* edge, unsigned int segment, double cache);

            Current* current;
            double _time;
            Edge** edges;
            
            void findMRCA(unsigned int segmentIndex);
            void computeTotalLength();
            void computeSegmentLength(unsigned int segmentIndex);

            unsigned int* numberOfEdgesPerSegment;
            Edge** _MRCA;
            
            EdgePool edgePool;
    };

}

#endif
