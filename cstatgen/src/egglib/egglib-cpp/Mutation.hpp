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

#ifndef EGGLIB_MUTATION_HPP
#define EGGLIB_MUTATION_HPP

#include <vector>
#include "Edge.hpp"

namespace egglib {

   /** \brief Very simple container of some information relative to a mutation
    * 
    * \ingroup coalesce
    *
    */
    class Mutation {

        public:
    
            /// Default constructor
            Mutation();
            
            /// Age
            //double age;
            
            /// Mutation index (for finding in Edge)
            unsigned int actualSiteIndex;
            
            /// Position
            double position;
            
            /// Segment index
            unsigned int segmentIndex;
            
            /// Pointer to edge
            //const Edge* edge;
            
        private:
        
            void init();

    };

}

#endif
