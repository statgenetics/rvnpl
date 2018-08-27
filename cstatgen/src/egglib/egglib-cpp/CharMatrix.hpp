/*
    Copyright 2009 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_CHARMATRIX_HPP
#define EGGLIB_CHARMATRIX_HPP


namespace egglib {

    /** \brief Interface for classes usable as a square matrix of characters
    *
    * \ingroup core
    *
    */
    class CharMatrix {
    
    public:
    
        /** \brief Gets number of rows or sequences
         * 
         */
        virtual unsigned int numberOfSequences() const = 0;


        /** \brief Gets number of columns or sites
         * 
         */
        virtual unsigned int numberOfSites() const = 0;
        
        
        /** \brief Gets the character at a given position
         * 
         * The accessor should be "fast" and does not guarantee to
         * perform out-of-bounds checks
         * 
         */
        virtual char character(unsigned int sequence, unsigned int site) const = 0;
        
        
       /** \brief Gets population index
        * 
        */
        virtual unsigned int populationLabel(unsigned int row) const = 0;
        
        
       /** \brief Get site position
        * 
        */
        virtual double sitePosition(unsigned int column) const = 0;

    };
}

#endif
