/*
    Copyright 2008-2009 Stéphane De Mita, Mathieu Siol

    This file is part of EggLib.

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

#ifndef EGGLIB_STADEN_HPP
#define EGGLIB_STADEN_HPP

#include <string>
#include <istream>
#include "Align.hpp"

namespace egglib {

    /** \brief Parser of Staden output format
     *
     * \ingroup core
     *
     * The parser is available as a static method. It takes either a
     * stream or a string containing data formatted by the program GAP4
     * of the Staden package (command 'dump contig to file').
     * 
     */
    class Staden {
        
        public:


           /** \brief Parses a string
            * 
            * \param string a string containing an alignment formatted
            * by the program GAP4 of the Staden package.
            * 
            * \param deleteConsensus if true, the sequence named
            * "CONSENSUS" is deleted from the file (if it is present).
            * 
            * \return An Align instance containing the data found in
            * the Staden while, after recoding the character following
            * the standard codes.
            *
            * This method opens a stream to the string and calls the
            * overloaded method.
            * 
            * The character replacement rules assume Staden default
            * convention, as follows:
            *    - "-" codes for an unknown base and is replaced by "N".
            *    - "*" codes for an alignment gap and is replaced by "-".
            *    - A white space represents missing data and is replaced
            * by "?".
            * 
            */
            static Align parse(const std::string& string, bool deleteConsensus=true);
            

           /** \brief Parses an open stream
            * 
            * \param stream the open containing an alignment formatted
            * by the program GAP4 of the Staden package.
            * 
            * \param deleteConsensus if true, the sequence named
            * "CONSENSUS" is deleted from the file (if it is present).
            * 
            * \return An Align instance containing the data found in
            * the Staden while, after recoding the character following
            * the standard codes.
            *
            * The character replacement rules assume Staden default
            * convention, as follows:
            *    - "-" codes for an unknown base and is replaced by "N".
            *    - "*" codes for an alignment gap and is replaced by "-".
            *    - A white space represents missing data and is replaced
            * by "?".
            * 
            */
            static Align parse(std::istream& stream, bool deleteConsensus=true);


        private:
        
            /// Not allowed to instantiate this class
            Staden() { }
            
            /// Not allowed to instantiate this class
            Staden(const Staden& source) { }
            
            /// Not allowed to instantiate this class
            ~Staden() { }


           /* Gets the start position of sequences
            *
            * The functions gives total number of characters before the start of sequences
            * and reads through until the next backspace (ignores the first line).
            */
            static void getShift();
 
            // Translates according to the Staden format
            static char transforme(char);
            
            // Imports one sequence
            static bool readOneSequence();
            
            // Imports and concatenates one sequence
            static bool readAppendOneSequence();
            
            // Replaces dots by the matching character from CONSENSUS
            static void undot(bool delete_consensus=true);

            // The number of characters before the start of sequences
            static int shift;
            
            // The dynamically filled container (will result in an aligment)
            static Container container;
            
            // The current position
            static int currpos;
            
            // The reading stream
            static std::istream* stream;
            
            // Stores unique 8 characters discriminating readings
            static std::vector<std::string> ID;
    };
}
    
#endif
