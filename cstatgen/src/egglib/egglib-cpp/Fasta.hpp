/*
    Copyright 2008-2009 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_FASTA_HPP
#define EGGLIB_FASTA_HPP

#include <istream>
#include <iostream>
#include <string>
#include "Container.hpp"

namespace egglib {

   /** \brief Fasta parser/formatter
    *
    * \ingroup core
    *
    * Reads a multifasta sequence file from a string, a stream or a file
    * and returns a Container. See the description of the format below.
    * Formats a fasta string from a sequence container object and places
    * it in a string, a stream of a file. All methods are static and the
    * class cannot be instantiated. The methods parsef and formatf will
    * open the file for you while the others will read/write directly
    * in a string.
    * 
    * Specifications of the fasta format:
    * 
    *    - The number of sequences is not limited.
    * 
    *    - Each sequence is preceded by a header limited to a single
    *      line and starting by a ">" character.
    * 
    *    - The header length is not limited and all characters are
    *      allowed but white spaces and special characters are
    *      discouraged.
    * 
    *    - Group indices are specified by \@0, \@1, \@2...  strings
    *      appearing at the end of the header string (just before the
    *      carriage return). Note that group labels are ignored by
    *      default.
    * 
    *    - Group indices are ignored unless specifically specified in a
    *      parser's options.
    * 
    *    - The sequence itself continues on following lines until the
    *      next ">" character or the end of the file.
    * 
    *    - White spaces, tab and carriage returns are allowed at any
    *      position There is no limitation in length and different
    *      sequences can have different lengths.
    * 
    *    - Although the standard is lower case characters, Fasta
    *      assumes upper case characters and only supports lower case
    *      characters (and converts them to upper case characters).
    *      Information coded by change in case is lost.
    *
    */
    class Fasta {

      public:

       /** \brief Imports a fasta file
         *
         * Imports the content of the file as is. Calls the method
         * pase(std::istream*, bool) by creating its own istream.
         *
         * \param fname the name of a fasta file.
         * 
         * \param importGroupLabels if set to true, scan automatically
         * for groups. The format is @ followed by an integer, placed
         * at the end of the header string(sequences without labels
         * will be treated as \@0).
         * 
         * \return A Container object containing the sequences.
         * 
         */
        static Container parsef(const char* fname, bool importGroupLabels=false);


       /** \brief Imports a fasta file
         *
         * Imports the content of the file as is. Calls the method
         * pase(std::istream*, bool) by creating its own istream. This
         * method expects a reference to a Container to which the
         * sequences will be appended.
         *
         * \param fname the name of a fasta file.
         * 
         * \param container a Container instance, empty or not.
         * 
         * \param importGroupLabels if set to true, scan automatically
         * for groups. The format is @ followed by an integer, placed
         * at the end of the header string(sequences without labels
         * will be treated as \@0).
         * 
         * \return Nothings: the new sequences are appended to the
         * Container passed as argument.
         * 
         */
        static void parsef(const char* fname, Container& container, bool importGroupLabels=false);


       /** \brief Imports a fasta file
         *
         * Imports the content of the file as is. Calls the method
         * pase(std::istream*, bool) by creating its own istream.
         *
         * \param str a string containing the data.
         * 
         * \param importGroupLabels if set to true, scan automatically
         * for groups. The format is @ followed by an integer, placed
         * at the end of the header string(sequences without labels
         * will be treated as \@0).
         * 
         * \return A Container object containing the sequences.
         * 
         */
        static Container parse(const std::string& str, bool importGroupLabels=false);


       /** \brief Imports a fasta file
         *
         * Imports the content of the file as is. Calls the method
         * pase(std::istream*, bool) by creating its own istream. This
         * method expects a reference to a Container to which the
         * sequences will be appended.
         *
         * \param str a string containing the data.
         * 
         * \param container a Container instance, empty or not.
         * 
         * \param importGroupLabels if set to true, scan automatically
         * for groups. The format is @ followed by an integer, placed
         * at the end of the header string(sequences without labels
         * will be treated as \@0).
         * 
         * \return Nothing: new sequences are appended to the Container
         * passed as argument.
         * 
         */
        static void parse(const std::string& str, Container& container, bool importGroupLabels=false);


       /** \brief Imports a fasta file from an open stream
         *
         * Imports the content of the file as is.
         *
         * \param stream an open stream (file or string) containing the
         * data.
         * 
         * \param importGroupLabels if set to true, scan automatically
         * for groups. The format is @ followed by an integer, placed
         * at the end of the header string(sequences without labels
         * will be treated as \@0).
         * 
         * \return A Container object containing the sequences.
         * 
         */
        static Container parse(std::istream& stream, bool importGroupLabels=false);


       /** \brief Imports a fasta file from an open stream
         *
         * Imports the content of the file as is. This
         * method expects a reference to a Container to which the
         * sequences will be appended.
         *
         * \param stream an open stream (file or string) containing the
         * data.
         * 
         * \param container a Container instance, empty or not.
         * 
         * \param importGroupLabels if set to true, scan automatically
         * for groups. The format is @ followed by an integer, placed
         * at the end of the header string(sequences without labels
         * will be treated as \@0).
         * 
         * \return Nothing: the new sequences are appended to the
         * Container passed as argument.
         * 
         */
        static void parse(std::istream& stream, Container& container, bool importGroupLabels=false);
        
        
       /** \brief Export sequences as fasta
        *
        * \param fname the name of the file where to place the result.
        * 
        * \param container Container object to export.
        * 
        * \param exportGroupLabels if set to true, exports group
        * indices as a \@x at the end of the sequence name, where x is
        * the group index. Otherwise, this information is discarded.
        * 
        * \param lineLength the number of characters to place on a
        * single line. If zero, no newlines are inserted within
        * sequences.
        * 
        */
        static void formatf(const char* fname, const Container& container, bool exportGroupLabels=false, unsigned int lineLength=50);


       /** \brief Export sequences as fasta
        *
        * \param file an open stream.
        * 
        * \param container Container object to export.
        * 
        * \param exportGroupLabels if set to true, exports group
        * indices as a \@x at the end of the sequence name, where x is
        * the group index. Otherwise, this information is discarded.
        * 
        * \param lineLength the number of characters to place on a
        * single line. If zero, no newlines are inserted within
        * sequences.
        * 
        */
        static void format(std::ostream& file, const Container& container, bool exportGroupLabels=false, unsigned int lineLength=50);


       /** \brief Export sequences as fasta
        * 
        * This medod creates internally an ostringstream, calls the
        * method format(ostream, container, bool) and returns the
        * resulting string.
        *
        * \param container Container object to export.
        * 
        * \param exportGroupLabels if set to true, exports group
        * indices as a \@x at the end of the sequence name, where x is
        * the group index. Otherwise, this information is discarded.
        * 
        * \param lineLength the number of characters to place on a
        * single line. If zero, no newlines are inserted within
        * sequences.
        * 
        * \return The formatted string.
        * 
        */
        static std::string format(const Container& container, bool exportGroupLabels=false, unsigned int lineLength=50);

          
          
      protected:
        
        /// This class cannot be instantiated
        Fasta() { }
        
        /// This class cannot be instantiated
        Fasta(const Fasta& source) { }
        
        /// This class cannot be or copied
        Fasta& operator=(const Fasta& source) { return *this; }
        
        /// This class cannot be instantiated
        virtual ~Fasta() { }

        
    };
}

#endif
