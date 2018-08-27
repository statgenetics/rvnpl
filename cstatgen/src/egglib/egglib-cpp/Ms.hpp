/*
    Copyright 2008,2009,2011 Stéphane De Mita and Mathieu Siol

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

#ifndef EGGLIB_GMS_HPP
#define EGGLIB_GMS_HPP

#include "DataMatrix.hpp"
#include <string>
#include <istream>

namespace egglib {

    /** \brief ms-like sequence format parser
     * 
     * The class provides parsing (input) and formatting (output)
     * operations in ms format, that is the format used by Richard
     * Hudson's program ms for outputting genotypes and by the
     * associated program samplestat for reading them. Both types of
     * operations are available through static methods using either
     * a string or a stream (which can be a stream to or from a file
     * or a string). In either case, types from the STL are used.
     * Although ms deals only with data coded with 0 and 1, the class Ms
     * offers the possibility of both importing and exporting data coded
     * with by integer. All methods have an option named "separated". If
     * this option is true, the parser or formatter introduces a slight
     * modification of the format: genotypes individual data are
     * separated by a white space ("1 0 1 1" instead of "1011", allowing
     * genotype values larger than 9: "1 0 11 1").
     *
     * \ingroup core
     *
     */
     class Ms {

      public:
         
       /** \brief Imports a sequence alignment
        * 
        * Creates a istringstream from the string and calls the
        * overloaded method.
        * 
        * \param str the string to parse.
        * \param ns the expected number of sequences.
        * \param separated true if a white space separator is placed
        * between genotype at each site.
        *
        * \return A sequence alignment as a data matrix.
        */
        static DataMatrix get(std::string, unsigned int ns, bool separated=false);


       /** \brief Imports a sequence alignment
        * 
        * Attemps to generate a DataMatrix object from the stream.
        * Reads only one simulation and throws a SeqlibFormatError
        * exception in case of format error.
        * 
        * Allows any number of white lines before the //, but no other
        * data. Supports \r at the end of lines (before the \n).
        * Accepted symbols are all integers (0-9).
        *
        * \param stream the stream to parse.
        * \param ns the expected number of sequences.
        * \param separated true if a white space separator is placed
        * between genotype at each site.
        * 
        * \return A sequence alignment as a data matrix.
        */
        static DataMatrix get(std::istream& stream, unsigned int ns, bool separated=false);


       /** \brief Exports a sequence alignment
        * 
        * Internally creates a stringstream, calls the overloaded method
        * and returns the outcome.
        *
        * \param dataMatrix the alignment object to write.
        * \param separated true if a white space separator must be placed
        * between the genotype at each site.
        * 
        */
        static std::string format(DataMatrix& dataMatrix, bool separated=false);
        
        
       /** \brief Exports a sequence alignment
        * 
        * Writes the formatted string to the stream 'on the fly'. The
        * formatted string is guaranteed to starts with a // line and
        * ends with an empty line. The client is expected to take care
        * of writing any header and add an additional white line between
        * simulations if needed. The method throws a SeqlibRuntimeError
        * if the stream is not writable. The data matrix should contain
        * only data within range 0-9 if separated is false (default) and
        * any positive (>=0) integer if separated is true. Note that
        * output generated with separated=true is never compatible with
        * the original ms format, and that output generated with
        * separator=false is compatible with the original ms format only
        * if all alleles are 0 or 1 (which is not checked by this
        * formatted).
        * 
        * \param stream the stream (file or string stream) where to
        * write the output.
        * \param dataMatrix the alignment object to write.
        * \param separated true if a white space separator must be placed
        * between the genotype at each site.
        * 
        */
        static void format(std::ostream& stream, DataMatrix& dataMatrix, bool separated=false);


       /** \brief Returns the last tMRCA read by any Ms instance
        * 
        * If a tMRCA value was present in the last simulation read by
        * any Ms instance, it will be returned by this method. A value
        * of -1. is returned if no simulation was read, or if the last
        * simulation didn't contain a tMRCA value or if the last
        * simulation provoked an exception before reaching the tMRCA
        * line.
        * 
        */
        static double tMRCA();


       /** \brief Returns the last "prob" read by any Ms instance
        * 
        * "prob" is returned by ms when a fixed number of segregating
        * sites is used in conjunction with a theta value. If a "prob"
        * value was present in the last simulation read by any Ms
        * instance, it will be returned by this method. A value of -1
        * is returned if no simulation was read, or if the last
        * simulation didn't contain a "prob" value or if the last
        * simulation provoked an exception before reaching the "prob"
        * line.
        * 
        */
        static double prob();
    

       /** \brief Returns the tree string found in the last simulation read by any Ms instance
        * 
        * If one or more trees were present in the last simulation read
        * by any Ms instance, they will be returned as a unique string
        * by this method. An empty string is returned if no simulation
        * was read, or if the last simulation, or if the last simulation
        * didn't contain any tree value or if the last simulation
        * provoked an exception before reaching the tree line.
        * 
        * Note: the trees are returned as a single line.
        * 
        */
        static std::string trees();

         
      private:
        // Line parser (the last \n is extracted and discarded - no error upon EOF)
        std::string next_line(std::istream& stream);
        
        /// tMRCA (-1 if not found in ms output)
        static double _tMRCA;
        
        /// probability (-1 if not found in ms output)
        static double _prob;
        
        /// tree string (maybe contain several trees) (empty string if not found in ms output)
        static std::string _trees;

        
        /// No instantiation allowed
        Ms() { }
        
        /// A fortiori no destruction allowed
        ~Ms() { }

        /// No copy allowed
        Ms(const Ms&) { }

        /// No copy allowed
        Ms& operator=(const Ms&) { return *this; }
                
    };
}
    
#endif
