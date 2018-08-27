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


#ifndef EGGLIB_ALIGN_HPP
#define EGGLIB_ALIGN_HPP

#include "Container.hpp"
#include "CharMatrix.hpp"
#include <vector>

/** \mainpage Summary
 * 
 * This is the automatically-generated reference manual of the C++
 * egglib-cpp library. The library is presented as several modules, but
 * note that they are only used to structure the documentation.
 * 
 * There is a single namespace (egglib) in which all classes are
 * defined. See an example of programming with egglib-cpp in the
 * EggLib package main documentation. Use "Modules" or "Classes" above
 * to navigate in the library reference manual.
 * 
 */


/** \defgroup core core
 *
 * \brief Central core of the C++ library of Egglib
 *
 * Data storage classes, parsers/formatters and tools, plus exception
 * types.
 * 
 */

namespace egglib {


   /** \brief Handles a sequence alignment
    *
    * \ingroup core
    * 
    * Creation from a file or string stream should be performed using
    * the class Fasta. Align objects can be created by deep copy from
    * both Align and Container type. In the latter case, the length are
    * artificially equalized by "?" characters. Align objects can be
    * created from a DataMatrix object (and all the way arround) using
    * the specific class DMAConverter.
    *
    * Sequences are represented by two strings (name and sequence) and
    * an integer (group) that can be accessed or modified by index.The
    * order of sequences is guaranteed to be conserved, as if Align was
    * a list of triplets (name, sequence, group).
    *
    * The data matrix is implemented as continuous array (char**) and
    * allows efficient access and modification of data. For very large
    * data matrices you might claim immediately the required memory
    * using the constructor Align(unsigned int, char**).
    * 
    */
    class Align : public Container, public CharMatrix {
        public:
        
           /** \brief Creates an empty alignment
            * 
            */
            Align();


           /** \brief Creates an alignment from a data matrix.
            * 
            * Allows you to create an object from data stored in a char*
            * array. The array's dimensions must be passed to the
            * constructor, and as a result there is not need to
            * terminate each sequence by a NULL character.
            * 
            * \param number_of_sequences the number of sequences (the
            * length of the first dimension of the array).
            * 
            * \param alignment_length the length of sequences (the
            * length of all lines of the array).
            * 
            * \param cstring_array the pointer to the data matrix.
            * 
            */
            Align(unsigned int number_of_sequences, unsigned int alignment_length, char const * const * const cstring_array);


            /** \brief Creates an alignment with given dimensions
             * 
             * Allows you to allocate directly a data matrix of a given
             * size. Names are empty strings, groups 0, and all
             * characters are ?.
             * 
            * \param number_of_sequences the number of sequences (the
            * length of the first dimension of the array).
            * 
            * \param alignment_length the length of sequences (the
            * length of all lines of the array).
            * 
            */
            Align(unsigned int number_of_sequences, unsigned int alignment_length);


           /** \brief Copy constructor
            * 
            */
            Align(const Align& align);


           /** \brief Copy constructor accepting a Container object
            * 
            * All but the longest sequences are padded with ? to match
            * the longest sequence's length.
            * 
            */
            Align(const Container& container);


           /** \brief Copy operator
            * 
            */
            Align& operator=(const Align& align);


           /** \brief Copy operator accepting a Container object
            * 
            * All but the longest sequences are padded with ? to match
            * the longest sequence's length.
            * 
            */
            Align& operator=(const Container& container);


           /** \brief Destructor
            * 
            */
            virtual ~Align();


           /** \brief Adds a sequence
            *
            * If the object already contains at least one sequence, the
            * new sequence must have the same length. Otherwise, a
            * EggUnalignedError is raised.
            *
            * \param name the name of the sequence.
            * \param sequence the sequence string.
            * \param group the group index of the sequence.
            * \return The new number of sequences.
            * 
            */
            virtual unsigned int append(const char* name, const char* sequence, unsigned int group=0);


           /** \brief Removes a position (column) of the alignment
            *
            * \param pos the position to remove in the alignment.
            * \return The new length of the alignment.
            *
            */
            virtual unsigned int removePosition(unsigned int pos);


           /** \brief Removes a sequence from the alignment
            *
            * \param pos the index of the sequence to remove.
            * \return The new number of sequences.
            * 
            */
            virtual unsigned int remove(unsigned int pos);


           /** \brief Replace a sequence string
            * 
            * The new sequence must have the same length than the
            * alignment. Otherwise, a EggUnalignedError is raised.
            * 
            * \param seq the index of the sequence to change.
            * \param sequence the new sequence.
            * 
            */
            virtual void sequence(unsigned int seq, const char* sequence);


           /** \brief Gets the name of a given sequence
            * 
            * \param pos the index of the sequence.
            * 
            * \return The sequence string for that particular sequence.
            * 
            */
            virtual inline const char* sequence(unsigned int pos) const { return Container::sequence(pos); }
            
            
           /** \brief Alignment length
            * 
            * Returns 0 if the alignment is empty.
            * 
            */
            virtual  unsigned int ls() const;


           /** \brief Length of a given sequence
            * 
            * Calling this function is exactly the same as calling ls()
            * (without arguments), regardless of the index provided,
            * except that an exception is thrown if the index is out of
            * bounds. Provided for compatibility with Container.
            * 
            * \param pos the index of the sequence.
            * \return the length of the alignment.
            * 
            */
            virtual unsigned int ls(unsigned int pos) const;


           /** \brief Fast and unsecure accessor
            * 
            * This accessor doesn't perform out-of-bound checking!
            * 
            * \param s the index of the sequence (line).
            * \param p the position in the alignment (column).
            * \return The character at the given position.
            * 
            */
            inline char character(unsigned int s, unsigned int p) const { return sequences[s][p]; }


           /** \brief Gets a nucleotide
            * 
            * This modifier does perform out-of-bound checking.
            * The specified position must exist.
            * 
            * \param sequence the index of the sequence (line).
            * \param position the position in the alignment (column).
            * \return the character at the given position.
            * 
            */
            virtual char get(unsigned int sequence, unsigned int position) const;


           /** \brief Sets a matrix position to a new character
            * 
            * This modifier does perform out-of-bound checking.
            * The specified position must exist.
            * 
            * \param sequence the index of the sequence (line).
            * \param position the position in the alignment (column).
            * \param ch the new character value.
            */
            virtual void set(unsigned int sequence, unsigned position, char ch);


           /** \brief Reverse a given column in binary data
            *
            * The specified column must contain only "0" ans "1" characters.
            * "0" is replaced by "1" and all the way around
            * 
            */
            void binSwitch(unsigned int pos);


           /** \brief Extracts specified positions (columns) of the alignment
            *
            * All the specified sites are extracted in the specified
            * order. This function is suitable for bootstrap (resample
            * allowing redrawing the same site) and permutations.
            * 
            * This function doesn't perform out-of-bound checking.
            * 
            * \param list_of_sites a vector containing alignment
            * positions.
            * 
            * \return A copy of the object containing the specified
            * set of positions.
            * 
            */
            Align vslice(std::vector<unsigned int> list_of_sites);


           /** \brief Extracts a range of positions (columns)
            * 
            * \param a the first position.
            * 
            * \param b the index immediately passed the last sequence to
            * extract.
            * 
            * \return A copy of the object containing the specified
            * range of sequences.
            * 
            * Positions a to b-1 are extracted, provided that the
            * indices fit in the current length of sequences. To extract
            * all sequences, use align.vslice(0, align.ls()).
            * 
            * Note: invalid ranges will be silently supported. If
            * a>=ls or b<=a, an empty object is returned. If b>ns,
            * ls will be substituted to a.
            */
            Align vslice(unsigned int a, unsigned int b);


           /** \brief Deletes all the content of the object
            * 
            */
            virtual void clear();


           /** \brief Same as ns()
            * 
            */
            inline unsigned int numberOfSequences() const {
                return _ns;
            }


           /** \brief Same as ls()
            * 
            */
            inline unsigned int numberOfSites() const {
                return _ls;
            }


           /** \brief Gets a group label (insecure)
            * 
            */
            inline unsigned int populationLabel(unsigned int sequenceIndex) const {
                return groups[sequenceIndex];
            }
            
            
           /** \brief Just return the passed value
            *
            */
            inline double sitePosition(unsigned int position) const {
                return (double) position;
            }


        protected:
        
            /// This function is not available for alignments
            virtual void appendSequence(unsigned int pos, const char* sequence) {}

            // Initializer (creates a valid empty alignment)
            virtual void init();
        
            // Makes a deep copy of the specified data matrix - if cstring_array is NULL, then ignores it and pads with ?'s
            virtual void setFromSource(unsigned int number_of_sequences, unsigned int alignment_length, const char* const * const cstring_array);

            // Copies from a Container
            virtual void copyObject(const Container&);
            
            // Copies from an Align
            virtual void copyObject(const Align&);
            
            // Alignment length
            unsigned int _ls;
    };
}

#endif
