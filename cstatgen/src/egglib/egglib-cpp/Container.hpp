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


#ifndef EGGLIB_CONTAINER_HPP
#define EGGLIB_CONTAINER_HPP


namespace egglib {

   /** \brief Handles a set of sequence alignment (aligned or not)
    *
    * \ingroup core
    * 
    * Creation from a file or string stream should be performed using
    * the class Fasta.
    *
    * Sequences are represented by two strings (name and sequence) and
    * an integer (group) that can be accessed or modified by index.The
    * order of sequences is guaranteed to be conserved, as if Container
    * was a list of triplets (name, sequence, group).
    *
    * The data matrix is implemented as continuous arrays (char**) and
    * allows efficient access and modification of data. For very large
    * data matrices you might claim immediately the required memory
    * using the constructor Container(unsigned int, char**).
    *
    */
    class Container {
    
        public:
      
          /** \brief Creates an empty object
            * 
            */
            Container();
        
        
           /** \brief Copy constructor
            * 
            */
            Container(const Container& source);
        
        
           /** \brief Assignment operator
            * 
            */
            Container& operator= (const Container& source);


           /** \brief Creates an object from a data matrix
            * 
            * Allows you to create an object from data stored in a char*
            * array. The array's size must be passed to the constructor.
            * Since sequences can have different lengths, you need to
            * terminate each sequence by a NULL character. This constructor
            * is dedicated to very performance-critical tasks. For usual
            * tasks, using the default constructor and subsequently adding
            * sequences with addSeq should be enough.
            * 
            * \param number_of_sequences the number of sequences (the length
            * of the first dimension of the array).
            * 
            * \param cstring_array the pointer to the data matrix.
            * 
            */
            Container(unsigned int number_of_sequences, char const* const* const cstring_array);


           /** \brief Destructor
            * 
            */
            virtual ~Container();

        
           /** \brief Clears all content of the object
            * 
            */
            virtual void clear();


           /** \brief Adds a sequence to the object
            *
            * \param name the name of the sequence, as a c-string.
            * \param sequence the sequence string, as a c-string.
            * \param group the group index of the sequence.
            *
            * \return The new number of sequences.
            * 
            */
            virtual unsigned int append(const char* name, const char* sequence, unsigned int group=0);
    
    
           /** \brief Removes a sequence from the object
            *
            * \param pos the index of the sequence to remove.
            * 
            * \return The new number of sequences.
            */
            virtual unsigned int remove(unsigned int pos);


           /** \brief Changes the name of a given sequence
            * 
            * \param pos the sequence index.
            * \param name the new name as a C-like string.
            * 
            */
            virtual void name(unsigned int pos, const char* name);


           /** \brief Changes the sequence string of a given sequence
            * 
            * \param pos the sequence index.
            * \param sequence the new sequence as a C-like string.
            * 
            */
            virtual void sequence(unsigned int pos, const char* sequence);


           /** \brief Appends a string to the a given sequence
            * 
            * \param pos the sequence index.
            * \param sequence the sequence to append at the end of the
            * current one.
            * 
            */
            virtual void appendSequence(unsigned int pos, const char* sequence);


           /** \brief Changes a character
            * 
            * \param sequence the sequence index.
            * \param position the character index.
            * \param ch the new character value.
            * 
            * The positions must fit in the current ranges.
            * 
            */
            virtual void set(unsigned int sequence, unsigned position, char ch);


           /** \brief Gets a given character
            * 
            * \param s the sequence index.
            * \param p the character index.
            * 
            * \return the character value.
            * 
            * The positions must fit in the current ranges.
            * 
            */
            virtual char get(unsigned int s, unsigned int p) const;


           /** \brief Changes the group index of a given sequence
            * 
            * \param pos the sequence index.
            * \param group the new group index value.
            * 
            */
            virtual void group(unsigned int pos, unsigned int group);

        
           /** \brief Extracts a range of sequences
            * 
            * \param a the index of the first sequence.
            * 
            * \param b the index immediately passed the last sequence to
            * extract.
            * 
            * \return A copy of the object containing the specified
            * range of sequences.
            * 
            * Sequences a to b-1 are extracted, provided that the
            * indices fit in the current number of sequences. To extract
            * all sequences, use container.hslice(0, container.ns()).
            * 
            * Note: invalid ranges will be silently supported. If
            * a>=ls or b<=a, an empty object is returned. If b>ns,
            * ls will be substituted to a.
            * 
            */
            Container hslice(unsigned int a, unsigned int b) const;


           /** \brief Gets the number of sequences
            * 
            */
            unsigned int ns() const;
        
        
           /** \brief Gets the length of a given sequence
            * 
            * \param pos the index of the sequence.
            * 
            * \return The length of that particular sequence.
            * 
            */
            virtual unsigned int ls(unsigned int pos) const ;
        
        
           /** \brief Gets the name of the a given sequence
            * 
            * \param pos the index of the sequence.
            * 
            * \return The name of that particular sequence.
            * 
            */
            virtual const char* name(unsigned int pos) const;

        
           /** \brief Gets the name of a given sequence
            * 
            * \param pos the index of the sequence.
            * 
            * \return The sequence string for that particular sequence.
            * 
            */
            virtual const char* sequence(unsigned int pos) const;



           /** \brief Gets the group index of a given sequence
            * 
            * \param pos the index of the sequence.
            * 
            * \return The group index of that particular sequence.
            * 
            */
            virtual unsigned int group(unsigned int pos) const;
        
        
           /** \brief Checks if all lengths are equal
            * 
            * Returns true if the length of all sequences are equal or
            * if there is less thant two sequences.
            * 
            */
            bool isEqual() const;


           /** \brief Equalizes sequence lengths
            *
            * Extends sequences as need to ensure that all sequences
            * have the same length.
            *
            * \param ch the character to use for padding.
            * 
            * \return The final length obtained, which is the length of
            * the longest sequence before the operation.
            * 
            */
            unsigned int equalize(char ch='?');

        
           /** \brief Finds a sequence by its name
            * 
            * Gets the position of the first sequence with the specified
            * name.
            * 
            * \param string a sequence name.
            * 
            * \param strict if true, seeks an exact match. If false,
            * compares only until the end of the requested name (for
            * example: ATCFF will match ATCFF_01 if strict is false).
            * 
            * \return The lowest index where the name matches, -1 if no
            * sequence has such name.
            * 
            */
            int find(const char* string, bool strict=true) const;


        protected:
            // The number of sequences
            unsigned int _ns;
        
            // The array of name lengths
            unsigned int* lnames;
    
            // The array of names
            char** names;
        
            // The array of sequences (as c-strings)
            char** sequences;
        
            // The array of groups
            unsigned int* groups;
        
            // Imports an array of c-strings
            virtual void setFromSource(unsigned int number_of_sequences, const char* const* const cstring_array);
        
            // Constructor helper
            virtual void copyObject(const Container&);
        
            // Constructor partial helper
            virtual void getNamesAndGroups(const Container&);
        
        private:
       
            // The array of sequence lengths
            unsigned int* lsequences;
        
            // Setup a valid empty object
            virtual void init();
    };
}
    
#endif
