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

#ifndef EGGLIB_DATAMATRIX_HPP
#define EGGLIB_DATAMATRIX_HPP


#include "CharMatrix.hpp"


namespace egglib {

    /** \brief Data table
    *
    * \ingroup core
    *
    * Holds a data matrix representing genotype data from simulations.
    * Data are stored as integers, to each site is associated a
    * position, and to each sequence is associated a group index (any
    * integer labelling, for example, a subpopulation). Supports the
    * CharMatrix interface with the condition that allele genotype
    * datum is >=0 and <=9.
    * 
    */
    class DataMatrix : public CharMatrix {
    
    public:
    
        /** \brief Default constructor
         *
         * The data table default dimensions are {0,0}
         *
         */
        DataMatrix();


        /** \brief Standard constructor
         *
         * The data table dimensions must be given.
         * Each cell default default is 0, and each site position is 0..
         *
         * \param numberOfSequences number of lines of the data table.
         * \param numberOfSites number of columns of the data table.
         *
         */
        DataMatrix(unsigned int numberOfSequences, unsigned int numberOfSites);


        /** \brief Copy constructor
         * 
         */
        DataMatrix(const DataMatrix&);


        /** \brief Copy constructor
         * 
         */
        DataMatrix(const CharMatrix&);

        
        /** \brief Assignment operator
         * 
         */
        virtual DataMatrix& operator=(const DataMatrix&);
        
        
        /** \brief Assignment operator
         * 
         */
        virtual DataMatrix& operator=(const CharMatrix&);


        /** \brief Destructor
         * 
         */
        virtual ~DataMatrix();
        
        
        /** \brief Gets number of sites
         * 
         */
        unsigned int numberOfSites() const;


        /** \brief Gets number of sequences
         * 
         */
        unsigned int numberOfSequences() const;
        
        
        /** \brief Sets a value of the data table
         * 
         */
        void set(unsigned int sequence, unsigned int site, int value);


        /** \brief Gets a value from the data table
         * 
         */
        int get(unsigned int sequence, unsigned int site) const;
        
        
        /** \brief Faster and unsecure version of get
         * 
         */
        inline int fget(unsigned int sequence, unsigned int site) const {
            return dataMatrix[sequence][site];
        }


        /** \brief Sets the position of a site
         * 
         */
        void sitePosition(unsigned int site, double value);


        /** \brief Gets the position of a site
         * 
         */
        double sitePosition(unsigned int site) const;


        /** \brief Sets the group label of a sequence
         * 
         */
        void populationLabel(unsigned int sequence, unsigned int value);


        /** \brief Gets the group label of a sequence
         * 
         */
        unsigned int populationLabel(unsigned int sequence) const;


        /** \brief Removes all information from the object
         * 
         */
        void clear();


        /** \brief Resizes the data matrix
         *
         * \param newNumberOfSequences number of sequences (rows)
         * \param newNumberOfSites number of sites (columns)
         *
         * If new values are larger, data already set is left unchanged.
         * New data are set to zero.
         *
         */
        void resize(unsigned int newNumberOfSequences, unsigned int newNumberOfSites);


       /** \brief Shifts allele value
        * 
        * \param minimum the minimum allele value.
        * 
        * Shifts all alleles at all sites to ensure that alleles alleles
        * are equal to or larger than minimum. The shifting is specific
        * to each site.
        * 
        */
        void shift(int minimum);

        /** \brief Gets the character at a given position
         * 
         * An exception is generated if the allele value at this
         * position is not >=0 and <=9. Not out-of-bound check is
         * performed.
         * 
         */
        char character(unsigned int sequence, unsigned int site) const;
        


    private:
        
        // Initializes to default values (for empty object)
        void init();
        
        // Copies from a source object
        virtual void copy(const CharMatrix&);

        // Copies from a source object
        virtual void copy(const DataMatrix&);
        
        // Number of lines of the data matrix
        unsigned int _numberOfSequences;
        
        // Number of columns of the data matrix
        unsigned int _numberOfSites;
        
        // Data matrix
        int **dataMatrix;
        
        // Vector of site positions
        double *positions;
        
        // Vector of group indices
        unsigned int *groups;
    };
}

#endif
