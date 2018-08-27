/*
    Copyright 2009,2013 Stéphane De Mita, Mathieu Siol

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


#ifndef EGGLIB_CONVERT_HPP
#define EGGLIB_CONVERT_HPP

#include "DataMatrix.hpp"
#include "Align.hpp"
#include "EggException.hpp"
#include "Random.hpp"
#include <string>

#include "config.h"

#ifdef HAVE_LIBBPP_SEQ
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Alphabet/RNA.h>
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>
#include <Bpp/Seq/Alphabet/StandardCodonAlphabet.h>
#include <Bpp/Seq/Alphabet/VertebrateMitochondrialCodonAlphabet.h>
#include <Bpp/Seq/Alphabet/InvertebrateMitochondrialCodonAlphabet.h>
#include <Bpp/Seq/Alphabet/EchinodermMitochondrialCodonAlphabet.h>
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#endif

namespace egglib {


   /** \brief Performs conversion between sequence holder types
    *
    * \ingroup core
    * 
    * Static methods of this class allows conversion between sequence
    * holder types implying parametrizable modifications.
    * 
    */
    class Convert {

        public:
        
           /** \brief DataMatrix to Align conversion
            * 
            * By defaut, this method generates an Align instance
            * containing only the polymorphic sites. The integers of
            * the DataMatrix will be converted as follow: 0 to A, 1 to
            * C, 2 to G and 3 to T. This behaviour can be largely
            * modified using options.
            * 
            * \param dataMatrix DataMatrix instance.
            * 
            * \param length length of the desired alignment. Non-varying
            * stretches of data will be introduced to reach the
            * specified length. By default the positions of segregating
            * sites will be determined from the positions given by the
            * DataMatrix object. Those positions are expressed in a
            * continuous range, and will be discretized. Mutations
            * falling on the same site will be moved of one position
            * left or right (always preserving the order of mutation
            * sites). If positions are all zero (the default of the
            * DataMatrix class) and if length is larger than the number
            * of segregating sites, then all segregating sites will
            * cluster on the left-hand side of the alignment.
            * 
            * \param random the address to a Random object allowing to 
            * draw random numbers (for randomizing positions and/or
            * non-varying states). If an address is provided but no
            * random numbers are required, it is ignored. If no address
            * if provided and random numbers are required, a Random
            * instance is built internally.
            * 
            * \param randomizePositions if true, the positions specified
            * in the DataMatrix objects are ignored and the positions of
            * mutations are drawn randomly along the interval (only if
            * the specified length is larger than the number of
            * segregating sites). If randomizePositions and false and
            * positions are not
            * 
            * \param enforceLength specify whether a
            * EggRuntimeError should be thrown when the number of
            * polymorphic sites is larger than the specified length. If
            * false (the default) and in cases where the specified
            * length is too short to harbor all polymorphic sites, the
            * alignment length will be increased as needed.
            * 
            * \param randomizeNonVaryingStates if true, the stretches of
            * conserved positions (between segregating sites) will be
            * randomly drawn from the current symbol mapping. Otherwise,
            * the symbol given by fixed will be used.
            * 
            * \param randomizeAlleles if true, alleles will be drawn
            * randomly from the mapped characters. Note that if a
            * genotype value is larger than the size of the mapping, it
            * will be replaced by the character given by unknown,
            * without randomization. In other words, with the mapping
            * "ACGT", alleles 0, 1, 2 and 3 will be randomly assigned
            * to these four characters, but larger and negative alleles
            * will be assigned to the unknown character.
            * 
            * \param mapping a string given the character to assign to
            * different character values read from the DataMatrix. If
            * the read value is 0, the first character of the string
            * will used, the the value is 1, the second character will
            * be used, and so on. If the integer read is out of range
            * (in particular, for any negative value), then the
            * character given by unknown will be used. An empty string
            * will always lead to alignments containing only the
            * character given by unknown. The string "01" is suitable
            * for binary data.
            * 
            * \param unknown the character to use if an integer genotype
            * value is not mapped in the mapping string (that is, if
            * the mapping string is too short).
            * 
            * \param nonVaryingState character to use for conserved
            * stretches of data. It doesn't have to be included in the
            * mapping. If randomizeNonVaryingState is true, this
            * argument is ignored.
            * 
            * \return The resulting Align object.
            * 
            */
            static Align align(
                DataMatrix& dataMatrix,
                unsigned int length=0,
                Random* random=NULL,
                bool randomizePositions=false,
                bool randomizeNonVaryingStates=false,
                bool randomizeAlleles=false,
                bool enforceLength=false,
                std::string mapping="ACGT",
                char unknown='?',
                char nonVaryingState='A'
            );


#ifdef HAVE_LIBBPP_SEQ

            /** \brief Converts an alignment to the equivalent Bio++ type
            *
            * During conversion, name information is lost (arbitrary
            * names are generated in order toprevent duplicate names).
            * The object is attached to an alphabet matching the passed
            * integer. The names are bare rank integers (starting at the
            * value giving by *offset*).
            *
            * \param align the source alignment object.
            * 
            * \param alphabetID an integer indicating which alphabet to
            * use:
            *       - 1 for DNA
            *       - 2 for RNA
            *       - 3 for proteins
            *       - 4 for standard codon
            *       - 5 for vertebrate mitochondrial codon
            *       - 6 for invertebrate mitochondrial codon
            *       - 7 for echinoderm mitochondrial codon
            *       .
            * Other values will result in an exception.
            * 
            * \param outgroupFlag an integer indicating whether to
            * include outgroup sequences:
            *       - 0 use all sequences
            *       - 1 use only sequences without 999 label (ingroup)
            *       - 2 use only sequences with 999 label (outgroup)
            *       .
            * Other values will result in an exception.
            * 
            * \param offset enter an integer to shift the names of the
            * resulting alignment (useful to merge alignment and ensure
            * that names are not duplicated).
            * 
            * \return A Bio++ alignment.
            * 
            */
            static bpp::AlignedSequenceContainer egglib2bpp(Align& align, unsigned int alphabetID, unsigned int outgroupFlag, unsigned int offset=0);

#endif



        protected:

           /** \brief This class cannot be instantiated
            * 
            */
            Convert() { }


           /** \brief This class cannot be instantiated
            * 
            */
            Convert(const Convert& source) { }


           /** \brief This class cannot be instantiated
            * 
            */
            Convert& operator=(const Convert& source) { return *this; }


           /** \brief This class cannot be instantiated
            * 
            */
            virtual ~Convert() { }

#ifdef HAVE_LIBBPP_SEQ
            static bpp::DNA dnaAlphabet;
            static bpp::RNA rnaAlphabet;
            static bpp::ProteicAlphabet proteicAlphabet;
            static bpp::StandardCodonAlphabet standardCodonAlphabet;
            static bpp::VertebrateMitochondrialCodonAlphabet vertebrateMitochondrialCodonAlphabet;
            static bpp::InvertebrateMitochondrialCodonAlphabet invertebrateMitochondrialCodonAlphabet;
            static bpp::EchinodermMitochondrialCodonAlphabet echinodermMitochondrialCodonAlphabet;
#endif

    };
}

#endif
