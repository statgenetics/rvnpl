/*
    Copyright 2008-2009,2013 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_BPPDIVERSITY_HPP
#define EGGLIB_BPPDIVERSITY_HPP

#include <vector>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Alphabet/RNA.h>
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>
#include <Bpp/Seq/Alphabet/StandardCodonAlphabet.h>
#include <Bpp/Seq/Alphabet/VertebrateMitochondrialCodonAlphabet.h>
#include <Bpp/Seq/Alphabet/InvertebrateMitochondrialCodonAlphabet.h>
#include <Bpp/Seq/Alphabet/EchinodermMitochondrialCodonAlphabet.h>
#include <Bpp/Seq/GeneticCode/StandardGeneticCode.h>
#include <Bpp/Seq/GeneticCode/VertebrateMitochondrialGeneticCode.h>
#include <Bpp/Seq/GeneticCode/InvertebrateMitochondrialGeneticCode.h>
#include <Bpp/Seq/GeneticCode/EchinodermMitochondrialGeneticCode.h>
#include <Bpp/PopGen/SequenceStatistics.h>
#include "Align.hpp"

namespace egglib {

   /** \brief Computes diversity statistics using third-party library Bio++
    *
    * \ingroup polymorphism
    *
    * Uses the PopGen library of <A href="http://kimura.univ-montp2.fr/BioPP/">Bio++</A>
    * and computes statistics of diversity.
    * 
    */
    class BppDiversity {
        public:

           /** \brief Builds an object
            * 
            */
            BppDiversity();


           /** \brief Destroys an object
            * 
            */
            virtual ~BppDiversity();


           /** \brief Performs polymorphism analysis using Bio++ tools
            * 
            * Automatically computes statistics. Some statistics
            * require the presence of outgroup sequences (with group
            * labels 999) in the passed object, and some statistics
            * are computed only if the data type is >=4. (Some both.)
            * Calling this method clears all previous data in the
            * instance (including if less statistics can be computed
            * with the new alignment compared with the previous one).
            * 
            * \param align an alignment object (Align class).
            * The presence of outgroup sequences will be automatically
            * detected based on the group member of the passed object.
            * Only label 999 will be considered and use as outgroup for
            * statistics were it is useful. If several outgroups are
            * used, they will all be passed to C++ tools.
            * 
            * \param dataType an integer what kind of data must be
            * analyzed:
            *       - 1 for DNA
            *       - 2 for RNA
            *       - 3 for protein sequences
            *       - 4 for standard codons
            *       - 5 for vertebrate mitochondrial codons
            *       - 6 for invertebrate mitochondrial codons
            *       - 7 for echinoderm mitochondrial codons
            *       .
            * Other values will result in an exception.
            * 
            */
            void load(Align& align, unsigned int dataType=1);


            /// true if an ougroup was available
            bool hasOutgroup() const;
            
            /// Number of polymorphic sites
            unsigned int S() const;

            /// Number of parsimony informative sites
            unsigned int Sinf() const;

            /// Number of singleton sites
            unsigned int Ssin() const;

            /// Minimal number of mutations
            unsigned int eta() const;

            /// Mutations on external branches (if an outgroup was used)
            unsigned int Sext() const;

            /// Heterozygosity
            double He() const;

            /// Squared heterozygosity
            double He2() const;

            ///// Average GC content
            //double GC() const;

            /// Watterson theta
            double tW() const;

            /// Tajima theta
            double T83() const;

            /// Number of haplotypes
            unsigned int K() const;

            /// Haplotypic diversity
            double H() const;

            /// Number of transitions
            unsigned int Ti() const;

            /// Number of transversions
            unsigned int Tv() const;

            /// Transition/transversion ratio.
            double TiTv() const;

            /// Number of codon sites with a stop codon (if coding)
            unsigned int nstop() const;

            /// Number of codon sites with one change (if coding)
            unsigned int ncodon1mut() const;

            /// Number of codon sites with a synonymous change (if coding)
            unsigned int nsyn() const;

            /// Synonymous Watterson theta (if coding)
            double tWS() const;

            /// Non-synonymous Watterson theta (if coding)
            double tWNS() const;

            /// Synonymous Pi (if coding)
            double PiS() const;

            /// Non-synonymous Pi (if coding)
            double PiNS() const;

            /// Average number of synonymous sites (if coding)
            double Ssites() const;

            /// Average number of non-synonymous sites (if coding)
            double NSsites() const;

            /// Aynonymous polymorphic sites (if coding)
            unsigned int SS() const;

            /// Non-synonymous polymorphic sites (if coding)
            unsigned int SNS() const;

            /// McDonald-Kreitman test table, as a vector (if coding + outgroup)
            const std::vector<unsigned int>& MK() const;

            /// Neutrality index (if coding + outgroup)
            double NI() const;

            /// Tajima's D
            double D() const;

            /// Tajima's D computed with eta instead of S
            double Deta() const;

            /// Fu and Li's D (if an outgroup was used)
            double Dfl() const;

            /// Fu and Li's D*
            double Dflstar() const;

            /// Fu's F (if an outgroup was used)
            double F() const;

            /// Fu's F*
            double Fstar() const;

            /// Hudson's estimator of rho (recombination rate)
            double rhoH() const;


        protected:

            // Initializes to default values
            void reset();

           /* \brief Compute within-population statistics
            * 
            * \param align the sequence alignment to analyze.
            * 
            */
            void computeSingle(bpp::AlignedSequenceContainer& align);


           /* \brief Compute between-population statistics
            * 
            * \param align1 the ingroup sequence alignment.
            * \param align2 the outgroup sequence alignment.
            * 
            */
            void computeDouble(bpp::AlignedSequenceContainer& ingroup, bpp::AlignedSequenceContainer& outgroup);


            // Bio++ stuff
            bpp::DNA dna;
            bpp::GeneticCode* p_code;
            bpp::StandardGeneticCode* p_code4;
            bpp::VertebrateMitochondrialGeneticCode* p_code5;
            bpp::InvertebrateMitochondrialGeneticCode* p_code6;
            bpp::EchinodermMitochondrialGeneticCode* p_code7;

            unsigned int v_Ss;
            unsigned int v_Sinf;
            unsigned int v_Ssin;
            unsigned int v_eta;
            unsigned int v_Sext;
            double v_He;
            double v_He2;
            //double v_GC;
            double v_tW;
            double v_T83; 
            unsigned int v_K;
            double v_H;
            unsigned int v_Ti;
            unsigned int v_Tv;
            double v_TiTv;
            unsigned int v_nstop;
            unsigned int v_ncodon1mut;
            unsigned int v_nsyn;
            double v_tWS;
            double v_tWNS;
            double v_PiS;
            double v_PiNS;
            double v_Ssites;
            double v_NSsites;
            unsigned int v_SS;
            unsigned int v_SNS;
            std::vector<unsigned int> v_MK;
            double v_NI;
            double v_D;
            double v_Deta;
            double v_Dfl;
            double v_Dflstar;
            double v_F;
            double v_Fstar;
            double v_rhoH;
            
            // true if data loaded
            bool b_loaded;
            bool b_outgroup;  // implies loaded
            bool b_coding;  // implies coding
            
            
        private:
        
            /// No copy allowed for this class
            BppDiversity(const BppDiversity& source) { }
            
            /// No copy allowed for this class
            BppDiversity& operator=(const BppDiversity& source) { return *this; }
            
    };
}
    
#endif
