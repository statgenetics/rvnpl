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

#ifndef EGGLIB_CONSENSUS_HPP
#define EGGLIB_CONSENSUS_HPP

#include "Align.hpp"
#include <sstream>
#include <string>
#include <vector>

namespace egglib {

   /** \brief Generates consensus sequences
    *
    * \ingroup polymorphism
    * 
    *
    * A consensus is generated when two sequences have the same name, 
    * ignoring everything after the first separator character (by
    * default, "_"). Hence, the names "foo", "foo_goo" and "foo_third"
    * will be treated as identical and the root will be "foo". The root
    * will be used to name the resulting sequence. Note that the
    * class works only for DNA sequences.
    *
    * Symbol convention:
    *     - A: adenosine
    *     - C: cytosine
    *     - G: guanine
    *     - T: thymine
    *     - M: A or C
    *     - R: A or G
    *     - W: A or T (weak)
    *     - S: C or G (strong)
    *     - Y: C or T
    *     - K: G or T
    *     - B: C or G or T(not A)
    *     - D: A or G or T (not C)
    *     - H: A or C or T (not G)
    *     - V: A or C or G (not T)
    *     - N: A or C or G or T
    *     - ?: nonsequenced position
    * 
    * Other symbols will be treated as ? (lowercase are supported).
    * 
    * Rigorous (alias liberal or strong) mode:
    *    - If two characters are the same, it is retained whatever it is
    * (A + A = A)
    * - Otherwise:
    *        - If one is the missing character (?) the other is retained
    * whatever it is (A + ? = A).
    *             - If characters are consistent, that is one contains
    * more information, that one is retained (A + M = A).
    *             - If characters are not consistent, the closest 
    * generic symbol is retained (A + C = M).
    *    .
    *     Note that the feedback of inconsistent characters in the
    * outcome is not garanteed.
    *     In fact, (A + A + G) will result in R (as expected) but (A +
    * G + A) will result in A, masking the problem.
    *     However, the position will indeed be counted as inconsistent.
    *  
    * Not rigorous (conservative/weak) mode:
    *     - If two characters are the same, it is retained whatever it
    * is (A + A = A).
    *     - Otherwise:
    *        - If one is ? the other is retained whatever it is (A + ?
    * = A).
    *        - Otherwise an inconsistent character (by default, Z) is
    * retained (A + C = Z).
    * 
    * Iterative process of consensus:
    *     - Each sequence is taken in turn.
    *     - Each pair involving the focus sequence is processed and a
    * consensus is generated.
    *     - When all pair have been processsed, the consensus already
    * generated are themselves iteratively processed until only one
    * remains.
    *     - Note that at each time the last two are taken first.
    * 
    * A transparent interface gives access to the data for all steps of
    * the consensus process, as vectors that covers all pairs (including
    * intermediate steps of the iterative procedure described above) as
    * well as singleton sequences. For the latter, the second name is
    * not filled and all counts are set to 0. Note also that the name of
    * such singleton sequence is shortened to the separator as well.
    * 
    */
    class Consensus {
       public:
         /** \brief Constructor
          * 
          */
          Consensus();
          
         /** \brief Destructor
          * 
          */
          virtual ~Consensus() {}
          
          /// Sets the character interpreted as missing (default: ?)
          void setMissing(char);
          
          /// Sets the character used to point to disagreements (default: Z)
          void setDisagreement(char);
          
          /// Checks all the characters
          bool check_sequences(Align& align);
          
         /** \brief Reduces the sequence alignment by making consensus sequences
          * 
          * \param align the original alignment.
          * 
          * \param separator the character used to separated the root
          * name of sequences to the variable part, as in (for the
          * default value: "sequence_read1".
          * 
          * \param rigorous consensus mode.
          * 
          * \return An Align instance with duplicated sequences consensed.
          * 
          */
          Align consensus(Align& align, char separator='_', bool rigorous=true);

          /// First name of consensed pairs
          const std::vector<std::string>& firstSequenceNames();

          /// Second names of consensed pairs
          const std::vector<std::string>& secondSequenceNames();
          
          /// Root names of consensed pairs
          const std::vector<std::string>& roots();
          
          /// Number of consistent positions for all consensed pairs
          const std::vector<int>& consistentPositions();

          /// Number of complementary positions for all consensed pairs
          const std::vector<int>& complementaryPositions();

          /// Number of uninformative positions for all consensed pairs
          const std::vector<int>& uninformativePositions();

          /// Number of ambiguous positions for all consensed pairs
          const std::vector<int>& ambiguousPositions();

          /// Number of at least partially resolved ambiguities for all consensed pairs
          const std::vector<int>& atLeastPartiallyResolvedAmbiguities();

          /// Vector of inconsistent positions ofr all consensed pairs
          const std::vector<std::vector<int> >& inconsistentPositions();


       private:

           /** \brief Copying this class is not allowed
            * 
            */
            Consensus(const Consensus& source) { }
            
           /** \brief Copying this class is not allowed
            * 
            */
            Consensus& operator=(const Consensus& source) { return *this; }
            
          // A private helper
          class PairwiseConsensus;
         
          // report data
          std::vector<std::string> t_firstSequenceNames;
          std::vector<std::string> t_secondSequenceNames;
          std::vector<std::string> t_roots;
          std::vector<int> t_consistentPositions;
          std::vector<int> t_complementaryPositions;
          std::vector<int> t_uninformativePositions;
          std::vector<int> t_ambiguousPositions;
          std::vector<int> t_atLeastPartiallyResolvedAmbiguities;
          std::vector<std::vector<int> > t_inconsistentPositions;

          // Code for missing data (usually ?)
          char MISSING;

          // Code for disgrement
          char DISAGREEMENT;

          // Helper class managing a single pair
          class PairwiseConsensus {
              public:
                 // Default object creation
                 PairwiseConsensus();
                 
                 // Object destruction
                 virtual ~PairwiseConsensus() {}
                 
                 // Usual object creation
                 PairwiseConsensus(std::string, std::string);
                 
                 // Fills an object created with the default constructor
                 void load(std::string,std::string);
                 
                 // Changes the MISSING character
                 void setUndeterminedCharacter(char);
                 
                 // Changes the DISAGREEMENT character
                 void setDisagreementCharacter(char);
                 
                 /* Uses the conservative mode of consensus
                  * 
                  * Tries to avoid to make decisions, and adds the
                  * character set by DISAGREEMENT upon inconsistencies
                  * 
                  * return The number of inconsistencies.
                  */
                 int generateSoftConsensus();
                 
                 /* Strict mode of consensus
                  * 
                  * The number of inconsistencies.
                  */
                 int generateHardConsensus();
                 
                 // Two fully resolved (including gap) and identical characters
                 int getConsistentPositions();

                 // One informative (including gap) and one missing
                 int getComplementaryPositions();

                 // None missing, but different and incompatible
                 int getInconsistentPositions();
                 
                 // Both are missing
                 int getUninformativePositions();
                 
                 // Both identical or one missing, but not fully resolved
                 int getAmbiguousPositions();
                 
                 // Different, not missing, complementary.
                 int getAtLeastPartiallyResolvedAmbiguities();
                 
                 // Accessor
                 int getThisInconsistentPosition(unsigned int);
                 
                 // Generates the consensus sequence
                 std::string getConsensus();

              private:
              
                inline bool isValid(char c) {
                    switch (c) {
                        case 'A': case 'a':
                        case 'C': case 'c':
                        case 'G': case 'g':
                        case 'T': case 't': 
                            return true;
                        default:
                            return false;
                    }
                }
              
                 // This initiates a series of embedded objects
                 void setCharacterContainers();
                 
                 // The first sequence
                 std::string seqA;
                 
                 // The second sequence
                 std::string seqB;
                 
                 // The resulting consensus
                 std::string cons;
                 
                 // The vecotr storing the inconsistent positions
                 std::vector<int> posIncons;
                 
                 // The length of the sequences
                 unsigned int ls;
                 
                 // Counter
                 int cntConsistentPositions;
                 
                 // Counter
                 int cntComplementaryPositions;
                 
                 // Counter
                 int cntAmbiguousPositions;
                 
                 // Counter
                 int cntInconsistentPositions;
                 
                 // Counter
                 int cntUninformativePositions;
                 
                 // Counter
                 int cntAtLeastPartiallyResolvedAmbiguities;

                 // Code for missing data (usually ?)
                 char MISSING;

                 // Code for disgrement
                 char DISAGREEMENT;
     
            public:
                // This class manages relationships different symbols
                class CharacterContainer {
                      public:
                         // Default value: @
                         CharacterContainer();
                         
                         // Initiates to a given symbol
                         CharacterContainer(const char&);
                         
                         // Assignment operator
                         CharacterContainer& operator=(const char&);
                         
                         // Sets the symbol
                         void setValue(char);
                         
                         // Set the descendants
                         void setSons(std::vector<CharacterContainer>);
                         
                         // Tests whether the symbol is the same
                         bool is(CharacterContainer);
                         
                         // Tests if the query is contained amongst the sons
                         bool has(CharacterContainer);
                         
                         // Tests if the query is contained amongst the sons
                         bool has(char);
                         
                         /* Tests whether the left character has the left one
                          * Should be called on the N object only.
                          */
                         char lhas(CharacterContainer,CharacterContainer);
                         
                         /* Creates the object with the proper sons
                          * Should be called on the N object only.
                          */
                         CharacterContainer init(char);
                         
                         // The symbol
                         char value;
                         
                         // The descendants
                         std::vector<CharacterContainer> sons;
                  };

             private:
                  // Symbol ?
                  CharacterContainer ccQ;
                  
                  // Symbol A
                  CharacterContainer ccA;
                  
                  // Symbol C
                  CharacterContainer ccC;
                  
                  // Symbol G
                  CharacterContainer ccG;
                  
                  // Symbol T
                  CharacterContainer ccT;
                  
                  // Symbol U
                  CharacterContainer ccU;
                  
                  // Symbol M
                  CharacterContainer ccM;
                  
                  // Symbol R
                  CharacterContainer ccR;
                  
                  // Symbol W
                  CharacterContainer ccW;
                  
                  // Symbol S
                  CharacterContainer ccS;
                  
                  // Symbol Y
                  CharacterContainer ccY;
                  
                  // Symbol K
                  CharacterContainer ccK;
                  
                  // Symbol B
                  CharacterContainer ccB;
                  
                  // Symbol D
                  CharacterContainer ccD;
                  
                  // Symbol H
                  CharacterContainer ccH;
                  
                  // Symbol V
                  CharacterContainer ccV;
                  
                  // Symbol N
                  CharacterContainer ccN;
                  
                  // Symbol -
                  CharacterContainer ccGAP;
          };
    };
}

#endif

