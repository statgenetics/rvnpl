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

#ifndef EGGLIB_BASEDIVERSITY_HPP
#define EGGLIB_BASEDIVERSITY_HPP

#include "CharMatrix.hpp"
#include "SitePolymorphism.hpp"
#include <string>

/** \defgroup polymorphism polymorphism
 *
 * \brief Diversity analyses
 *
 * Two classes are contained in this module: NucleotideDiversity, that
 * performs site-centered polymorphism analyses, and HaplotypeDiversity,
 * that performs haplotype-centered analyses. The detection of
 * polymorphic sites is common to both, through the base class
 * BaseDiversity. However this phase must be repeated when stats from
 * the two classes are needed. To reduce the computational burden, the
 * function reserve() can be use, that directly allocates needed memory
 * when the eventual number of polymorphic sites is known prior to
 * analysis (even if not precisely). For both classes, a set of
 * statistics are computed immediately upon load of a data set. For
 * NucleotideDiversity, additional statistics are computed per group
 * upon use of the corresponding accessors. This number of operations
 * performed several times is strictly limited. This is particularly
 * useful when different statistics are needed for a given alignment.
 * However, this system allows not computing unnecessary statistics to
 * a certain extend.
 * 
 */

namespace egglib {

    /** \brief Base class of diversity classes
    *
    * Mutualizes the analysis of polymorphic sites through the method
    * importSites() and related accessors.
    * 
    * \ingroup polymorphism
    *
    */
    class BaseDiversity {
    
        public:
    
           /** \brief Constructor
            * 
            */ 
            BaseDiversity();
            
           /** \brief Destructor
            * 
            */ 
            virtual ~BaseDiversity();
            
           /** \brief Reserve sufficient memory for a given number of
            * polymorphic sites.
            * 
            * This method makes importSite function faster when you
            * already know how many polymorphic sites to expect, since
            * the necessary memory will be allocated prior the screening
            * of data. It is possible to use reserve() even if with a
            * number of sites that is not matching what importSites()
            * will find.
            * 
            * \param numberOfSites a strictly positive integer.
            * 
            */
            virtual void reserve(unsigned int numberOfSites);

            /// Gets a site
            const SitePolymorphism* get_site(unsigned int index) const;

            /// Gets a site position
            unsigned int get_position(unsigned int index) const;

           /** \brief Predefined mapping string for DNA data
            * 
            */
            static const std::string dnaMapping;


           /** \brief Predefined mapping string for RNA data
            * 
            */
            static const std::string rnaMapping;


           /** \brief Predefined mapping string for amino acid data
            * 
            */
            static const std::string aaMapping;


            /// Clears and re-initializes object
            virtual void reset();


        protected:
    
            virtual void init();
            virtual void clear();
    
            // 
            void importSites(CharMatrix& data, bool allowMultipleMutations,
                double minimumExploitableData, unsigned int ignoreFrequency,
                std::string characterMapping, bool useZeroAsAncestral,
                bool ignoreOutgroup);

            // 
            void analyzeSite(CharMatrix& data, unsigned int index, double maxMissingData, bool ignoreOutgroup); // analyzes a site, adds a Site to the Site container if the site is polymorphic
            unsigned int getPopIndex(unsigned int label) const;  // returns v_npop if not found
            
            SitePolymorphism** v_sites;  // holder of polymorphic site addresses
            bool* v_orientables;         // stores whether the sites are orientable or not
            unsigned int* v_sitePositions;   // stores position of sites

            unsigned int  v_reserved;
            unsigned int  v_ns;       // maximum number of sequences analyzed (max of sites' ns)
            unsigned int  v_S;        // number of polymorphic sites
            unsigned int  v_So;       // number of orientable sites
            unsigned int  v_eta;      // number of mutation (whatever multiple)
            double        v_nseff;    // average number of analyzed sequence
            unsigned int  v_lseff;    // number of analyzed sites
            double        v_nseffo;   // average number of analyzed sequences for analyzes with outgroup
            unsigned int  v_lseffo;   // number of analyzed sites for analyzes with outgroup
            unsigned int  v_npop;     // number of populations
            unsigned int *v_popLabel; // label of each pop

            // options
            bool          p_allowMultipleMutations;
            double        p_minimumExploitableData;
            std::string   p_characterMapping;
            unsigned int  p_pos_sep_mapping;
            bool          p_useZeroAsAncestral;
            unsigned int  p_ignoreFrequency;


        private:

            BaseDiversity(const BaseDiversity& source) { }

            BaseDiversity& operator=(const BaseDiversity& source) {
                return *this;
            }

    };
}

#endif
