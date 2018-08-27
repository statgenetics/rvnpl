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

#ifndef EGGLIB_LINKAGEDISEQUILIBRUM_HPP
#define EGGLIB_LINKAGEDISEQUILIBRUM_HPP


#include "BaseDiversity.hpp"
#include "EggException.hpp"


namespace egglib {

   /** \brief Analyzes linkage disequilibrium per pair of polymorphic sites
    *
    * \ingroup polymorphism
    * 
    * The class considers an alignment and detects polymorphic sites
    * using the BaseDiversity functionality (shared with other classes
    * of the module). Only sites with exactly two alleles are
    * considered. Statistics of pairwise linkage disequilibrium can
    * be accessed by pair index (note that out-of-range errors are not
    * checked). Population labels are ignored (but outgroups are
    * excluded from the analysis).
    *
    */
    class LinkageDisequilibrium : public BaseDiversity {

      public:

        /// Default constructor
        LinkageDisequilibrium();

        /// Destructor
        virtual ~LinkageDisequilibrium();

       /** \brief Analyzes polymorphic sites of an alignment
        * 
        * \param data an alignment object (subclass of CharMatrix).
        * The presence of outgroup or of different populations will
        * be detected based on the populationLabel members of the
        * passed object. The populationLabel 999 will be interpreted
        * as outgroups. If several outgroups are passed, sites were
        * the outgroups are not consistent will be treated as "non-
        * orientable".
        * 
        * \param minimumExploitableData site where the non-missing
        * data (as defined by characterMapping) are at a frequency
        * larger than this value will be removed from the analysis.
        * Use 1. to take only 'complete' sites into account and 0.
        * to use all sites.
        * 
        * \param ignoreFrequency removes sites that are polymorphic
        * because of an allele at absolute frequency smaller than or
        * equal to this value. If ignoreFrequency=1, no sites are
        * removed, if ignoreFrequency=1, singleton sites are
        * ignored. Such sites are completely removed from the
        * analysis (not counted in lseff). Note that if more than
        * one mutation is allowed, the site is removed only if all
        * the alleles but one are smaller than or equal to this
        * value. For example, an alignment column AAAAAAGAAT is
        * ignored with an ignoreFrequency of 1, but AAAAAAGGAT is
        * conserved (including the third allele T which is a
        * singleton).
        * 
        * \param characterMapping a string giving the list of
        * characters that should be considered as valid data. If a
        * space is present in the string, the characters left of the
        * space will be treated as valid data and the characters
        * right of the space will be treated as missing data, that
        * is tolerated but ignored. All characters not in the string
        * will cause an EggInvalidCharacterError to be raised.
        */
        void load(CharMatrix& data,
                double minimumExploitableData=1.,
                unsigned int ignoreFrequency=0,
                std::string characterMapping=dnaMapping);
    
    
        /// Number of pairs contained in the instance
        unsigned int numberOfPairs() const;

        /// Alignment distance between a given pair
        int d(unsigned int pair_index);

        /// D statistic for a given pair
        double D(unsigned int pair_index);

        /// D' statistic for a given pair
        double Dp(unsigned int pair_index);

        /// r statistic for a given pair
        double r(unsigned int pair_index);

        /// r2 statistic for a given pair
        double r2(unsigned int pair_index);

        /// position of the first site for a given pair
        unsigned int site1(unsigned int pair_index);

        /// position of the second site for a given pair
        unsigned int site2(unsigned int pair_index);

        /// correlation coefficient between r2 and distance
        double correl() const;
        
       /** \brief Computes the minimal number of recombination events
        * 
        * The computation is performed as described in Hudson, RR and
        * NL Kaplan. 1985. Statistical properties of the number of
        * recombination events in the history of a sample of DNA
        * sequences. Genetics 111: 147-164. The returned parameter is
        * the minimal number of recombination events, given by the
        * number of non-overlapping pairs of segregating sites violating
        * the rule of the four gamete. Only sites with two alleles are
        * considered. Note that homoplasy (multiple mutations) mimicks
        * recombination. The result of this function is not stored
        * in this instance, and re-computed at each call.
        * 
        * \param data the same CharMatrix instance as passed to the load
        * method. The instance must not have been modified.
        * 
        */
        unsigned int Rmin(CharMatrix& data) const;



      protected:
      
        // adds a pair of polymorphic sites
        // assume position2>position1,
        //  sites are polymorphic with exactly 2 alleles
        void add(CharMatrix& data, unsigned int position1, unsigned int position2);

        // Constructor help
        void init();
                
        // Destructor helper
        void clear();
        
        // Resizes arrays
        void reset();
        
        // Small helper
        inline double min(double a, double b) { return (a>b)?a:b;}

        // Small helper
        inline double max(double a, double b) { return (a>b)?b:a;}

        // Small helper
        inline void check(unsigned int pos) {  if (pos>=_n) throw EggRuntimeError("tried to access an invalid index"); }

       /* Performs correlation
        *
        * This function works independently from the rest of the class.
        *
        * \param n length of data arrays.
        * \param x first data vector.
        * \param y second data vector.
        * \param r variable to receive the correlation coefficient.
        * \param a variable to receive the regression slope.
        */
        static void _correl(unsigned int n, const int* x, const double* y, double& r, double& a);

        // Distance between pairs
        int* _d;
        
        // D (classical) measure of LD
        double *_D;
        
        // D'
        double *_Dp;
        
        // r, correlation coefficient
        double *_r;
        
        // square r
        double *_r2;
        
        // Data array (not managed by the instance)
        unsigned int *_site1;

        // Data array (not managed by the instance)
        unsigned int *_site2;
        
        // Number of pairs
        unsigned int _n;

     private:
     
        /// Copy constructor not available
        LinkageDisequilibrium(const LinkageDisequilibrium&) { }

        /// Assignment operator not available
        LinkageDisequilibrium& operator=(const LinkageDisequilibrium&) {
            return *this;
        }


        class Interval {
            public:
                Interval(unsigned int, unsigned int);
                unsigned int a() const;
                unsigned int b() const;
                bool good() const;
                void set_false();
            private:
                unsigned int _a;
                unsigned int _b;
                unsigned int _good;
        };


  };
}

#endif
