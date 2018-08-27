/*
    Copyright 2008-2010 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_MICROSATELLITEDIVERSITY_HPP
#define EGGLIB_MICROSATELLITEDIVERSITY_HPP

#include "DataMatrix.hpp"
#include <cstdlib>

namespace egglib {

   /** \brief Analyzes microsatellite data
    *
    * \ingroup polymorphism
    * 
    * Use the load() method to analyze data. All sites will be analyzed
    * and accessors allow to access the value of a given statistics for
    * a given site. There is no out-of-bound checking implemented in
    * accessors.
    * 
    */
    class MicrosatelliteDiversity {
        
          public:
            
           /** \brief Creates an object
            * 
            */
            MicrosatelliteDiversity();
            
            
           /** \brief Destroys an object
            * 
            */
            virtual ~MicrosatelliteDiversity();


           /** \brief Performs the analysis
            *
            * \param dataMatrix the object to analyze.
            * 
            * \param missingData the integer identifying missing data.
            * 
            * \param noMissingData if true, no allele will be
            * excluded (including the one identified by the argument
            * missingData).
            * 
            */
            void load(const DataMatrix& dataMatrix,
                    int missingData=999, bool noMissingData=false);
            
            
            /// Number of sites (or markers)
            unsigned int numberOfSites() const;
            
            /// Heterozygosity
            double He(unsigned int siteIndex) const;
            
            /// Number of alleles
            unsigned int numberOfAlleles(unsigned int siteIndex) const;
            
            /// Variance of allele size
            double sizeVariance(unsigned int siteIndex) const;
            
            /// IAM-based estimator of theta
            double thetaAssumingIAM(unsigned int siteIndex) const;
            
            /// SMM-based estimator of theta, calculated from He
            double thetaAssumingSMMfromHe(unsigned int siteIndex) const;

            /// SMM-based estimator of theta, calculated from VarSize
            double thetaAssumingSMMfromSizeVariance(unsigned int siteIndex) const;
            
            
        protected:
        
            unsigned int  v_numberOfSites;
            double       *v_He;
            unsigned int *v_numberOfAlleles;
            double       *v_sizeVariance;
            double       *v_thetaAssumingIAM;
            double       *v_thetaAssumingSMMfromHe;
            double       *v_thetaAssumingSMMfromSizeVariance;
            
            void init();
            void clear();
        
        
        private:
        
        
            /// No copy allowed
            MicrosatelliteDiversity(const MicrosatelliteDiversity& source) {
            }
            
            /// No copy allowed
            MicrosatelliteDiversity& operator=(const MicrosatelliteDiversity& source) {
                return *this;
            }
        
    };
}

#endif
