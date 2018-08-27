/*
    Copyright 2011 Stéphane De Mita
    
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


#ifndef EGGLIB_ABC_HPP
#define EGGLIB_ABC_HPP

#include <vector>
#include <string>


namespace egglib {

   /** \brief Model estimation by Approximate Bayesian Computation
    *
    * \ingroup abc
    * 
    * It is required to set the number of statistics and at least one
    * input file name before performing analysis. The analysis itself
    * consists in several steps. (1) Computation of the threshold, which
    * requires to read through all files and imports statistics. In the
    * process, the standard deviation of all statisticcs will be
    * calculated and will be available. (2) Computation of Euclidean
    * distances and weights and generation of a second sample file with
    * weights and (non-standardized) statistics with only non-null
    * weights exported. This step requires that the observed statistics
    * have been set (between steps 1 and 2). (3) Local-linear regression
    * using the fit method. While in the previous steps several models
    * can be mixed, in this step a single model can be processed at a
    * time. The output is a simple file with adjusted parameters only.
    *
    */
    class ABC {
        
        public:
        
        
            /// \brief Modes for parameter transformation
			enum TransformMode { NONE, LOG, TAN };


            /// Constructor
            ABC();

            
            /// Destructor
            ~ABC();

          
           /** \brief Sets number of statistics
            * 
            * If data was already present in the instance, it will all
            * be cleared.
            * 
            */
            void number_of_statistics(unsigned int ns);


            /// Adds a file name
            void add_fname(const char* fname, unsigned int number_of_params);

            
            /// Gets number of imported data samples
            unsigned int number_of_samples() const;


           /** \brief Gets regression threshold
            * 
            * If several files are loaded, the data will be aggregated
            * (note that they must all contain the same number of
            * statistics). At least one file must have been set, and the
            * number of statistics must have been set as well.
            * 
            * \param tolerance rejection threshold (proportion of points
            * in the local region.
            * 
            */
            void get_threshold(double tolerance);


           /** \brief Gets a standard deviation
            * 
            * The threshold() method must have been called, and the
            * index must not be out of bound.
            * 
            */
            double sd(unsigned int index) const;


           /** \brief Sets a summary statistics
            * 
            * The number of statistics must have been set, and the index
            * must not be out of bound.
            * 
            */
            void obs(unsigned int index, double value);

                        
           /** \brief Performs rejection step
            * 
            * The observed value must have been entered (otherwise the
            * results will be meaningless), and the threshold must have
            * been computed.
            * 
            * \param outfname the name of the intermediary file.
            * \param exportlabels if true: exports a tag at the
            *        beginning of each line to identify the file or
            *        origin of each accepted sample (starting from 1).
            * \return the number of points in the local region.
            * 
            */
            unsigned int rejection(const char* outfname, bool exportlabels=false);
            
           /** \brief Performs regression step
            * 
            * \param infname input file name (generated using rejection)
            * \param outfname output file name (final posterior)
            * \param mode transformation mode
            * \param header outfile file header (name of parameters; if
            *        empty string, no header is printed)
            *
            * \return Number of data point processed
            */
            unsigned int regression(const char* infname, 
                                    const char* outfname,
                                    TransformMode mode,
                                    const char* header = "");
            
            
           /// Gets Euclidean threshold
           double threshold() const; 
            
        private:
        
            ABC(const ABC& src) {}
            ABC& operator=(const ABC& src) {return *this;}
            
            void init();
            void clear();
            
            double *_sd;
            double *_obs;
            unsigned int _nstats;
            unsigned int _nsam;
            double _threshold;
            std::vector<std::string> _fnames;
            std::vector<unsigned int> _nparams;

    };

}

#endif
