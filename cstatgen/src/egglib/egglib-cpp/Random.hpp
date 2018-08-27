/*
    Copyright 2008,2009,2012 Stéphane De Mita, Mathieu Siol
    Adapted from MStrat, developed by Charles-Edouard Coste,
    Thomas M. Bataillon, Mathieu Cotisson, Guy Decoux, Chistophe Rozale,
    Daniel J. Schoen and Jacques L. David.
    
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

#ifndef EGGLIB_RANDOM_HPP
#define EGGLIB_RANDOM_HPP

namespace egglib {

    /** \brief Pseudo-random number generator
     *
     * \ingroup core
     *
     * Random is a pseudo-random number generator, adapted from a part of MStrat,
     * developed by Charles-Edouard Coste, Thomas M. Bataillon, Mathieu Cotisson,
     * Guy Decoux, Chistophe Rozale, Daniel J. Schoen and Jacques L. David.
     *
     * It uses two different seeds. By default, they are initialized to available
     * arbitrary values. However, a given sequence can be repeated by passing the
     * same two seeds.
     * 
     */
    class Random {
       public:
            /** \brief Initializes using default seeds
             *
             * Uses the current system time and the memory address of the object as an attempt to generate unique sequences.
             */
            Random();

            /** \brief Initializes using given seeds
             *
             * This constructor can be used to reproduce a given sequence.
             */
            Random(double seed1, double seed2);
            
           /** \brief Draws a number from an exponential distribution
            * 
            * \param expectation the distribution mean (also 1/lambda
            * where lambda is the rate parameter).
            *
            */
            double erand(double expectation);

           /** \brief Draws an integer from a uniform distribution bound by 0 and max (max is not included)
            * 
            * max is not included.
            * 
            */
            unsigned int irand(unsigned int max);

            /** \brief Draws an integer from a Poisson distribution with parameter p
             *
             * The Poisson transformation algorithm was taken from (in French)
             * http://www.u-picardie.fr/~cochard/IEM/demos/C107/C107_3.htm.
             */
            unsigned int prand(double p);

           /** \brief Draws a number from a normal distribution of expectation 0 and variance 1
            * 
            * The algorithm used is the polar form of the Box-Muller
            * algorithm. \todo use the Ziggurat algorithm for the
            * nrand() method of Random.
            * 
            */
            double nrand();
            
           /** \brief Draws a number from a geometric law
            * 
            * \param param the parameter of the law
            * 
            */
            unsigned int grand(double);

           /** \brief Draws a number from a uniform distribution between 0 and 1
            * 
            */
            double uniform();
            
           /** \brief Gets the current value of the first seed
            * 
            */
            double seed1() const;

           /** \brief Gets the current value of the second seed
            * 
            */
            double seed2() const;
            
           /** \brief Sets the current value of the first seed
            * 
            */
            void seed1(double);

           /** \brief Sets the current value of the second seed
            * 
            */
            void seed2(double);

        private:
            // First seed
            double _seed1;
            
            // Second seed
            double _seed2;
            
            /* since the normal random generator draws two numbers at
             * a time, one is cached and returned at any subsequent call
             */
            bool b_ncached;
            double v_ncached;
            
    };
}

#endif
