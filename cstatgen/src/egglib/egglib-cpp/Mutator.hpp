/*
    Copyright 2009, 2010, 2012 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_MUTATOR_HPP
#define EGGLIB_MUTATOR_HPP


#include "DataMatrix.hpp"
#include "Random.hpp"
#include "Arg.hpp"
#include "Mutation.hpp"


namespace egglib {
    

   /** \brief Implements mutation models
    *
    * \ingroup coalesce
    * 
    * Works with a previously built Ancestral Reconbination Graph. The
    * user must sets options using the setter-based interface. After
    * that he or she can call the method mute() that will generates
    * a DataMatrix object.
    * 
    * Genotype data are represented by integer numbers. Regardless of
    * the mutation model, the ancestral state is always 0. The user can
    * set the rate of mutation (or, alternatively, fix the number of
    * mutations that occurred - which is the number of segregating sites
    * only with an infinite site model).
    * 
    * Other options fall into two separate groups: the positions of the
    * mutated sites and the process of mutation (how new alleles are
    * generated).
    * 
    * Concerning allele generation, several mutation models are available
    * (coded by single letters):
    *   - F: fixed number of alleles. Among other markers, this model is
    *        appropriate for simulating nucleotides. The user is able
    *        to choose the number of alleles (where 2 is the standard
    *        for an infinite site model and 4 for a finite site model).
    *        Mutator allows assigning independent weights between all
    *        different transition types and can draw randomly the
    *        ancestral states, providing a way to emulate evolution of
    *        nucleotides with multiple mutations at the same site and
    *        reversion.
    *   - I: infinite number of alleles. At a given site, each mutation
    *        raises a new allele. The value of the alleles is therefore
    *        irrelevant (it only denotes its order of appearance). This
    *        model does not permit homoplasy.
    *   - S: stepwise mutation model. In this model the value of the
    *        alleles are interpreted as a size (typically for simulating
    *        a microsatellite marker). Each mutation either increases
    *        or decreases the allele size by an increment of one.
    *   - T: two-phase mutation model. This model is a generalization
    *        of the stepwise mutation model (S). For a mutation, the
    *        increment (either increase or decrease) is 1 with the
    *        probability given by the parameter (1-TPMproba). With
    *        probability TPMproba, the increment is drawn from a
    *        geometric distribution of parameter given by the other
    *        parameter (TPMparam).
    * 
    * By default, the program will assume an infinite site model (ISM).
    * Each mutation will occur to a new position drawn from the interval
    * [0,1]. It is possible to set any mutation model with an ISM 
    * (including microsatellite-like models I, S and T). Alternatively,
    * the user can specify a finite number of sites available for
    * mutation. For a microsatellite marker, the user will want to
    * specify a single site. It is possible to set a finite number of
    * sites for all mutation models. In all cases, the mutations will
    * be forced to target these sites. It is possible to apply weights
    * independently to all sites. The higher the weight value
    * (comparatively to the other sites), the higher the probability
    * that this site mutes. The weights needs not to be relative. In
    * addition, the user can set the positions of the different sites.
    * Nothings forces him or her to place them in order. Note that this
    * does not affect the mutation process, but on the amount of
    * recombination that will be allowed between sites.
    *
    */
    class Mutator {

        public:
        
           /** \brief Initializes with default values
            * 
            * List of default values:
            *   - theta = 0
            *   - fixedNumberOfMutations = 0
            *   - model = F (fixed number of alleles)
            *   - fixed number of alleles = 2
            *   - infinite site model
            *   - TPM parameters are both preset to 0.5
            * 
            */
            Mutator();

            
           /** \brief Destroys the object
            * 
            */
            ~Mutator();
            
            
           /** \brief Copy constructor
            * 
            */
            Mutator(const Mutator&);
        
        
           /** \brief Assignement operator
            * 
            */
            Mutator& operator=(const Mutator&);


           /** \brief Restores default values of all parameters
            * 
            */
            void reset();
            
            
           /** \brief Gets the fixed number of mutations
            * 
            */
            unsigned int fixedNumberOfMutations() const;
        

           /** \brief Sets the fixed number of mutations
            * 
            * The value can be 0. It is not allowed to set both the
            * fixed number of mutations and the mutation rate at
            * non-zero value
            * 
            */
            void fixedNumberOfMutations(unsigned int);


           /** \brief Gets the mutation rate
            * 
            */
            double mutationRate() const;
        
        
           /** \brief Sets the mutation rate
            * 
            * The value cannot be negative. The value can be 0.  It is
            * not allowed to set both the fixed number of mutations and
            * the mutation rate at non-zero value
            * 
            */
            void mutationRate(double);


           /** \brief Gets the mutation model
            * 
            * See the class documentation for the signification of the
            * different one-letter codes.
            * 
            */
            char mutationModel() const;
            
            
           /** \brief Sets the mutation model
            * 
            * The passed character must be one of F, I, S and T. See the
            * class documentation for their signification.
            * 
            */
            void mutationModel(char);
            
            
           /** \brief Gets the fixed number of possible alleles
            * 
            */
            unsigned int numberOfAlleles() const;
           
            
           /** \brief Sets the fixed number of possible alleles
            * 
            * The value must be larger than 1. This parameter is
            * meaningful only for the fixed number allele model of
            * mutation, and ignored otherwise.
            *
            */
            void numberOfAlleles(unsigned int);
           
           
           /** \brief Sets a transition weight
            * 
            * \param i row (previous allele index).
            * \param j column (next allele index).
            * \param value weight to apply.
            * 
            * Indices i and j must be different. Weights can be any
            * strictly positive value.
            * 
            */
            void transitionWeight(unsigned int i, unsigned int j, double value);


           /** \brief Gets a transition weight
            * 
            * \param i row (previous allele index).
            * \param j column (next allele index).
            * 
            * Indices i and j must be different.
            * 
            */
            double transitionWeight(unsigned int i, unsigned int j);
            
            
           /** \brief Set to true to draw ancestral alleles randomly
            * 
            * By default, the ancestral allele is always 0. If this
            * variable is set to true, the ancestral allele will be
            * randomly drawn from the defined number of alleles. This
            * option is always ignored unless in combination with the
            * Fixed Allele Number model.
            * 
            */
            void randomAncestralAllele(bool flag);

            
            /** \brief true if ancestral alleles must be drawn randomly
             * 
             */
            bool randomAncestralAllele() const;

           
           /** \brief Gets the TPM probability parameter
            * 
            */
            double TPMproba() const;
            
            
           /** \brief Sets the TPM probability parameter
            * 
            * This parameter is considered only if the mutation model
            * is T (two-phase mutation model). It gives the probability
            * that a mutation step is not fixed to be 1. If TPMproba is
            * 0, the mutation model is SMM.
            * 
            * The value must be >=0. and <=1. 
            * 
            */
            void TPMproba(double value);
           
            
           /** \brief Gets the TPM distribution parameter
            * 
            */
            double TPMparam() const;
            
            
           /** \brief Sets the TPM distribution parameter
            * 
            * This parameter is considered only if the mutation model
            * is T (two-phase mutation model). It gives the parameter
            * of the geometric distribution which is used to generate
            * the mutation step (if it is not one).
            * 
            * The value must be >=0. and <=1. 
            * 
            */
            void TPMparam(double value);


           /** \brief Gets the number of mutable sites
            * 
            * A value a zero must be interpreted as the infinite site
            * model. Note that after all calls all data from the tables
            * sitePositions and siteWeights will be reset.
            * 
            */
            unsigned int numberOfSites() const;
           
            
           /** \brief Sets the number of mutable sites
            * 
            * The value of zero is accepted and imposed the infinite
            * site model.
            * 
            */
            void numberOfSites(unsigned int);
            
            
           /** \brief Gets the position of a given site
            * 
            */
            double sitePosition(unsigned int siteIndex) const;

            
           /** \brief Set the position of a given site
            * 
            * The position must be >=0 and <=1
            * 
            */
            void sitePosition(unsigned int siteIndex, double position);


           /** \brief Gets the mutation weight of a given site
            * 
            */
            double siteWeight(unsigned int siteIndex) const;

            
           /** \brief Set the site weight of a given site
            * 
            * The weight must be strictly positive.
            * 
            */
            void siteWeight(unsigned int siteIndex, double weight);


           /** \brief Performs mutation
            * 
            * \param arg Ancestral recombination graph instance. If the
            * ARG is partially built or not a all, or improperly so,
            * the behaviour of this method is not defined.
            * 
            * \param random The address of a Random instance to be
            * used for generating random numbers.
            * 
            * \return A DataMatrix instance containing simulated data.
            * 
            */
            DataMatrix mute(Arg* arg, Random* random);


          /** \brief Gets the last number of mutations
            *
            * Returns the number of mutations of the last call of mute( ).
            * By default, this method returns 0.
            *
            */
            unsigned int numberOfMutations() const; 


        private:
        
            void clear();
            void init();
            void copy(const Mutator&);

            //int nextAllele(int allele, Random* random);
            int TPMstep(double inTPMproba, Random* random);
            void apply_mutation(unsigned int matrixIndex,
                                unsigned int actualSite, DataMatrix& data,
                                const Edge* edge, int allele,
                                unsigned int segment, Random* random);

        
            char _model;
            double _mutationRate;
            unsigned int _fixedNumberOfMutations;
            unsigned int _numberOfAlleles;
            double** _transitionWeights;
            bool _randomAncestralAllele;
            unsigned int _numberOfSites;
            double* _sitePositions;
            double* _siteWeights;
            double _TPMproba;
            double _TPMparam;
            int maxAllele;
            unsigned int _numberOfMutations;
            std::vector<Mutation> _cache_mutations;
            unsigned int _cache_mutations_reserved;

    };


    bool compare(Mutation mutation1, Mutation mutation2); // returns True if mutation1 is older

}




#endif

