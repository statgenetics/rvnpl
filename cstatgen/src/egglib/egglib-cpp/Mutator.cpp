/*
    Copyright 2009,2010,2012 Stéphane De Mita, Mathieu Siol

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

#include "Mutator.hpp"
#include "Mutation.hpp"
#include "EggException.hpp"
#include <cstdlib>
#include <algorithm>


namespace egglib {

    Mutator::Mutator() {
        init();
    }


    Mutator::~Mutator() {
        clear();
    }


    Mutator::Mutator(const Mutator& src) {
        init();
        copy(src);
    }
        

    Mutator& Mutator::operator=(const Mutator& src) {
        clear();
        init();
        copy(src);
        return *this;
    }
                

    void Mutator::reset() {
        clear();
        init();
    }


    void Mutator::init() {
        _mutationRate = 0.;
        _fixedNumberOfMutations = 0;
        _numberOfSites = 0;
        _sitePositions = NULL;
        _siteWeights = NULL;
        _model = 'F';
        _numberOfAlleles = 0;
        _transitionWeights = NULL;
        _randomAncestralAllele = false;
        numberOfAlleles(2);
        _TPMproba = 0.5;
        _TPMparam = 0.5;
        maxAllele = 0;
        _numberOfMutations = 0;
        _cache_mutations_reserved = 0;
    }
    
    
    void Mutator::clear() {
        if (_sitePositions) free(_sitePositions);
        if (_siteWeights) free(_siteWeights);
        if (_transitionWeights) {
            for (unsigned int i=0; i<_numberOfAlleles; i++) {
                if (_transitionWeights[i]) free(_transitionWeights[i]);
            }
            free(_transitionWeights);
        }
        _cache_mutations.clear();
        _cache_mutations.reserve(0);
    }


    void Mutator::copy(const Mutator& src) {
 
        fixedNumberOfMutations(src.fixedNumberOfMutations());
        mutationRate(src.mutationRate());
        mutationModel(src.mutationModel());
        numberOfAlleles(src.numberOfAlleles());
        for (unsigned int i=0; i<src.numberOfAlleles(); i++) {
            for (unsigned int j=0; j<src.numberOfAlleles(); j++) {
                if (i==j) continue;
                transitionWeight(i, j, src._transitionWeights[i][j]);
            }
        }
        randomAncestralAllele(src.randomAncestralAllele());
        TPMproba(src.TPMproba());
        TPMparam(src.TPMparam());
        numberOfSites(src.numberOfSites());
        for (unsigned int i=0; i<src.numberOfSites(); i++) {
            sitePosition(i, src.sitePosition(i));
            siteWeight(i, src.siteWeight(i));
        }
        
        _numberOfMutations = src.numberOfMutations();

    }

                
                
    unsigned int Mutator::fixedNumberOfMutations() const {
        return _fixedNumberOfMutations;
    }


    unsigned int Mutator::numberOfMutations() const {
        return _numberOfMutations;
    }        


    void Mutator::fixedNumberOfMutations(unsigned int value) {
        // safety check

        if (value != 0 && _mutationRate != 0) throw EggRuntimeError("at least one of the mutation rate and the fixed number of mutation must be zero at all time");
        _fixedNumberOfMutations = value;
    }


    double Mutator::mutationRate() const {
        return _mutationRate;
    }


    void Mutator::mutationRate(double value) {
        if (value<0.) throw EggArgumentValueError("the mutation rate cannot be negative");
        if (value != 0 && _fixedNumberOfMutations != 0) throw EggRuntimeError("at least one of the mutation rate and the fixed number of mutation must be zero at all time");
        _mutationRate = value;
    }
            

    char Mutator::mutationModel() const {
        return _model;
    }


    void Mutator::mutationModel(char value) {
        if (value!='F' && value!= 'I' && value!= 'S' && value!= 'T') throw EggArgumentValueError("unknown mutation model");
        _model = value;
    }
                

    unsigned int Mutator::numberOfAlleles() const {
        return _numberOfAlleles;
    }
               

    void Mutator::numberOfAlleles(unsigned int value) {
        if (value<2) throw EggArgumentValueError("the fixed number of alleles must be at least 2");

        // frees the previous matrix
        if (_numberOfAlleles) {
            for (unsigned int i=0; i<_numberOfAlleles; i++) {
                free(_transitionWeights[i]);
            }
            free(_transitionWeights);
        }
        else if (_transitionWeights) throw EggRuntimeError("memory management error (1) in Mutator");
        
        // allocates (and pre-fills with 1.'s) the new matrix
        _transitionWeights = (double**) malloc(value * sizeof(double*));
        if (!_transitionWeights) throw EggMemoryError();
        for (unsigned int i=0; i<value; i++) {
            _transitionWeights[i] = (double*) malloc(value * sizeof(double));
            if (!_transitionWeights[i]) throw EggMemoryError();
            for (unsigned int j=0; j<value; j++) {
                if (i==j) _transitionWeights[i][j] = -1.;
                else _transitionWeights[i][j] = 1.;
            }
        }

        // stores the new value
        _numberOfAlleles = value;
    }           
                

    double Mutator::transitionWeight(unsigned int i, unsigned int j) {
        if (i>=_numberOfAlleles || j>=_numberOfAlleles || i==j) throw EggArgumentValueError("error while accessing transition weight - invalid index");
        return _transitionWeights[i][j];
    }


    void Mutator::transitionWeight(unsigned int i, unsigned int j, double value) {
        if (i>=_numberOfAlleles || j>=_numberOfAlleles || i==j) throw EggArgumentValueError("error while setting transition weight - invalid index");
        if (value<=0.) throw EggArgumentValueError("error while setting transition weight - invalid weight value");
        _transitionWeights[i][j] = value;
    }


    bool Mutator::randomAncestralAllele() const {
        return _randomAncestralAllele;
    }
    
    
    void Mutator::randomAncestralAllele(bool flag) {
        _randomAncestralAllele = flag;
    }


    unsigned int Mutator::numberOfSites() const {
        return _numberOfSites;
    }
               

    void Mutator::numberOfSites(unsigned int value) {
        // if numberOfSites remains 0, nothing to do (_numberOfSites is already = value)
        if (_numberOfSites==0 && value==0);

        // if 0, only needs to free the tables
        if (value==0) {
            if (!_sitePositions || !_siteWeights) {
                if (!_sitePositions && !_siteWeights) {
                    free(_sitePositions);
                    free(_siteWeights);
                }
                else {  // one of the table is already free
                    throw EggRuntimeError("memory management (2) error in Mutator");
                }
            }
            // if both tables are free beforehand, do nothing
        }
        
        // then we need to allocate and fill the arrays
        else {
            
            // memory allocation
            _sitePositions = (double*) realloc(_sitePositions, value*sizeof(double));
            if (!_sitePositions) {
                throw EggMemoryError();
            }
            _siteWeights = (double*) realloc(_siteWeights, value*sizeof(double));
            if (!_siteWeights)  {
                throw EggMemoryError();
            }
            
            // filling - edge case
            if (value==1) {
                _sitePositions[0] = 0.5;
                _siteWeights[0] = 1.;
            }
            
            // filling - general case
            else {
                _sitePositions[0] = 0.;
                _siteWeights[0] = 1.;
                double incr = 1./(value-1);
                double cur = incr;
                for (unsigned int c=1; c<value; c++) {
                    _sitePositions[c]=cur;
                    _siteWeights[c]=1.;
                    cur+=incr;
                }
            }
        }
        
        // this is done in the end
        _numberOfSites = value;
    }
                

    double Mutator::sitePosition(unsigned int siteIndex) const {
        if (siteIndex>=_numberOfSites) throw EggArgumentValueError("invalid mutable site index");
        return _sitePositions[siteIndex];
    }


    void Mutator::sitePosition(unsigned int siteIndex, double value) {
        if (siteIndex>=_numberOfSites) throw EggArgumentValueError("invalid mutable site index");
        if (value<0. || value>1.) throw EggArgumentValueError("site positions must be between 0 and 1");
        _sitePositions[siteIndex] = value;
    }


    double Mutator::siteWeight(unsigned int siteIndex) const {
        if (siteIndex>=_numberOfSites) throw EggArgumentValueError("invalid mutable site index");
        return _siteWeights[siteIndex];
    }


    void Mutator::siteWeight(unsigned int siteIndex, double value) {
        if (siteIndex>=_numberOfSites) throw EggArgumentValueError("invalid mutable site index");
        if (value<=0.) throw EggArgumentValueError("site weights must be strictly positive");
        _siteWeights[siteIndex] = value;
    }


    /******************************************************************/

    DataMatrix Mutator::mute(Arg* arg, Random* random) {

        double slen = 1./arg->numberOfSegments;

        // computes the tree length
        double L = 0;
        double* Li = NULL; // the array of tree lengths (at each of the sites)
        unsigned int* segmentIndices=NULL;
        double sumWeights=0.;   // the last  variable supposed to be used only with finite site model 
        
        if (_numberOfSites==0) {
            L = arg->totalLength / arg->numberOfSegments;
        }
        else {
            Li = (double*) malloc(_numberOfSites*sizeof(double));
            if (!Li) throw EggMemoryError();
            segmentIndices = (unsigned int*) malloc(_numberOfSites*sizeof(unsigned int));
            if (!segmentIndices) {
                if (Li) free(Li);
                throw EggMemoryError();
            }
            unsigned int segmentIndex=0;
            for (unsigned int i=0; i<_numberOfSites; i++) {
             
                while(1) {

                    if ((_sitePositions[i]-slen*(segmentIndex+1))<0.000000000001) {
                        break;
                    }

                    segmentIndex++;

                    if (segmentIndex>=arg->numberOfSegments) {
                        if (Li) free(Li);
                        if (segmentIndices) free(segmentIndices);
                        throw EggRuntimeError("cannot find mutable site's segment!");
                    }
                }
                
                segmentIndices[i] = segmentIndex;
                
                Li[i] = arg->segmentLengths[segmentIndex] * _siteWeights[i];
                L+=Li[i];
                sumWeights+=_siteWeights[i];
            }
            if (!sumWeights) {
                if (Li) free(Li);
                if (segmentIndices) free(segmentIndices);
                throw EggRuntimeError("sum of site weights is zero: cannot perform mutation");
            }
            L/=sumWeights;
        }

        // determines the number of mutations
                
        if (L==0 || (!_fixedNumberOfMutations && !_mutationRate)) {
            _numberOfMutations = 0;
        }
        else {
            if (_fixedNumberOfMutations) {
                _numberOfMutations = _fixedNumberOfMutations;
            }
            else {
                if (L==0.) {
                    _numberOfMutations = 0;
                }
                else {
                    _numberOfMutations = random->prand(_mutationRate * L);
                }
            }
        }


        unsigned int actualNumberOfSites;
        
        if (_numberOfSites==0) {
            actualNumberOfSites = _numberOfMutations;
            //segmentIndices = (unsigned int*) malloc(actualNumberOfSites*sizeof(unsigned int));
            //if (!segmentIndices) {
            //    if (Li) free(Li);
            //    throw EggMemoryError();
            //}
        }
        else { 
            actualNumberOfSites = _numberOfSites;
            // segmentIndices already allocated and set
        }

        arg->set_actualNumberOfSites(actualNumberOfSites);

        //double *actualSitePositions = (double*) malloc(actualNumberOfSites*sizeof(double));
        //if (!actualSitePositions) throw EggMemoryError();
        
        //if (_numberOfSites) {
        //    for (unsigned int i=0; i<actualNumberOfSites; i++) {
        //        actualSitePositions[i] = _sitePositions[i];
        //    }
        //}

        // increases the Mutation table if needed
        
        _cache_mutations.resize(_numberOfMutations);
        
        if (_numberOfMutations > _cache_mutations_reserved) {
            _cache_mutations.reserve(_numberOfMutations);
            _cache_mutations_reserved = _numberOfMutations;
        }

        // creates the mutations
        
//        std::map<double, std::vector<Mutation> > mutation_mapping;
            
        for (unsigned int i=0; i<_numberOfMutations; i++) {
            
            unsigned int actualSiteIndex;
            double treeLengthAtThatPosition;
            unsigned int segmentIndex;
            
            // determines which site is muted
            if (!_numberOfSites) {

                actualSiteIndex = i;

                // picks a segment
                double X = random->uniform() * L * arg->numberOfSegments;
                double acc=0;
                bool flag=false;

                for (unsigned int j=0; j<arg->numberOfSegments; j++) {

                    acc+=arg->segmentLengths[j];
                    if (acc>=X) {

                        flag  = true;
                        segmentIndex = j;
                        //segmentIndices[i] = j;
                        _cache_mutations[i].segmentIndex = j;
                        
                        treeLengthAtThatPosition = arg->segmentLengths[j];
                        _cache_mutations[i].position = (double)j/arg->numberOfSegments + 
                                (X-(acc-arg->segmentLengths[j])) /
                                (arg->segmentLengths[j]*arg->numberOfSegments);
                        break;
                    }
                }
                
                if (!flag) {
                    if (Li) free(Li);
                    if (segmentIndices) free(segmentIndices);
                    throw EggRuntimeError("could not draw site (the computed tree length was incorrect)");
                }

            }
            else {
                
                double X = random->uniform() * L  * sumWeights;
                double acc = 0.;
                bool flag = false;  // used for safety check
                
                for (unsigned int j=0; j<_numberOfSites; j++) {
                    acc += Li[j];
                    if (acc>=X) {
                        flag = true;
                        actualSiteIndex = j;
                        segmentIndex = segmentIndices[j];
                        treeLengthAtThatPosition = Li[j] / _siteWeights[j];
                        _cache_mutations[i].segmentIndex = segmentIndices[j];
                        _cache_mutations[i].position = _sitePositions[j];
                        break;
                    }
                }
                
                if (!flag) {
                    if (Li) free(Li);
                    if (segmentIndices) free(segmentIndices);
                    throw EggRuntimeError("could not draw site (the computed weighted tree length was incorrect)");
                }
            }
            
            // determines where (which branch) the mutation occurs

            double X = random->uniform() * treeLengthAtThatPosition;

            Edge *mutedEdge = arg->mute(segmentIndex, X);
            mutedEdge->numberOfMutationsPerActualSite[actualSiteIndex]++;
            _cache_mutations[i].actualSiteIndex = actualSiteIndex;


            // adds the mutation to the holder
            //Mutation mutation = arg->mute(segmentIndex, X);
            //mutation.position = mutedPosition;
            //std::vector<Mutation> empty;
            //if (mutation_mapping.count(mutedPosition)==0) {
            //    mutation_mapping[mutedPosition] = empty; // deep-copied
            //}
            //mutation_mapping[mutedPosition].push_back(mutation);
        }

        // position sorting is performed automatically when using the multimap

        // sorting operations are unstable (as of now, I think stability is not an issue)

        std::sort(_cache_mutations.begin(), _cache_mutations.end(), compare);

        // takes mutations from left to right

        DataMatrix data(arg->numberOfSamples, actualNumberOfSites);

        unsigned int siteIndex;
        unsigned int segmentIndex;


        for (unsigned int i=0; i<actualNumberOfSites; i++) {

            // sets the position
            if (_numberOfSites) data.sitePosition(i, _sitePositions[i]);
            else  data.sitePosition(i, _cache_mutations[i].position);


            // defines the ancestral state
            int allele = 0;
            if (_randomAncestralAllele && _model=='F') {
                    allele= (int) random->irand(_numberOfAlleles);
            }

            // gets mutation location
            if (_numberOfSites) {
                siteIndex = i;
                segmentIndex = segmentIndices[i];
            }
            else {
                siteIndex = _cache_mutations[i].actualSiteIndex;
                segmentIndex = _cache_mutations[i].segmentIndex;
            }

            maxAllele = 0;
            apply_mutation(i, siteIndex, data, arg->MRCA(segmentIndex),
                                        allele, segmentIndex, random);
           
        }


        //unsigned int mutation_index = 0;
        //std::map<double, std::vector<Mutation> >::iterator it;



   /*     for (it = mutation_mapping.begin(); it != mutation_mapping.end(); ++it) {
           
            // sort mutations at that position (descendingly in age)

            std::vector<Mutation>& mutations = it->second;   // only a reference

            if (_numberOfMutations) {  // not needed if infinitely many sites
                std::sort(mutations.begin(), mutations.end(), compare);
            }

            // defines the ancestral state
            int allele = 0;
            if (_randomAncestralAllele && _model=='F') {
                    allele= (int) random->irand(_numberOfAlleles);
            }
            
            // operates mutations
            maxAllele=0;
            for (std::vector<Mutation>::iterator mut=mutations.begin(); mut!=mutations.end(); ++mut) {
                allele = nextAllele(allele, random); 
                apply_mutation(data, mut->edge, mutation_index, mut->segment, allele);
                data.sitePosition(mutation_index, it->first);  // sets the position
            }
            mutation_index++;
        }
    */

        if (Li) free(Li);
        if (segmentIndices) free(segmentIndices);

        // returns the datamatrix
        return data;
    }

    /******************************************************************/

    void Mutator::apply_mutation(unsigned int matrixIndex,
                                 unsigned int actualSite, DataMatrix& data,
                                 const Edge* edge, int allele, 
                                 unsigned int segment, Random* random) {

        // changes allele as needed
        
        if (edge->numberOfMutationsPerActualSite[actualSite] > 0) {
                // warning: not all changes implemented (only flat matrices - with all equal rates - allowed)
            
            switch (_model) {
                // fixed number of alleles model
                case 'F':
                    {
                        for (unsigned int n=0; n<edge->numberOfMutationsPerActualSite[actualSite]; n++) {
                            
                            bool flag=false;
                            double sumWeights = 0.;
                            for (unsigned int i=0; i<_numberOfAlleles; i++) {
                                if ((int)i!=allele) sumWeights+= _transitionWeights[allele][i];
                            }
                            double X = random->uniform() * sumWeights;
                            double c=0.;
                            for (unsigned int i=0; i<_numberOfAlleles; i++) {
                                if ((int)i==allele) continue;
                                c+=_transitionWeights[allele][i];
                                if (c>X) {
                                    flag = true;
                                    allele = i;
                                    break;
                                }
                            }
                            if (!flag) throw EggRuntimeError("hole in Mutator::apply_mutation( )");
                        }
                        break;
                    }

                // infinite allele model (only once)
                case 'I':
                
                    maxAllele++;   // pre-increment to avoid putting the initial allele


                    allele= maxAllele; 
                    break;
                
                // stepwise mutation model (incrementation performed in apply_mutation)
                case 'S':
                    for (unsigned int n=0; n<edge->numberOfMutationsPerActualSite[actualSite]; n++) {
                        allele += TPMstep(0., random);
                    }
                    break;

                // two-phase mutation model (incrementation performed in apply_mutation)
                case 'T':
                    for (unsigned int n=0; n<edge->numberOfMutationsPerActualSite[actualSite]; n++) {
                        allele += TPMstep(_TPMproba, random);
                    }
                    break;

                // this is only destined to balance the code
                default:
                    throw EggRuntimeError("unknown mutation model code");        
            }
        }
        
        // descends or applies allele


        switch (edge->numberOfSons) {
            
            case 0:
                data.set(edge->label(), matrixIndex, allele);
                break;
                
            
            case 2:
                if (edge->son2->segment(segment)) {
                    apply_mutation(matrixIndex, actualSite, data, edge->son2, allele, segment, random);
                }

            case 1:
                if (edge->son1->segment(segment)) {
                    apply_mutation(matrixIndex, actualSite, data, edge->son1, allele, segment, random);
                }
                break;
            
            default:
                throw EggRuntimeError("Edge with more than 2 sons!");
            
        }

    }

    double Mutator::TPMproba() const {
        return _TPMproba;
    }


    double Mutator::TPMparam() const {
        return _TPMparam;
    }


    void Mutator::TPMproba(double value) {
        if (value<0. || value>1.) throw EggArgumentValueError("value for TPM probability parameter out of range");
        _TPMproba = value;
    }


    void Mutator::TPMparam(double value) {
        if (value<0. || value>1.) throw EggArgumentValueError("value for TPM distribution parameter out of range");
        _TPMparam = value;
    }


    int Mutator::TPMstep(double inTPMproba, Random* random) {
        unsigned int step = 1; // by default: SMM
        if (inTPMproba!=0. && random->uniform()<inTPMproba) {
            step = random->grand(_TPMparam);
        }
        int sign = (random->uniform()>0.5)?1:-1;
        return sign * step;
    }


  /*  int Mutator::nextAllele(int allele, Random* random) {
        switch (_model) {
            // fixed number of alleles model
            case 'F':
                {
                    double sumWeights = 0.;
                    for (unsigned int i=0; i<_numberOfAlleles; i++) {
                        if ((int)i!=allele) sumWeights+= _transitionWeights[allele][i];
                    }
                    double X = random->uniform() * sumWeights;
                    double c=0.;
                    for (unsigned int i=0; i<_numberOfAlleles; i++) {
                        if ((int)i==allele) continue;
                        c+=_transitionWeights[allele][i];
                        if (c>X) return i;
                    }
                    throw EggRuntimeError("hole in Mutator::nextAllele( )");
                }

            // infinite allele model
            case 'I': return ++maxAllele;  // pre-increment to avoid putting the initial allele
            
            // stepwise mutation model (incrementation performed in apply_mutation)
            case 'S': return / *allele+* /TPMstep(0., random);

            // two-phase mutation model (incrementation performed in apply_mutation)
            case 'T': return / *allele+* /TPMstep(_TPMproba, random);

            // this is only destined to balance the code
            default:
                throw EggRuntimeError("unknown mutation model code");
        }
        return 0;
    }  */



    bool compare(Mutation mutation1, Mutation mutation2) {
        return mutation1.segmentIndex < mutation2.segmentIndex;
    }

}
