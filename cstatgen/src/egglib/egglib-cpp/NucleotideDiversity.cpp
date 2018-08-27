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


#include "NucleotideDiversity.hpp"
#include "EggException.hpp"
#include <sstream>
#include <cstdlib>
#include <cmath>


namespace egglib {
  
    
//***** Constructors ***************************************************
    NucleotideDiversity::NucleotideDiversity() {
        init();
    }


    NucleotideDiversity::~NucleotideDiversity() {
        clear();
    }
        
        
//***** Helpers ********************************************************
    void NucleotideDiversity::init() {
        
        BaseDiversity::init();

        b_analysisSites = false;

        b_diversity = false;
        v_Pi = 0.;
        v_thetaW = 0.;
        v_average_Pi = 0.;
        v_pop_Pi = NULL;
        v_D = 0.;

        b_outgroupDiversity = false;
        v_thetaH = 0.;
        v_thetaL = 0.;
        v_H = 0.;
        v_Z = 0.;
        v_E = 0.;
            
        b_differentiation = false;
            
        v_pairwiseFixedDifferences = NULL;
        v_pairwiseCommonAlleles = NULL;
        v_pairwiseSharedAlleles = NULL;
        v_popPolymorphic = NULL;
        v_popSpecific = NULL;
        v_popSpecificDerived = NULL;
        v_countFixedDifferences = 0;
        v_countCommonAlleles = 0;
        v_countSharedAlleles = 0;
        v_countSpecificAlleles = 0;
        v_countSpecificDerivedAlleles = 0;

        b_triConfigurations = false;
        
        v_triConfigurations = NULL;

    }
        
        
    void NucleotideDiversity::clear() {

        BaseDiversity::clear();
        
        if (b_diversity) {
            free(v_pop_Pi);
        }
        
        if (b_differentiation) {

            free(v_pairwiseFixedDifferences);
            free(v_pairwiseCommonAlleles);
            free(v_pairwiseSharedAlleles);
            free(v_popPolymorphic);
            free(v_popSpecific);
            free(v_popSpecificDerived);
        }
        
        if (b_triConfigurations) {
            free(v_triConfigurations);
        }

    }
        




//***** FIRST STEP OF ANALYSIS *****************************************

    void NucleotideDiversity::load(
        CharMatrix& data,
        bool allowMultipleMutations,
        double minimumExploitableData,
        unsigned int ignoreFrequency,
        std::string characterMapping,
        bool useZeroAsAncestral)
    {
    
        clear();
        init();
        
        importSites(data, allowMultipleMutations, minimumExploitableData, ignoreFrequency,
                        characterMapping, useZeroAsAncestral, false);

        b_analysisSites = true;
    }
    
    
    
//***** SECOND STEP OF ANALYSIS ****************************************

    void NucleotideDiversity::diversity() {

        if (!b_analysisSites) throw EggRuntimeError("Diversity: data not available (data not loaded)");
        if (b_diversity) throw EggRuntimeError("Diversity: cannot re-compute data (implementation error)");
        
        // prepares and initializes the population Pi vector
        
        
        v_pop_Pi = (double*) realloc(v_pop_Pi, v_npop*sizeof(double));
        if (!v_pop_Pi) throw EggMemoryError();
        for (unsigned int i=0; i<v_npop; i++) {
            v_pop_Pi[i] = 0.;
        }

        // escapes immediatly if not Diversity
        b_diversity=true; 
        if (!v_S) return;

        // safety check
        if (v_npop==0) {
            throw EggRuntimeError("error while computing diversity: data but no pop");
        }

        // sums diversity
        for (unsigned int site=0; site<v_S; site++) {

            unsigned int ns = v_sites[site]->ns();

            // total diversity
            double sumSquareFreqs = 0.;
            for (unsigned int allele=0; allele<v_sites[site]->numberOfAlleles(); allele++) {
                double p = 1. * v_sites[site]->alleleFrequency(allele) / ns;
                sumSquareFreqs+= p*p ;
            }
            double pi = (1-sumSquareFreqs) * ns / (ns-1);
            v_Pi += pi;
            
            // per pop diversities
            for (unsigned int pop=0; pop<v_npop; pop++) {

                unsigned int ns_pop = v_sites[site]->ns(pop);

                sumSquareFreqs = 0.;
                for (unsigned int allele=0; allele<v_sites[site]->numberOfAlleles(); allele++) {
                    double p = 1.* v_sites[site]->alleleFrequency(pop,
                                                        allele) / ns_pop;

                    sumSquareFreqs+= p*p;
                }
                v_pop_Pi[pop] += (1-sumSquareFreqs) * ns_pop / (ns_pop-1);
            }
        }
        
        // average diversity (over pops)
        for (unsigned int pop=0; pop<v_npop; pop++) {
            v_average_Pi += v_pop_Pi[pop];
        }
        v_average_Pi/=v_npop;
        
        // Tajima's D
        double a1= 0., a2= 0.;
        for (int i=1; i<round(v_nseff); i++) {
            a1+= 1./i;
            a2+= 1./(i*i);
        }
        double b1= (v_nseff+1)/(3.*(v_nseff-1));
        double b2= 2.*(v_nseff*v_nseff+v_nseff+3)/(9.*v_nseff*(v_nseff-1));
        double c1= b1 - 1./a1;
        double c2= b2 - (v_nseff+2)/(a1*v_nseff) + a2/(a1*a1);
        double e1= c1/a1;
        double e2= c2/(a1*a1 + a2);
        
        v_thetaW = v_S/a1;

        double V= e1*v_S + e2*v_S*(v_S-1);
        V= (V<0.00000000001)?0:V;
        v_D = (v_Pi-v_thetaW)/sqrt(V);

        // express statistics per site

        if (v_lseff) {
            v_Pi /= v_lseff;
            v_thetaW /= v_lseff;
            for (unsigned int pop=0; pop<v_npop; pop++) {
                v_pop_Pi[pop] /= v_lseff;
            }
            v_average_Pi /= v_lseff;
        }

    }
    

//***** THIRD STEP OF ANALYSIS *****************************************

    void NucleotideDiversity::outgroupDiversity() {
        
        if (!b_analysisSites) throw EggRuntimeError("Diversity: data not available (data not loaded)");
        if (!b_diversity) diversity();
        if (b_outgroupDiversity) throw EggRuntimeError("Diversity: cannot re-compute data (implementation error)");
        
        // escapes immediatly if no diversity
        b_outgroupDiversity=true; 
        if (!v_S) return;

        unsigned int *esse;
        esse= (unsigned int*) malloc(v_ns*sizeof(unsigned int));
        if (v_ns && !esse) throw EggMemoryError();

        for (unsigned int i=0; i<v_ns; i++) {
            esse[i]= 0;
        }

        for (unsigned int i=0; i<v_S; i++) {
            if (v_orientables[i]) {
                esse[v_sites[i]->derivedAlleleFrequency()]++;
            }
        }

        v_thetaH = 0.;
        for (unsigned int i=1; i<v_ns; i++) {
            v_thetaH += i*i*esse[i];
        }
        v_thetaH *= 2./(v_ns*(v_ns-1));
        v_H = v_Pi*v_lseff - v_thetaH;

        v_thetaL = 0.;
        for (unsigned int i=1; i<v_ns; i++) v_thetaL += i*esse[i];
        v_thetaL/= (v_ns-1);

        double a1= 0.;
        for (int i=1; i<round(v_nseffo); i++) a1+= 1./i;
        
        double bn= 0.;
        for (int i=1; i<round(v_nseffo); i++) bn+= 1./(i*i);
        
        double bnp1= 0.;
        bnp1 = bn + 1./(round(v_nseffo)*round(v_nseffo));

        double theta= v_S / a1;
        double theta2=  v_S*(v_S-1.)/(a1*a1+bn);

        double VarZ = theta*(v_nseffo-2)/(6*(v_nseffo-1)) +
                theta2*(18*v_nseffo*v_nseffo*(3*v_nseffo+2)*bnp1-(88*v_nseffo*v_nseffo*v_nseffo+9*v_nseffo*v_nseffo-13*v_nseffo+6))/
                    (9*v_nseffo*(v_nseffo-1)*(v_nseffo-1));

        v_Z = (v_Pi*v_lseff - v_thetaL) / sqrt(VarZ);

        double VarE = theta *(v_nseffo/(2*(v_nseffo-1))-1/a1)
             +  theta2*(
                   bn/(a1*a1) + 2*bn*(v_nseffo/(v_nseffo-1))*(v_nseffo/(v_nseffo-1))
                 - 2*(v_nseffo*bn-v_nseffo+1)/((v_nseffo-1)*a1)
                 - (3*v_nseffo+1)/(v_nseffo-1));

        v_E = (v_thetaL - v_thetaW*v_lseff) / sqrt(VarE);

        if (esse) free(esse);
        
        if (v_lseff) {
            v_thetaH/=v_lseff;
            v_thetaL/=v_lseff;
        }

    }



//***** FOURTH STEP OF ANALYSIS ****************************************

    void NucleotideDiversity::differentiation() {

        if (!b_analysisSites) throw EggRuntimeError("Diversity: data not available (data not loaded)");
        if (b_differentiation) throw EggRuntimeError("Diversity: cannot re-compute data (implementation error)");
        
        // allocates arrays
        v_pairwiseFixedDifferences = (unsigned int*) malloc((v_npop*(v_npop-1)/2)*sizeof(unsigned int));
        if (!v_pairwiseFixedDifferences) throw EggMemoryError();
        for (unsigned int i=0; i<(v_npop*(v_npop-1)/2); i++) {
            v_pairwiseFixedDifferences[i] = 0;
        }
        v_pairwiseCommonAlleles = (unsigned int*) malloc((v_npop*(v_npop-1)/2)*sizeof(unsigned int));
        if (!v_pairwiseCommonAlleles) throw EggMemoryError();
        for (unsigned int i=0; i<(v_npop*(v_npop-1)/2); i++) {
            v_pairwiseCommonAlleles[i] = 0;
        }
        v_pairwiseSharedAlleles = (unsigned int*) malloc((v_npop*(v_npop-1)/2)*sizeof(unsigned int));
        if (!v_pairwiseSharedAlleles) throw EggMemoryError();
        for (unsigned int i=0; i<(v_npop*(v_npop-1)/2); i++) {
            v_pairwiseSharedAlleles[i] = 0;
        }
        v_popPolymorphic = (unsigned int*) malloc(v_npop*sizeof(unsigned int));
        if (!v_popPolymorphic) throw EggMemoryError();
        for (unsigned int i=0; i<v_npop; i++) {
            v_popPolymorphic[i] = 0;
        }
        v_popSpecific = (unsigned int*) malloc(v_npop*sizeof(unsigned int));
        if (!v_popSpecific) throw EggMemoryError();
        for (unsigned int i=0; i<v_npop; i++) {
            v_popSpecific[i] = 0;
        }
        v_popSpecificDerived = (unsigned int*) malloc(v_npop*sizeof(unsigned int));
        if (!v_popSpecificDerived) throw EggMemoryError();
        for (unsigned int i=0; i<v_npop; i++) {
            v_popSpecificDerived[i] = 0;
        }
        // escapes immediatly if no diversity
        b_differentiation=true; 
        if (!v_S) return;

        // big loop to count site patterns
        unsigned int c,site,pop1,pop2;
        bool iP, hSA, hSAD, FD, CA, SA;
        bool flag1,flag2,flag3,flag4,flag5;
        for (site=0; site<v_S; site++) {
            c=0;
            flag1=flag2=flag3=flag4=flag5=true;
            for (pop1=0; pop1<v_npop; pop1++) {
                // properties of the population
                iP = v_sites[site]->isPolymorphic(pop1);
                hSA = v_sites[site]->hasSpecificAllele(pop1, false);
                if (v_orientables[site]) {
                    hSAD = v_sites[site]->hasSpecificAllele(pop1, true);
                }
                else hSAD = false;
                if (iP) v_popPolymorphic[pop1]++;
                if (hSA)  v_popSpecific[pop1]++;
                if (v_orientables[site] && hSAD) v_popSpecificDerived[pop1]++;

                if (flag4 && hSA) {
                    flag4=0; // there may be several "specific polymorphism" per site in case of ISM violation - in that case we had only one 
                    v_countSpecificAlleles++;
                }
                if (flag5 && hSAD) {
                    flag5=0; // there may be several "specific polymorphism" per site in case of ISM violation - in that case we had only one 
                    v_countSpecificDerivedAlleles++;
                }

                // pairwise patterns
                for (pop2=pop1+1; pop2<v_npop; pop2++) {


                    FD=v_sites[site]->haveFixedDifference(pop1,pop2);
                    CA=v_sites[site]->haveCommonAllele(pop1,pop2);
                    SA=v_sites[site]->haveSharedAllele(pop1,pop2);
                    v_pairwiseFixedDifferences[c]+= FD?1:0;
                    v_pairwiseCommonAlleles[c]+= CA?1:0;
                    v_pairwiseSharedAlleles[c]+= SA?1:0;
                    if (flag1 && FD) {
                        v_countFixedDifferences++;
                        flag1=0;
                    }
                    if (flag2 && CA) {
                        v_countCommonAlleles++;
                        flag2=0;
                    }
                    if (flag3 && SA) {
                        v_countSharedAlleles++;
                        flag3=0;
                    }
                    c++;
                }
            }
        }

    }


//***** FIFTH STEP OF ANALYSIS *****************************************

    void NucleotideDiversity::triConfigurations() {

        if (!b_analysisSites) throw EggRuntimeError("Diversity: data not available (data not loaded)");
        if (b_triConfigurations) throw EggRuntimeError("Diversity: cannot re-compute data (implementation error)");
        
        // allocates and initializes array
        v_triConfigurations = (unsigned int*) malloc(13*sizeof(unsigned int));
        if (!v_triConfigurations) throw EggMemoryError();
        for (unsigned int i=0; i<13; i++) {
            v_triConfigurations[i] = 0;
        }
        
        // escapes immediatly if no diversity or if the number of populations is incorrect
        b_triConfigurations=true; 
        if (!v_S) return;
        if (v_npop!=3) return;

        // sums up the patterns
        for (unsigned int site=0; site<v_S; site++) {
            if (v_sites[site]->numberOfAlleles()!=2) continue;
            if (v_sites[site]->ns(0)==0 || v_sites[site]->ns(1)==0 || v_sites[site]->ns(2)==0) continue;
            v_triConfigurations[0]+=  (v_sites[site]->isPolymorphic(0)         && !v_sites[site]->haveFixedDifference(1,2) && !v_sites[site]->isPolymorphic(1) && !v_sites[site]->isPolymorphic(2))?1:0;
            v_triConfigurations[1]+=  (v_sites[site]->isPolymorphic(0)         &&  v_sites[site]->haveFixedDifference(1,2))?1:0;
            v_triConfigurations[2]+=  (v_sites[site]->isPolymorphic(1)         && !v_sites[site]->haveFixedDifference(0,2) && !v_sites[site]->isPolymorphic(0) && !v_sites[site]->isPolymorphic(2))?1:0;
            v_triConfigurations[3]+=  (v_sites[site]->isPolymorphic(1)         &&  v_sites[site]->haveFixedDifference(0,2))?1:0;
            v_triConfigurations[4]+=  (v_sites[site]->isPolymorphic(2)         && !v_sites[site]->haveFixedDifference(0,1) && !v_sites[site]->isPolymorphic(0) && !v_sites[site]->isPolymorphic(1))?1:0;
            v_triConfigurations[5]+=  (v_sites[site]->isPolymorphic(2)         &&  v_sites[site]->haveFixedDifference(0,1))?1:0;
            v_triConfigurations[6]+=  (v_sites[site]->haveSharedAllele(0,1)    && !v_sites[site]->isPolymorphic(2)        )?1:0;
            v_triConfigurations[7]+=  (v_sites[site]->haveSharedAllele(0,2)    && !v_sites[site]->isPolymorphic(1)        )?1:0;
            v_triConfigurations[8]+=  (v_sites[site]->haveSharedAllele(1,2)    && !v_sites[site]->isPolymorphic(0)        )?1:0;
            v_triConfigurations[9]+=  (v_sites[site]->haveSharedAllele(0,1)    &&  v_sites[site]->haveSharedAllele(0,2)   )?1:0;
            v_triConfigurations[10]+= (v_sites[site]->haveFixedDifference(0,1) &&  v_sites[site]->haveFixedDifference(0,2))?1:0;
            v_triConfigurations[11]+= (v_sites[site]->haveFixedDifference(0,1) &&  v_sites[site]->haveFixedDifference(1,2))?1:0;
            v_triConfigurations[12]+= (v_sites[site]->haveFixedDifference(0,2) &&  v_sites[site]->haveFixedDifference(1,2))?1:0;
        }

    }
    
    





//***** siteAnalysis accessors **************************************


    unsigned int NucleotideDiversity::S() const {
        if (!b_analysisSites) throw EggRuntimeError("NucleotideDiversity: data not available (data not loaded)");
        return v_S;
    }

    unsigned int NucleotideDiversity::So() const {
        if (!b_analysisSites) throw EggRuntimeError("NucleotideDiversity: data not available (data not loaded)");
        return v_So;
    }
    
    unsigned int NucleotideDiversity::eta() const {
        if (!b_analysisSites) throw EggRuntimeError("NucleotideDiversity: data not available (data not loaded)");
        return v_eta;
    }
    
    double NucleotideDiversity::nseff() const {
        if (!b_analysisSites) throw EggRuntimeError("NucleotideDiversity: data not available (data not loaded)");
        return v_nseff;
    }
    
    unsigned int NucleotideDiversity::lseff() const {
        if (!b_analysisSites) throw EggRuntimeError("NucleotideDiversity: data not available (data not loaded)");
        return v_lseff;
    }
    
    double NucleotideDiversity::nseffo() const {
        if (!b_analysisSites) throw EggRuntimeError("NucleotideDiversity: data not available (data not loaded)");
        return v_nseffo;
    }
    
    unsigned int NucleotideDiversity::lseffo() const {
        if (!b_analysisSites) throw EggRuntimeError("NucleotideDiversity: data not available (data not loaded)");
        return v_lseffo;
    }
    
    unsigned int NucleotideDiversity::npop() const  {
        if (!b_analysisSites) throw EggRuntimeError("NucleotideDiversity: data not available (data not loaded)");
        return v_npop;
    }
    
    unsigned int NucleotideDiversity::popLabel(unsigned int popIndex) const {
        if (!b_analysisSites) throw EggRuntimeError("NucleotideDiversity: data not available (data not loaded)");
        return v_popLabel[popIndex];
    }




//***** diversity accessors ********************************************

    double NucleotideDiversity::Pi() {
        if (!b_diversity) diversity();
        return v_Pi;
    }

    double NucleotideDiversity::thetaW() {
        if (!b_diversity) diversity();
        return v_thetaW;
    }

    double NucleotideDiversity::average_Pi() {
        if (!b_diversity) diversity();
        return v_average_Pi;
    }

    double NucleotideDiversity::pop_Pi(unsigned int popIndex) {
        if (!b_diversity) diversity();
        return v_pop_Pi[popIndex];
    }

    double NucleotideDiversity::D() {
        if (!b_diversity) diversity();
        return v_D;
    }


//***** outgroup diversity accessor ***********************************/

    double NucleotideDiversity::H() {
        if (!b_outgroupDiversity) outgroupDiversity();
        return v_H;
    }

    double NucleotideDiversity::Z() {
        if (!b_outgroupDiversity) outgroupDiversity();
        return v_Z;
    }

    double NucleotideDiversity::E() {
        if (!b_outgroupDiversity) outgroupDiversity();
        return v_E;
    }

    double NucleotideDiversity::thetaH() {
        if (!b_outgroupDiversity) outgroupDiversity();
        return v_thetaH;
    }

    double NucleotideDiversity::thetaL() {
        if (!b_outgroupDiversity) outgroupDiversity();
        return v_thetaL;
    }


//***** differentiation accessors *************************************/

    unsigned int NucleotideDiversity::FixedDifferences() {
        if (!b_differentiation) differentiation();
        return v_countFixedDifferences;
    }
    
    unsigned int NucleotideDiversity::CommonAlleles() {
        if (!b_differentiation) differentiation();
        return v_countCommonAlleles;
    }
    
    unsigned int NucleotideDiversity::SharedAlleles() {
        if (!b_differentiation) differentiation();
        return v_countSharedAlleles;
    }
    
    unsigned int NucleotideDiversity::SpecificAlleles() {
        if (!b_differentiation) differentiation();
        return v_countSpecificAlleles;
    }
    
    unsigned int NucleotideDiversity::SpecificDerivedAlleles() {
        if (!b_differentiation) differentiation();
        return v_countSpecificDerivedAlleles;
    }
    
    unsigned int NucleotideDiversity::Polymorphisms(unsigned int pop) {
        if (!b_differentiation) differentiation();
        if (pop>=v_npop) throw EggArgumentValueError("invalid population index while accessing polymorphism data");
        return v_popPolymorphic[pop];
    }
    
    unsigned int NucleotideDiversity::SpecificAlleles(unsigned int pop) {
        if (!b_differentiation) differentiation();
        if (pop>=v_npop) throw EggArgumentValueError("invalid population index while accessing polymorphism data");
        return v_popSpecific[pop];
    }
    
    unsigned int NucleotideDiversity::SpecificDerivedAlleles(unsigned int pop) {
        if (!b_differentiation) differentiation();
        if (pop>=v_npop) throw EggArgumentValueError("invalid population index while accessing polymorphism data");
        return v_popSpecificDerived[pop];
    }
    
    unsigned int NucleotideDiversity::FixedDifferences(unsigned int pop1, unsigned int pop2) {
        if (!b_differentiation) differentiation();
        if (pop1>=v_npop || pop2>=v_npop || pop1>=pop2) {
            throw EggArgumentValueError("invalid population indices while accessing polymorphism data");
        }
        unsigned int c=0;
        for (unsigned int i=0; i<v_npop; i++) {
            for (unsigned int j=i+1; j<v_npop; j++) {
                if (i==pop1 && j==pop2) break;
                c++;
            }
        }
        return v_pairwiseFixedDifferences[c];
    }
    
    unsigned int NucleotideDiversity::CommonAlleles(unsigned int pop1, unsigned int pop2) {
        if (!b_differentiation) differentiation();
        if (pop1>=v_npop || pop2>=v_npop || pop1>=pop2) {
            throw EggArgumentValueError("invalid population indices while accessing polymorphism data");
        }
        unsigned int c=0;
        for (unsigned int i=0; i<v_npop; i++) {
            for (unsigned int j=i+1; j<v_npop; j++) {
                if (i==pop1 && j==pop2) break;
                c++;
            }
        }
        return v_pairwiseCommonAlleles[c];
    }
    
    unsigned int NucleotideDiversity::SharedAlleles(unsigned int pop1, unsigned int pop2) {
        if (!b_differentiation) differentiation();
        if (pop1>=v_npop || pop2>=v_npop || pop1>=pop2) {
            throw EggArgumentValueError("invalid population indices while accessing polymorphism data");
        }
        unsigned int c=0;
        for (unsigned int i=0; i<v_npop; i++) {
            for (unsigned int j=i+1; j<v_npop; j++) {
                if (i==pop1 && j==pop2) break;
                c++;
            }
        }
        return v_pairwiseSharedAlleles[c];
    }


// ***** tri configuration accessor ************************************

    unsigned int NucleotideDiversity::triConfiguration(unsigned int index) {
        if (!b_triConfigurations) triConfigurations();
        return v_triConfigurations[index];
    }


// ***** site position accessors ***************************************

        std::vector<unsigned int> NucleotideDiversity::polymorphic_positions() const {
            std::vector<unsigned int> v(v_sitePositions, v_sitePositions+v_S);
            return v;
        }

        std::vector<unsigned int> NucleotideDiversity::singleton_positions() const {
            std::vector<unsigned int> v;
            
            for (unsigned int i=0; i<v_S; i++) {
                if (v_sites[i]->numberOfAlleles() == 2 &&
                    (v_sites[i]->alleleFrequency(0)==1 ||
                     v_sites[i]->alleleFrequency(1)==1))
                {
                         v.push_back(v_sitePositions[i]);
                }
            }
            
            return v;
        }

}
