/*
    Copyright 2009,2011 Stéphane De Mita, Mathieu Siol

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


#include <vector>
#include <deque>
#include <algorithm>

#include "Convert.hpp"
#include "EggException.hpp"
#include "config.h"

namespace egglib {


    Align Convert::align(
        DataMatrix& dataMatrix,
        unsigned int length,
        Random* random,
        bool randomizePositions,
        bool randomizeNonVaryingStates,
        bool randomizeAlleles,
        bool enforceLength,
        std::string mapping,
        char unknown,
        char nonVaryingState
        )
    {
    
        // check random
        bool killRandom = false;
        if (random==NULL && (randomizePositions || randomizeNonVaryingStates)) {
            random = new(std::nothrow) Random;
            if (!random) throw EggMemoryError();
            killRandom=true;
        }
        
        // determines matrix dimensions
        unsigned int ns = dataMatrix.numberOfSequences();
        unsigned int ls = 0;
        unsigned int S = dataMatrix.numberOfSites();

        if (length!=0 && enforceLength && length<S) {
            throw EggRuntimeError("cannot convert DataMatrix to Align: too many segregating sites for the specified alignment length");
        }
        if (length<S) ls = S;
        else ls = length;
        
        // ##### DETERMINES POSITIONS #####
        std::vector<unsigned int> positions(S,0);
        unsigned int c;

        // don't worry of anything if S=0
        if (S) {
                        
            // if S=ls then all sites are contiguous
            if (length<=S) {
                for (unsigned int i=0; i<S; i++) {
                    positions[i] = i;
                }
            }
                    
            else  {
                        
                // places the sites as logically as possible (but shift to the left when needed, and back to the right if the extreme left end is jammed)
                std::vector<double> vector(S); // a queue before 2.0.3

                // uses positions from the DataMatrix object
                if (!randomizePositions) {
                    // queue of positions x length
                    for (unsigned int i=0; i<S; i++) {
                        vector[i] = dataMatrix.sitePosition(i)*ls;
                    }
                    std::sort(vector.begin(), vector.end()); // added in 2.0.3
                }

                // or, randomizes them
                else {
                    std::vector<double> randoms;
                    for (unsigned int i=0; i<S; i++) randoms.push_back( random->uniform() );
                    std::sort(randoms.begin(), randoms.end());
                    c=0;
                    for (std::vector<double>::iterator it = randoms.begin(); it!=randoms.end(); it++) {
                        vector[c++] = *it*ls;
                    }
                    
                }
                
                // adds the sites as they come
                c=0;
                unsigned int vector_pos = 0;
                for (unsigned int i=0; i<ls; i++) {
                    if (vector_pos==S) break;
                    if (i+1>vector[vector_pos]) {
                        positions[c++] = i;
                        vector_pos++;
                    }
                }

                // if there is still sites left to place, shift all to the left
                if (vector_pos<S) {
                 
                    if (positions[S-1] != 0) {
                        throw EggRuntimeError("an unexpected error occurred while placing segregating sites [REF:1]");
                    }
                    positions[S-1] = ls-1;
                    unsigned int i=S-1;
                    while (true) {
                        if (i--==0) throw EggRuntimeError("an unexpected error occurred while placing segregating sites [REF:2]");
                        if (positions[i]>=positions[i+1] || (positions[i]==0 && i!=0)) {
                            positions[i]=positions[i+1]-1;
                        }
                        else break;
                    }
                }

                // sanity checks
                if (positions[S-1]>=ls) {
                    throw EggRuntimeError("an unexpected error occurred while placing segregating sites [REF:3]");
                }
                for (unsigned int i=0; i<S-1; i++) {
                    if (positions[i]>=positions[i+1]) throw EggRuntimeError("an unexpected error occurred while placing segregating sites [REF:4]");
                }
            } 

        }



        // ##### PRODUCES MUTATIONS #####

        // instantiates the alignment object
        Align align(ns, ls);
            
        // iterates of the sites (variable or not)
        unsigned int segregatingSite=0;
        for (unsigned int site=0; site<ls; site++) {

            // for a fixed site
            if (segregatingSite==S || site!=positions[segregatingSite]) {

                char c=nonVaryingState;
                if (randomizeNonVaryingStates) {
                    if (!mapping.size()) throw EggRuntimeError("cannot convert to align: cannot draw non-varying state (no characters were provided)");
                    c = mapping[ random->irand(mapping.size()) ];
                }
                for (unsigned int i=0; i<ns; i++) {
                    align.set(i, site, c);
                }
            }
            
            // for a polymorphic site
            else {
                
                // if needed , reshuffles the mapping
                std::string rmapping;
                if (randomizeAlleles) {
                    std::deque<char> set(mapping.begin(), mapping.end());
                    for (unsigned int i=0; i<mapping.size()-1; i++) {
                        unsigned int X = random->irand(set.size());
                        rmapping.push_back(set[X]);
                        set.erase(set.begin()+X);
                    }
                    rmapping.push_back(set.front());
                }
                
                // sets the alleles
                for (unsigned int i=0; i<ns; i++) {
                    char c = unknown;
                    int a = dataMatrix.get( i, segregatingSite );
                    if (a>=0 && a<(int)mapping.size()) {
                        if (!randomizeAlleles) c = mapping[a];
                        else c = rmapping[a];
                    }
                    align.set(i, site, c);

                }


                // points next segregating site
                segregatingSite++;
            }
        }
        
        // sanity check
        if (segregatingSite!=S) throw EggRuntimeError("");


        // ##### COPIES GROUP LABELS #####
        for (unsigned int i=0; i<ns; i++) align.group( i, dataMatrix.populationLabel(i) );

        
        // the end
        if (killRandom) delete random;
        return align;
    }


#ifdef HAVE_LIBBPP_SEQ

    bpp::AlignedSequenceContainer Convert::egglib2bpp(Align& align, unsigned int alphabetID, unsigned int outgroupFlag, unsigned int offset) {

        bpp::Alphabet* alph = NULL;
        switch (alphabetID) {
            case 1:
                alph = &dnaAlphabet;
                break;
            case 2:
                alph = &rnaAlphabet;
                break;
            case 3:
                alph = &proteicAlphabet;
                break;
            case 4:
                alph = &standardCodonAlphabet;
                break;
            case 5:
                alph = &vertebrateMitochondrialCodonAlphabet;
                break;
            case 6:
                alph = &invertebrateMitochondrialCodonAlphabet;
                break;
            case 7:
                alph = &echinodermMitochondrialCodonAlphabet;
                break;
            default:
                throw EggArgumentValueError("invalid alphabetID value");
        }
        
        if (outgroupFlag>2) throw EggArgumentValueError("invalid outgroupFlag value");
        
        bpp::AlignedSequenceContainer bppAlign(alph);
        unsigned c=offset;
        for (unsigned int i=0; i<align.ns(); i++) {

            if (outgroupFlag==1 && align.group(i)==999) continue;
            if (outgroupFlag==2 && align.group(i)!=999) continue;

            std::ostringstream str;
            str << c++;
            bpp::BasicSequence sequence(str.str(), align.sequence(i), alph);
            bppAlign.addSequence(sequence);

        }
        return bppAlign;
        
    }


bpp::DNA Convert::dnaAlphabet;
bpp::RNA Convert::rnaAlphabet;
bpp::ProteicAlphabet Convert::proteicAlphabet;
bpp::StandardCodonAlphabet Convert::standardCodonAlphabet(&Convert::dnaAlphabet);
bpp::VertebrateMitochondrialCodonAlphabet Convert::vertebrateMitochondrialCodonAlphabet(&Convert::dnaAlphabet);
bpp::InvertebrateMitochondrialCodonAlphabet Convert::invertebrateMitochondrialCodonAlphabet(&Convert::dnaAlphabet);
bpp::EchinodermMitochondrialCodonAlphabet Convert::echinodermMitochondrialCodonAlphabet(&Convert::dnaAlphabet);

#endif


}
