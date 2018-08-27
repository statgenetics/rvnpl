/*
    Copyright 2008-2009,2013 Stéphane De Mita, Mathieu Siol

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

#include "BppDiversity.hpp"
#include "EggException.hpp"
#include "Convert.hpp"
#include <Bpp/Exceptions.h>

namespace egglib {

    BppDiversity::~BppDiversity() {
        if (p_code4) delete p_code4;
        if (p_code5) delete p_code5;
        if (p_code6) delete p_code6;
        if (p_code7) delete p_code7;
    }


    BppDiversity::BppDiversity() {
        p_code4 = new(std::nothrow) bpp::StandardGeneticCode(&dna);
        if (!p_code4) throw EggMemoryError();

        p_code5 = new(std::nothrow) bpp::VertebrateMitochondrialGeneticCode(&dna);
        if (!p_code5) throw EggMemoryError();

        p_code6 = new(std::nothrow) bpp::InvertebrateMitochondrialGeneticCode(&dna);
        if (!p_code6) throw EggMemoryError();

        p_code7 = new(std::nothrow) bpp::EchinodermMitochondrialGeneticCode(&dna);
        if (!p_code7) throw EggMemoryError();
    }
    
    
/* ****************************************************************** */

    
    void BppDiversity::reset() {
             b_loaded = false;
             b_coding = false;
           b_outgroup = false;
                 v_Ss = 0;
               v_Sinf = 0;
               v_Ssin = 0;
                v_eta = 0;
               v_Sext = 0;
                 v_He = 0.;
                v_He2 = 0.;
                 //v_GC = 0.;
                 v_tW = 0.;
                v_T83 = 0.;
                  v_K = 0;
                  v_H = 0.;
                 v_Ti = 0;
                 v_Tv = 0;
               v_TiTv = 0.;
              v_nstop = 0;
         v_ncodon1mut = 0;
               v_nsyn = 0;
                v_tWS = 0.;
               v_tWNS = 0.;
                v_PiS = 0.;
               v_PiNS = 0.;
             v_Ssites = 0.;
            v_NSsites = 0.;
                 v_SS = 0;
                v_SNS = 0;
                  v_MK.clear();
                 v_NI = 0.;
                  v_D = 0.;
               v_Deta = 0.;
                v_Dfl = 0.;
            v_Dflstar = 0.;
                  v_F = 0.;
              v_Fstar = 0.;
               v_rhoH = 0.;
    }


/* ****************************************************************** */

    void BppDiversity::load(Align& align, unsigned int dataType) {

        reset();
        b_loaded = true;
        b_coding = (dataType>=4 && dataType<=7);
        if (dataType == 4) p_code = p_code4;
        if (dataType == 5) p_code = p_code5;
        if (dataType == 6) p_code = p_code6;
        if (dataType == 7) p_code = p_code7;

        try {

            bpp::AlignedSequenceContainer ingroup = Convert::egglib2bpp(align, dataType, 1);

            // intra pop stats
            computeSingle(ingroup);
            bpp::AlignedSequenceContainer outgroup = Convert::egglib2bpp(align, dataType, 2, ingroup.getNumberOfSequences());

            if (outgroup.getNumberOfSequences()>0) {
                b_outgroup = true;

                // inter pop stats
                computeDouble(ingroup, outgroup);
            }
        }
        catch (bpp::Exception bppExcept) {
            reset();
            throw EggRuntimeError(bppExcept.what());
        }
    }
    
    
/* ****************************************************************** */


    void BppDiversity::computeSingle(bpp::AlignedSequenceContainer& align) {

                      v_Ss = bpp::SequenceStatistics::polymorphicSiteNumber(align);
                    v_Sinf = bpp::SequenceStatistics::parsimonyInformativeSiteNumber(align);
                    v_Ssin = bpp::SequenceStatistics::countSingleton(align);
                     v_eta = bpp::SequenceStatistics::totNumberMutations(align);
                      v_He = bpp::SequenceStatistics::heterozygosity(align);
                     v_He2 = bpp::SequenceStatistics::squaredHeterozygosity(align);
                      //v_GC = bpp::SequenceStatistics::gcContent(align);
                      v_tW = bpp::SequenceStatistics::watterson75(align);
                     v_T83 = bpp::SequenceStatistics::tajima83(align);
                       v_K = bpp::SequenceStatistics::DVK(align);
                       v_H = bpp::SequenceStatistics::DVH(align);
                      v_Ti = bpp::SequenceStatistics::getNumberOfTransitions(align);
                      v_Tv = bpp::SequenceStatistics::getNumberOfTransversions(align);
   if (v_Tv!=0)     v_TiTv = bpp::SequenceStatistics::getTransitionsTransversionsRatio(align);
   else             v_TiTv = 0.;
   if (v_Ss!=0)        v_D = bpp::SequenceStatistics::tajimaDSS(align);
   else                v_D = 0.;
   if (v_eta!=0)    v_Deta = bpp::SequenceStatistics::tajimaDTNM(align);
   else             v_Deta = 0.;
   if (v_eta!=0) v_Dflstar = bpp::SequenceStatistics::fuliDstar(align);
   else          v_Dflstar = 0.;
   if (v_eta!=0)   v_Fstar = bpp::SequenceStatistics::fuliFstar(align);
   else            v_Fstar = 0.;
                    v_rhoH = bpp::SequenceStatistics::hudson87(align);

        if (b_coding) {
            v_ncodon1mut = bpp::SequenceStatistics::monoSitePolymorphicCodonNumber(align);
                 v_nstop = bpp::SequenceStatistics::stopCodonSiteNumber(align);
                  v_nsyn = bpp::SequenceStatistics::synonymousPolymorphicCodonNumber(align, *p_code);
                   v_tWS = bpp::SequenceStatistics::watterson75Synonymous(align, *p_code);
                  v_tWNS = bpp::SequenceStatistics::watterson75NonSynonymous(align, *p_code);
                   v_PiS = bpp::SequenceStatistics::piSynonymous(align, *p_code);
                  v_PiNS = bpp::SequenceStatistics::piNonSynonymous(align, *p_code);
                v_Ssites = bpp::SequenceStatistics::meanSynonymousSitesNumber(align, *p_code);
               v_NSsites = bpp::SequenceStatistics::meanNonSynonymousSitesNumber(align, *p_code);
                    v_SS = bpp::SequenceStatistics::synonymousSubstitutionsNumber(align, *p_code);
                   v_SNS = bpp::SequenceStatistics::nonSynonymousSubstitutionsNumber(align, *p_code);
        }
    }


    void BppDiversity::computeDouble(bpp::AlignedSequenceContainer& ingroup, bpp::AlignedSequenceContainer& outgroup) {
            v_Sext = bpp::SequenceStatistics::totMutationsExternalBranchs(ingroup, outgroup);
    if (v_eta!=0) v_Dfl = bpp::SequenceStatistics::fuliD(ingroup, outgroup);
                else v_Dfl = 0.;
    if (v_eta!=0) v_F = bpp::SequenceStatistics::fuliF(ingroup, outgroup);
                else v_F = 0.;

            if (b_coding) {
                std::vector<size_t> temp = bpp::SequenceStatistics::MKtable(ingroup, outgroup, *p_code);
                v_MK.resize(4);
                    v_MK[0] = (unsigned int) temp[0];
                        v_MK[1] = (unsigned int) temp[1];
                            v_MK[2] = (unsigned int) temp[2];
                                v_MK[3] = (unsigned int) temp[3];
                v_NI = bpp::SequenceStatistics::neutralityIndex(ingroup, outgroup, *p_code);
            }
    }


/* ****************************************************************** */

    bool BppDiversity::hasOutgroup() const {
        return b_outgroup;
    }

    unsigned int BppDiversity::S() const {
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_Ss;           
    }
    
    unsigned int BppDiversity::Sinf() const {
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_Sinf;
    }

    unsigned int BppDiversity::Ssin() const {
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_Ssin;         
    }

    unsigned int BppDiversity::eta() const {
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_eta;          
    }

    unsigned int BppDiversity::Sext() const {
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_outgroup) throw EggRuntimeError("BppDiversity: cannot access data (no outgroup)");
        return v_Sext;         
    }

    double BppDiversity::He() const {
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_He;           
    }

    double BppDiversity::He2() const {
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_He2;          
    }

//    double BppDiversity::GC() const { 
//        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
//        return v_GC;           
//    }

    double BppDiversity::tW() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_tW;           
    }

    double BppDiversity::T83() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_T83;          
    }

    unsigned int BppDiversity::K() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_K;            
    }

    double BppDiversity::H() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_H;            
    }

    unsigned int BppDiversity::Ti() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_Ti;           
    }

    unsigned int BppDiversity::Tv() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_Tv;           
    }

    double BppDiversity::TiTv() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_TiTv;         
    }

    unsigned int BppDiversity::nstop() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_nstop;        
    }

    unsigned int BppDiversity::ncodon1mut() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_ncodon1mut;   
    }

    unsigned int BppDiversity::nsyn() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_nsyn;         
    }

    double BppDiversity::tWS() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_tWS;          
    }

    double BppDiversity::tWNS() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_tWNS;         
    }

    double BppDiversity::PiS() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_PiS;          
    }

    double BppDiversity::PiNS() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_PiNS;         
    }

    double BppDiversity::Ssites() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_Ssites;       
    }

    double BppDiversity::NSsites() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_NSsites;      
    }

    unsigned int BppDiversity::SS() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_SS;           
    }

    unsigned int BppDiversity::SNS() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        return v_SNS;          
    }

    const std::vector<unsigned int>& BppDiversity::MK() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_outgroup) throw EggRuntimeError("BppDiversity: cannot access data (no outgroup)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
       return v_MK;           
    }

    double BppDiversity::NI() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_coding) throw EggRuntimeError("BppDiversity: cannot access data (not coding sequences)");
        if (!b_outgroup) throw EggRuntimeError("BppDiversity: cannot access data (no outgroup)");
        return v_NI;           
    }

    double BppDiversity::D() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_D;            
    }

    double BppDiversity::Deta() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_Deta;         
    }

    double BppDiversity::Dfl() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_outgroup) throw EggRuntimeError("BppDiversity: cannot access data (no outgroup)");
        return v_Dfl;          
    }

    double BppDiversity::Dflstar() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_Dflstar;      
    }

    double BppDiversity::F() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        if (!b_outgroup) throw EggRuntimeError("BppDiversity: cannot access data (no outgroup)");
        return v_F;            
    }

    double BppDiversity::Fstar() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_Fstar;        
    }

    double BppDiversity::rhoH() const { 
        if (!b_loaded) throw EggRuntimeError("BppDiversity: cannot access data (not computed)");
        return v_rhoH;         
    }

}
