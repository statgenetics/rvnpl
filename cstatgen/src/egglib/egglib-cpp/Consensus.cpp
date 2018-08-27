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

#include "Consensus.hpp"
#include "EggException.hpp"

namespace egglib {

    Consensus::Consensus() {
        MISSING = '?';
        DISAGREEMENT = 'Z';
    }

    void Consensus::setMissing(char c) {
        MISSING = c;
    }
    
    void Consensus::setDisagreement(char c) {
        DISAGREEMENT = c;
    }

    Align Consensus::consensus(Align& align, char separator, bool rigorous) {
        t_firstSequenceNames.clear();
        t_secondSequenceNames.clear();
        t_roots.clear();
        t_consistentPositions.clear();
        t_complementaryPositions.clear();
        t_uninformativePositions.clear();
        t_ambiguousPositions.clear();
        t_atLeastPartiallyResolvedAmbiguities.clear();
        t_inconsistentPositions.clear();
        
        Align consensus;
        
        std::string tname1, tname2, name1, name2;
        std::vector< std::vector<unsigned int> > complements(align.ns());
        std::vector<bool> isComplement(align.ns(),false);
        PairwiseConsensus pw;
        pw.setUndeterminedCharacter(MISSING);
        pw.setDisagreementCharacter(DISAGREEMENT);

        for (unsigned int i=0; i<align.ns(); i++) {
            if (isComplement[i]) continue;
            for (unsigned int j=i+1; j<align.ns(); j++) {
                tname1 = align.name(i);
                name1 = tname1.substr(0,tname1.find(separator));
                tname2 = align.name(j);
                name2 = tname2.substr(0,tname2.find(separator));
                if (name1==name2) {
                    complements[i].push_back(j);
                    complements[j].push_back(i);
                    isComplement[j] = true;

                    pw.load(align.sequence(i),align.sequence(j));
                    if (!rigorous) pw.generateSoftConsensus();
                    else pw.generateHardConsensus();
                    consensus.append(name1.c_str(),pw.getConsensus().c_str());
                    
                    t_firstSequenceNames.push_back(tname1);
                    t_secondSequenceNames.push_back(tname2);
                    t_roots.push_back(name1);
                    t_consistentPositions.push_back( pw.getConsistentPositions() );
                    t_complementaryPositions.push_back( pw.getComplementaryPositions() );
                    t_uninformativePositions.push_back( pw.getUninformativePositions() );
                    t_ambiguousPositions.push_back( pw.getAmbiguousPositions() );
                    t_atLeastPartiallyResolvedAmbiguities.push_back( pw.getAtLeastPartiallyResolvedAmbiguities() );
                    t_inconsistentPositions.push_back( std::vector<int>(pw.getInconsistentPositions()) );
                    for (int i=0; i<pw.getInconsistentPositions(); i++) {
                        t_inconsistentPositions.back()[i] = pw.getThisInconsistentPosition(i);
                    }
                }
            }
            if (complements[i].size()==0) {
                std::string name = align.name(i);
                std::string root = name.substr(0, name.find(separator));
                consensus.append(root.c_str(), align.sequence(i));  // even singletons names are truncated

                t_firstSequenceNames.push_back(name);
                t_secondSequenceNames.push_back("");
                t_roots.push_back(root);
                t_consistentPositions.push_back( 0 );
                t_complementaryPositions.push_back( 0 );
                t_uninformativePositions.push_back( 0 );
                t_ambiguousPositions.push_back( 0 );
                t_atLeastPartiallyResolvedAmbiguities.push_back( 0 );
                t_inconsistentPositions.push_back( std::vector<int>() );

            }
            else {
                while(complements[i].size()>1) {
                    pw.load(consensus.sequence(consensus.ns()-2),consensus.sequence(consensus.ns()-1));
                    if (!rigorous) pw.generateSoftConsensus();
                    else pw.generateHardConsensus();
                    consensus.remove(consensus.ns()-1);
                    consensus.remove(consensus.ns()-1);
                    std::string name = align.name(i);
                    std::string root = name.substr(0, name.find(separator));
                    consensus.append(root.c_str(), pw.getConsensus().c_str());
                    complements[i].pop_back();

                    t_firstSequenceNames.push_back(name);
                    t_secondSequenceNames.push_back(name);
                    t_roots.push_back(root);
                    t_consistentPositions.push_back( pw.getConsistentPositions() );
                    t_complementaryPositions.push_back( pw.getComplementaryPositions() );
                    t_uninformativePositions.push_back( pw.getUninformativePositions() );
                    t_ambiguousPositions.push_back( pw.getAmbiguousPositions() );
                    t_atLeastPartiallyResolvedAmbiguities.push_back( pw.getAtLeastPartiallyResolvedAmbiguities() );
                    t_inconsistentPositions.push_back( std::vector<int>(pw.getInconsistentPositions()) );
                    for (int i=0; i<pw.getInconsistentPositions(); i++) {
                        t_inconsistentPositions.back()[i] = pw.getThisInconsistentPosition(i);
                    }
                }
            }
        }

        return consensus;
    }



    // implementation of the private classes Consensus::PairwiseConsensus

    Consensus::PairwiseConsensus::PairwiseConsensus() {
        setCharacterContainers();
        MISSING= '?';
        DISAGREEMENT= 'Z';
    };

    Consensus::PairwiseConsensus::PairwiseConsensus(std::string seqA, std::string seqB) {
        setCharacterContainers();
        MISSING= '?';
        DISAGREEMENT= 'Z';
        load(seqA,seqB);
    }

    void Consensus::PairwiseConsensus::load(std::string seqA, std::string seqB) {
        ls= (seqA.size()>seqB.size())?seqA.size():seqB.size();
        this->seqA = seqA;
        this->seqB = seqB;

        if (seqA.size()<ls) seqA.resize(ls,'?');
        if (seqB.size()<ls) seqB.resize(ls,'?');

        cons.resize(ls,'Z');
    }

    void Consensus::PairwiseConsensus::setUndeterminedCharacter(char c) { MISSING=c; }
    void Consensus::PairwiseConsensus::setDisagreementCharacter(char c) { DISAGREEMENT=c; }

    int Consensus::PairwiseConsensus::generateHardConsensus() {
        cntConsistentPositions= 0;
        cntComplementaryPositions= 0;
        cntUninformativePositions= 0;
        cntAmbiguousPositions= 0;
        cntAtLeastPartiallyResolvedAmbiguities= 0;
        cntInconsistentPositions= 0;
        posIncons.clear();

        CharacterContainer cc1, cc2;

        for (unsigned int i=0;i<ls;i++) {
            cc1 = ccN.init(seqA[i]);
            cc2 = ccN.init(seqB[i]);
            if(cc1.is(cc2)) {
                cons[i]= seqA[i];
                if (isValid(seqA[i])) cntConsistentPositions++;
                else if (!ccN.has(cc1)) cntUninformativePositions++;
                else cntAmbiguousPositions++;
                continue;
            }
            if (ccN.has(cc1) && !ccN.has(cc2)) {
                cons[i] = seqA[i];
                cntComplementaryPositions++;
                continue;
            }
            if (ccN.has(cc2) && !ccN.has(cc1)) {
                cons[i] = seqB[i];
                cntComplementaryPositions++;
                continue;
            }
            if (cc1.has(cc2)) {
                cons[i] = seqB[i];
                cntAtLeastPartiallyResolvedAmbiguities++;
                continue;
            }
            if (cc2.has(cc1)) {
                cons[i] = seqA[i];
                cntAtLeastPartiallyResolvedAmbiguities++;
                continue;
            }
            cons[i] = ccN.lhas(cc1,cc2);
            cntInconsistentPositions++;
            posIncons.push_back(i);
        }
        if (ls!=(unsigned int)cntConsistentPositions+cntComplementaryPositions+
                 cntUninformativePositions+cntAmbiguousPositions+
                 cntInconsistentPositions+cntAtLeastPartiallyResolvedAmbiguities)  throw EggRuntimeError("Consensus: unexpected situation!");
        return cntInconsistentPositions;
    }

    int Consensus::PairwiseConsensus::generateSoftConsensus() {
        cntConsistentPositions= 0;
        cntComplementaryPositions= 0;
        cntUninformativePositions= 0;
        cntAmbiguousPositions= 0;
        cntAtLeastPartiallyResolvedAmbiguities= -1;
        cntInconsistentPositions= 0;
        posIncons.clear();

        for (unsigned int i=0;i<ls;i++) {
            if(seqA[i]==seqB[i]) {
                cons[i]= seqA[i];
                if (isValid(seqA[i])) cntConsistentPositions++;
                else if (seqA[i]==MISSING) cntUninformativePositions++;
                else cntAmbiguousPositions++;
            }
            else if (isValid(seqA[i])) {
                 if (seqB[i]==MISSING) {
                      cons[i]= seqA[i];
                      cntComplementaryPositions++;
                 }
                 else {
                      cons[i]= DISAGREEMENT;
                      cntInconsistentPositions++;
                      posIncons.push_back(i);
                 }
            }
            else if (isValid(seqB[i])) {
                 if (seqA[i]==MISSING) {
                      cons[i]= seqB[i];
                      cntComplementaryPositions++;
                 }
                 else {
                      cons[i]= DISAGREEMENT;
                      cntInconsistentPositions++;
                      posIncons.push_back(i);
                 }
            }
            else if (seqA[i]==MISSING) {
                 cntAmbiguousPositions++;
                 cons[i]= seqB[i];
            }
            else if (seqB[i]==MISSING) {
                 cntAmbiguousPositions++;
                 cons[i]= seqA[i];
            }
            else {
                 cntInconsistentPositions++;
                 cons[i]= DISAGREEMENT;
            }
        }

        if (ls!=(unsigned int)cntConsistentPositions+cntComplementaryPositions+
                 cntUninformativePositions+cntAmbiguousPositions+
                 cntInconsistentPositions) throw EggRuntimeError("Consensus: unexpected situation!");
        return cntInconsistentPositions;
    }

    int Consensus::PairwiseConsensus::getConsistentPositions() { return cntConsistentPositions; }
    int Consensus::PairwiseConsensus::getInconsistentPositions() { return cntInconsistentPositions; }
    int Consensus::PairwiseConsensus::getComplementaryPositions() { return cntComplementaryPositions; }
    int Consensus::PairwiseConsensus::getAmbiguousPositions() { return cntAmbiguousPositions;}
    int Consensus::PairwiseConsensus::getUninformativePositions() { return cntUninformativePositions; }
    int Consensus::PairwiseConsensus::getAtLeastPartiallyResolvedAmbiguities() { return cntAtLeastPartiallyResolvedAmbiguities; }

    int Consensus::PairwiseConsensus::getThisInconsistentPosition( unsigned int i) {return posIncons[i]; }

    std::string Consensus::PairwiseConsensus::getConsensus() { return cons; }

    void Consensus::PairwiseConsensus::setCharacterContainers() {
        std::vector<CharacterContainer> buff(0);
        ccA.setValue('A');
        ccC.setValue('C');
        ccG.setValue('G');
        ccT.setValue('T');
        ccGAP.setValue('-');

        buff.push_back(ccA);
        buff.push_back(ccC);
        ccM.setValue('M'); ccM.setSons(buff);

        buff.clear();
        buff.push_back(ccA);
        buff.push_back(ccG);
        ccR.setValue('R'); ccR.setSons(buff);

        buff.clear();
        buff.push_back(ccA);
        buff.push_back(ccT);
        ccW.setValue('W'); ccW.setSons(buff);

        buff.clear();
        buff.push_back(ccC);
        buff.push_back(ccG);
        ccS.setValue('S'); ccS.setSons(buff);

        buff.clear();
        buff.push_back(ccC);
        buff.push_back(ccT);
        ccY.setValue('Y'); ccY.setSons(buff);

        buff.clear();
        buff.push_back(ccG);
        buff.push_back(ccT);
        ccK.setValue('K'); ccK.setSons(buff);

        buff.clear();
        buff.push_back(ccS);
        buff.push_back(ccY);
        buff.push_back(ccK);
        buff.push_back(ccC);
        buff.push_back(ccG);
        buff.push_back(ccT);
        ccB.setValue('B'); ccB.setSons(buff);

        buff.clear();
        buff.push_back(ccR);
        buff.push_back(ccW);
        buff.push_back(ccK);
        buff.push_back(ccA);
        buff.push_back(ccG);
        buff.push_back(ccT);
        ccD.setValue('D'); ccD.setSons(buff);

        buff.clear();
        buff.push_back(ccM);
        buff.push_back(ccW);
        buff.push_back(ccY);
        buff.push_back(ccA);
        buff.push_back(ccC);
        buff.push_back(ccT);
        ccH.setValue('H'); ccH.setSons(buff);

        buff.clear();
        buff.push_back(ccM);
        buff.push_back(ccR);
        buff.push_back(ccS);
        buff.push_back(ccA);
        buff.push_back(ccC);
        buff.push_back(ccG);
        ccV.setValue('V'); ccV.setSons(buff);
        buff.clear();

        buff.clear();
        buff.push_back(ccB);
        buff.push_back(ccD);
        buff.push_back(ccH);
        buff.push_back(ccV);
        buff.push_back(ccM);
        buff.push_back(ccR);
        buff.push_back(ccW);
        buff.push_back(ccS);
        buff.push_back(ccY);
        buff.push_back(ccK);
        buff.push_back(ccA);
        buff.push_back(ccC);
        buff.push_back(ccG);
        buff.push_back(ccT);
        buff.push_back(ccGAP);
        ccN.setValue('N'); ccN.setSons(buff);

        buff.push_back('?');
        buff.push_back(ccN);
        ccQ.setValue('?'); ccQ.setSons(buff);
    }

    Consensus::PairwiseConsensus::CharacterContainer::CharacterContainer() { value = '@'; }
    Consensus::PairwiseConsensus::CharacterContainer::CharacterContainer(const char& c)  { value = c; }

    Consensus::PairwiseConsensus::CharacterContainer& Consensus::PairwiseConsensus::CharacterContainer::operator=(const char& c) {
        value = c;
        sons.clear();
        return *this;
    }

    void Consensus::PairwiseConsensus::CharacterContainer::setValue(char c) { value = c; }
    void Consensus::PairwiseConsensus::CharacterContainer::setSons(std::vector<CharacterContainer> inSons) {
         sons.clear();
         sons.resize(inSons.size());
         for (unsigned i=0; i<inSons.size(); i++) sons[i] = inSons[i];
    }

    bool Consensus::PairwiseConsensus::CharacterContainer::is(CharacterContainer c) { return value == c.value; }

    bool Consensus::PairwiseConsensus::CharacterContainer::has(CharacterContainer c) {
         if (value==c.value) return true;   // s'assurer de l'innocuite de ce test: est-ce que l'identite implique la possession?
         for (unsigned int i=0; i<sons.size(); i++) if (sons[i].is(c)) return true;
         return false;
    }

    bool Consensus::PairwiseConsensus::CharacterContainer::has(char c) {
         for (unsigned int i=0; i<sons.size(); i++) if (sons[i].value == c) return true;
         return false;
    }

    char Consensus::PairwiseConsensus::CharacterContainer::lhas(CharacterContainer c1, CharacterContainer c2) {
         if (c1.is(c2)) return c1.value;
         unsigned int C = 0;
         for (unsigned int i=0; i<sons.size(); i++) {
             if (sons[i].is(c1) || sons[i].is(c2)) C++;
             if (sons[i].has(c1) && sons[i].has(c2)) return sons[i].lhas(c1,c2);
         }
         if (C==2) return value;
         else return '@';
    }

    Consensus::PairwiseConsensus::CharacterContainer Consensus::PairwiseConsensus::CharacterContainer::init(char c) {
         CharacterContainer output('@');
         if (is(c)) {
            output.setValue(value);
            output.setSons(sons);
         }
         for (unsigned int i=0; i<sons.size(); i++) if (sons[i].value == c) {
            output.setValue(sons[i].value);
            output.setSons(sons[i].sons);
         }
         return output;
    }

    bool Consensus::check_sequences(Align& align) {
        PairwiseConsensus::CharacterContainer cc;
        std::vector<PairwiseConsensus::CharacterContainer> buff;
        buff.push_back('B'); buff.push_back('b');
        buff.push_back('D'); buff.push_back('d');
        buff.push_back('H'); buff.push_back('h');
        buff.push_back('V'); buff.push_back('v');
        buff.push_back('M'); buff.push_back('m');
        buff.push_back('R'); buff.push_back('r');
        buff.push_back('W'); buff.push_back('w');
        buff.push_back('S'); buff.push_back('s');
        buff.push_back('Y'); buff.push_back('y');
        buff.push_back('K'); buff.push_back('k');
        buff.push_back('A'); buff.push_back('a');
        buff.push_back('C'); buff.push_back('c');
        buff.push_back('G'); buff.push_back('g');
        buff.push_back('T'); buff.push_back('t');
        buff.push_back('?'); buff.push_back('?');
        buff.push_back('N'); buff.push_back('n');
        buff.push_back('-');
        cc.setSons(buff);

         for (unsigned int i=0; i<align.ns(); i++) for (unsigned int j=0; j<align.ls(); j++) if (!cc.has(align.sequence(i)[j])) throw EggRuntimeError("Consensus: strange character");
         return true;
    }

    const std::vector<std::string>& Consensus::firstSequenceNames() {
        return t_firstSequenceNames;
    }

    const std::vector<std::string>& Consensus::secondSequenceNames() {
        return t_secondSequenceNames;
    }
          
    const std::vector<std::string>& Consensus::roots() {
        return t_roots;
    }
          
    const std::vector<int>& Consensus::consistentPositions() {
        return t_consistentPositions;
    }

    const std::vector<int>& Consensus::complementaryPositions() {
        return t_complementaryPositions;
    }

    const std::vector<int>& Consensus::uninformativePositions() {
        return t_uninformativePositions;
    }
          
    const std::vector<int>& Consensus::ambiguousPositions() {
        return t_ambiguousPositions;
    }

    const std::vector<int>& Consensus::atLeastPartiallyResolvedAmbiguities() {
        return t_atLeastPartiallyResolvedAmbiguities;
    }

    const std::vector<std::vector<int> >& Consensus::inconsistentPositions() {
        return t_inconsistentPositions;
    }





}
