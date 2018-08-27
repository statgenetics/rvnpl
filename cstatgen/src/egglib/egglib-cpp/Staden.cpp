/*
    Copyright 2008,2009,2011 Stéphane De Mita, Mathieu Siol

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

#include "Staden.hpp"
#include "EggException.hpp"
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>

namespace egglib {
    
    int Staden::shift;
    Container Staden::container;
    int Staden::currpos;
    std::istream* Staden::stream;
    std::vector<std::string> Staden::ID;

/* HELPERS ************************************************************/

    void Staden::getShift() {
        // goes until the consensus sequence
        int i;
        std::streampos pos = stream->tellg();
        char c[] = "01234567890123456";
        while (strcmp(c,"        CONSENSUS")) {
            for(i=0;i<16;i++) c[i] = c[i+1];
            c[16] = stream->get();
            if (stream->eof()) {
                shift = -1;
                return;
            }
        }
        // goes until the last space
        shift = 17;
        while (c[0]==' ') {
            shift++;
            c[0] = stream->get();
        }
        shift--;

        // gets back to the initial place
        stream->seekg(pos);

        // forward until the next line (first sequence)
        while(c[0]!='\n') c[0] = stream->get();
    }


    char Staden::transforme(char c) {
        if (c==' ') return '?';
        if (c=='-') return 'N';
        if (c=='*') return '-';
        return c;
    }


    bool Staden::readOneSequence() {
        // is there a sequence here?
        if (stream->peek()=='\r') {
            stream->get();
            return false;
        } // I actually assume no \r except in \r\n
        
        if (stream->peek()=='\n') {
            return false;
        }
        
        // skip 8 chars -> stores IDs
        ID.resize(ID.size()+1);
        for (int i=0; i<8; i++) ID[ID.size()-1].push_back(stream->get());

        // read a taxon name
        char *name  = (char*) malloc(sizeof(char)* (shift-7));   // one more in reserve to avoid name[0] if shift is 8
        if (shift>7 && !name) throw EggMemoryError();
        stream->read(name,shift-8);
        name[shift-8] = '\0';
        std::string tax(name);
        while (tax[tax.length()-1]==' ') tax.erase(tax.length()-1);

        // read a line of sequence
        std::string seq;
        char c;
        while ((c=stream->peek())!='\n') {
            if (c=='\r') {
                stream->get();
                continue;
            }
            if (c==-1) EggFormatError("string", "Staden GAP4 dump file", "invalid format");
            seq.append(1,transforme(stream->get()));
        }
        stream->get();  // actually reads the line-ending backspace
        container.append(tax.c_str(),seq.c_str());
        if (name) free(name);
        return true;
    }


    bool Staden::readAppendOneSequence() {
        // is there a sequence here?
        if (stream->peek()=='\r') {
            stream->get();
            return false;
        } // I actually assume no \r except in \r\n
        
        if (stream->peek()=='\n' || stream->eof()) {
            return false;
        }
        
        // skip 8 chars, read ID
        std::string sID;
        for (int i=0; i<8; i++) {
            sID.push_back(stream->get());
        }
        
        // read a taxon name
        char *name  = (char*) malloc(sizeof(char)* (shift-7));   // one more in reserve to avoid name[0] if shift is 8
        if (shift>7 && !name) EggMemoryError();
        stream->read(name,shift-8);
        name[shift-8] = '\0';
        std::string tax(name);
        while (tax[tax.length()-1]==' ') tax.erase(tax.length()-1);

        // check if the taxon is in list
        unsigned int c = 0;
        while (c<ID.size()) if (ID[c]==sID) break; else c++;
        int rank = c;
        
        // if not in the list
        if (c==ID.size()) {
            std::string tp(currpos,'?');
            container.append(tax.c_str(),tp.c_str());
            ID.push_back(sID);
        }
        
        // add the string of sequence
        std::string seq;
        char ch;
        while ((ch=stream->peek())!='\n') {
            if (ch=='\r') {
                stream->get();
                continue;
            }
            if (ch==-1) EggFormatError("string", "Staden GAP4 dump file", "invalid format");
            seq.append(1,transforme(stream->get()));
        }

        stream->get();

        container.appendSequence(rank, seq.c_str());

        if (name) free(name);
        return true;
    }


    void Staden::undot(bool delete_consensus) {
        int pos = container.find("CONSENSUS");
        if (pos<0) return;
        for (unsigned int i=0; i<container.ns(); i++) {
            for (unsigned int j=0; j<strlen(container.sequence(i)); j++) {
                if (container.sequence(i)[j]=='.') container.set(i, j, container.sequence(pos)[j]);
            }
        }
        if (delete_consensus) container.remove(pos);
    }


/* PARSERS ************************************************************/

    Align Staden::parse(const std::string& string, bool deleteConsensus) {
        std::istringstream localStream(string.c_str());
        return parse(localStream, deleteConsensus);
    }


    Align Staden::parse(std::istream& inStream, bool delete_consensus) {
        stream = &inStream;
        container.clear();
        ID.clear();
        
        if (!stream->good()) {
            throw EggFormatError("stream", "Staden GAP4 dump file", "cannot read in stream");        
        }
        
        /* getShift gives the total number of characters before the start of sequences
                  and read until the next backspace */
        getShift(); // result stored in shift
        if (shift<0) {
            throw EggFormatError("string", "Staden GAP4 dump file", "strange file");
        }
        currpos = 0;

        // read sequences from first block
        while(readOneSequence());

        // readOneSequence() returns false when read a second \n avec the end of sequence line
        // skip indices line

        while(!stream->eof()) {
            currpos = strlen(container.sequence(0));
            stream->get();  // skip the blank line between blocks
            while (stream->get()!='\n') if (stream->eof()) break; // skip line of indices

            while(readAppendOneSequence());
//            container.equalize();
        }
        
        undot(delete_consensus);
//        Align align;
//        for(unsigned int i=0; i<container.ns(); i++) align.append(container.name(i), container.sequence(i), container.group(i));
        return container;
    }

}
