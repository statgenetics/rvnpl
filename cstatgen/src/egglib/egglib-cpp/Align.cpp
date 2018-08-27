/*
    Copyright 2008-2009 Stéphane De Mita, Mathieu Siol
    
    This file is part of EggLib.

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


#include "Align.hpp"
#include <cstdlib>
#include <cstring>
#include "EggException.hpp"



namespace egglib {


    // Constructors / destructor / assignment operators

    Align::Align() {
        init();
    }


    Align::Align(const Align& s) { 
        init();
        copyObject(s);    
    }


    Align::Align(const Container& s) { 
        init();
        copyObject(s);
    }


    Align::Align(unsigned int number_of_sequences, unsigned int alignment_length, const char* const * const cstring_array) {
        init();
        setFromSource(number_of_sequences, alignment_length, cstring_array);
    }


    Align::Align(unsigned int number_of_sequences, unsigned int alignment_length) {
        init();
        setFromSource(number_of_sequences, alignment_length, NULL);
    }
    

    Align::~Align() { clear(); }


    Align& Align::operator=(const Container& s) {
        init();
        copyObject(s);
        return *this;
    }


    Align& Align::operator=(const Align& s) {
        init();
        copyObject(s);
        return *this;
    }



    // Constructors/destructor helps

    void Align::init(){
        _ns= 0;
        _ls= 0;
        lnames= NULL;
        names= NULL;
        sequences= NULL;
        groups= NULL; 
    } 


    void Align::copyObject(const Container& s) {
        if (_ns) clear();
        _ns= s.ns();
        if (!(lnames= (unsigned int*) realloc(lnames, _ns* sizeof(unsigned int)))) throw EggMemoryError();
        if (!(groups= (unsigned int*) realloc(groups, _ns* sizeof(unsigned int)))) throw EggMemoryError();
        if (!(sequences= (char**) realloc(sequences, _ns* sizeof(char*)))) throw EggMemoryError();
        if (!(names= (char**) realloc(names, _ns* sizeof(char*)))) throw EggMemoryError();
        for (unsigned int i=0; i<s.ns(); i++) if ((s.ls(i))>_ls) _ls= s.ls(i);
        for (unsigned int i=0; i<_ns; i++) {
            if (!(sequences[i]= (char*) malloc((_ls+1)* sizeof(char)))) throw EggMemoryError();
            strncpy(sequences[i], s.sequence(i), _ls);
            sequences[i][_ls]='\0';
            lnames[i]= strlen(s.name(i));
            if (!(names[i]= (char*) malloc((lnames[i]+1)* sizeof(char)))) throw EggMemoryError();
            strcpy(names[i], s.name(i));
        }
        getNamesAndGroups(s);
        for (unsigned int i=0; i<_ns; i++) for (unsigned int j=0; j<_ls; j++) if (!sequences[i][j]) sequences[i][j] = '?';
    }


    void Align::copyObject(const Align& s) {
        setFromSource(s._ns, s._ls, s.sequences);
        getNamesAndGroups(s);
    }


    void Align::setFromSource(unsigned int number_of_sequences, unsigned int alignment_length, const char* const * const cstring_array) {
        if (_ns) clear();
        _ns= number_of_sequences;
        _ls= alignment_length;
        if (!(lnames= (unsigned int*) realloc(lnames, _ns* sizeof(unsigned int)))) throw EggMemoryError();
        if (!(groups= (unsigned int*) realloc(groups, _ns* sizeof(unsigned int)))) throw EggMemoryError();
        if (!(names= (char**) realloc(names, _ns* sizeof(char*)))) throw EggMemoryError();
        if (!(sequences= (char**) realloc(sequences, _ns* sizeof(char*)))) throw EggMemoryError();

        for (unsigned int i=0; i<_ns; i++) {
            lnames[i]= 0;
            if (!(names[i]= (char*) malloc((lnames[i]+1)* sizeof(char)))) throw EggMemoryError();
            names[i][0]= '\0';
            groups[i]= 0;
            if (!(sequences[i]= (char*) malloc((_ls+1)* sizeof(char)))) throw EggMemoryError();
            if (cstring_array) strncpy(sequences[i], cstring_array[i], _ls);
            else for (unsigned int j=0; j<_ls; j++) sequences[i][j] = '?';
            sequences[i][_ls]='\0';
        }
    }


    void Align::clear() {
        for (unsigned int i=0; i<_ns; i++) {
            if (names[i]) free(names[i]);
            if (sequences[i]) free(sequences[i]);
        }
        if (lnames) free(lnames);
        if (names) free(names);
        if (groups) free(groups);
        if (sequences) free(sequences);
        _ns= 0;
        _ls= 0;
        lnames= NULL;
        names= NULL;
        sequences= NULL;
        groups= NULL;
    }



    // Accessors

    unsigned int Align::ls() const {
        return _ls;
    }


    unsigned int Align::ls(unsigned int pos) const {
        if (pos>=_ns) throw EggArgumentValueError("cannot access length of Align sequence: invalid index");
        return _ls;
    }


    char Align::get(unsigned int sequence, unsigned int position) const {

        if (sequence>=_ns) throw EggArgumentValueError("cannot access Align data: invalid sequence index");
        if (position>=_ls) throw EggArgumentValueError("cannot access Align data: invalid position index");
        return sequences[sequence][position];
    }



    // Basic modifiers (replacement)

    void Align::set(unsigned int sequence, unsigned position, char ch) {
        if (sequence>=_ns) throw EggArgumentValueError("cannot set Align data: invalid sequence index");
        if (position>=_ls) throw EggArgumentValueError("cannot set Align data: invalid position index");
        sequences[sequence][position] = ch;
    }


    void Align::sequence(unsigned int i, const char* sequence) {
        if (i>=_ns) throw EggArgumentValueError("cannot set Align sequence: invalid index");
        if (strlen(sequence)!=_ls) throw EggUnalignedError();
        strcpy(sequences[i], sequence);
    }


    void Align::binSwitch(unsigned int p) {
        if (p>=_ls) throw EggArgumentValueError("cannot switch Align column states: invalid position index");
        for (unsigned int i=0; i<_ns; i++) 
            switch (sequences[i][p]) {
                case '0':
                    sequences[i][p]= '1';
                    break;
                case '1':
                    sequences[i][p]= '0';
                    break;
                default:
                    throw EggRuntimeError("tried to switch non-binary data in Align");
            }
    }



    // Deep modifiers (change the size of the data matrix)

    unsigned int Align::append(const char* name, const char* sequence, unsigned int group) {
        if (_ls==0) _ls= strlen(sequence);
        else if (_ls!=strlen(sequence)) throw EggUnalignedError();

        _ns++;

        // (re)allocates the arrays
        if (!(lnames= (unsigned int*) realloc(lnames, _ns* sizeof(unsigned int)))) throw EggMemoryError();
        if (!(groups= (unsigned int*) realloc(groups, _ns* sizeof(unsigned int)))) throw EggMemoryError();
        if (!(names= (char**) realloc(names, _ns* sizeof(char*)))) throw EggMemoryError();
        if (!(sequences= (char**) realloc(sequences, _ns* sizeof(char*)))) throw EggMemoryError();

        // allocates the space for the first sequence
        if (!(names[_ns-1]= (char*) malloc((strlen(name)+1)* sizeof(char)))) throw EggMemoryError();
        if (!(sequences[_ns-1]= (char*) malloc((_ls+1)* sizeof(char)))) throw EggMemoryError();

        // fills the arrays
        lnames[_ns-1]= strlen(name);
        groups[_ns-1]= group;

        // name and sequence
        for (unsigned int i=0; i<lnames[_ns-1]; i++) names[_ns-1][i]= name[i];
        for (unsigned int i=0; i<_ls; i++) sequences[_ns-1][i]= sequence[i];

        names[_ns-1][lnames[_ns-1]]= '\0';
        sequences[_ns-1][_ls]= '\0';

        return _ns;
    }


    unsigned int Align::remove(unsigned int pos) {
        if (pos>=_ns) throw EggArgumentValueError("cannot remove Align sequence: invalid index");

        // free the two slots
        if (names[pos]) free(names[pos]);
        if (sequences[pos]) free(sequences[pos]);

        // shifts all sequences with larger index
        for (unsigned int i=pos; i<(_ns-1); i++) {
            names[i]= names[i+1];
            lnames[i]= lnames[i+1];
            sequences[i]= sequences[i+1];
            groups[i]= groups[i+1];
        }

        // reallocates (down) the arrays
        _ns--;

        if (!_ns) {
            _ls=0;
            free(lnames);       lnames= NULL;
            free(names);        names= NULL;
            free(sequences);    sequences= NULL;
            free(groups);       groups= NULL;
        }
        else {
            if (!(lnames= (unsigned int*) realloc(lnames, _ns* sizeof(unsigned int)))) throw EggMemoryError();
            if (!(groups= (unsigned int*) realloc(groups, _ns* sizeof(unsigned int)))) throw EggMemoryError();
            if (!(names= (char**) realloc(names, _ns* sizeof(char*)))) throw EggMemoryError();
            if (!(sequences= (char**) realloc(sequences, _ns* sizeof(char*)))) throw EggMemoryError();
        }

        return _ns;
    }

    unsigned int Align::removePosition(unsigned int pos) {
        if (pos>=_ls) throw EggArgumentValueError("cannot remove Align position: invalid index");
        _ls--;
        for (unsigned int i=0; i<_ns; i++) {
            for (unsigned j=pos; j<_ls; j++) sequences[i][j]= sequences[i][j+1];
            if (!(sequences[i]= (char*) realloc(sequences[i], (_ls+1)* sizeof(char)))) throw EggMemoryError();
            sequences[i][_ls] = '\0';
        }
        return _ls;
    }



    // Partial copies

    Align Align::vslice(std::vector<unsigned int> sites) {
        Align sa;
        for (unsigned int i=0; i<_ns; i++) {
            std::string str(sites.size(), '?');
            for (unsigned int j=0; j<sites.size(); j++) str[j] = sequences[i][j];
            sa.append(names[i], str.c_str(), groups[i]);
        }
        return sa;
    }


    Align Align::vslice(unsigned int a, unsigned int b) {
        Align sa;
        if (a>_ls) a=_ls;
        if (b<=a) b=a;
        if (b>_ls) b=_ls;
        for (unsigned int i=0; i<_ns; i++) {
            std::string str(b-a, '?');
            for (unsigned int j=a; j<b; j++) str[j-a] = sequences[i][j];
            sa.append(names[i], str.c_str(), groups[i]);
        }
        return sa;
    }

}
