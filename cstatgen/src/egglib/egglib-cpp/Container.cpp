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

#include "Container.hpp"
#include "EggException.hpp"
#include <cstdlib>
#include <cstring>



namespace egglib {


    // Constructors and affiliated

    Container::Container() {
        init();    
    }


    Container::Container(unsigned int number_of_sequences, char const* const* const cstring_array) {
        init();
        setFromSource(number_of_sequences, cstring_array);
    }


    Container& Container::operator=(const Container& s) {
        init();
        copyObject(s);
        return *this;
    }


    Container::Container(const Container& s) {
        init();
        copyObject(s); 
    }


    Container::~Container() {
        clear();
    }



    // Helpers

    void Container::init(){
         _ns= 0;
        lnames= NULL;
        names= NULL;
        lsequences= NULL;
        sequences= NULL;
        groups= NULL;
     }


    void Container::copyObject(const Container& s) {
        setFromSource(s._ns, s.sequences);
        getNamesAndGroups(s);
    }


    void Container::getNamesAndGroups(const Container& s) {
        for (unsigned int i=0; i<_ns; i++) {
            groups[i]= s.groups[i];
            lnames[i]= strlen(s.names[i]);
            if (!(names[i]= (char*) realloc(names[i], (lnames[i]+1)* sizeof(char)))) throw EggMemoryError();
            strcpy(names[i],s.names[i]);
        }
    }


    void Container::setFromSource(unsigned int number_of_sequences, const char* const* const cstring_array) {
        if (_ns) clear();
        _ns= number_of_sequences;
        if (!(lnames= (unsigned int*) realloc(lnames, _ns* sizeof(unsigned int)))) throw EggMemoryError();
        if (!(lsequences= (unsigned int*) realloc(lsequences, _ns* sizeof(unsigned int)))) throw EggMemoryError();
        if (!(groups= (unsigned int*) realloc(groups, _ns* sizeof(unsigned int)))) throw EggMemoryError();
        if (!(names= (char**) realloc(names, _ns* sizeof(char*)))) throw EggMemoryError();
        if (!(sequences= (char**) realloc(sequences, _ns* sizeof(char*)))) throw EggMemoryError();
            
        for (unsigned int i=0; i<_ns; i++) {
            lnames[i]= 0;
            if (!(names[i]= (char*) malloc((lnames[i]+1)* sizeof(char)))) throw EggMemoryError();
            names[i][0]= '\0';
            groups[i]= 0;        
            lsequences[i]= strlen(cstring_array[i]);
            if (!(sequences[i]= (char*) malloc((lsequences[i]+1)* sizeof(char)))) throw EggMemoryError();
            strcpy(sequences[i], cstring_array[i]);
        }
    }


    void Container::clear() {
        for (unsigned int i=0; i<_ns; i++) {
            if (names[i]) free(names[i]);
            if (sequences[i]) free(sequences[i]);
        }
        if (lnames) free(lnames);
        if (lsequences) free(lsequences);
        if (names) free(names);
        if (groups) free(groups);
        if (sequences) free(sequences);
        _ns= 0;
        lnames= NULL;
        names= NULL;
        lsequences= NULL;
        sequences= NULL;
        groups= NULL;
    }



    // Heavy modifiers

    unsigned int Container::append(const char* name, const char* sequence, unsigned int group) {
        _ns++;

        // (re)allocates the arrays
        if (!(lnames= (unsigned int*) realloc(lnames, _ns* sizeof(unsigned int)))) throw EggMemoryError();
        if (!(lsequences= (unsigned int*) realloc(lsequences, _ns* sizeof(unsigned int)))) throw EggMemoryError();
        if (!(groups= (unsigned int*) realloc(groups, _ns* sizeof(unsigned int)))) throw EggMemoryError();
        if (!(names= (char**) realloc(names, _ns* sizeof(char*)))) throw EggMemoryError();
        if (!(sequences= (char**) realloc(sequences, _ns* sizeof(char*)))) throw EggMemoryError();

        // allocates the space for the first sequence
        if (!(names[_ns-1]= (char*) malloc((strlen(name)+1)* sizeof(char)))) throw EggMemoryError();
        if (!(sequences[_ns-1]= (char*) malloc((strlen(sequence)+1)* sizeof(char)))) throw EggMemoryError();
            
        // fills the arrays
        lnames[_ns-1]= strlen(name);
        lsequences[_ns-1]= strlen(sequence);
        groups[_ns-1]= group;
        
        // name and sequence
        for (unsigned int i=0; i<lnames[_ns-1]; i++) names[_ns-1][i]= name[i];
        for (unsigned int i=0; i<lsequences[_ns-1]; i++) sequences[_ns-1][i]= sequence[i];
        names[_ns-1][strlen(name)]= '\0';
        sequences[_ns-1][strlen(sequence)]= '\0';
        
        return _ns;
    }


    unsigned int Container::remove(unsigned int pos) {
        if (pos>=_ns) throw EggArgumentValueError("cannot remove sequence from Container: invalid index");

        // free the two slots
        if (names[pos]) free(names[pos]);
        if (sequences[pos]) free(sequences[pos]);

        // shifts all sequences with larger index
        for (unsigned int i=pos; i<(_ns-1); i++) {
            names[i]= names[i+1];
            lnames[i]= lnames[i+1];
            sequences[i]= sequences[i+1];
            lsequences[i]= lsequences[i+1];
            groups[i]= groups[i+1];
        }
        
        // reallocates (down) the arrays
        _ns--;
        if (!_ns) {
            free(lnames);       lnames= NULL;
            free(names);        names= NULL;
            free(sequences);    sequences= NULL;
            free(groups);       groups= NULL;
        }
        else {
            if (!(lnames= (unsigned int*) realloc(lnames, _ns* sizeof(unsigned int)))) throw EggMemoryError();
            if (!(lsequences= (unsigned int*) realloc(lsequences, _ns* sizeof(unsigned int)))) throw EggMemoryError();
            if (!(groups= (unsigned int*) realloc(groups, _ns* sizeof(unsigned int)))) throw EggMemoryError();
            if (!(names= (char**) realloc(names, _ns* sizeof(char*)))) throw EggMemoryError();
            if (!(sequences= (char**) realloc(sequences, _ns* sizeof(char*)))) throw EggMemoryError();
        }
        return _ns;
    }


    void Container::appendSequence(unsigned int pos, const char* sequence) {
        if (pos>=_ns) throw EggArgumentValueError("cannot append sequence string to Container; invalid index");
        unsigned int old= lsequences[pos];
        lsequences[pos] += strlen(sequence);
        if (!(sequences[pos]= (char*) realloc(sequences[pos], (lsequences[pos]+1)* sizeof(char)))) throw EggMemoryError();
        strcpy(sequences[pos]+old, sequence);
    }



    // simple accessors

    unsigned int Container::ns() const {
        return _ns;
    }

    
    unsigned int Container::ls(unsigned int pop) const {
        if (pop>=_ns) throw EggArgumentValueError("cannot access length of Container's sequence: invalid index");
        return lsequences[pop];
    }


    int Container::find(const char* query, bool strict)  const {
        if (strict) {
            for (unsigned int i=0; i<_ns; i++) if (!strcmp( names[i], query )) return i;
        }
        else {
            for (unsigned int i=0; i<_ns; i++) if (!strncmp( names[i], query, strlen(query) )) return i;
        }
        return -1;
    }


    const char* Container::name(unsigned int i) const {
        if (i>=_ns) throw EggArgumentValueError("cannot access name of Container's sequence: invalid index");
        return names[i];
    }


    const char* Container::sequence(unsigned int i) const {
        if (i>=_ns) throw EggArgumentValueError("cannot access Container's sequence: invalid index");
        return sequences[i];
    }


    unsigned int Container::group(unsigned int i) const {
        if (i>=_ns) throw EggArgumentValueError("cannot access group index of Container's sequence: invalid index");
        return groups[i];
    }


    char Container::get(unsigned int s, unsigned int p) const {
        if (s>=_ns) throw EggArgumentValueError("cannot access Container data: invalid sequence index");
        if (p>=lsequences[s]) throw EggArgumentValueError("cannot access Container data: invalid position index");
        return sequences[s][p];
    }



    // setters

    void Container::name(unsigned int i, const char* name) {
        if (i>=_ns) throw EggArgumentValueError("cannot set name of Container's sequence: invalid index");
        lnames[i] = strlen(name);
        if (!(names[i]= (char*) realloc(names[i], (lnames[i]+1)* sizeof(char)))) throw EggMemoryError();
        strcpy(names[i], name);
    }


    void Container::group(unsigned int i, unsigned int group) {
        if (i>=_ns) throw EggArgumentValueError("cannot set group index of Container's sequence: invalid index");
        groups[i] = group;
    }


    void Container::sequence(unsigned int i, const char* sequence) {
        if (i>=_ns) throw EggArgumentValueError("cannot set Container sequence: invalid index");
        lsequences[i] = strlen(sequence);
        if (!(sequences[i]= (char*) realloc(sequences[i], (lsequences[i]+1)* sizeof(char)))) throw EggMemoryError();
        strcpy(sequences[i], sequence);
    }


    void Container::set(unsigned int sequence, unsigned position, char ch) {
        if (sequence>=_ns) throw EggArgumentValueError("cannot set Container data: invalid sequence index");
        if (position>=lsequences[sequence]) throw EggArgumentValueError("cannot set Container data: invalid position index");
        sequences[sequence][position] = ch;
    }


    // specials

    bool Container::isEqual() const {
        for (unsigned int i=1; i<_ns; i++) {
            if (lsequences[i]!=lsequences[0]) return false;
        }
        return true;
    }


    unsigned int Container::equalize(char ch) {
        unsigned int i,j,ls = 0;
        for (i=0; i<_ns; i++)    if (lsequences[i] > ls) ls= lsequences[i];
        for (i=0; i<_ns; i++) {
            if (!(sequences[i]= (char*) realloc(sequences[i], (ls+1)* sizeof(char)))) throw EggMemoryError();
            for (j=lsequences[i]; j<ls; j++) sequences[i][j]= ch; // should overwrite the previous \0
            sequences[i][ls]='\0';
            lsequences[i] = ls;
        }
        return ls;
    }


    Container Container::hslice(unsigned int a, unsigned int b) const {
        if (a>_ns) a=_ns;
        if (b<=a) b=a;
        if (b>_ns) b=_ns;

        Container sc(b-a, sequences+a);
        for (unsigned i=0; i<sc.ns(); i++) sc.group(i, groups[a+i]);
        return sc;
    }

}
