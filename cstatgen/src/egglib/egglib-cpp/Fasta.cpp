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

#include "Fasta.hpp"
#include "EggException.hpp"
#include <sstream>
#include <fstream>
#include <cstdlib>


namespace egglib {


    Container Fasta::parsef(const char* fname, bool importGroupLabels) {
        std::ifstream stream(fname);
        if (!stream.is_open()) {
            throw EggOpenFileError(fname);
        }
        return parse(stream, importGroupLabels);
    }


    void Fasta::parsef(const char* fname, Container& container, bool importGroupLabels) {
        std::ifstream stream(fname);
        if (!stream.is_open()) {
            throw EggOpenFileError(fname);
        }
        parse(stream, container, importGroupLabels);
    }


    Container Fasta::parse(const std::string& str, bool importGroupLabels) {
        std::istringstream stream(str);
        return parse(stream, importGroupLabels);
    }


    void Fasta::parse(const std::string& str, Container& container, bool importGroupLabels) {
        std::istringstream stream(str);
        parse(stream, container, importGroupLabels);
    }


    Container Fasta::parse(std::istream& stream, bool importGroupLabels) {
        Container container;
        parse(stream, container, importGroupLabels);
        return container;
    }


    void Fasta::parse(std::istream& stream, Container& container, bool importGroupLabels) {
        if (!stream.good()) throw EggRuntimeError("cannot import fasta data: invalid stream");


        // ignores everything until the first >
        char c=' ';
        while (c!='>') {
            c= stream.get();
            if (stream.eof()) {
                // accepting empty files
                return;
            }
        }

        // main loop over sequences (each starts after having read a >)
        while (stream.good()) {

            std::string name;
            std::string sequence;
            std::string label;
            bool Qlabel=false;

            // gets the name
            c= stream.get();
            while (c!='\n') {
                if (stream.eof()) throw EggFormatError("string", "fasta format", "interrupted string");

                bool ignore= false;
                if (importGroupLabels) {
                    if (Qlabel) {
                        if (c==' ') {
                            Qlabel=false;
                            ignore=true;
                        }
                        else label.push_back(c);
                    }
                    else if (c=='@') {
                        Qlabel=true;
                    }
                }
                if (c!='\r' && (!importGroupLabels || (!Qlabel && !ignore))) name.push_back(c);
                c= stream.get();
            }

            name.push_back('\0');

            // tries to get the group index string
            int group = 0;
            if (importGroupLabels) group= atoi(label.c_str());

            // gets the sequence
            c = stream.get();
            while (c!='>' && !stream.eof()) {
                char d = (c!='\n' && c!='\r' && c!='\r')?c:'\0';
                if (d) sequence.push_back(d);
                c= stream.get();
            }
            sequence.push_back('\0');

            container.append(name.c_str(), sequence.c_str(), group);
        }
    }


    void Fasta::formatf(const char* fname, const Container& container, bool exportGroupLabels, unsigned int lineLength) {
        std::ofstream stream(fname);
        if (!stream.is_open()) {
            throw EggOpenFileError(fname);
        }
        format(stream, container, exportGroupLabels, lineLength);
    }


    std::string Fasta::format(const Container& container, bool exportGroupLabels, unsigned int lineLength) {
        std::ostringstream stream;
        format(stream, container, exportGroupLabels, lineLength);
        return stream.str();
    }


    void Fasta::format(std::ostream& stream, const Container& container, bool exportGroupLabels, unsigned int lineLength) {
        for (unsigned i=0; i<container.ns(); i++) {
            stream << ">" << container.name(i);
            if (exportGroupLabels) stream << "@" << container.group(i);
            stream << std::endl;
            unsigned c = 0;
            for (unsigned int j=0; j<container.ls(i); j++) {
                stream << container.get(i, j);
                if (lineLength!=0) {
                    c++;
                    if (c==lineLength) {
                        c = 0;
                        stream << std::endl;
                    }
                }
            }
            if (c!=0 || lineLength==0) stream << std::endl;
        }
    }

}

