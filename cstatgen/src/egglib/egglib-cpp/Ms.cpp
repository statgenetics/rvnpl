/*
    Copyright 2008,2009,2011 Stéphane De Mita and Mathieu Siol

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

#include "Ms.hpp"
#include "EggException.hpp"
#include <sstream>
#include <cstdlib>


namespace egglib {

    double Ms::_tMRCA = -1.;
    double Ms::_prob = -1.;
    std::string Ms::_trees = "";
    

/* PARSER *************************************************************/

    DataMatrix Ms::get(std::string str, unsigned int ns, bool separated) {

        std::istringstream stream(str);
        return get(stream, ns, separated);
    }
    
    
    DataMatrix Ms::get(std::istream& stream, unsigned int ns, bool separated) {

        _tMRCA = -1.;
        _prob = -1.;
        _trees = "";

        // simulation start symbol
        std::string line;
        while ( (line=="" || line=="\r") && stream.good() ) {
            getline(stream, line);
        }
        
        if (!stream.good()) throw EggFormatError("string", "ms format", "cannot detect start of simulation (//)");

        if (line!="//" && line!="//\r") {
            std::string message =  "the following line was found where the start of a simulation (//) is expected: " + line;
            throw EggFormatError("string", "ms format", message.c_str());
        }

        getline(stream, line);

        std::istringstream sstream("");
        unsigned int S = 0;

        while (true) {

            // reads the tMRCA if there
            if ( ! line.compare(0, 6, "time:\t") ) {
                if (_tMRCA!=-1) throw EggFormatError("string", "ms format", "tMRCA found multiple times!");
                std::string times = line.substr(6);
                sstream.clear();
                sstream.str(times);
                sstream >> _tMRCA;
                getline(stream, line);
                continue;
            }

            // reads "prob" if there
            if ( ! line.compare(0, 6, "prob: ") ) {
                if (_prob!=-1) throw EggFormatError("string", "ms format", "prob found multiple times!");
                std::string prob = line.substr(6);
                sstream.clear();
                sstream.str(prob);
                sstream >> _prob;
                getline(stream, line);
                continue;
            }

            // reads the unique tree if there
            if (line.size()>3 && line[0]=='(' && line[line.size()-1]==';') {
                if (_trees!="") throw EggFormatError("string", "ms format", "tree found multiple times!");
                _trees = line;
                getline(stream, line);
                continue;
            }

            // reads a recombined tree (there may be several)
            if (line.size()>3 && line[0]=='[' && line[line.size()-1]==';') {
                unsigned int pos = line.find("]");
                if (pos >= line.size()) throw EggFormatError("string", "ms format", "invalid tree line!");
                if (line[pos+1]!='(') throw EggFormatError("string", "ms format", "invalid tree line!");
                _trees = _trees + line.substr(pos+1);
                getline(stream, line);
                continue;
            }

            // gets to the line with the number of segregating sites (has to be current line)
            if (line.compare(0, 10, "segsites: ")) {
                std::string message =  "the following line was found instead of segsites: " + line;
                throw EggFormatError("string", "ms format", message.c_str());
            }

            // extracts the number of segregating sites
            std::string ssegr = line.substr(10);
            sstream.clear();
            sstream.str(ssegr);
            sstream >> S;

            break;
        }

        // initiates the data matrix
        DataMatrix dataMatrix(ns, S);
        
        // from ms, if the last simulation is empty, lacks a white space
        // for this reason we just escape now if S = 0
        if (S==0) return dataMatrix;
        
        // gets positions
        getline(stream, line);
        if (line.compare(0, 11, "positions: ")) {
            std::string message =  "the following line was found instead of positions: " + line;
            throw EggFormatError("string", "ms format", message.c_str());
        }

        // gets positions
        sstream.clear();
        sstream.str(line.substr(11));
        double position = 0.;
        for (unsigned int i=0; i<S; i++) {
            sstream >> position;
            if (!sstream.good()) throw EggFormatError("string", "ms format", "error while parsing positions line");
            dataMatrix.sitePosition(i, position);
        }

        // gets sequences
        std::string dest= "@";
        std::string valid = "0123456789";
        unsigned int i = 99;
        for (unsigned int n=0; n<ns; n++) {

            if (!stream.good()) throw EggFormatError("string", "ms format", "error while parsing sequences, not enough genotypes");
            
            // prepares reading the line
            getline(stream, line);
            
            sstream.clear();
            sstream.str(line);
            
            // read each value
            for (unsigned int s=0; s<S; s++) {
                
                // reads normal data
                if (!separated) {
                    if (!sstream.good()) throw EggFormatError("string", "ms format", "error while parsing sequences, not enough sites");
                    dest[0] = sstream.get();
                    
                    if (valid.find(dest[0])==std::string::npos) {
                        std::string message = "invalid character found in sequences: " + dest;
                        throw EggFormatError("string", "ms format", message.c_str());
                    }
                    i = atoi(dest.c_str());
                    dataMatrix.set(n, s, i);
                }
                
                // reads data separated with a white space
                else {
                    if (!sstream.good()) throw EggFormatError("string", "ms format", "error while parsing sequences, not enough sites");
                    sstream >> dest;
                    for (std::string::iterator it=dest.begin(); it!=dest.end(); it++) {
                        if (valid.find(*it)==std::string::npos) {
                            std::string message = "invalid genotype: " + dest;
                            throw EggFormatError("string", "ms format", message.c_str());
                        }
                    }
                    i = atoi(dest.c_str());
                    dataMatrix.set(n, s, i);
                }
            }
            
            // checks that we have nothing left more than carriage returns
            dest[0] = sstream.get();
            if (sstream.good() && dest[0]!='\n' && dest[0]!='\r') {
                std::string message = "sequence line longer than expected: " + sstream.str();
                throw EggFormatError("string", "ms format", message.c_str());
            }
        }
        
        // checks trailing white line
        getline(stream, line);
        if (line!="" && line!="\r") {
            std::string message = "the following line where found instead of the white line expected after simulations: " + line;
            throw EggFormatError("string", "ms format", message.c_str());
        }
        
        // returns the data
        return dataMatrix;
        
    }



    /* FORMATTER **********************************************************/

    std::string Ms::format(DataMatrix& dataMatrix, bool separated) {
        std::ostringstream stream;
        format(stream, dataMatrix, separated);
        return stream.str();
    }
            
            

    void Ms::format(std::ostream& stream, DataMatrix& dataMatrix, bool separated) {

        // checks the stream
        if (!stream.good()) throw EggRuntimeError("ms formatter used with an invalid (unwritable) stream");

        // simulation header
        stream << "//" << std::endl;
        stream << "segsites: " << dataMatrix.numberOfSites() << std::endl;

        // data printed only if at least one site
        if (dataMatrix.numberOfSites()>0) {
            
            // prints positions lines
            stream << "positions:";
            for (unsigned int s=0; s<dataMatrix.numberOfSites(); s++) {
                stream << " " << dataMatrix.sitePosition(s);
            }
            stream << std::endl;
            
            // prints each sequence
            for (unsigned int n=0; n<dataMatrix.numberOfSequences(); n++) {
                for (unsigned int s=0; s<dataMatrix.numberOfSites(); s++) {

                    // checks the stream (again)
                    if (!stream.good()) throw EggRuntimeError("ms formatter used with an invalid (unwritable) stream");
                    
                    int value = dataMatrix.fget(n, s);
                    
                    // checks the data can be exported 
                    if (value<0 || (separated==false && dataMatrix.get(n,s)>9)) {
                        std::ostringstream errorStream;
                        errorStream << "this genotype value cannot be exported: " << value;
                        std::string errorMessage = errorStream.str();
                        throw EggFormatError("DataMatrix object", "data for ms format", errorMessage.c_str());
                    }
                    
                    // exports the data
                    stream << value;
                    
                    // add a separator if need
                    if (separated && s!=(dataMatrix.numberOfSites()-1)) {
                        stream << " ";
                    }
                }
                
                // finishing the genotype line
                stream << std::endl;
            }
        }
        
        // concludes the string with a white line (whatever S)
        stream << std::endl;
        
        // checks the stream (last time)
        if (!stream.good()) throw EggRuntimeError("ms formatter used with an invalid (unwritable) stream");
    }


    double Ms::tMRCA() {
        return _tMRCA;
    }

    double Ms::prob() {
        return _prob;
    }

    std::string Ms::trees() {
        return _trees;
    }
}

