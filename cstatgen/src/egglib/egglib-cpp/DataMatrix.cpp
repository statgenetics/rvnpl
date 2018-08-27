/*
    Copyright 2009-2010 Stéphane De Mita, Mathieu Siol

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

#include <cstdlib>
#include "DataMatrix.hpp"
#include "EggException.hpp"


namespace egglib {

    DataMatrix::DataMatrix() {
        init();
    }


    void DataMatrix::init() {
        _numberOfSequences = 0;
        _numberOfSites = 0;
        dataMatrix = NULL;
        positions = NULL;
        groups = NULL;
    }


    void DataMatrix::resize(unsigned int newNumberOfSequences, unsigned int newNumberOfSites) {
        // matrix dimensions
        if (!newNumberOfSequences) {
            clear();
            return;
        }
        unsigned int oldNumberOfSequences = _numberOfSequences;
        unsigned int oldNumberOfSites = _numberOfSites;
        
        // realloc number of sequences
        if (oldNumberOfSequences!=newNumberOfSequences) {
            _numberOfSequences = newNumberOfSequences;
            dataMatrix = (int**) realloc(dataMatrix, _numberOfSequences * sizeof(int*));
            if (dataMatrix==NULL) throw EggMemoryError();
            for (unsigned int i=oldNumberOfSequences; i<newNumberOfSequences; i++) {
                dataMatrix[i] = NULL;
            }
            groups = (unsigned int*) realloc(groups, _numberOfSequences * sizeof(unsigned int));
            if (groups==NULL) throw EggMemoryError();
        }

        // realloc number of sites
        if (oldNumberOfSites!=newNumberOfSites) {
            _numberOfSites = newNumberOfSites;
            if (!newNumberOfSites) {
                for (unsigned int i=0; i<_numberOfSequences; i++) {
                    if (dataMatrix[i]) free(dataMatrix[i]);
                    dataMatrix[i] = NULL;
                }
                if (positions) free(positions);
                positions=NULL;
            }
            else {
                for (unsigned int i=0; i<_numberOfSequences; i++) {
                    dataMatrix[i] = (int*) realloc(dataMatrix[i], _numberOfSites * sizeof(int));
                    if (dataMatrix[i]==NULL) throw EggMemoryError();
                }
                positions = (double*) realloc(positions, _numberOfSites * sizeof(double));
                if (positions==NULL) throw EggMemoryError();
            }
        }

        // fills the new cells
        for (unsigned int i=oldNumberOfSequences; i<newNumberOfSequences; i++) {
            for (unsigned int j=0; j<_numberOfSites; j++) {
                dataMatrix[i][j] = 0;
            }
            groups[i] = 0.;
        }
        for (unsigned int i=0; i<_numberOfSequences; i++) {
            for (unsigned int j=oldNumberOfSites; j<newNumberOfSites; j++) {
                dataMatrix[i][j] = 0;
            }
        }
        for (unsigned int i=oldNumberOfSites; i<newNumberOfSites; i++) {
            positions[i] = 0.;
        }
    }


    DataMatrix::DataMatrix(unsigned int numberOfSequences, unsigned int numberOfSites) {
        init();
        resize(numberOfSequences, numberOfSites);
    }


    DataMatrix::DataMatrix(const DataMatrix& source) {
        init();
        copy(source);
    }


    DataMatrix::DataMatrix(const CharMatrix& source) {
        init();
        copy(source);
    }


    DataMatrix& DataMatrix::operator=(const DataMatrix& source) {
        clear();
        copy(source);
        return *this;
    }


    DataMatrix& DataMatrix::operator=(const CharMatrix& source) {
        clear();
        copy(source);
        return *this;
    }


    DataMatrix::~DataMatrix() {
        clear();
    }


    void DataMatrix::copy(const CharMatrix& source) {
        resize(source.numberOfSequences(), source.numberOfSites());
        for (unsigned int i=0; i<_numberOfSequences; i++) {
            for (unsigned int j=0; j<_numberOfSites; j++) {
                set(i, j, source.character(i,j));
            }
            populationLabel(i, source.populationLabel(i));
        }
        for (unsigned int i=0; i<_numberOfSites; i++) {
            sitePosition(i, source.sitePosition(i));
        }
    }


    void DataMatrix::copy(const DataMatrix& source) {
        resize(source.numberOfSequences(), source.numberOfSites());
        for (unsigned int i=0; i<_numberOfSequences; i++) {
            for (unsigned int j=0; j<_numberOfSites; j++) {
                set(i, j, source.get(i,j));
            }
            populationLabel(i, source.populationLabel(i));
        }
        for (unsigned int i=0; i<_numberOfSites; i++) {
            sitePosition(i, source.sitePosition(i));
        }
    }
    
    
    unsigned int DataMatrix::numberOfSites() const {
        return _numberOfSites;
    }


    unsigned int DataMatrix::numberOfSequences() const {
        return _numberOfSequences;
    }

            
    void DataMatrix::set(unsigned int sequence, unsigned int site, int value) {
        if (sequence>=_numberOfSequences) throw EggArgumentValueError("invalid sequence index");
        if (site>=_numberOfSites) throw EggArgumentValueError("invalid site index");
        dataMatrix[sequence][site] = value;
    }


    int DataMatrix::get(unsigned int sequence, unsigned int site) const {
        if (sequence>=_numberOfSequences) throw EggArgumentValueError("invalid sequence index");
        if (site>=_numberOfSites) throw EggArgumentValueError("invalid site index");
        return dataMatrix[sequence][site];
    }


    void DataMatrix::sitePosition(unsigned int site, double value) {
        if (site>=_numberOfSites) throw EggArgumentValueError("invalid site index");
        positions[site] = value;
    }
        

    double DataMatrix::sitePosition(unsigned int site) const {
        if (site>=_numberOfSites) throw EggArgumentValueError("invalid site index");
        return positions[site];
    }


    void DataMatrix::populationLabel(unsigned int sequence, unsigned int value) {
        if (sequence>=_numberOfSequences) throw EggArgumentValueError("invalid sequence index");
        groups[sequence] = value;
    }
        

    unsigned int DataMatrix::populationLabel(unsigned int sequence) const {
        if (sequence>=_numberOfSequences) throw EggArgumentValueError("invalid sequence index");
        return groups[sequence];
    }


    void DataMatrix::clear() {
        if (positions) free(positions);
        if (groups) free(groups);
        if (dataMatrix) {
            for (unsigned int i=0; i<_numberOfSequences; i++) {
                if (dataMatrix[i]) free(dataMatrix[i]);
            }
            free(dataMatrix);
        }
        init();
    }



    void DataMatrix::shift(int minimum) {
        for (unsigned int s=0; s<_numberOfSites; s++) {
            
            // finds the offset between the smallest allele and the minimum
            int offset=0;
            for (unsigned int n=0; n<_numberOfSequences; n++) {
                if (dataMatrix[n][s]<minimum && (minimum-dataMatrix[n][s])>offset) {
                    offset = minimum-dataMatrix[n][s];
                }
            }
            
            // uses the offset
            for(unsigned int n=0; n<_numberOfSequences; n++) {
                dataMatrix[n][s]+=offset;
            }
        }
    }



    char DataMatrix::character(unsigned int sequence, unsigned int site) const {
        if (dataMatrix[sequence][site]<0 || dataMatrix[sequence][site]>9) {
            throw EggRuntimeError("DataMatrix cannot export genotype value as character");
        }
        return '0' + dataMatrix[sequence][site];
    }

}
